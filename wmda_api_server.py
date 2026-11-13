#!/usr/bin/env python3
"""
WMDA HLA Imputation API Server.

A dedicated REST service for the WMDA match-validation project. It serves ONE
graph ("wmda") built exclusively from the HFX frequency sets in WMDA_HFX_DIR
(default: /Users/mmaiers/src/wmda/match_validation/MV4/v2/HFX). Every
``<name>.hfx`` file in that directory (hfset-1 ... hfset-41) is one population
of the graph. No other frequency source is reachable through this service.

The graph is generated on first start (or with --rebuild) into
graphs_wmda/wmda/, then loaded once and reused across requests.

Swagger UI:  http://<host>:<port>/apidocs

Run:
    python3 wmda_api_server.py --port 5001
"""

import argparse
import logging
import os
from pathlib import Path

from flasgger import Swagger
from flask import Flask, jsonify, redirect, request
from flask_cors import CORS

from graph_factory import build_generation_config, generate_graph, hfx_loci_map
from graph_registry import GraphRegistry
from imputation_service import find_malformed_alleles

logger = logging.getLogger("wmda_api_server")

WMDA_HFX_DIR = Path(
    os.environ.get(
        "WMDA_HFX_DIR", "/Users/mmaiers/src/wmda/match_validation/MV4/v2/HFX"
    )
).expanduser()
GRAPHS_ROOT = os.environ.get("WMDA_GRAPHS_ROOT", "graphs_wmda")
GRAPH_NAME = "wmda"
# HFX bundles for this project report cohortSize 0, so each set counts as 1.0
# and the trim cutoff equals the threshold itself: haplotypes with freq below it
# are dropped. Across the 41 sets (2.26M haplotypes total) an untrimmed graph
# produced >11 GB of edges and was unusable; 1e-4 keeps ~59k haplotypes.
# Override with WMDA_FREQ_TRIM (e.g. 1e-5 keeps ~340k, slower/larger).
DEFAULT_TRIM = float(os.environ.get("WMDA_FREQ_TRIM", "1e-4"))

# Population priors: use the assigned pop1/pop2 STRICTLY. The prior matrix gets
# weight only at cell (pop1, pop2) via alpha; delta (same-pop diagonal), gamma
# (pop1-or-pop2 row/column), eta (uniform) and beta (identity) are all zero, so
# every haplotype pair outside the assigned populations is multiplied by 0.
# planb is turned off (see build_wmda_graph) because the engine's plan-B
# all-populations pass and plan C would otherwise discard this prior.
STRICT_PRIORITY = {"alpha": 1.0, "eta": 0, "beta": 0, "gamma": 0, "delta": 0}

app = Flask(__name__)
CORS(app)
app.config["SWAGGER"] = {"title": "WMDA HLA Imputation API", "uiversion": 3}
Swagger(
    app,
    template={
        "swagger": "2.0",
        "info": {
            "title": "WMDA HLA Imputation API",
            "description": (
                "Impute HLA genotypes from a GL string using only the WMDA "
                "HFX frequency sets (one graph, one population per set).\n\n"
                "GL format: loci '^', haplotypes '+', ambiguity '/'. Example: "
                "`A*01:01+A*02:01^B*07:02+B*08:01^C*07:01+C*07:02"
                "^DRB1*03:01+DRB1*15:01^DQB1*02:01+DQB1*06:02`"
            ),
            "version": "1.0.0",
        },
    },
)

registry = GraphRegistry()


def available_sets():
    """Names of the HFX frequency sets (populations) in WMDA_HFX_DIR."""

    def _key(n):
        head, _, tail = n.rpartition("-")
        return (head, int(tail)) if tail.isdigit() else (n, 0)

    if not WMDA_HFX_DIR.is_dir():
        return []
    return sorted((p.stem for p in WMDA_HFX_DIR.glob("*.hfx")), key=_key)


def config_path():
    return Path(GRAPHS_ROOT) / GRAPH_NAME / "config.json"


def graph_is_built():
    return (Path(GRAPHS_ROOT) / GRAPH_NAME / "csv" / "nodes.csv").is_file()


def build_wmda_graph():
    """Generate the single WMDA graph with every HFX set as a population."""
    sets = available_sets()
    if not sets:
        raise RuntimeError(f"No .hfx files found in {WMDA_HFX_DIR}")
    # The engine looks haplotypes up as "~".join(sorted(alleles)), i.e. loci in
    # alphabetical order. The HFX metadata lists A,B,C,DRB1,DQB1, so re-number
    # loci alphabetically or full-haplotype (plan A) lookups never match.
    loci = sorted(hfx_loci_map(WMDA_HFX_DIR / f"{sets[0]}.hfx"))
    loci_map = {locus: i + 1 for i, locus in enumerate(loci)}
    config, _graph_dir, _csv_dir = build_generation_config(
        name=GRAPH_NAME,
        populations=sets,
        loci_map=loci_map,
        freq_data_dir=str(WMDA_HFX_DIR),
        graphs_root=GRAPHS_ROOT,
        freq_trim_threshold=DEFAULT_TRIM,
        source_type="hfx",
    )
    config["priority"] = dict(STRICT_PRIORITY)
    config["UNK_priors"] = "SR"
    # Plan B (and plan C, which is gated on planb) relax/discard the prior
    # matrix when plan A finds nothing, so disable it to keep pop1/pop2 strict.
    config["planb"] = False
    logger.info("generating graph %r from %d HFX sets", GRAPH_NAME, len(sets))
    generate_graph(config, config_path())


def _service():
    return registry.get(GRAPH_NAME)


def _as_num(val, cast):
    """Parse a positive number; None/""/0/negative/garbage -> None (= default).

    Swagger UI pre-fills numeric fields with 0, which would otherwise mean
    "return zero results" and surface as a bogus "miss".
    """
    if val in (None, ""):
        return None
    try:
        num = cast(val)
    except (TypeError, ValueError):
        return None
    return num if num > 0 else None


def _as_expand(val):
    if isinstance(val, bool):
        return val
    s = str(val).strip().lower() if val is not None else ""
    if s in ("true", "1", "yes", "on"):
        return True
    if s in ("false", "0", "no", "off"):
        return False
    return "auto"


def _do_impute():
    data = (
        (request.get_json(silent=True) or {})
        if request.method == "POST" and request.is_json
        else request.args
    )
    glstring = data.get("glstring") or data.get("gl_string") or data.get("gl")
    if not glstring:
        return jsonify({"error": "Missing required parameter: glstring"}), 400
    bad = find_malformed_alleles(glstring)
    if bad:
        return (
            jsonify(
                {
                    "error": "Malformed allele(s) in glstring; expected LOCUS*FIELDS "
                    "like A*02:01",
                    "bad_alleles": bad,
                }
            ),
            400,
        )

    # Strict priors need an assigned population. pop2 defaults to pop1.
    sets = available_sets()
    pop1 = data.get("pop1")
    pop2 = data.get("pop2") or pop1
    if not pop1:
        return (
            jsonify(
                {
                    "error": "Missing required parameter: pop1 (a WMDA frequency "
                    "set, e.g. hfset-32). pop2 defaults to pop1.",
                    "available_sets": sets,
                }
            ),
            400,
        )
    for val in (pop1, pop2):
        if val not in sets:
            return (
                jsonify(
                    {"error": f"Unknown frequency set: {val}", "available_sets": sets}
                ),
                400,
            )

    try:
        result = _service().impute(
            gl_string=glstring,
            pop1=pop1,
            pop2=pop2,
            output=data.get("output", "both"),
            number_of_results=_as_num(data.get("number_of_results"), int),
            number_of_pop_results=_as_num(data.get("number_of_pop_results"), int),
            marginal_ratio_threshold=_as_num(
                data.get("marginal_ratio_threshold"), float
            ),
            expand=_as_expand(data.get("expand")),
        )
        result["graph"] = GRAPH_NAME
        return jsonify(result), 200
    except KeyError:
        return jsonify({"error": "WMDA graph not available"}), 503
    except ValueError as e:
        logger.warning("bad impute request (glstring=%r): %s", glstring, e)
        return jsonify({"error": str(e)}), 400
    except Exception as e:  # noqa: BLE001
        logger.exception("imputation failed (glstring=%r)", glstring)
        return jsonify({"error": f"Imputation failed: {e}"}), 500


@app.route("/impute", methods=["GET"])
def impute_get():
    """
    Impute a single GL string against the WMDA graph (query-string form).
    ---
    tags:
      - imputation
    parameters:
      - in: query
        name: glstring
        type: string
        required: true
        description: "GL string. Loci '^', haplotypes '+', ambiguity '/'."
        default: "A*01:01+A*02:01^B*07:02+B*08:01^C*07:01+C*07:02^DRB1*03:01+DRB1*15:01^DQB1*02:01+DQB1*06:02"
      - in: query
        name: pop1
        type: string
        required: true
        description: "WMDA frequency set the imputation is restricted to, e.g. hfset-32 (strict prior)."
      - in: query
        name: pop2
        type: string
        required: false
        description: "Second set for the other haplotype. Defaults to pop1."
      - in: query
        name: output
        type: string
        enum: [umug, pmug, both]
        required: false
      - in: query
        name: number_of_results
        type: integer
        required: false
      - in: query
        name: number_of_pop_results
        type: integer
        required: false
      - in: query
        name: marginal_ratio_threshold
        type: number
        required: false
      - in: query
        name: expand
        type: string
        enum: ["auto", "true", "false"]
        required: false
    responses:
      200:
        description: Ranked umug/pmug genotypes and populations (WMDA sets).
      400:
        description: Missing/invalid parameters or unknown frequency set.
      503:
        description: WMDA graph not available.
    """
    return _do_impute()


@app.route("/impute", methods=["POST"])
def impute_post():
    """
    Impute a single GL string against the WMDA graph (JSON body form).
    ---
    tags:
      - imputation
    parameters:
      - in: body
        name: body
        required: true
        schema:
          type: object
          required: [glstring, pop1]
          properties:
            glstring:
              type: string
              example: "A*01:01+A*02:01^B*07:02+B*08:01^C*07:01+C*07:02^DRB1*03:01+DRB1*15:01^DQB1*02:01+DQB1*06:02"
            pop1: {type: string, example: "hfset-32"}
            pop2: {type: string, example: "hfset-32"}
            output: {type: string, enum: [umug, pmug, both], example: "both"}
            number_of_results:
              type: integer
              example: 10
              description: "Omit or <= 0 for the server default (1000)."
            number_of_pop_results:
              type: integer
              example: 10
              description: "Omit or <= 0 for the server default (100)."
            marginal_ratio_threshold:
              type: number
              example: 0.001
              description: "Omit or <= 0 to disable truncation."
            expand: {type: string, enum: ["auto", "true", "false"], example: "auto"}
    responses:
      200:
        description: Ranked umug/pmug genotypes and populations (WMDA sets).
      400:
        description: Missing/invalid parameters or unknown frequency set.
      503:
        description: WMDA graph not available.
    """
    return _do_impute()


@app.route("/sets", methods=["GET"])
def list_sets():
    """
    List the WMDA HFX frequency sets (the populations of the graph).
    ---
    tags:
      - wmda
    responses:
      200:
        description: Frequency set names, source directory and graph state.
    """
    sets = available_sets()
    return (
        jsonify(
            {
                "graph": GRAPH_NAME,
                "freq_data_dir": str(WMDA_HFX_DIR),
                "count": len(sets),
                "sets": sets,
                "built": graph_is_built(),
                "loaded": registry.is_loaded(GRAPH_NAME),
            }
        ),
        200,
    )


@app.route("/health", methods=["GET"])
def health():
    """
    Health check.
    ---
    tags:
      - meta
    responses:
      200:
        description: Readiness of the WMDA graph.
    """
    return (
        jsonify(
            {
                "status": "healthy" if GRAPH_NAME in registry.names else "not_ready",
                "graph": GRAPH_NAME,
                "loaded": registry.is_loaded(GRAPH_NAME),
            }
        ),
        200,
    )


@app.route("/", methods=["GET"])
def index():
    return redirect("/apidocs")


def prepare(rebuild=False, eager=True):
    """Ensure the graph exists on disk, register it, optionally load it."""
    if rebuild or not graph_is_built():
        build_wmda_graph()
    registry.register(GRAPH_NAME, str(config_path()), make_default=True)
    if eager:
        svc = _service()
        logger.info(
            "graph %r ready: %d populations, loci %s",
            GRAPH_NAME,
            len(svc.populations),
            list(svc.config["loci_map"]),
        )


def main():
    parser = argparse.ArgumentParser(description="WMDA HLA Imputation API Server")
    parser.add_argument("--port", type=int, default=5001)
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Regenerate the WMDA graph from the HFX sets.",
    )
    parser.add_argument(
        "--lazy",
        action="store_true",
        help="Do not load the graph until the first request.",
    )
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    prepare(rebuild=args.rebuild, eager=not args.lazy)
    print(f"WMDA HLA Imputation API: {len(available_sets())} sets from {WMDA_HFX_DIR}")
    print(f"Swagger UI: http://{args.host}:{args.port}/apidocs")
    app.run(host=args.host, port=args.port, debug=False, use_reloader=False)


os.environ.setdefault("TQDM_DISABLE", "1")

if __name__ == "__main__":
    main()
