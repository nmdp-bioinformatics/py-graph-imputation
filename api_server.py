#!/usr/bin/env python3
"""
HLA Imputation API Server.

Serves single GL-string imputations backed by one or more frequency graphs,
and can GENERATE new graphs on demand for a chosen set of populations/loci
(the "factory"), auto-registering them for imputation.

Each graph is loaded ONCE (lazily on first use, or eagerly at startup) and
reused across requests. Populations and loci are auto-detected from the graph
files, so you only point the server at a config and it adapts.

Swagger UI:  http://<host>:<port>/apidocs

Run (single graph, eager load):
    python3 api_server.py --config conf/minimal-configuration.json --port 5000

Run (multiple graphs, lazy load):
    python3 api_server.py \
        --graph cau=conf/minimal-configuration.json \
        --graph wmda=conf/WMDA-minimal-configuration.json \
        --port 5000
"""

import argparse
import os
import logging
import threading

from flask import Flask, request, jsonify, redirect
from flask_cors import CORS
from flasgger import Swagger

from imputation_service import ImputationService
from graph_factory import GraphFactory

logger = logging.getLogger("api_server")

app = Flask(__name__)
CORS(app)

app.config["SWAGGER"] = {"title": "HLA Imputation API", "uiversion": 3}
swagger = Swagger(
    app,
    template={
        "swagger": "2.0",
        "info": {
            "title": "HLA Imputation API",
            "description": (
                "Impute phased/unphased HLA genotypes and populations from a GL "
                "string, and generate new frequency-backed graphs on demand.\n\n"
                "GL format: loci separated by '^', haplotypes by '+', allele "
                "ambiguity by '/'. Example: "
                "`A*01:01+A*02:01^B*07:02+B*08:01^DRB1*03:01+DRB1*15:01`"
            ),
            "version": "3.1.0",
        },
    },
)


class GraphRegistry:
    """Registers named graphs (name -> config file) and lazily builds services."""

    def __init__(self):
        self._configs = {}
        self._services = {}
        self._default = None
        # Guards _configs/_services membership checks and the per-name lock map.
        self._registry_lock = threading.Lock()
        # One build lock per graph name, so concurrent first-time requests for
        # the SAME graph build it once instead of racing to build it N times.
        self._build_locks = {}

    def register(self, name, config_file, make_default=False):
        with self._registry_lock:
            self._configs[name] = config_file
            if make_default or self._default is None:
                self._default = name

    @property
    def names(self):
        return list(self._configs.keys())

    @property
    def default_name(self):
        return self._default

    def is_loaded(self, name):
        return name in self._services

    def config_file(self, name):
        return self._configs.get(name)

    def _build_lock_for(self, name):
        """Return (creating once) the dedicated build lock for a graph name."""
        with self._registry_lock:
            lock = self._build_locks.get(name)
            if lock is None:
                lock = threading.Lock()
                self._build_locks[name] = lock
            return lock

    def get(self, name=None):
        if name in (None, ""):
            name = self._default
        if name not in self._configs:
            raise KeyError(name)
        # Fast path: already built (dict reads are atomic under the GIL).
        service = self._services.get(name)
        if service is not None:
            return service
        # Slow path: serialize builds of THIS graph on its own lock so only the
        # first thread builds it; the rest wait and reuse the cached instance.
        with self._build_lock_for(name):
            service = self._services.get(name)
            if service is None:
                service = ImputationService(self._configs[name])
                self._services[name] = service
            return service

    def load_all(self):
        for name in self._configs:
            self.get(name)

    def discover(self, graphs_root="graphs"):
        """
        Register any completed graph found under graphs_root/<name>/config.json
        (a graph is "completed" when its csv/nodes.csv exists). Does not build
        the graph; registration is lazy. Returns the list of newly registered
        names.
        """
        from pathlib import Path
        import json as _json

        newly = []
        root = Path(graphs_root)
        if not root.is_dir():
            return newly
        for cfg in sorted(root.glob("*/config.json")):
            name = cfg.parent.name
            if name in self._configs:
                continue
            # Only register if the graph actually finished generating.
            try:
                conf = _json.loads(cfg.read_text())
                csv_dir = Path(conf.get("graph_files_path", cfg.parent / "csv"))
                nodes = csv_dir / conf.get("node_csv_file", "nodes.csv")
            except Exception:  # noqa: BLE001
                nodes = cfg.parent / "csv" / "nodes.csv"
            if Path(nodes).is_file():
                self.register(name, str(cfg))
                newly.append(name)
        return newly


registry = GraphRegistry()
factory = GraphFactory(graphs_root="graphs")


def _register_generated(name, config_path):
    """Callback used by the factory once a graph finishes generating."""
    registry.register(name, config_path)


def _parse_request():
    """Pull imputation parameters from JSON body or query string."""
    if request.method == "POST" and request.is_json:
        data = request.get_json(silent=True) or {}
    else:
        data = request.args

    def _as_int(val):
        if val in (None, ""):
            return None
        try:
            return int(val)
        except (TypeError, ValueError):
            return None

    def _as_float(val):
        if val in (None, ""):
            return None
        try:
            return float(val)
        except (TypeError, ValueError):
            return None

    def _as_expand(val):
        # Accept: missing -> "auto"; true/false booleans or strings;
        # anything else falls back to "auto".
        if val is None or val == "":
            return "auto"
        if isinstance(val, bool):
            return val
        s = str(val).strip().lower()
        if s in ("true", "1", "yes", "on"):
            return True
        if s in ("false", "0", "no", "off"):
            return False
        if s == "auto":
            return "auto"
        return "auto"

    return {
        "glstring": data.get("glstring") or data.get("gl_string") or data.get("gl"),
        "pop1": data.get("pop1"),
        "pop2": data.get("pop2"),
        "output": data.get("output", "both"),
        "graph": data.get("graph"),
        "number_of_results": _as_int(data.get("number_of_results")),
        "number_of_pop_results": _as_int(data.get("number_of_pop_results")),
        "marginal_ratio_threshold": _as_float(data.get("marginal_ratio_threshold")),
        "expand": _as_expand(data.get("expand")),
    }


def _do_impute():
    """Shared imputation logic for both GET and POST /impute."""
    params = _parse_request()

    if not params["glstring"]:
        return jsonify({"error": "Missing required parameter: glstring"}), 400

    graph_name = params["graph"]
    if graph_name and graph_name not in registry.names:
        return (
            jsonify(
                {
                    "error": f"Unknown graph: {graph_name}",
                    "available_graphs": registry.names,
                }
            ),
            400,
        )

    try:
        service = registry.get(graph_name)
    except KeyError:
        logger.error("no graph available for request (graph=%r)", graph_name)
        return jsonify({"error": "No graph available"}), 503

    try:
        result = service.impute(
            gl_string=params["glstring"],
            pop1=params["pop1"],
            pop2=params["pop2"],
            output=params["output"],
            number_of_results=params["number_of_results"],
            number_of_pop_results=params["number_of_pop_results"],
            marginal_ratio_threshold=params["marginal_ratio_threshold"],
            expand=params["expand"],
        )
        result["graph"] = graph_name or registry.default_name
        return jsonify(result), 200
    except ValueError as e:
        logger.warning(
            "bad impute request (graph=%r, glstring=%r): %s",
            graph_name,
            params["glstring"],
            e,
        )
        return jsonify({"error": str(e)}), 400
    except Exception as e:  # noqa: BLE001
        logger.exception(
            "imputation failed (graph=%r, glstring=%r)",
            graph_name,
            params["glstring"],
        )
        return jsonify({"error": f"Imputation failed: {e}"}), 500


@app.route("/impute", methods=["GET"])
def impute_get():
    """
    Impute a single GL string (quick, query-string form).
    ---
    tags:
      - imputation
    parameters:
      - in: query
        name: glstring
        type: string
        required: true
        description: "GL string. Loci '^', haplotypes '+', ambiguity '/'."
        example: "A*01:01+A*02:01^B*07:02+B*08:01^DRB1*03:01+DRB1*15:01"
      - in: query
        name: graph
        type: string
        required: false
        description: Which registered graph to use (default = server default).
      - in: query
        name: pop1
        type: string
        required: false
      - in: query
        name: pop2
        type: string
        required: false
      - in: query
        name: output
        type: string
        enum: [umug, pmug, both]
        required: false
        description: "umug=unphased (internally MUUG), pmug=phased, both."
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
        description: >
          Optional (e.g. 0.001). Truncate ranked results at the first point where
          next_freq / cumulative_sum < threshold, then renormalize kept results
          into a probability that sums to 1.
      - in: query
        name: expand
        type: string
        enum: ["auto", "true", "false"]
        required: false
        description: >
          Expand ambiguous codes (XX, MAC, serology) to member alleles via
          py-ard before imputing. "auto" (default) expands only when such codes
          are detected; "true" always expands; "false" never expands. When an
          expansion happens the response includes "expanded": true and
          "input_gl_string" (the original), and "gl_string" holds the expanded
          form actually imputed.
    responses:
      200:
        description: Imputation result with ranked umug/pmug genotypes and populations.
      400:
        description: Missing/invalid parameters or unknown graph.
      503:
        description: No graph available.
    """
    return _do_impute()


@app.route("/impute", methods=["POST"])
def impute_post():
    """
    Impute a single GL string (JSON body form).
    ---
    tags:
      - imputation
    parameters:
      - in: body
        name: body
        required: true
        schema:
          type: object
          required: [glstring]
          properties:
            glstring:
              type: string
              example: "A*01:01+A*02:01^B*07:02+B*08:01^DRB1*03:01+DRB1*15:01"
            graph: {type: string}
            pop1: {type: string, example: "CAU"}
            pop2: {type: string, example: "CAU"}
            output: {type: string, enum: [umug, pmug, both]}
            number_of_results: {type: integer}
            number_of_pop_results: {type: integer}
            marginal_ratio_threshold: {type: number, example: 0.001}
            expand:
              description: >
                Expand ambiguous codes (XX, MAC, serology) via py-ard before
                imputing. "auto" (default) expands only when codes are detected;
                true always; false never.
              example: "auto"
    responses:
      200:
        description: Imputation result with ranked umug/pmug genotypes and populations.
      400:
        description: Missing/invalid parameters or unknown graph.
      503:
        description: No graph available.
    """
    return _do_impute()


@app.route("/graphs", methods=["GET"])
def list_graphs():
    """
    List registered graphs and their load state.
    ---
    tags:
      - graphs
    responses:
      200:
        description: Registered graphs with detected populations/loci when loaded.
    """
    graphs = []
    for name in registry.names:
        entry = {
            "name": name,
            "config_file": registry.config_file(name),
            "loaded": registry.is_loaded(name),
            "default": name == registry.default_name,
        }
        if registry.is_loaded(name):
            svc = registry.get(name)
            entry["populations"] = svc.populations
            entry["loci"] = list(svc.config["loci_map"].keys())
        graphs.append(entry)
    return jsonify({"graphs": graphs, "default": registry.default_name}), 200


@app.route("/graphs", methods=["POST"])
def create_graph():
    """
    Generate and register a new graph for the chosen populations/loci (the factory).
    ---
    tags:
      - graphs
    parameters:
      - in: body
        name: body
        required: true
        schema:
          type: object
          required: [populations]
          properties:
            name:
              type: string
              example: "cau_6loc"
              description: >
                Unique graph name (also used as its directory). Optional — if
                omitted, a name is auto-suggested from the populations and loci
                (e.g. ["CAU"] + 6 loci -> "cau_6loc").
            populations:
              type: array
              items: {type: string}
              example: ["CAU"]
            source:
              type: string
              enum: [hfx4, hfx9]
              description: >
                HFX frequency source. "hfx4" = 4-locus (A,B,C,DRB1), the
                default. "hfx9" = 9-locus (A,C,B,DRB3/4/5,DRB1,DQA1,DQB1,DPA1,
                DPB1) — large & slower to build.
            loci:
              type: integer
              enum: [4, 9]
              description: "Alternative to 'source': pick an HFX source by loci count."
            freq_trim_threshold:
              type: number
              description: >
                Expected-count cutoff (frequency x cohort size): 1 keeps
                haplotypes expected at least once in the sample, 5 keeps those
                expected >= 5x, etc. Higher = smaller/faster graph. Default 1.
    responses:
      202:
        description: Generation started in the background. Poll GET /graphs/{name}.
      400:
        description: Invalid request (bad populations, loci, or duplicate in progress).
    """
    if not request.is_json:
        return jsonify({"error": "Expected JSON body"}), 400
    body = request.get_json(silent=True) or {}

    name = body.get("name")
    populations = body.get("populations")
    source = body.get("source")  # e.g. "hfx4", "hfx9", "freqs5", "freqs6"
    loci = body.get("loci")  # optional; used to pick a source if source omitted
    freq_trim = body.get("freq_trim_threshold")  # None -> source-aware default

    if not populations or not isinstance(populations, list):
        return jsonify({"error": "populations (non-empty array) is required"}), 400

    # Determine the loci count for name suggestion (best-effort).
    try:
        loci_for_name = int(loci) if loci is not None else None
    except (TypeError, ValueError):
        loci_for_name = None

    # Auto-suggest a name when none is provided.
    if not name or not isinstance(name, str):
        name = factory.suggest_name(
            populations, loci_for_name, source=source, taken=registry.names
        )

    if name in registry.names:
        return jsonify({"error": f"Graph '{name}' already registered"}), 400

    try:
        job = factory.start(
            name=name,
            populations=populations,
            loci=loci_for_name,
            source=source,
            freq_trim_threshold=(float(freq_trim) if freq_trim is not None else None),
            on_ready=_register_generated,
        )
        return jsonify({"message": "Generation started", "job": job}), 202
    except ValueError as e:
        return jsonify({"error": str(e)}), 400


@app.route("/graphs/<name>", methods=["GET"])
def graph_status(name):
    """
    Get the status of a graph (generation job and/or registration).
    ---
    tags:
      - graphs
    parameters:
      - in: path
        name: name
        type: string
        required: true
    responses:
      200:
        description: Job status (pending|running|ready|error) and/or registration info.
      404:
        description: Unknown graph.
    """
    job = factory.job_status(name)
    registered = name in registry.names
    if job is None and not registered:
        return jsonify({"error": f"Unknown graph: {name}"}), 404
    return (
        jsonify(
            {
                "name": name,
                "registered": registered,
                "loaded": registry.is_loaded(name) if registered else False,
                "job": job,
            }
        ),
        200,
    )


@app.route("/graphs/register", methods=["POST"])
def register_graph():
    """
    Register an already-generated graph with the running server (no restart).
    ---
    tags:
      - graphs
    parameters:
      - in: body
        name: body
        required: true
        schema:
          type: object
          properties:
            name:
              type: string
              description: >
                Graph name. If config_path is omitted, resolves to
                graphs/<name>/config.json.
            config_path:
              type: string
              description: Path to the graph's config.json (optional if name given).
    responses:
      200:
        description: Graph registered (lazy-loaded on first impute).
      400:
        description: Missing/invalid name or config, or already registered.
      404:
        description: Config file not found.
    """
    from pathlib import Path

    body = request.get_json(silent=True) or {}
    name = body.get("name")
    config_path = body.get("config_path")

    if not name and not config_path:
        return jsonify({"error": "Provide 'name' and/or 'config_path'"}), 400

    if config_path and not name:
        name = Path(config_path).parent.name
    if name and not config_path:
        config_path = f"graphs/{name}/config.json"

    if name in registry.names:
        return jsonify({"error": f"Graph '{name}' already registered"}), 400
    if not Path(config_path).is_file():
        return jsonify({"error": f"Config not found: {config_path}"}), 404

    registry.register(name, config_path)
    return (
        jsonify(
            {"message": "Graph registered", "name": name, "config_path": config_path}
        ),
        200,
    )


@app.route("/graphs/rescan", methods=["POST"])
def rescan_graphs():
    """
    Scan the graphs directory and register any completed graphs not yet known.
    ---
    tags:
      - graphs
    responses:
      200:
        description: List of newly registered graph names.
    """
    newly = registry.discover(factory.graphs_root)
    return jsonify({"registered": newly, "graphs": registry.names}), 200


@app.route("/factory/sources", methods=["GET"])
def factory_sources():
    """
    List available frequency sources and the populations each offers.
    ---
    tags:
      - graphs
    responses:
      200:
        description: Frequency sources by loci count with available populations.
    """
    return jsonify({"sources": factory.available_sources()}), 200


@app.route("/factory/suggest-name", methods=["GET"])
def factory_suggest_name():
    """
    Suggest a graph name for a set of populations and loci.
    ---
    tags:
      - graphs
    parameters:
      - in: query
        name: populations
        type: string
        required: true
        description: Comma-separated population codes (e.g. "CAU,HIS").
      - in: query
        name: source
        type: string
        enum: [hfx4, hfx9]
        required: false
      - in: query
        name: loci
        type: integer
        enum: [4, 9]
        required: false
    responses:
      200:
        description: A suggested unique name.
      400:
        description: Missing populations.
    """
    raw = request.args.get("populations", "")
    populations = [p.strip() for p in raw.split(",") if p.strip()]
    if not populations:
        return jsonify({"error": "populations (comma-separated) is required"}), 400
    source = request.args.get("source")
    loci_arg = request.args.get("loci")
    try:
        loci = int(loci_arg) if loci_arg not in (None, "") else None
    except (TypeError, ValueError):
        loci = None
    name = factory.suggest_name(populations, loci, source=source, taken=registry.names)
    return (
        jsonify(
            {
                "suggested_name": name,
                "populations": populations,
                "source": source,
                "loci": loci,
            }
        ),
        200,
    )


@app.route("/health", methods=["GET"])
def health_check():
    """
    Health check.
    ---
    tags:
      - meta
    responses:
      200:
        description: Server readiness and loaded graphs.
    """
    return (
        jsonify(
            {
                "status": "healthy" if registry.names else "not_ready",
                "graphs": registry.names,
                "loaded": [n for n in registry.names if registry.is_loaded(n)],
                "default": registry.default_name,
            }
        ),
        200,
    )


@app.route("/", methods=["GET"])
def index():
    """Redirect root to the Swagger UI."""
    return redirect("/apidocs")


def _parse_graph_arg(value):
    if "=" not in value:
        raise argparse.ArgumentTypeError(
            f"--graph expects NAME=CONFIG_PATH, got: {value!r}"
        )
    name, config_file = value.split("=", 1)
    name, config_file = name.strip(), config_file.strip()
    if not name or not config_file:
        raise argparse.ArgumentTypeError(f"Invalid --graph value: {value!r}")
    return name, config_file


def configure_graphs(
    graph_specs=None,
    config=None,
    graphs_dir="graphs",
    discover=True,
    eager=False,
    eager_names=None,
    eager_limit=None,
):
    """Register graphs (named + default), auto-discover on-disk graphs, and
    optionally eager-load them. Shared by the CLI (`main`) and the module-import
    bootstrap used under gunicorn ``--preload`` so both paths behave identically.

    Eager loading here (in the gunicorn *parent*, before workers fork) builds
    each graph once; forked workers then share that memory copy-on-write instead
    of each rebuilding — and no request pays the cold-start cost.

    Limiting which graphs eager-load (the rest stay lazy — built on first use):
      - ``eager_names``: explicit list of graph names to load (most precise).
      - ``eager_limit``: cap the count; loads the first N registered graphs.
      - ``eager=True`` with neither: load ALL registered graphs (default).
    """
    factory.graphs_root = graphs_dir

    for gval in graph_specs or []:
        name, config_file = _parse_graph_arg(gval)
        registry.register(name, config_file)

    if config:
        registry.register("default", config, make_default=not graph_specs)

    if discover:
        discovered = registry.discover(graphs_dir)
        if discovered:
            logger.info("Discovered existing graphs: %s", discovered)

    if not registry.names:
        registry.register("default", "conf/minimal-configuration.json")

    # Decide which graphs to eager-load. Eager triggers if eager=True OR an
    # explicit subset was requested.
    if eager or eager_names or eager_limit:
        all_names = registry.names
        if eager_names:
            unknown = [n for n in eager_names if n not in all_names]
            if unknown:
                logger.warning("eager_names not registered (ignored): %s", unknown)
            to_load = [n for n in all_names if n in set(eager_names)]
        else:
            to_load = list(all_names)
        if eager_limit is not None:
            to_load = to_load[: max(0, int(eager_limit))]

        skipped = [n for n in all_names if n not in to_load]
        logger.info(
            "Eager loading %d of %d graph(s) before serving: %s",
            len(to_load),
            len(all_names),
            to_load,
        )
        if skipped:
            logger.info("Lazy (load on first request): %s", skipped)

        for name in to_load:
            svc = registry.get(name)
            logger.info(
                "graph %r ready: %d pops, %d loci",
                name,
                len(svc.populations),
                len(svc.config["loci_map"]),
            )


def _bootstrap_from_env():
    """Configure graphs at *module import* time from environment variables.

    This runs in the gunicorn master process when launched with ``--preload``
    (gunicorn imports the app module but never calls ``main()``), so setting
    ``GRIM_EAGER=1`` builds graphs once in the parent before workers fork.

    Environment variables:
      GRIM_EAGER      truthy -> run bootstrap (and eager-load graphs).
      GRIM_GRAPHS     comma-separated NAME=CONFIG specs (like repeated --graph).
      GRIM_CONFIG     single config registered as graph "default".
      GRIM_GRAPHS_DIR directory to auto-discover graphs from (default: graphs).
      GRIM_NO_DISCOVER truthy -> skip auto-discovery.
      GRIM_EAGER_GRAPHS comma-separated graph NAMES to eager-load (rest lazy).
      GRIM_EAGER_LIMIT  eager-load at most N graphs.
      GRIM_LOG_LEVEL  logging level (default: INFO).
    """
    if not os.environ.get("GRIM_EAGER"):
        return

    logging.basicConfig(
        level=getattr(
            logging, os.environ.get("GRIM_LOG_LEVEL", "INFO").upper(), logging.INFO
        ),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    graphs_env = os.environ.get("GRIM_GRAPHS", "").strip()
    graph_specs = [g for g in graphs_env.split(",") if g.strip()] if graphs_env else []

    eager_graphs_env = os.environ.get("GRIM_EAGER_GRAPHS", "").strip()
    eager_names = (
        [n.strip() for n in eager_graphs_env.split(",") if n.strip()]
        if eager_graphs_env
        else None
    )
    eager_limit_env = os.environ.get("GRIM_EAGER_LIMIT", "").strip()
    eager_limit = int(eager_limit_env) if eager_limit_env else None

    configure_graphs(
        graph_specs=graph_specs,
        config=os.environ.get("GRIM_CONFIG") or None,
        graphs_dir=os.environ.get("GRIM_GRAPHS_DIR", "graphs"),
        discover=not os.environ.get("GRIM_NO_DISCOVER"),
        eager=True,
        eager_names=eager_names,
        eager_limit=eager_limit,
    )


# Silence tqdm progress bars under any entry point (CLI or gunicorn import).
os.environ.setdefault("TQDM_DISABLE", "1")

# Run the import-time bootstrap (no-op unless GRIM_EAGER is set). This lets
# `gunicorn --preload` build graphs in the master before forking workers.
_bootstrap_from_env()


def main():
    parser = argparse.ArgumentParser(description="HLA Imputation API Server")
    parser.add_argument(
        "--config",
        default=None,
        help="Single config JSON (registered as graph 'default').",
    )
    parser.add_argument(
        "--graph",
        action="append",
        default=[],
        metavar="NAME=CONFIG",
        help="Register a named graph. Repeatable.",
    )
    parser.add_argument(
        "--eager",
        action="store_true",
        help="Load all registered graphs at startup (default: lazy).",
    )
    parser.add_argument(
        "--eager-graphs",
        default=None,
        metavar="NAME[,NAME...]",
        help="Comma-separated graph names to eager-load (the rest stay lazy). "
        "Implies eager for just those graphs.",
    )
    parser.add_argument(
        "--eager-limit",
        type=int,
        default=None,
        metavar="N",
        help="Eager-load at most N graphs (the first N registered). "
        "Combine with --eager or --eager-graphs.",
    )
    parser.add_argument("--port", type=int, default=5000)
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument(
        "--graphs-dir",
        default="graphs",
        help="Directory scanned for existing graphs to auto-register (default: graphs)",
    )
    parser.add_argument(
        "--no-discover",
        action="store_true",
        help="Do not auto-register completed graphs found in --graphs-dir.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level for the server/service (default: INFO). "
        "One of DEBUG, INFO, WARNING, ERROR, CRITICAL.",
    )
    args = parser.parse_args()

    # Configure application logging. Goes to stderr by default, so it is
    # captured alongside Flask's access log when you redirect the process
    # output to a file (e.g. `... > /tmp/api_server.log 2>&1`).
    logging.basicConfig(
        level=getattr(logging, str(args.log_level).upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    factory.graphs_root = args.graphs_dir

    eager_names = (
        [n.strip() for n in args.eager_graphs.split(",") if n.strip()]
        if args.eager_graphs
        else None
    )
    configure_graphs(
        graph_specs=args.graph,
        config=args.config,
        graphs_dir=args.graphs_dir,
        discover=not args.no_discover,
        eager=args.eager,
        eager_names=eager_names,
        eager_limit=args.eager_limit,
    )

    any_eager = args.eager or eager_names or args.eager_limit
    print("=" * 80)
    print("HLA Imputation API Server")
    for name in registry.names:
        marker = " (default)" if name == registry.default_name else ""
        loaded = " [loaded]" if registry.is_loaded(name) else " [lazy]"
        print(f"  graph: {name}{marker}{loaded} -> {registry.config_file(name)}")
    if not any_eager:
        print("Lazy loading: each graph builds on its first /impute request.")

    print(f"Swagger UI: http://{args.host}:{args.port}/apidocs")
    print(f"Serving on {args.host}:{args.port}")
    print("=" * 80)

    app.run(host=args.host, port=args.port, debug=args.debug, use_reloader=False)


if __name__ == "__main__":
    main()
