#!/usr/bin/env python3
"""
Graph factory (HFX-only).

Given a set of populations and an HFX frequency source, this builds a
self-contained imputation graph into its own output directory and writes a
matching config JSON. The config can then be registered with the API server so
imputation can be served against the new graph.

The generation pipeline:
    1. produce_hpf_from_hfx -> reads <POP>.hfx bundles via hfx_tools, writes
                               hpf.csv + population counts (real cohort sizes)
    2. generate_graph       -> reads hpf.csv, writes nodes/edges/top_links/info_node

Frequency sources (see GraphFactory.FREQ_SOURCES), all HFX bundles:
    - hfx4 : 4-locus (A,B,C,DRB1) at ~/src/HF/HFX/4loc            [default]
    - hfx9 : 9-locus (A,C,B,DRB3/4/5,DRB1,DQA1,DQB1,DPA1,DPB1) at ~/src/HF/HFX/9loc

Loci order/names and sample (cohort) sizes are read directly from each bundle's
metadata, so graphs always match the data and freq_trim_threshold behaves as an
expected-count cutoff (freq * cohortSize).

Everything for one graph lives under:  <graphs_root>/<name>/
    <name>/
      config.json
      csv/
        hpf.csv
        pop_counts_file.txt
        nodes.csv, edges.csv, top_links.csv, info_node.csv
"""

import csv
import json
import threading
import time
import traceback
from pathlib import Path

from graph_generation import generate_neo4j_multi_hpf

try:
    from hfx_tools.io import read_hfx_bundle

    HFX_AVAILABLE = True
except ImportError:  # pragma: no cover - hfx_tools optional at import time
    read_hfx_bundle = None
    HFX_AVAILABLE = False


def hfx_loci_map(hfx_path):
    """
    Read an HFX bundle's metadata and return a positional loci_map
    ({locus: 1-based index}) matching the order loci appear in its haplotypes.

    Prefers metadata.outputResolution ordering; falls back to parsing the first
    haplotype's '~'-separated alleles.
    """
    if not HFX_AVAILABLE:
        raise RuntimeError("hfx_tools is not installed; cannot read HFX bundles")
    bundle = read_hfx_bundle(Path(hfx_path))
    loci = [
        entry.get("locus")
        for entry in bundle.metadata.get("outputResolution", [])
        if entry.get("locus")
    ]
    if not loci and bundle.rows:
        loci = [a.split("*")[0] for a in bundle.rows[0][0].split("~")]
    if not loci:
        raise ValueError(f"Could not determine loci from HFX bundle: {hfx_path}")
    return {locus: i + 1 for i, locus in enumerate(loci)}


def _hfx_cohort_size(bundle):
    """
    Extract the sample size from an HFX bundle's metadata.

    Prefers metadata.cohortDescription.cohortSize; falls back to a
    population entry's populationSize. Returns None if unavailable.
    """
    cohort = bundle.metadata.get("cohortDescription", {}) or {}
    size = cohort.get("cohortSize")
    if size:
        return float(size)
    pop = cohort.get("population")
    if isinstance(pop, dict) and pop.get("populationSize"):
        return float(pop["populationSize"])
    if isinstance(pop, list):
        for entry in pop:
            if isinstance(entry, dict) and entry.get("populationSize"):
                return float(entry["populationSize"])
    return None


def produce_hpf_from_hfx(config):
    """
    HFX-aware replacement for generate_hpf.produce_hpf.

    Reads each population's <freq_data_dir>/<POP>.hfx bundle via hfx_tools and
    writes the same artifacts produce_hpf would: the hpf.csv (hap,pop,freq) at
    config["freq_file"] and the population counts file at
    config["pops_count_file"].

    HFX bundles carry the sample size in ``metadata.cohortDescription``
    (``cohortSize`` / per-population ``populationSize``). The graph generator
    applies the trim as ``freq < freq_trim_threshold / pop_count``, so using the
    real cohort size makes ``freq_trim_threshold`` an EXPECTED-COUNT cutoff, as
    originally designed: e.g. threshold ``1`` keeps haplotypes expected to occur
    at least once in the sample, ``5`` keeps those expected >= 5 times, etc.

    The third column of the counts file (the relative population ratio, used
    later as an imputation prior) is derived from the cohort sizes.
    """
    if not HFX_AVAILABLE:
        raise RuntimeError("hfx_tools is not installed; cannot read HFX bundles")

    pops = config["populations"]
    freq_data_dir = Path(config["freq_data_dir"])
    freq_file = Path(config["freq_file"])
    pops_count_file = Path(config["pops_count_file"])

    freq_file.parent.mkdir(parents=True, exist_ok=True)
    pops_count_file.parent.mkdir(parents=True, exist_ok=True)

    print("*" * 100)
    print("Conversion to HPF from HFX bundles:")
    print(f"\tPopulations: {pops}")
    print(f"\tHFX Source Dir: {freq_data_dir}")
    print(f"\tOutput File: {freq_file}")
    print("*" * 100)

    per_pop_rows = {}
    pop_count = {}
    for pop in pops:
        hfx_path = freq_data_dir / f"{pop}.hfx"
        if not hfx_path.is_file():
            raise FileNotFoundError(f"HFX bundle not found: {hfx_path}")
        print(f"Reading HFX bundle:\t {hfx_path}")
        bundle = read_hfx_bundle(hfx_path)
        rows = [(hap, freq) for hap, freq in bundle.rows if freq and freq > 0.0]
        per_pop_rows[pop] = rows
        pop_count[pop] = _hfx_cohort_size(bundle) or 1.0

    total_count = sum(pop_count.values()) or 1.0
    # count column = real cohort size so the generator's (freq_trim / count)
    # yields an expected-count cutoff. ratio column feeds imputation priors.
    with open(pops_count_file, "w") as pcf:
        for pop in pops:
            ratio = pop_count[pop] / total_count
            pcf.write(f"{pop},{pop_count[pop]},{ratio}\n")

    print(f"Writing hpf File:\t {freq_file}")
    with open(freq_file, "w", newline="") as csv_file:
        writer = csv.writer(csv_file, delimiter=",", quoting=csv.QUOTE_NONE)
        writer.writerow(["hap", "pop", "freq"])
        for pop in pops:
            for hap, freq in per_pop_rows[pop]:
                writer.writerow([hap, pop, freq])


def _plan_b_matrix(n_loci):
    """Plan B fallback matrix for n_loci (whole -> progressive splits -> singles)."""
    idx = list(range(1, n_loci + 1))
    matrix = [[idx[:]]]
    if n_loci >= 2:
        head, tail = idx[:-2], idx[-2:]
        matrix.append(([head] if head else []) + [tail])
        if len(head) >= 2:
            matrix.append([[head[0]], head[1:], tail])
        matrix.append(([head] if head else []) + [[idx[-2]], [idx[-1]]])
        if len(head) >= 2:
            matrix.append([[head[0]], head[1:], [idx[-2]], [idx[-1]]])
    matrix.append([[i] for i in idx])
    return matrix


def build_generation_config(
    name,
    populations,
    loci_map,
    freq_data_dir,
    graphs_root="graphs",
    freq_trim_threshold=1e-5,
    number_of_results=1000,
    number_of_pop_results=100,
    epsilon=1e-3,
    source_type="freqs",
):
    """
    Build (and return) the config dict for a generated graph, plus the paths
    it will occupy. Does not run generation.
    """
    graph_dir = Path(graphs_root) / name
    csv_dir = graph_dir / "csv"

    n_loci = len(loci_map)
    config = {
        "populations": list(populations),
        "freq_trim_threshold": freq_trim_threshold,
        "priority": {
            "alpha": 0.4999999,
            "eta": 0,
            "beta": 1e-7,
            "gamma": 1e-7,
            "delta": 0.4999999,
        },
        "UNK_priors": "SR",
        "loci_map": dict(loci_map),
        "factor_missing_data": 0.0001,
        "Plan_B_Matrix": _plan_b_matrix(n_loci),
        "planb": True,
        "number_of_options_threshold": 100000,
        "epsilon": epsilon,
        "number_of_results": number_of_results,
        "number_of_pop_results": number_of_pop_results,
        "output_MUUG": True,
        "output_haplotypes": True,
        # Frequency source: "freqs" (<POP>.freqs.gz) or "hfx" (<POP>.hfx bundle)
        "freq_source_type": source_type,
        "freq_data_dir": str(freq_data_dir),
        # All generated artifacts live inside this graph's own csv dir
        "freq_file": str(csv_dir / "hpf.csv"),
        "graph_files_path": str(csv_dir) + "/",
        "pops_count_file": str(csv_dir / "pop_counts_file.txt"),
        "node_csv_file": "nodes.csv",
        "edges_csv_file": "edges.csv",
        "info_node_csv_file": "info_node.csv",
        "top_links_csv_file": "top_links.csv",
        "max_haplotypes_number_in_phase": 100,
    }
    return config, graph_dir, csv_dir


def generate_graph(config, config_path):
    """
    Run the two-step generation pipeline for the given config dict.
    Writes the config JSON to config_path first, then produces hpf + graph.
    """
    config_path = Path(config_path)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    Path(config["graph_files_path"]).mkdir(parents=True, exist_ok=True)

    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)

    # Step 1: HFX bundles -> hpf.csv + pop_counts_file.txt
    produce_hpf_from_hfx(config)

    # Step 2: hpf.csv -> nodes/edges/top_links/info_node
    # NOTE: generate_neo4j_multi_hpf.generate_graph() calls argparse internally,
    # which would otherwise parse THIS process's sys.argv (breaking CLI callers).
    # Temporarily blank argv so it uses the config_file default we pass in.
    import sys

    saved_argv = sys.argv
    try:
        sys.argv = [saved_argv[0] if saved_argv else "generate_graph"]
        generate_neo4j_multi_hpf.generate_graph(
            config_file=str(config_path),
            em_pop=None,
            em=False,
            use_default_path=False,
        )
    finally:
        sys.argv = saved_argv

    return str(config_path)


class GraphFactory:
    """
    Generates graphs in the background and tracks their status.

    Status states: pending -> running -> ready | error
    """

    # Known HFX frequency sources, keyed by a short source id. All sources are
    # HFX bundles; loci order/names are auto-detected from each bundle's
    # metadata (loci_map=None here). Add more directories here as needed.
    HFX_4LOC_DIR = str(Path("~/src/HF/HFX/4loc").expanduser())
    HFX_9LOC_DIR = str(Path("~/src/HF/HFX/9loc").expanduser())

    FREQ_SOURCES = {
        "hfx4": {
            "dir": HFX_4LOC_DIR,
            "ext": ".hfx",
            "loci": 4,  # A,B,C,DRB1 (detected)
        },
        "hfx9": {
            "dir": HFX_9LOC_DIR,
            "ext": ".hfx",
            "loci": 9,  # A,C,B,DRB3/4/5,DRB1,DQA1,DQB1,DPA1,DPB1 (detected)
        },
    }

    # Default source used when a request does not specify one.
    DEFAULT_SOURCE = "hfx4"

    # Convenience mapping from a bare loci count to a source id.
    _LOCI_TO_SOURCE = {4: "hfx4", 9: "hfx9"}

    # Default expected-count trim (freq * cohortSize). 1.0 keeps haplotypes
    # expected to occur at least once in the sample. Higher = smaller/faster.
    DEFAULT_TRIM = 1.0

    def __init__(self, graphs_root="graphs"):
        self.graphs_root = graphs_root
        self._jobs = {}  # name -> job dict
        self._lock = threading.Lock()

    def _resolve_source(self, source=None, loci=None):
        """Resolve a source id from an explicit source or a loci count."""
        if source:
            if source not in self.FREQ_SOURCES:
                raise ValueError(
                    f"Unknown source '{source}'. Options: {sorted(self.FREQ_SOURCES)}"
                )
            return source
        if loci is not None:
            if loci not in self._LOCI_TO_SOURCE:
                raise ValueError(
                    f"No source for loci={loci}. Options: {sorted(self._LOCI_TO_SOURCE)}"
                )
            return self._LOCI_TO_SOURCE[loci]
        return self.DEFAULT_SOURCE

    def _pops_in_source(self, spec):
        d = Path(spec["dir"])
        ext = spec["ext"]
        if not d.is_dir():
            return []
        return sorted(fp.name[: -len(ext)] for fp in d.glob(f"*{ext}"))

    # -- introspection -----------------------------------------------------
    def available_sources(self):
        """List HFX frequency sources and the populations available in each."""
        sources = []
        for source_id, spec in self.FREQ_SOURCES.items():
            pops = self._pops_in_source(spec)
            loci_names = None
            if pops and HFX_AVAILABLE:
                try:
                    loci_names = list(
                        hfx_loci_map(
                            Path(spec["dir"]) / f"{pops[0]}{spec['ext']}"
                        ).keys()
                    )
                except Exception:  # noqa: BLE001
                    loci_names = None
            sources.append(
                {
                    "source": source_id,
                    "type": "hfx",
                    "loci": spec["loci"],
                    "freq_data_dir": spec["dir"],
                    "loci_map": loci_names,
                    "default": source_id == self.DEFAULT_SOURCE,
                    "available_populations": pops,
                }
            )
        return sources

    def job_status(self, name):
        with self._lock:
            job = self._jobs.get(name)
            return dict(job) if job else None

    def all_jobs(self):
        with self._lock:
            return [dict(j) for j in self._jobs.values()]

    def suggest_name(self, populations, loci=None, source=None, taken=None):
        """
        Suggest a graph name from the populations and source/loci.

        Examples:
            (["CAU"], source="hfx4")   -> "cau_hfx4"
            (["CAU","HIS"], loci=9)    -> "cau-his_hfx9"
            many populations           -> "<first>-plus<N>_<source>"

        Ensures uniqueness against existing jobs and `taken` (e.g. registered
        graph names) by appending -2, -3, ... when needed.
        """
        # The suffix is the resolved HFX source id (e.g. "hfx4").
        try:
            suffix = self._resolve_source(source=source, loci=loci)
        except ValueError:
            suffix = self.DEFAULT_SOURCE

        pops = [p.strip().lower() for p in (populations or []) if p and p.strip()]
        if not pops:
            base = f"graph_{suffix}"
        elif len(pops) <= 3:
            base = "-".join(pops) + f"_{suffix}"
        else:
            base = f"{pops[0]}-plus{len(pops) - 1}_{suffix}"

        taken_set = set(taken or [])
        with self._lock:
            taken_set.update(self._jobs.keys())

        if base not in taken_set:
            return base
        i = 2
        while f"{base}-{i}" in taken_set:
            i += 1
        return f"{base}-{i}"

    # -- generation --------------------------------------------------------
    def start(
        self,
        name,
        populations,
        loci=None,
        source=None,
        freq_trim_threshold=None,
        on_ready=None,
    ):
        """
        Kick off generation in a background thread. Returns the initial job dict.

        Args:
            name: unique graph name.
            populations: list of population codes present in the source.
            loci: optional loci count (4|9) used to pick a source when
                  `source` is not given.
            source: optional explicit HFX source id ("hfx4", "hfx9").
                    Defaults to the factory DEFAULT_SOURCE (hfx4).
            freq_trim_threshold: expected-count cutoff (freq * cohortSize).
                  Defaults to 1.0 (keep haplotypes expected >= 1x in sample).
            on_ready: optional callback(name, config_path) invoked on success.
        """
        if not HFX_AVAILABLE:
            raise ValueError("hfx_tools is not installed; cannot read HFX bundles.")

        source_id = self._resolve_source(source=source, loci=loci)
        spec = self.FREQ_SOURCES[source_id]
        freq_dir = Path(spec["dir"])
        ext = spec["ext"]

        if freq_trim_threshold is None:
            freq_trim_threshold = self.DEFAULT_TRIM

        # Validate populations exist in the source.
        missing = [p for p in populations if not (freq_dir / f"{p}{ext}").is_file()]
        if missing:
            raise ValueError(
                f"Populations not found in {freq_dir}: {missing}. "
                f"Available: {self._pops_in_source(spec)}"
            )

        # Loci order/names are detected from the HFX bundle metadata.
        loci_map = hfx_loci_map(freq_dir / f"{populations[0]}{ext}")

        with self._lock:
            existing = self._jobs.get(name)
            if existing and existing["status"] in ("pending", "running"):
                raise ValueError(f"Graph '{name}' is already being generated")

        config, graph_dir, csv_dir = build_generation_config(
            name=name,
            populations=populations,
            loci_map=loci_map,
            freq_data_dir=str(freq_dir),
            graphs_root=self.graphs_root,
            freq_trim_threshold=freq_trim_threshold,
            source_type="hfx",
        )
        config_path = graph_dir / "config.json"

        job = {
            "name": name,
            "status": "pending",
            "populations": list(populations),
            "source": source_id,
            "type": "hfx",
            "loci": spec["loci"],
            "loci_map": list(loci_map.keys()),
            "freq_trim_threshold": freq_trim_threshold,
            "config_path": str(config_path),
            "graph_dir": str(graph_dir),
            "started_at": time.time(),
            "finished_at": None,
            "error": None,
        }
        with self._lock:
            self._jobs[name] = job

        thread = threading.Thread(
            target=self._run,
            args=(name, config, config_path, on_ready),
            daemon=True,
        )
        thread.start()
        return dict(job)

    def _set(self, name, **fields):
        with self._lock:
            if name in self._jobs:
                self._jobs[name].update(fields)

    def _run(self, name, config, config_path, on_ready):
        self._set(name, status="running")
        try:
            generate_graph(config, config_path)
            self._set(name, status="ready", finished_at=time.time())
            if on_ready:
                try:
                    on_ready(name, str(config_path))
                except Exception:  # noqa: BLE001 - registration must not crash job
                    self._set(
                        name,
                        error="generated but registration failed:\n"
                        + traceback.format_exc(),
                    )
        except Exception:  # noqa: BLE001 - capture full trace for the status endpoint
            self._set(
                name,
                status="error",
                finished_at=time.time(),
                error=traceback.format_exc(),
            )
