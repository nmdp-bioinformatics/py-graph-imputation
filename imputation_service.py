#!/usr/bin/env python3
"""
In-memory HLA imputation service.

Loads a configuration and builds the NetworkX graph ONCE, then imputes
individual GL strings on demand without touching the filesystem for input
or output. This is the reusable core behind the web API.

The heavy config wiring is intentionally kept here (rather than in the API
layer) so the same service can be used from a notebook, a script, or a server.
"""

import json
import logging
import os
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

from grim.imputation.impute import Imputation, clean_up_gl  # noqa: F401
from grim.imputation.networkx_graph import Graph


# Module logger. The application (e.g. api_server) is responsible for
# configuring handlers/levels; by default this logger inherits the root config.
logger = logging.getLogger("imputation_service")


# Every allele in a GL string must look like LOCUS*FIELDS (e.g. A*02:01,
# DRB1*11:XX, B*44:AB). Used by the API layers to reject malformed input with a
# clear 400 instead of a generic "allele '' not in graph" problem.
_ALLELE_RE = re.compile(r"^[A-Z][A-Z0-9]*\*[A-Za-z0-9:]+$")


def find_malformed_alleles(gl_string):
    """Return the tokens of a GL string that are not syntactically alleles."""
    tokens = re.split(r"[\^+/~|]", (gl_string or "").strip())
    return [t for t in tokens if not _ALLELE_RE.match(t.strip())]


# ---------------------------------------------------------------------------
# py-ard integration (ambiguity expansion: XX codes, MAC codes, serology, ...)
# ---------------------------------------------------------------------------
#
# The imputation engine only understands allele-level ambiguity via the '/'
# separator; it does not know how to expand an XX code (e.g. "A*01:XX") or a
# MAC code (e.g. "A*01:AB") into its member alleles. py-ard does exactly that.
#
# We reduce an ambiguous GL string with the "lgx" method, which expands codes
# and reduces every allele to the 2-field ARD level (e.g. "A*01:01") -- the same
# resolution the graphs are built at -- producing a plain '/'-ambiguous GL
# string the engine can impute over.
#
# py-ard's first init() downloads a large IPD-IMGT/HLA database, so it is
# initialized lazily and cached: a plain 2-field GL string never pays the cost.

_ARD_CACHE = {}

# Tokens that indicate a GL string carries ambiguity py-ard must expand.
# - XX code:      A*01:XX
# - MAC code:     A*01:AB   (uppercase-letter second field)
# - serology:     A1, Bw4, Cw7, DR15 ... (locus letters + digits, no '*')
_XX_RE = re.compile(r"\*\d+:XX\b", re.IGNORECASE)
_MAC_RE = re.compile(r"\*\d+:[A-Z]{2,}\b")
_SEROLOGY_RE = re.compile(r"(?:^|[/+^~])(?:Cw|Bw|DR|DQ|DP|A|B|C)\d+(?:$|[/+^~])")


def _get_ard(imgt_version=None):
    """Lazily initialize and cache a py-ard ARD object.

    Returns None (and never raises) if py-ard is not installed or fails to
    initialize, so callers can degrade gracefully.
    """
    key = imgt_version or "Latest"
    if key in _ARD_CACHE:
        return _ARD_CACHE[key]
    try:
        import pyard

        ard = pyard.init(imgt_version) if imgt_version else pyard.init()
    except Exception as e:  # noqa: BLE001 - any failure => no expansion available
        logger.warning(
            "py-ard init failed (version=%s); ambiguity codes will not be "
            "expanded: %s",
            key,
            e,
        )
        ard = None
    _ARD_CACHE[key] = ard
    return ard


def _has_ambiguous_codes(gl_string):
    """Heuristic: does this GL string contain XX/MAC/serology codes?

    A plain '/'-ambiguous or unambiguous 2-field GL string returns False and
    is imputed directly (no py-ard cost).
    """
    if not gl_string:
        return False
    if _XX_RE.search(gl_string) or _MAC_RE.search(gl_string):
        return True
    if _SEROLOGY_RE.search(gl_string):
        return True
    return False


def expand_gl_string(gl_string, imgt_version=None):
    """Expand XX/MAC/serology codes in a GL string via py-ard's 'lgx' redux.

    Returns (expanded_gl_string, note). ``note`` is None on a successful
    expansion, otherwise a short human-readable reason the original string was
    returned unchanged (py-ard unavailable, invalid input, etc.).
    """
    ard = _get_ard(imgt_version)
    if ard is None:
        return gl_string, "py-ard not available; codes not expanded"
    try:
        expanded = ard.redux(gl_string, "lgx")
        return expanded, None
    except Exception as e:  # noqa: BLE001 - invalid allele/code, etc.
        return gl_string, f"py-ard could not expand input: {e}"


def normalize_homozygous_loci(gl_string):
    """Duplicate single-allele loci into a homozygous genotype ('X' -> 'X+X').

    The imputation engine's genotype parser requires every locus to be a full
    genotype with two haplotypes separated by '+'. A locus written as a single
    allele (or a single '/'-ambiguity list, e.g. an expanded 'A*01:XX') is
    otherwise rejected, yielding an empty candidate set. Treating a lone locus
    allele as homozygous is the standard interpretation and lets inputs like
    "A*01:01^B*08:01" (or "A*01:XX^B*08:XX") impute.

    Returns (normalized_gl_string, changed).
    """
    if not gl_string:
        return gl_string, False
    changed = False
    loci = gl_string.split("^")
    for i, locus in enumerate(loci):
        # Mirror the engine: a leading '+' means an empty first haplotype.
        if locus.startswith("+"):
            locus = locus[1:]
        if locus == "":
            continue
        if "+" not in locus:
            loci[i] = f"{locus}+{locus}"
            changed = True
        else:
            loci[i] = locus
    return "^".join(loci), changed


# ---------------------------------------------------------------------------
# Config building
# ---------------------------------------------------------------------------

DEFAULT_PLAN_B_MATRIX = [
    [[1, 2, 3, 4, 5]],
    [[1, 2, 3], [4, 5]],
    [[1], [2, 3], [4, 5]],
    [[1, 2, 3], [4], [5]],
    [[1], [2, 3], [4], [5]],
    [[1], [2], [3], [4], [5]],
]


def _read_graph_loci_map(node_file):
    """
    Derive the loci_map ({locus_name: 1-based index}) from a graph's nodes.csv.

    Node names look like "A*01:01~C*07:01~B*08:01~DRB1*03:01~DQB1*02:01~..."
    The locus order in a FULL-haplotype node defines the index each locus gets
    (A=1, C=2, B=3, ...), which must match how the graph frequency vectors and
    node labels were built.

    We pick the node whose name has the most loci (a full haplotype) to read
    the complete ordering. Returns None if the file can't be read.
    """
    try:
        path = Path(node_file)
        if not path.is_file():
            return None
        best_loci = []
        with open(path) as f:
            header = f.readline()  # skip header
            # Scan a bounded number of rows; the first full-width node is enough.
            for count, line in enumerate(f):
                parts = line.rstrip("\n").split(",")
                if len(parts) < 2:
                    continue
                name = parts[1]
                alleles = name.split("~")
                loci = [a.split("*")[0] for a in alleles if "*" in a]
                if len(loci) > len(best_loci):
                    best_loci = loci
                # Stop once we've seen enough rows or found a wide node.
                if count > 5000:
                    break
        if not best_loci:
            return None
        return {locus: i + 1 for i, locus in enumerate(best_loci)}
    except OSError:
        return None


def _default_plan_b_matrix(n_loci):
    """
    Build a Plan B fallback matrix for n_loci, following the built-in pattern:
    start with the whole haplotype, then progressively split loci apart.

    For n_loci == 5 this reproduces DEFAULT_PLAN_B_MATRIX.
    """
    idx = list(range(1, n_loci + 1))
    matrix = [
        [idx[:]],  # whole haplotype
    ]
    if n_loci >= 2:
        head = idx[:-2]
        tail = idx[-2:]
        # [1..n-2],[n-1,n]
        first = ([head] if head else []) + [tail]
        matrix.append(first)
        # [1],[2..n-2],[n-1,n]  (only meaningful when head has >=2 entries)
        if len(head) >= 2:
            matrix.append([[head[0]], head[1:], tail])
        # [1..n-2],[n-1],[n]
        matrix.append(([head] if head else []) + [[idx[-2]], [idx[-1]]])
        if len(head) >= 2:
            matrix.append([[head[0]], head[1:], [idx[-2]], [idx[-1]]])
    # fully split: [1],[2],...,[n]
    matrix.append([[i] for i in idx])
    return matrix


def _read_graph_populations(info_node_path):
    """
    Read the ordered population list from a graph's info_node.csv.

    Format (Neo4j-style header):
        INFO_NODE_ID:ID(INFO_NODE),populations:STRING[],INFO_NODE:LABEL
        1,CAU;HIS;NAM;AFA;API,INFO_NODE

    Returns a list like ["CAU", "HIS", ...] or None if unavailable.
    """
    try:
        path = Path(info_node_path)
        if not path.is_file():
            return None
        with open(path) as f:
            header = f.readline().strip().split(",")
            try:
                pop_col = next(
                    i for i, h in enumerate(header) if h.startswith("populations")
                )
            except StopIteration:
                pop_col = 1
            for line in f:
                parts = line.strip().split(",")
                if len(parts) > pop_col and parts[pop_col]:
                    return parts[pop_col].split(";")
    except OSError:
        return None
    return None


def build_config(conf_file):
    """
    Build the internal `config` dict the Imputation/Graph classes expect.

    This mirrors grim.run_impute_def.run_impute but strips out every input /
    output *file* path, because the service never reads subjects from disk nor
    writes results to disk. Only the graph CSV paths matter here.
    """
    conf_path = Path(conf_file)
    if not conf_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {conf_path}")

    with open(conf_path) as f:
        json_conf = json.load(f)

    # The graph CSVs are resolved relative to the current working directory
    # (same behaviour as the original CLI when run from the repo root).
    graph_files_path = json_conf.get("graph_files_path", "output/csv/")
    if not graph_files_path.endswith("/"):
        graph_files_path += "/"

    loci_map = json_conf.get("loci_map", {"A": 1, "B": 2, "C": 3, "DQB1": 4, "DRB1": 5})

    # The graph nodes encode a fixed number of loci in a fixed order (e.g.
    # A~C~B~DRB1~DQB1~DPB1). The config's loci_map must match that order/count
    # exactly or lookups silently break. Derive it from the graph itself so the
    # service adapts to whatever graph is loaded (5-locus, 6-locus, etc.).
    node_file = graph_files_path + json_conf.get("node_csv_file", "nodes.csv")
    graph_loci_map = _read_graph_loci_map(node_file)
    if graph_loci_map:
        loci_map = graph_loci_map

    # The graph's frequency vectors have one entry PER population, in the order
    # the graph was generated. The authoritative population list lives in
    # info_node.csv (populations:STRING[] column). We trust that over the
    # config's "populations" field so the service always aligns with whatever
    # graph is actually loaded (avoids index-out-of-bounds on mismatch).
    populations = json_conf.get("populations")
    info_node_name = json_conf.get("info_node_csv_file", "info_node.csv")
    info_node_path = Path(graph_files_path) / info_node_name
    graph_pops = _read_graph_populations(info_node_path)
    if graph_pops:
        populations = graph_pops

    config = {
        "planb": json_conf.get("planb", True),
        "pops": populations,
        "priority": json_conf.get("priority"),
        "epsilon": json_conf.get("epsilon", 1e-3),
        "number_of_results": json_conf.get("number_of_results", 1000),
        "number_of_pop_results": json_conf.get("number_of_pop_results", 100),
        "output_MUUG": json_conf.get("output_MUUG", True),
        "output_haplotypes": json_conf.get("output_haplotypes", True),
        "node_file": graph_files_path + json_conf.get("node_csv_file", "nodes.csv"),
        "top_links_file": graph_files_path
        + json_conf.get("top_links_csv_file", "top_links.csv"),
        "edges_file": graph_files_path + json_conf.get("edges_csv_file", "edges.csv"),
        "factor_missing_data": json_conf.get("factor_missing_data", 0.01),
        "loci_map": loci_map,
        "matrix_planb": json_conf.get("Plan_B_Matrix", DEFAULT_PLAN_B_MATRIX),
        "pops_count_file": json_conf.get("pops_count_file", ""),
        "use_pops_count_file": bool(json_conf.get("pops_count_file", False)),
        "number_of_options_threshold": json_conf.get(
            "number_of_options_threshold", 100000
        ),
        "max_haplotypes_number_in_phase": json_conf.get(
            "max_haplotypes_number_in_phase", 100
        ),
        "nodes_for_plan_A": json_conf.get("Plan_A_Matrix", []),
        "save_mode": json_conf.get("save_space_mode", False),
        "UNK_priors": json_conf.get("UNK_priors", "MR"),
    }

    # full_loci is derived from the loci_map values (as in run_impute_def)
    all_loci_set = {str(v) for v in loci_map.values()}
    config["full_loci"] = "".join(sorted(all_loci_set))

    # Ensure the Plan B matrix covers every locus index. The stock config
    # matrix only references indices 1..5; on a 6+ locus graph that would omit
    # loci from the fallback. Regenerate a matching matrix when needed.
    n_loci = len(loci_map)
    max_index_in_matrix = 0
    for grouping in config["matrix_planb"]:
        for group in grouping:
            if group:
                max_index_in_matrix = max(max_index_in_matrix, max(group))
    if max_index_in_matrix != n_loci:
        config["matrix_planb"] = _default_plan_b_matrix(n_loci)

    return config


def resolve_count_by_prob(config):
    """
    Build the population-count weight vector aligned to config["pops"].

    The original Imputation constructor blindly reads every line of the
    pops_count_file into a vector sized to len(populations), which breaks
    when the counts file lists MORE populations than the graph was built
    for (e.g. a 5-population counts file with a single-population CAU graph).

    Here we match counts to the configured populations by name, defaulting
    to 1.0 for any population missing from the file. Returns None when no
    usable counts file is configured, letting Imputation default to ones.
    """
    pops = config.get("pops") or []
    counts_path = config.get("pops_count_file", "")

    if not config.get("use_pops_count_file") or not counts_path:
        return None
    if not os.path.isfile(counts_path):
        return None

    file_counts = {}
    with open(counts_path) as f_count:
        for line in f_count:
            parts = line.strip().split(",")
            if len(parts) >= 3:
                try:
                    file_counts[parts[0]] = float(parts[2])
                except ValueError:
                    continue

    weights = np.ones(len(pops))
    for i, pop in enumerate(pops):
        if pop in file_counts:
            weights[i] = file_counts[pop]
    return weights


# ---------------------------------------------------------------------------
# Result formatting
# ---------------------------------------------------------------------------


def _truncate_and_normalize(sorted_items, limit, marginal_ratio_threshold):
    """
    Given (key, frequency) pairs already sorted descending by frequency, apply
    the count cap then the marginal-frequency-ratio truncation, and compute a
    renormalized probability for each kept item.

    Marginal ratio rule: walk results in rank order accumulating cum_sum. Before
    admitting result k (k >= 1), compute freq_k / cum_sum_before_k. If that ratio
    is below marginal_ratio_threshold, stop (result k and beyond are dropped).
    The first result is always kept. Probabilities are freq / sum(kept freqs).

    Returns a list of (key, frequency, probability).
    """
    capped = sorted_items[:limit] if limit is not None else list(sorted_items)

    if marginal_ratio_threshold and marginal_ratio_threshold > 0:
        kept = []
        cum_sum = 0.0
        for _, (key, freq) in enumerate(capped):
            if kept and cum_sum > 0:
                ratio = freq / cum_sum
                if ratio < marginal_ratio_threshold:
                    break
            kept.append((key, freq))
            cum_sum += freq
    else:
        kept = list(capped)

    total = sum(freq for _, freq in kept)
    out = []
    for key, freq in kept:
        prob = (freq / total) if total > 0 else 0.0
        out.append((key, freq, prob))
    return out


def _rank_muug(res, limit, marginal_ratio_threshold=None):
    """
    Format an unphased (UMUG) result.

    res["Haps"] is a dict {genotype_string: probability}
    res["Pops"] is a dict {pop_string:      probability}
    """
    if not res or not isinstance(res.get("Haps"), dict):
        return []
    items = sorted(res["Haps"].items(), key=lambda kv: kv[1], reverse=True)
    ranked = _truncate_and_normalize(items, limit, marginal_ratio_threshold)
    return [
        {"rank": r, "genotype": g, "frequency": float(f), "probability": float(p)}
        for r, (g, f, p) in enumerate(ranked)
    ]


def _rank_muug_pops(res, limit, marginal_ratio_threshold=None):
    if not res or not isinstance(res.get("Pops"), dict):
        return []
    items = sorted(res["Pops"].items(), key=lambda kv: kv[1], reverse=True)
    ranked = _truncate_and_normalize(items, limit, marginal_ratio_threshold)
    return [
        {"rank": r, "pop": pop, "frequency": float(f), "probability": float(p)}
        for r, (pop, f, p) in enumerate(ranked)
    ]


def _rank_haps(res, limit, marginal_ratio_threshold=None):
    """
    Format a phased (PMUG) result.

    res["Haps"]  is a list of [hap1, hap2]
    res["Probs"] is a list of probabilities
    Duplicate genotypes (hap1+hap2) are summed, matching write_best_prob.
    """
    if not res or not isinstance(res.get("Haps"), list):
        return []
    haps = res["Haps"]
    probs = res["Probs"]

    summed = defaultdict(float)
    for i in range(len(haps)):
        key = haps[i][0] + "+" + haps[i][1]
        alt = haps[i][1] + "+" + haps[i][0]
        if alt in summed:
            summed[alt] += probs[i]
        else:
            summed[key] += probs[i]

    items = sorted(summed.items(), key=lambda kv: kv[1], reverse=True)
    ranked = _truncate_and_normalize(items, limit, marginal_ratio_threshold)
    return [
        {"rank": r, "genotype": g, "frequency": float(f), "probability": float(p)}
        for r, (g, f, p) in enumerate(ranked)
    ]


def _rank_haps_pops(res, limit, marginal_ratio_threshold=None):
    """Phased population pairs, summed by pop pair (matches write_best_prob)."""
    if not res or not isinstance(res.get("Pops"), list):
        return []
    pops = res["Pops"]
    probs = res["Probs"]

    summed = defaultdict(float)
    for i in range(len(pops)):
        key = pops[i][0] + "," + pops[i][1]
        alt = pops[i][1] + "," + pops[i][0]
        if alt in summed:
            summed[alt] += probs[i]
        else:
            summed[key] += probs[i]

    items = sorted(summed.items(), key=lambda kv: kv[1], reverse=True)
    ranked = _truncate_and_normalize(items, limit, marginal_ratio_threshold)
    return [
        {"rank": r, "pops": pp, "frequency": float(f), "probability": float(p)}
        for r, (pp, f, p) in enumerate(ranked)
    ]


# ---------------------------------------------------------------------------
# Service
# ---------------------------------------------------------------------------


class ImputationService:
    """
    Loads config + graph once, then imputes single GL strings in memory.

    Usage:
        svc = ImputationService("conf/minimal-configuration.json")
        result = svc.impute("A*01:01+A*02:01^B*07:02+B*08:01^DRB1*03:01+DRB1*15:01")
    """

    def __init__(self, conf_file):
        self.conf_file = conf_file
        self.config = build_config(conf_file)
        self.populations = self.config["pops"]

        # Build the graph a single time — this is the expensive step.
        logger.info(
            "building graph from config=%r (populations=%s)",
            conf_file,
            self.populations,
        )
        try:
            self.graph = Graph(self.config)
            self.graph.build_graph(
                self.config["node_file"],
                self.config["top_links_file"],
                self.config["edges_file"],
            )
        except Exception:
            logger.exception("failed to build graph from config=%r", conf_file)
            raise

        # One Imputation instance reused across requests.
        # Resolve population weights ourselves (name-matched) and pass them in,
        # so Imputation.__init__ does not re-read the counts file with a
        # possibly-mismatched population count.
        count_by_prob = resolve_count_by_prob(self.config)
        self.config["use_pops_count_file"] = False
        self.imputation = Imputation(self.graph, self.config, count_by_prob)
        logger.info("graph ready from config=%r", conf_file)

    @staticmethod
    def _attach_expansion(
        response,
        original_gl_string,
        gl_expanded,
        expand_note,
        gl_homozygous_filled=False,
    ):
        """Record *how* (if at all) the input GL string was rewritten.

        Only boolean flags and any warning note are surfaced; the GL string
        values themselves (original and expanded) are intentionally NOT included
        in the response — an expanded XX code can be enormous. The strings are
        logged server-side instead. ``expanded`` marks a py-ard code expansion;
        ``homozygous_filled`` marks single-allele loci filled to X+X.
        """
        if gl_expanded:
            response["expanded"] = True
        if gl_homozygous_filled:
            response["homozygous_filled"] = True
        if expand_note:
            response["expand_note"] = expand_note

    def impute(
        self,
        gl_string,
        pop1=None,
        pop2=None,
        output="both",
        number_of_results=None,
        number_of_pop_results=None,
        marginal_ratio_threshold=None,
        expand="auto",
    ):
        """
        Impute a single GL string.

        Args:
            gl_string: GL string, loci separated by '^', haplotypes by '+',
                       allele ambiguity by '/'. Example:
                       "A*01:01+A*02:01^B*07:02+B*08:01^DRB1*03:01+DRB1*15:01"
            pop1, pop2: optional population priors (must exist in the graph's
                        configured populations, otherwise ignored).
            output: "umug" (unphased; internally MUUG), "pmug" (phased), or
                    "both". Aliases accepted: muug/unphased -> umug;
                    haplotype/phased -> pmug.
            number_of_results / number_of_pop_results: optional overrides for
                        how many ranked genotype / population results to return.
            marginal_ratio_threshold: optional float (e.g. 0.001). When set,
                        the ranked results are truncated at the first point
                        where the next result's frequency divided by the
                        cumulative frequency-sum so far falls below the
                        threshold, and the kept results are renormalized into a
                        "probability" that sums to 1. Applied AFTER the
                        number_of_results cap.

        Returns a dict with the requested sections keyed "umug" and/or "pmug".
        """
        if not gl_string or not gl_string.strip():
            raise ValueError("gl_string is required")

        # --- Normalize alleles to 2-field via py-ard (lgx) ------------------
        # expand:
        #   "auto"  -> reduce EVERY input through py-ard lgx (default). This
        #              normalizes 3/4-field, G/P groups, serology, MAC and XX
        #              codes to the 2-field alleles the engine requires;
        #              already-2-field alleles pass through unchanged.
        #   True    -> same as auto (always run lgx redux).
        #   False   -> never touch the input; pass the GL string through as-is
        #              (only safe if the caller guarantees 2-field alleles).
        original_gl_string = gl_string
        expand_note = None
        gl_expanded = False
        # The imputation engine only understands 2-field alleles. It cannot
        # handle 3/4-field alleles (A*03:01:01), G/P group names (A*03:01:01G,
        # A*03:01P), serology, MAC, or XX codes. py-ard's "lgx" redux normalizes
        # ALL of these down to the 2-field level the graph is built at, and
        # leaves already-2-field alleles unchanged. So "auto" reduces EVERY
        # input through py-ard — not just ones that look ambiguous.
        want_expand = expand is True or expand == "auto" or expand is None
        if want_expand:
            new_gl, expand_note = expand_gl_string(
                gl_string, self.config.get("imgt_version")
            )
            if expand_note is None and new_gl and new_gl != gl_string:
                gl_string = new_gl
                gl_expanded = True

        # Treat any single-allele locus as homozygous (X -> X+X); the engine's
        # genotype parser requires two haplotypes per locus. Applied after
        # expansion so an expanded single-code locus is normalized too.
        gl_string, gl_homozygous_filled = normalize_homozygous_loci(gl_string)

        output = (output or "both").lower()
        want_umug = output in ("umug", "muug", "unphased", "both")
        want_pmug = output in (
            "pmug",
            "haplotype",
            "haplotypes",
            "phased",
            "both",
        )
        if not (want_umug or want_pmug):
            raise ValueError(
                f"output must be one of: umug, pmug, both (got {output!r})"
            )

        n_results = (
            number_of_results
            if number_of_results is not None
            else self.config["number_of_results"]
        )
        n_pop_results = (
            number_of_pop_results
            if number_of_pop_results is not None
            else self.config["number_of_pop_results"]
        )
        if marginal_ratio_threshold is None:
            marginal_ratio_threshold = self.config.get("marginal_ratio_threshold")

        priority = self.config["priority"]
        epsilon = self.config["epsilon"]
        n = 1000
        planb = self.config["planb"]

        # binary vector: 1 per locus gap (len(full_loci) - 1), same as impute_file
        subject_bin = [1] * (len(self.imputation.full_loci) - 1)

        # Reset per-request plan bookkeeping (as impute_file does).
        self.imputation.plan = "a"
        self.imputation.option_1 = 0
        self.imputation.option_2 = 0

        try:
            _, res_muugs, res_haps = self.imputation.impute_one(
                "query",
                gl_string,
                subject_bin,
                pop1,
                pop2,
                priority,
                epsilon,
                n,
                want_umug,
                want_pmug,
                planb,
                False,  # em
            )
        except KeyError as e:
            # The engine raises KeyError when an allele in the (expanded) input
            # is not a node in this graph — e.g. an unexpanded serology/MAC/XX
            # token (py-ard unavailable or an allele the graph doesn't carry).
            # Convert to a graceful 'problem' instead of a 500.
            problem = {
                "status": "problem",
                "message": (
                    f"Imputation failed: allele {e} is not present in graph "
                    f"'{self.conf_file}' (it may be an unexpanded code or an "
                    f"allele the graph does not carry)."
                ),
            }
            self._attach_expansion(
                problem,
                original_gl_string,
                gl_expanded,
                expand_note,
                gl_homozygous_filled,
            )
            logger.warning(
                "impute problem: allele %s not in graph. input=%r imputed=%r",
                e,
                original_gl_string,
                gl_string,
            )
            return problem

        if res_muugs is None and res_haps is None:
            problem = {
                "status": "problem",
                "message": "Imputation failed (invalid input or empty candidate set).",
            }
            self._attach_expansion(
                problem,
                original_gl_string,
                gl_expanded,
                expand_note,
                gl_homozygous_filled,
            )
            logger.warning(
                "impute problem: no candidates. input=%r imputed=%r pop1=%r pop2=%r",
                original_gl_string,
                gl_string,
                pop1,
                pop2,
            )
            return problem

        response = {
            "populations": self.populations,
            "status": "ok",
        }
        self._attach_expansion(
            response,
            original_gl_string,
            gl_expanded,
            expand_note,
            gl_homozygous_filled,
        )
        if pop1 or pop2:
            response["input_pops"] = [pop1, pop2]

        if want_umug:
            response["umug"] = {
                "genotypes": _rank_muug(res_muugs, n_results, marginal_ratio_threshold),
                "populations": _rank_muug_pops(
                    res_muugs, n_pop_results, marginal_ratio_threshold
                ),
            }
        if want_pmug:
            response["pmug"] = {
                "genotypes": _rank_haps(res_haps, n_results, marginal_ratio_threshold),
                "populations": _rank_haps_pops(
                    res_haps, n_pop_results, marginal_ratio_threshold
                ),
            }

        # Mark a miss if nothing came back at all.
        empty = True
        if want_umug and response["umug"]["genotypes"]:
            empty = False
        if want_pmug and response["pmug"]["genotypes"]:
            empty = False
        if empty:
            response["status"] = "miss"
            response["message"] = "No imputation output for the given input."
            logger.info(
                "impute miss: no output. input=%r imputed=%r pop1=%r pop2=%r",
                original_gl_string,
                gl_string,
                pop1,
                pop2,
            )

        if expand_note:
            logger.warning(
                "expansion note for input=%r: %s", original_gl_string, expand_note
            )

        return response
