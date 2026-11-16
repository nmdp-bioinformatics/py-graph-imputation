#!/usr/bin/env python3
"""
Command-line HFX graph generator.

Runs the same pipeline as the API's POST /graphs, but synchronously in your
terminal with timing for each step, so you can experiment with populations,
loci source, and freq_trim_threshold and see how long generation takes.

Examples:
    # Single population, aggressive trim (expected-count cutoff = 100)
    python3 generate_graph_cli.py --name cau_test --populations CAU \
        --source hfx4 --trim 100

    # Five populations, default trim (=1)
    python3 generate_graph_cli.py --name us5 \
        --populations AFA,CAU,HIS,NAM,API --source hfx4

    # 9-locus (large)
    python3 generate_graph_cli.py --name cau9 --populations CAU --source hfx9 --trim 10

Trim is an EXPECTED-COUNT cutoff for HFX (frequency x cohort size):
    --trim 1     keep haplotypes expected >= 1x in the sample
    --trim 100   keep haplotypes expected >= 100x  (much smaller / faster)
"""

import argparse
import time
from pathlib import Path

import graph_factory as gf


def main():
    parser = argparse.ArgumentParser(description="Generate an HFX imputation graph")
    parser.add_argument("--name", required=True, help="Graph name (output dir)")
    parser.add_argument(
        "--populations",
        required=True,
        help="Comma-separated population codes, e.g. CAU or AFA,CAU,HIS,NAM,API",
    )
    parser.add_argument(
        "--source",
        default="hfx4",
        choices=sorted(gf.GraphFactory.FREQ_SOURCES.keys()),
        help="HFX source id (default: hfx4)",
    )
    parser.add_argument(
        "--trim",
        type=float,
        default=None,
        help="Expected-count trim threshold (default: 1.0). Higher = smaller/faster.",
    )
    parser.add_argument(
        "--graphs-root",
        default="graphs",
        help="Root directory for generated graphs (default: graphs)",
    )
    args = parser.parse_args()

    populations = [p.strip() for p in args.populations.split(",") if p.strip()]
    spec = gf.GraphFactory.FREQ_SOURCES[args.source]
    freq_dir = Path(spec["dir"])
    ext = spec["ext"]

    # Validate populations exist.
    missing = [p for p in populations if not (freq_dir / f"{p}{ext}").is_file()]
    if missing:
        avail = sorted(fp.name[: -len(ext)] for fp in freq_dir.glob(f"*{ext}"))
        parser.error(
            f"Populations not found in {freq_dir}: {missing}\nAvailable: {avail}"
        )

    trim = args.trim if args.trim is not None else gf.GraphFactory.DEFAULT_TRIM
    loci_map = gf.hfx_loci_map(freq_dir / f"{populations[0]}{ext}")

    print("=" * 80)
    print("HFX graph generation")
    print(f"  name:        {args.name}")
    print(f"  populations: {populations}")
    print(f"  source:      {args.source}  ({freq_dir})")
    print(f"  loci:        {list(loci_map.keys())}")
    print(f"  trim:        {trim}  (expected-count cutoff = freq * cohortSize)")
    print("=" * 80)

    config, graph_dir, csv_dir = gf.build_generation_config(
        name=args.name,
        populations=populations,
        loci_map=loci_map,
        freq_data_dir=str(freq_dir),
        graphs_root=args.graphs_root,
        freq_trim_threshold=trim,
        source_type="hfx",
    )
    config_path = graph_dir / "config.json"

    t0 = time.time()
    gf.generate_graph(config, config_path)
    dt = time.time() - t0

    # Report resulting graph size.
    nodes_file = csv_dir / "nodes.csv"
    n_nodes = (sum(1 for _ in open(nodes_file)) - 1) if nodes_file.is_file() else 0
    print("=" * 80)
    print(f"DONE in {dt:.1f}s")
    print(f"  config: {config_path}")
    print(f"  nodes:  {n_nodes}")
    print(f"  files:  {sorted(p.name for p in csv_dir.glob('*.csv'))}")
    print("=" * 80)
    print("Register it with a running server by pointing --graph at the config,")
    print(
        f"or restart: python3 api_server.py --graph {args.name}={config_path} --eager"
    )


if __name__ == "__main__":
    main()
