# -*- coding: utf-8 -*-
"""Tests for HFX (Haplotype Frequency Exchange) support in graph generation."""

import csv
import gzip
import json
import zipfile
from types import SimpleNamespace

import pytest

pytest.importorskip("hfx_tools")

from graph_generation.generate_hpf import produce_hpf  # noqa: E402
from graph_generation.generate_hpf_from_hfx import (  # noqa: E402
    hfx_cohort_size,
    hfx_loci_map,
    produce_hpf_from_hfx,
    sort_haplotype,
)

LOCI = ["A", "B", "C", "DRB1", "DQB1"]  # deliberately not alphabetical
ROWS = [
    ("A*02:01g~B*07:02g~C*07:02g~DRB1*15:01g~DQB1*06:02g", 0.6),
    ("A*01:01g~B*08:01g~C*07:01g~DRB1*03:01g~DQB1*02:01g", 0.4),
    ("A*03:01g~B*35:01g~C*04:01g~DRB1*01:01g~DQB1*05:01g", 0.0),  # dropped
]


def make_hfx(directory, name, rows=ROWS, loci=LOCI, cohort_size=0):
    """Write a minimal <name>.hfx (metadata.json + frequencies.csv)."""
    metadata = {
        "version": "0.1.1",
        "metadata": {
            "hfeMethod": {"method": "unknown"},
            "cohortDescription": {
                "species": "Homo sapiens",
                "population": [{"name": name}],
                "cohortSize": cohort_size,
                "dataSource": "unknown",
            },
            "nomenclatureUsed": {"database": "IPD-IMGT/HLA", "version": "unknown"},
            "frequencyLocation": "file://frequencies.csv",
        },
    }
    if loci:
        metadata["metadata"]["outputResolution"] = [
            {"locus": locus, "resolution": "g"} for locus in loci
        ]
    path = directory / f"{name}.hfx"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("metadata.json", json.dumps(metadata))
        z.writestr(
            "frequencies.csv",
            "haplotype,frequency\n" + "".join(f"{h},{f}\n" for h, f in rows),
        )
    return path


def conf_for(tmp_path, pops, **extra):
    out = tmp_path / "out"
    conf = {
        "populations": pops,
        "freq_source_type": "hfx",
        "freq_data_dir": str(tmp_path),
        "graph_files_path": str(out),
        "pops_count_file": str(out / "pop_counts_file.txt"),
        "freq_file": str(out / "hpf.csv"),
    }
    conf.update(extra)
    return conf


def read_hpf(tmp_path):
    with open(tmp_path / "out" / "hpf.csv", newline="") as f:
        return list(csv.DictReader(f))


def read_counts(tmp_path):
    with open(tmp_path / "out" / "pop_counts_file.txt") as f:
        return [line.strip().split(",") for line in f]


def test_sort_haplotype_orders_alleles_by_locus():
    assert (
        sort_haplotype("A*01:01g~DRB1*03:01g~B*08:01g~DQB1*02:01g~C*07:01g")
        == "A*01:01g~B*08:01g~C*07:01g~DQB1*02:01g~DRB1*03:01g"
    )


def test_loci_map_is_alphabetical_whatever_the_metadata_order(tmp_path):
    path = make_hfx(tmp_path, "pop")
    assert hfx_loci_map(path) == {"A": 1, "B": 2, "C": 3, "DQB1": 4, "DRB1": 5}


def test_loci_map_falls_back_to_the_first_haplotype(tmp_path):
    path = make_hfx(tmp_path, "pop", loci=None)
    assert hfx_loci_map(path) == {"A": 1, "B": 2, "C": 3, "DQB1": 4, "DRB1": 5}


@pytest.mark.parametrize(
    "cohort, expected",
    [
        ({"cohortSize": 250}, 250.0),
        ({"cohortSize": 0}, None),
        ({}, None),
        ({"cohortSize": 0, "population": [{"name": "x", "populationSize": 80}]}, 80.0),
        ({"population": {"populationSize": 12}}, 12.0),
    ],
)
def test_cohort_size(cohort, expected):
    bundle = SimpleNamespace(metadata={"cohortDescription": cohort})
    assert hfx_cohort_size(bundle) == expected


def test_hpf_has_sorted_haplotypes_and_skips_zero_frequencies(tmp_path):
    make_hfx(tmp_path, "popA")
    produce_hpf_from_hfx(conf_for(tmp_path, ["popA"]))
    rows = read_hpf(tmp_path)
    assert [r["pop"] for r in rows] == ["popA", "popA"]
    assert rows[0]["hap"] == "A*02:01g~B*07:02g~C*07:02g~DQB1*06:02g~DRB1*15:01g"
    assert [float(r["freq"]) for r in rows] == [0.6, 0.4]


def test_population_names_may_contain_hyphens(tmp_path):
    make_hfx(tmp_path, "hfset-10")
    produce_hpf_from_hfx(conf_for(tmp_path, ["hfset-10"]))
    assert {r["pop"] for r in read_hpf(tmp_path)} == {"hfset-10"}


def test_counts_use_cohort_size_and_default_to_one(tmp_path):
    make_hfx(tmp_path, "big", cohort_size=300)
    make_hfx(tmp_path, "unknown", cohort_size=0)
    produce_hpf_from_hfx(conf_for(tmp_path, ["big", "unknown"]))
    (big, unknown) = read_counts(tmp_path)
    assert big[:2] == ["big", "300.0"]
    assert unknown[:2] == ["unknown", "1.0"]
    assert float(big[2]) == pytest.approx(300 / 301)
    assert float(big[2]) + float(unknown[2]) == pytest.approx(1.0)


def test_several_populations_share_one_hpf(tmp_path):
    make_hfx(tmp_path, "one")
    make_hfx(tmp_path, "two", rows=ROWS[:1])
    produce_hpf_from_hfx(conf_for(tmp_path, ["one", "two"]))
    assert [r["pop"] for r in read_hpf(tmp_path)] == ["one", "one", "two"]


def test_missing_bundle_is_reported_by_path(tmp_path):
    with pytest.raises(FileNotFoundError, match="nope.hfx"):
        produce_hpf_from_hfx(conf_for(tmp_path, ["nope"]))


def test_produce_hpf_uses_hfx_when_the_config_says_so(tmp_path):
    make_hfx(tmp_path, "popA")
    conf_file = tmp_path / "conf.json"
    conf_file.write_text(json.dumps(conf_for(tmp_path, ["popA"])))
    produce_hpf(str(conf_file))
    assert len(read_hpf(tmp_path)) == 2


def test_produce_hpf_still_reads_freqs_gz_by_default(tmp_path):
    out = tmp_path / "out"
    with gzip.open(tmp_path / "CAU.freqs.gz", "wt") as f:
        f.write("Haplo,Count,Freq\nA*01:01~B*08:01,10,0.75\nA*02:01~B*07:02,5,0.25\n")
    conf = conf_for(tmp_path, ["CAU"])
    del conf["freq_source_type"]  # the default is the existing .freqs.gz behaviour
    conf_file = tmp_path / "conf.json"
    conf_file.write_text(json.dumps(conf))
    produce_hpf(str(conf_file))
    assert [(r["hap"], r["pop"]) for r in read_hpf(tmp_path)] == [
        ("A*01:01~B*08:01", "CAU"),
        ("A*02:01~B*07:02", "CAU"),
    ]
    assert read_counts(tmp_path)[0][:2] == ["CAU", "15.0"]
    assert out.is_dir()


def test_unknown_freq_source_type_is_rejected(tmp_path):
    conf_file = tmp_path / "conf.json"
    conf_file.write_text(json.dumps(conf_for(tmp_path, ["x"], freq_source_type="csv")))
    with pytest.raises(ValueError, match="freq_source_type"):
        produce_hpf(str(conf_file))
