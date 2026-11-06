# -*- coding: utf-8 -*-
"""End-to-end test of ImputationService on a graph built from an HFX bundle."""

import json
import pathlib
import sys
import zipfile

import pytest

pytest.importorskip("hfx_tools")
pytest.importorskip("grim.imputation.cutils")  # compiled extension

from test_hfx import make_hfx  # noqa: E402

from graph_generation.generate_hpf import produce_hpf  # noqa: E402
from grim.grim import graph_freqs  # noqa: E402
from imputation_service import ImputationService  # noqa: E402

POP = "popA"
LOCI = ["A", "B", "C", "DQB1", "DRB1"]
H1 = "A*01:01g~B*08:01g~C*07:01g~DQB1*02:01g~DRB1*03:01g"
H2 = "A*02:01g~B*07:02g~C*07:02g~DQB1*06:02g~DRB1*15:01g"
H3 = "A*03:01g~B*35:01g~C*04:01g~DQB1*03:01g~DRB1*11:01g"
GENOTYPE = (
    "A*01:01+A*02:01^B*07:02+B*08:01^C*07:01+C*07:02"
    "^DQB1*02:01+DQB1*06:02^DRB1*03:01+DRB1*15:01"
)


@pytest.fixture
def service(tmp_path, monkeypatch):
    """Build a graph from a 3-haplotype HFX and load it into an ImputationService."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "data" / "hfx").mkdir(parents=True)
    make_hfx(
        tmp_path / "data" / "hfx",
        POP,
        rows=[(H1, 0.5), (H2, 0.3), (H3, 0.2)],
        loci=LOCI,
    )
    with zipfile.ZipFile(
        pathlib.Path(__file__).parents[2] / "example-conf-data.zip"
    ) as z:
        conf = json.loads(z.read("conf/minimal-configuration.json"))
    conf.update(
        populations=[POP],
        freq_source_type="hfx",
        freq_data_dir="data/hfx",
        freq_trim_threshold=1e-6,
        planb=False,
        priority={"alpha": 1.0, "eta": 0, "beta": 0, "gamma": 0, "delta": 0},
    )
    (tmp_path / "conf").mkdir()
    (tmp_path / "conf" / "hfx.json").write_text(json.dumps(conf))
    produce_hpf(conf_file="conf/hfx.json")
    # graph generation parses sys.argv with argparse; hide pytest's own arguments.
    monkeypatch.setattr(sys, "argv", ["graph_freqs"])
    graph_freqs(conf_file="conf/hfx.json")
    return ImputationService("conf/hfx.json")


def test_service_imputes_from_an_hfx_graph(service):
    result = service.impute(GENOTYPE, pop1=POP, pop2=POP, output="umug", expand=False)
    assert result["status"] == "ok"
    genotypes = result["umug"]["genotypes"]
    assert [g["genotype"] for g in genotypes] == [
        "A*01:01+A*02:01^B*07:02+B*08:01^C*07:01+C*07:02"
        "^DQB1*02:01+DQB1*06:02^DRB1*03:01+DRB1*15:01"
    ]
    assert genotypes[0]["probability"] == pytest.approx(1.0)


def test_genotype_frequency_is_twice_the_haplotype_product(service):
    result = service.impute(GENOTYPE, pop1=POP, pop2=POP, output="umug", expand=False)
    # H1 (0.5) and H2 (0.3) are different haplotypes: 2 * 0.5 * 0.3.
    assert result["umug"]["genotypes"][0]["frequency"] == pytest.approx(0.3)


def test_missing_locus_is_imputed_from_the_frequencies(service):
    without_dqb1 = (
        "A*01:01+A*02:01^B*07:02+B*08:01^C*07:01+C*07:02^DRB1*03:01+DRB1*15:01"
    )
    result = service.impute(
        without_dqb1, pop1=POP, pop2=POP, output="umug", expand=False
    )
    assert result["status"] == "ok"
    assert "DQB1*02:01+DQB1*06:02" in result["umug"]["genotypes"][0]["genotype"]


def test_genotype_not_in_the_frequencies_is_a_miss(service):
    unseen = (
        "A*24:02+A*68:01^B*07:02+B*08:01^C*07:01+C*07:02"
        "^DQB1*02:01+DQB1*06:02^DRB1*03:01+DRB1*15:01"
    )
    result = service.impute(unseen, pop1=POP, pop2=POP, output="umug", expand=False)
    assert result["status"] == "miss"
    assert result["umug"]["genotypes"] == []


def test_empty_input_is_rejected(service):
    with pytest.raises(ValueError, match="gl_string is required"):
        service.impute("  ")
