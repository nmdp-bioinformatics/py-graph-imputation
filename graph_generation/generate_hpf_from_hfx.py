"""Produce the HPF csv and population counts file from HFX bundles.

HFX (Haplotype Frequency Exchange) is a bundle format for haplotype frequencies
(https://github.com/societyforimmunepolymorphism/hfx). This module is the HFX
counterpart of ``generate_hpf.produce_hpf`` (which reads ``<POP>.freqs.gz``): it
reads ``<freq_data_dir>/<POP>.hfx`` for every population and writes the same two
artifacts, so the rest of the graph generation is unchanged:

* the HPF csv at ``freq_file`` with the columns ``hap,pop,freq``; and
* the population counts file at ``pops_count_file`` (``pop,count,ratio``).

Select it in the configuration file with ``"freq_source_type": "hfx"``
(``"freqs"``, the default, keeps the existing behaviour)::

    {
      "populations": ["hfset-1"],
      "freq_source_type": "hfx",
      "freq_data_dir": "data/hfx",
      ...
    }

Reading bundles needs the optional ``hfx-tools`` package
(``pip install "py-graph-imputation[hfx]"``).
"""

import csv
import pathlib


def _read_hfx_bundle(path):
    try:
        from hfx_tools.io import read_hfx_bundle
    except ImportError as e:
        raise ImportError(
            "Reading HFX bundles needs hfx-tools: "
            'pip install "py-graph-imputation[hfx]"'
        ) from e
    return read_hfx_bundle(pathlib.Path(path))


def sort_haplotype(haplotype):
    """Return ``haplotype`` with its alleles in alphabetical (locus) order.

    The imputation engine looks haplotypes up by ``"~".join(sorted(alleles))``, so
    haplotypes must be written in that order whatever order the bundle uses.
    """
    return "~".join(sorted(haplotype.split("~")))


def hfx_loci_map(hfx_path):
    """Return ``{locus: 1-based position}`` for the loci of an HFX bundle.

    Positions follow alphabetical locus order (``A, B, C, DQB1, DRB1``) to match
    :func:`sort_haplotype`, which is what the engine's lookups use. Loci come from
    ``metadata.outputResolution``, falling back to the first haplotype.
    """
    bundle = _read_hfx_bundle(hfx_path)
    loci = [
        entry.get("locus")
        for entry in bundle.metadata.get("outputResolution", [])
        if entry.get("locus")
    ]
    if not loci and bundle.rows:
        loci = [allele.split("*")[0] for allele in bundle.rows[0][0].split("~")]
    if not loci:
        raise ValueError(f"Could not determine loci from HFX bundle: {hfx_path}")
    return {locus: i + 1 for i, locus in enumerate(sorted(loci))}


def hfx_cohort_size(bundle):
    """Sample size from the bundle metadata, or ``None`` if it is not given.

    Prefers ``cohortDescription.cohortSize``, then a population's
    ``populationSize``. A size of 0 means "unknown".
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


def produce_hpf_from_hfx(conf, project_dir=""):
    """Convert each population's ``<POP>.hfx`` into the HPF csv and counts file.

    ``conf`` is the loaded configuration (a dict) and uses the same keys as
    ``produce_hpf``: ``populations``, ``freq_data_dir``, ``graph_files_path``,
    ``pops_count_file`` and ``freq_file``; ``project_dir`` is prefixed to them.

    The count written for a population is its cohort size, so the generator's
    ``freq_trim_threshold / count`` cutoff is an expected-count cutoff. A bundle
    with no cohort size (``cohortSize`` 0) counts as 1.0, which makes the trim an
    absolute frequency cutoff. Unlike ``produce_hpf``, population names may contain
    ``-``.
    """
    pops = conf.get("populations")
    freq_data_dir = pathlib.Path(project_dir + conf.get("freq_data_dir"))
    output_dir = pathlib.Path(project_dir + conf.get("graph_files_path"))
    pops_count_file = pathlib.Path(project_dir + conf.get("pops_count_file"))
    freq_file = pathlib.Path(project_dir + conf.get("freq_file"))

    output_dir.mkdir(parents=True, exist_ok=True)
    freq_file.parent.mkdir(parents=True, exist_ok=True)
    pops_count_file.parent.mkdir(parents=True, exist_ok=True)

    print("*" * 100)
    print("Conversion to HPF file from HFX bundles:")
    print("\tPopulation: {}".format(pops))
    print("\tHFX Directory: {}".format(freq_data_dir))
    print("\tOutput File: {}".format(freq_file))
    print("*" * 100)

    rows_by_pop = {}
    count_by_pop = {}
    for pop in pops:
        hfx_path = freq_data_dir / f"{pop}.hfx"
        if not hfx_path.is_file():
            raise FileNotFoundError(f"HFX bundle not found: {hfx_path}")
        print("Reading HFX Bundle:\t {}".format(hfx_path))
        bundle = _read_hfx_bundle(hfx_path)
        # Ignore rows with 0 frequency, as produce_hpf does.
        rows_by_pop[pop] = [
            (sort_haplotype(hap), freq)
            for hap, freq in bundle.rows
            if freq and freq > 0.0
        ]
        count_by_pop[pop] = hfx_cohort_size(bundle) or 1.0

    total = sum(count_by_pop.values())
    with open(pops_count_file, "w") as pop_ratio_file:
        for pop in pops:
            pop_ratio_file.write(
                "{},{},{}\n".format(pop, count_by_pop[pop], count_by_pop[pop] / total)
            )

    print("Writing hpf File:\t {}".format(freq_file))
    with open(freq_file, mode="w", newline="") as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=",", quoting=csv.QUOTE_NONE)
        csv_writer.writerow(["hap", "pop", "freq"])
        for pop in pops:
            for haplotype, freq in rows_by_pop[pop]:
                csv_writer.writerow([haplotype, pop, freq])
