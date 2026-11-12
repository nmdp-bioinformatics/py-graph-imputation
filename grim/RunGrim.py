import json
import os
import pandas as pd

from graph_generation.generate_hpf import produce_hpf
from grim.grim import graph_freqs
from grim.grim import impute
from filter_top_3 import change_donor_file
from filter_by_rest import change_output_by_extra_gl
from .imputation.impute import Imputation


def remove_empty_rows(file_path):
    df = pd.read_csv(file_path)

    df_cleaned = df.dropna(how="all")

    df_cleaned.to_csv(file_path, index=False)


def run_original_grim(
    path_configuration,
    hap_pop_pair=True,
    Producehpf=False,
    dominant3=True,
):
    """Run the whole GRIM pipeline off `path_configuration`.

    `num_processes` in the configuration is how many subjects to impute at a
    time - 1 (the default) for a single process. Whatever the count, the
    workers all impute against the same graph, which is built once here.

    With `dominant3`, the loci `change_donor_file` holds aside are checked
    inside the imputation, against the full candidate set. Setting
    `filter_extra_gl_before_truncation` to false in the configuration puts that
    check back where it was - a pass over `don.pmug` after the imputation, by
    which point the candidates have already been cut to `number_of_results`.
    """
    with open(path_configuration, "r") as f:
        config = json.load(f)

    filter_before_truncation = config.get("filter_extra_gl_before_truncation", True)

    # first step in py-graph-imputation
    if Producehpf:

        produce_hpf(conf_file=path_configuration)

        path_hpf = config["freq_file"]
        # remove empty rows from hpf otherwise doesnt work
        remove_empty_rows(path_hpf)

        # second step in py-graph-imputation
        graph_freqs(conf_file=path_configuration)

    # changing donor file to 3 most imporatnt gls and returning short_gl,extra_gl for each row in donor
    extra_gl_by_id = None
    if dominant3:
        path_donor = config["imputation_in_file"]

        gls, lines = change_donor_file(path_donor)  # change so wont change donor file

        if filter_before_truncation:
            # Hand the held-aside loci to the imputation, which checks its own
            # results against them before it truncates and so needs no pass
            # afterwards.
            extra_gl_by_id = {
                str(subject_id): gls["extra_gl"][idx]
                for idx, subject_id in enumerate(gls["subject_id"])
            }

    # imputation
    impute(
        conf_file=path_configuration,
        hap_pop_pair=hap_pop_pair,
        extra_gl_by_id=extra_gl_by_id,
    )

    # change the output and filter by the extra_gl
    if dominant3:
        path_pmug = os.path.join(
            config["imputation_out_path"], config["imputation_out_hap_freq_filename"]
        )
        path_umug = os.path.join(
            config["imputation_out_path"], config["imputation_out_umug_freq_filename"]
        )
        path_umug_pops = os.path.join(
            config["imputation_out_path"], config["imputation_out_umug_pops_filename"]
        )
        path_pmug_pops = os.path.join(
            config["imputation_out_path"], config["imputation_out_hap_pops_filename"]
        )
        path_miss = os.path.join(
            config["imputation_out_path"], config["imputation_out_miss_filename"]
        )

        if not filter_before_truncation:
            change_output_by_extra_gl(
                config,
                gls,
                path_pmug,
                path_umug,
                path_umug_pops,
                path_pmug_pops,
                path_miss,
            )  # filter reasults in our origianl file, add miss to existing miss

        # changing to original donor file
        with open(path_donor, "w") as file:
            for line in lines:
                file.write(line)
        file.close()


class Impute_for_em(object):
    def __init__(
        self,
        config=None,
        graph=None,
        count_by_prob=None,
    ):
        self.imputation = Imputation(graph, config, count_by_prob)

    def impute_for_em(self, config, planb, em_mr, em=True, dominant3=True):
        if dominant3:
            path_donor = config["imputation_input_file"]
            gls, lines = change_donor_file(
                path_donor
            )  # change so wont change donor file

        # imputation
        self.imputation.impute_file(config, planb, True, em)

        # change the output and filter by the extra_gl
        if dominant3:
            path_pmug = os.path.join(config["imputation_out_hap_freq_file"])
            path_pmug_pops = os.path.join(config["imputation_out_hap_pops_file"])
            path_miss = os.path.join(config["imputation_out_miss_file"])

            change_output_by_extra_gl(
                config, gls, path_pmug, None, None, path_pmug_pops, path_miss
            )  # filter results in our original file, add miss to existing miss

            # changing to original donor file
            with open(path_donor, "w") as file:
                for line in lines:
                    file.write(line)
            file.close()


if __name__ == "__main__":
    conf_file = "conf/minimal-configuration.json"
    run_original_grim(conf_file, True, True, True)
