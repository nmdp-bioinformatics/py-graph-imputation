import json
import os
import pandas as pd

from graph_generation.generate_hpf import produce_hpf
from grim.grim import graph_freqs
from grim.grim import impute
from filter_top_3 import change_donor_file
from filter_by_rest import change_output_by_extra_gl


def remove_empty_rows(file_path):
    df = pd.read_csv(file_path)

    df_cleaned = df.dropna(how="all")

    df_cleaned.to_csv(file_path, index=False)


def run_original_grim(
    path_configuration, hap_pop_pair=True, Producehpf=False, dominant3=True
):
    with open(path_configuration, "r") as f:
        config = json.load(f)

    # first step in py-graph-imputation
    if Producehpf:

        produce_hpf(conf_file=path_configuration)

        path_hpf = config["freq_file"]
        # remove empty rows from hpf otherwise doesnt work
        remove_empty_rows(path_hpf)

        # second step in py-graph-imputation
        graph_freqs(conf_file=path_configuration)

    # changing donor file to 3 most imporatnt gls and returning short_gl,extra_gl for each row in donor
    if dominant3:
        path_donor = config["imputation_in_file"]

        gls, lines = change_donor_file(path_donor)  # change so wont change donor file

    # imputation
    impute(conf_file=path_configuration, hap_pop_pair=hap_pop_pair)

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

        change_output_by_extra_gl(
            config, gls, path_pmug, path_umug, path_umug_pops, path_pmug_pops, path_miss
        )  # filter reasults in our origianl file, add miss to existing miss

        # changing to original donor file
        with open(path_donor, "w") as file:
            for line in lines:
                file.write(line)
        file.close()


if __name__ == "__main__":
    conf_file = "conf/minimal-configuration.json"
    run_original_grim(conf_file, True, True, True)
