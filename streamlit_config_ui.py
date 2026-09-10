#!/usr/bin/env python3
"""
Streamlit UI for configuring HLA imputation parameters and launching the API server.
"""

import json
import os
import subprocess
import streamlit as st
from pathlib import Path
from datetime import datetime

# Configuration
BASE_DIR = Path(__file__).parent
CONF_DIR = BASE_DIR / "conf"
DATA_DIR = BASE_DIR / "data"
FREQS_DIR = DATA_DIR / "freqs"

# Ensure directories exist
CONF_DIR.mkdir(exist_ok=True)


def get_available_populations():
    """Scan data/freqs directory for available population frequency files."""
    if not FREQS_DIR.exists():
        return []

    populations = set()

    # Scan top-level items
    for item in FREQS_DIR.iterdir():
        if item.is_file() and item.suffix in [".gz", ".freqs"]:
            pop_name = item.stem.replace(".freqs", "").upper()
            populations.add(pop_name)
        elif item.is_dir() and not item.name.startswith("."):
            populations.add(item.name.upper())
            # Also check subdirectories for frequency files
            for subitem in item.rglob("*.freqs*"):
                if subitem.is_file():
                    subpop = subitem.stem.replace(".freqs", "").upper()
                    populations.add(subpop)

    return sorted(populations)


def load_template_config():
    """Load the minimal configuration as a template."""
    template_path = CONF_DIR / "minimal-configuration.json"
    if template_path.exists():
        with open(template_path) as f:
            return json.load(f)
    return {}


def save_config(config_dict, config_name):
    """Save configuration to a JSON file in conf/ directory."""
    config_path = CONF_DIR / f"{config_name}.json"
    with open(config_path, "w") as f:
        json.dump(config_dict, f, indent=2)
    return config_path


def main():
    st.set_page_config(page_title="HLA Imputation Config", layout="wide")
    st.title("HLA Imputation Configuration UI")

    # Sidebar for general settings
    with st.sidebar:
        st.header("Configuration Settings")
        config_name = st.text_input(
            "Configuration Name",
            value=f"config_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
            help="Name for this configuration (will be saved as JSON in conf/)",
        )

        # Show available populations for debugging
        with st.expander("ℹ️ Available Populations"):
            avail_pops = get_available_populations()
            if avail_pops:
                st.write(f"Found {len(avail_pops)} population(s):")
                for pop in avail_pops:
                    st.code(pop)
            else:
                st.warning("No populations found in data/freqs/")

    # Load template
    template = load_template_config()

    # Main configuration form
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Data Selection")

        # Population selection
        available_pops = get_available_populations()
        if not available_pops:
            st.error("❌ No populations found in data/freqs/")

        selected_populations = st.multiselect(
            "Populations (select one or more)",
            options=available_pops,
            default=available_pops[:1] if available_pops else [],
            help="Select one or more populations from available frequency files. Use Ctrl/Cmd+Click to multi-select.",
        )

        if not selected_populations:
            st.warning("⚠️ Please select at least one population")

        # Loci configuration
        st.subheader("Loci Configuration")
        all_loci = template.get("loci_map", {})
        loci_options = list(all_loci.keys())

        col_loci1, col_loci2 = st.columns(2)
        with col_loci1:
            input_loci = st.multiselect(
                "Input Loci",
                options=loci_options,
                default=loci_options,
                help="Loci to be used as input",
            )

        with col_loci2:
            output_loci = st.multiselect(
                "Output Loci",
                options=loci_options,
                default=loci_options,
                help="Loci to be imputed",
            )

        # Phasing
        phasing = st.radio(
            "Phasing",
            options=["phased", "unphased"],
            help="Whether input is phased or unphased",
        )

    with col2:
        st.subheader("Algorithm Parameters")

        num_results = st.number_input(
            "Number of Results",
            value=template.get("number_of_results", 1000),
            min_value=1,
            help="Maximum number of imputation results to return",
        )

        num_pop_results = st.number_input(
            "Number of Population Results",
            value=template.get("number_of_pop_results", 100),
            min_value=1,
            help="Maximum population results per combination",
        )

        epsilon = st.number_input(
            "Epsilon (ε)",
            value=template.get("epsilon", 1e-3),
            format="%.6f",
            help="Convergence threshold for the algorithm",
        )

        frequency_trim = st.number_input(
            "Frequency Trim Threshold",
            value=template.get("freq_trim_threshold", 0.1),
            min_value=0.0,
            max_value=1.0,
            step=0.01,
            help="Trim frequencies below this threshold",
        )

        st.subheader("Plan B Strategy")
        use_plan_b = st.checkbox(
            "Enable Plan B",
            value=template.get("planb", True),
            help="Use plan B fallback when optimization fails",
        )

        factor_missing = st.number_input(
            "Factor Missing Data",
            value=template.get("factor_missing_data", 0.0001),
            format="%.6f",
            help="Weight factor for missing data in Plan B",
        )

        st.subheader("Priority Settings")
        st.write("Priority weights (should sum to 1):")

        col_p1, col_p2 = st.columns(2)
        with col_p1:
            alpha = st.number_input("Alpha (α)", value=0.4999999, format="%.7f")
            beta = st.number_input("Beta (β)", value=1e-7, format="%.7f")
        with col_p2:
            delta = st.number_input("Delta (δ)", value=0.4999999, format="%.7f")
            gamma = st.number_input("Gamma (γ)", value=1e-7, format="%.7f")

        eta = st.number_input("Eta (η)", value=0.0, format="%.7f")

    # Display preview
    st.divider()
    st.subheader("Configuration Preview")

    # Build the configuration object
    config_obj = {
        "populations": selected_populations,
        "freq_trim_threshold": frequency_trim,
        "priority": {
            "alpha": alpha,
            "eta": eta,
            "beta": beta,
            "gamma": gamma,
            "delta": delta,
        },
        "UNK_priors": template.get("UNK_priors", "SR"),
        "FULL_LOCI": "".join(sorted(set(input_loci + output_loci))),
        "loci_map": {
            k: v for k, v in all_loci.items() if k in (input_loci + output_loci)
        },
        "factor_missing_data": factor_missing,
        "Plan_B_Matrix": template.get("Plan_B_Matrix", []),
        "planb": use_plan_b,
        "number_of_options_threshold": template.get(
            "number_of_options_threshold", 100000
        ),
        "epsilon": epsilon,
        "number_of_results": num_results,
        "number_of_pop_results": num_pop_results,
        "output_MUUG": True,
        "output_haplotypes": True,
        "freq_data_dir": "data/freqs",
        "freq_file": "output/hpf.csv",
        "graph_files_path": "output/csv/",
        "node_csv_file": "nodes.csv",
        "edges_csv_file": "edges.csv",
        "info_node_csv_file": "info_node.csv",
        "top_links_csv_file": "top_links.csv",
        "imputation_in_file": "data/subjects/donor.csv",
        "imputation_out_umug_freq_filename": "imputation_output/don.umug",
        "imputation_out_umug_pops_filename": "imputation_output/don.umug.pops",
        "imputation_out_hap_freq_filename": "imputation_output/don.pmug",
        "imputation_out_hap_pops_filename": "imputation_output/don.pmug.pops",
        "imputation_out_umug_mc_filename": "imputation_output/don.umug.mc",
        "imputation_out_hap_mc_filename": "imputation_output/don.pmug.mc",
        "imputation_out_umug_trs_filename": "imputation_output/don.umug.trs",
        "imputation_out_hap_trs_filename": "imputation_output/don.pmug.trs",
        "phasing": phasing,
    }

    st.json(config_obj)

    # Action buttons
    st.divider()
    col_btn1, col_btn2, col_btn3 = st.columns(3)

    with col_btn1:
        if st.button("💾 Save Configuration", use_container_width=True):
            try:
                config_path = save_config(config_obj, config_name)
                st.success(f"✅ Configuration saved to {config_path}")
                st.session_state.last_config = config_name
                st.session_state.last_config_path = str(config_path)
            except Exception as e:
                st.error(f"❌ Error saving configuration: {e}")

    with col_btn2:
        if st.button("🚀 Launch API Server", use_container_width=True):
            # First save the config
            try:
                config_path = save_config(config_obj, config_name)
                st.success(f"✅ Configuration saved to {config_path}")

                st.info(
                    "📋 **To launch the API server, run this command in your terminal:**"
                )
                st.code(
                    f"python3 api_server.py --config {config_name}.json --port 5000",
                    language="bash",
                )
                st.info("The server will be available at `http://localhost:5000`")

            except Exception as e:
                st.error(f"❌ Error: {e}")

    with col_btn3:
        if st.button("📋 View Saved Configs", use_container_width=True):
            configs = list(CONF_DIR.glob("*.json"))
            if configs:
                st.write("**Available configurations:**")
                for cfg in sorted(configs):
                    st.code(cfg.name)
            else:
                st.info("No saved configurations found")


if __name__ == "__main__":
    main()
