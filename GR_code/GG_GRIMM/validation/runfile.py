import argparse
import cProfile
import json
import pickle
import pathlib
import sys
import os

sys.path.insert(0, os.path.join(".."))

from GR_code.GG_GRIMM.imputation.imputegl import Imputation
from GR_code.GG_GRIMM.imputation.imputegl.networkx_graph import Graph


def run_GRIMM(res_100, sim=False):

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config",
                        required=False,
                        default="GR_code/GG_GRIMM/conf/haplogic/5loci-don-configuration.json",
                        help="Configuration JSON file",
                        type=str)
    args = parser.parse_args()
    configuration_file = args.config
    project_dir = "GR_code/GG_GRIMM/"
    output_dir = "GR_code/GG_GRIMM/validation/output/"

    # Read configuration file and load properties
    with open(configuration_file) as f:
        json_conf = json.load(f)

    config = {
        "planb": json_conf.get('planb', True),
        "pops": json_conf.get('populations'),
        "priority": json_conf.get('priority'),
        "epsilon": json_conf.get('epsilon', 1e-3),
        "number_of_results": json_conf.get('number_of_results', 1000),
        "number_of_pop_results": json_conf.get('number_of_pop_results', 100),
        "output_MUUG": json_conf.get("output_MUUG", True),
        "output_haplotypes": json_conf.get("output_haplotypes", False),
        "node_file": project_dir + json_conf.get("node_csv_file"),
        "top_links_file": project_dir + json_conf.get("top_links_csv_file"),
        "edges_file": project_dir + json_conf.get("edges_csv_file"),
        "imputation_input_file": project_dir + json_conf.get("imputation_in_file"),
        "bin_imputation_input_file": project_dir + json_conf.get("bin_imputation_in_file", "None"),
        "imputation_out_umug_freq_file": output_dir + json_conf.get("imputation_out_umug_freq_filename"),
        "imputation_out_umug_pops_file": output_dir + json_conf.get("imputation_out_umug_pops_filename"),
        "imputation_out_hap_freq_file": output_dir + json_conf.get("imputation_out_hap_freq_filename"),
        "imputation_out_hap_pops_file": output_dir + json_conf.get("imputation_out_hap_pops_filename"),
        "imputation_out_miss_file": output_dir + json_conf.get("imputation_out_miss_filename"),
        "imputation_out_problem_file": output_dir + json_conf.get("imputation_out_problem_filename"),
        "factor_missing_data": json_conf.get("factor_missing_data", 0.01),
        "loci_map": json_conf.get("loci_map", {"A": 1, "B": 3, "C": 2, "DQB1": 4, "DRB1": 5}),
        "matrix_planb": json_conf.get("Plan_B_Matrix", [
            [[1, 2, 3, 4, 5]],
            [[1, 2, 3], [4, 5]],
            [[1], [2, 3], [4, 5]],
            [[1, 2, 3], [4], [5]],
            [[1], [2, 3], [4], [5]],
            [[1], [2], [3], [4], [5]]
        ]),
        "pops_count_file": project_dir + json_conf.get("pops_count_file", ''),
        "use_pops_count_file": json_conf.get("pops_count_file", False),
        "number_of_options_threshold": json_conf.get("number_of_options_threshold", 100000),
        "max_haplotypes_number_in_phase": json_conf.get("max_haplotypes_number_in_phase", 100),
        "nodes_for_plan_A": json_conf.get("Plan_A_Matrix", []),
        "save_mode": json_conf.get("save_space_mode", False)

    }

    if res_100:  # 100 results from grimm (only if 10 results failed_count with inconsistent)
        config["number_of_results"] = 100

    # Display the configurations we are using
    print('****************************************************************************************************')
    print("Performing imputation based on:")
    print("\tPopulation: {}".format(config["pops"]))
    print("\tPriority: {}".format(config["priority"]))
    print("\tEpsilon: {}".format(config["epsilon"]))
    print("\tPlan B: {}".format(config["planb"]))
    print("\tNumber of Results: {}".format(config["number_of_results"]))
    print("\tNumber of Population Results: {}".format(config["number_of_pop_results"]))
    print("\tNodes File: {}".format(config["node_file"]))
    print("\tTop Links File: {}".format(config["edges_file"]))
    print("\tInput File: {}".format(config["imputation_input_file"]))
    print("\tOutput UMUG Format: {}".format(config["output_MUUG"]))
    print("\tOutput UMUG Freq Filename: {}".format(config["imputation_out_umug_freq_file"]))
    print("\tOutput UMUG Pops Filename: {}".format(config["imputation_out_umug_pops_file"]))
    print("\tOutput Haplotype Format: {}".format(config["output_haplotypes"]))
    print("\tOutput HAP Freq Filename: {}".format(config["imputation_out_hap_freq_file"]))
    print("\tOutput HAP Pops Filename: {}".format(config["imputation_out_hap_pops_file"]))
    print("\tOutput Miss Filename: {}".format(config["imputation_out_miss_file"]))
    print("\tOutput Problem Filename: {}".format(config["imputation_out_problem_file"]))
    print("\tFactor Missing Data: {}".format(config["factor_missing_data"]))
    print("\tLoci Map: {}".format(config["loci_map"]))
    print("\tPlan B Matrix: {}".format(config["matrix_planb"]))
    print("\tPops Count File: {}".format(config["pops_count_file"]))
    print("\tUse Pops Count File: {}".format(config["use_pops_count_file"]))
    print("\tNumber of Options Threshold: {}".format(config["number_of_options_threshold"]))
    print("\tMax Number of haplotypes in phase: {}".format(config["max_haplotypes_number_in_phase"]))
    if config["nodes_for_plan_A"]:
        print("\tNodes in plan A: {}".format(config["nodes_for_plan_A"]))
    print("\tSave space mode: {}".format(config["save_mode"]))
    print('****************************************************************************************************')

    all_loci_set = set()
    for _, val in config["loci_map"].items():
        all_loci_set.add(str(val))

    config["full_loci"] = ''.join(sorted(all_loci_set))

    """ 
    at the first time you run the code, need to run the part of creating graph. 
    then, can save it to pickle and in the next time to load the pickle (to save time) 
    """
    # Perform imputation
    graph = Graph(config)
    graph.build_graph(config["node_file"], config["top_links_file"], config["edges_file"])
    # pickle.dump(graph, open('../graph.pkl', 'wb'))

    # graph = pickle.load(open(project_dir + 'graph.pkl', 'rb'))
    imputation = Imputation(graph, config)

    # Create output directory if it doesn't exist
    pathlib.Path(output_dir).mkdir(parents=False, exist_ok=True)

    # Write out the results from imputation
    imputation.impute_file(config)

