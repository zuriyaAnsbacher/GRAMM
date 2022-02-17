import os
from shutil import move
from GR_code.GG_GRAMM.code.runfile_update import run_GRAMM
from GR_code.GG_GRIMM.validation.runfile import run_GRIMM
from GR_code.GG_Post.code.create_results_updated import run_Post_GRIMM
from GR_code.GG_Post.code.visualization_updated import visualize


def run_all(input2GRAMM, alleles_names, files_address, res_100, is_serology, race_dict, open_ambiguity_sim):
    """
    called by server, run GRAMM, GRIMM, Post-GRAMM
    :param input2GRAMM: file with input from user (after processing)
    :param alleles_names: ['A', 'B', 'C', 'DRB1', 'DQB1']
    :param files_address: path to directory (that stores the output files)
    :param res_100: default: False. if the process do not find results to the family in first iteration, change to True
        and then GRIMM receives more options
    :param is_serology: if data is serology, TRue. otherwise, False
    :param race_dict: dict that contains the races of the family that the user inserted
    :param open_ambiguity_sim: in some simulations and Israeli-Australian file data, it keep the data of the ambiguity
        of the families, that have been removed before, to use it in GRAMM
    :return: path to results file, path to errors file, flag that=True if the process did not find any valid results
    """
    errors_in_families, aux_tools = run_GRAMM(input2GRAMM, files_address, alleles_names, race_dict,
                                              open_ambiguity_sim, is_serology)
    grimm_path = os.path.abspath("GR_code/GG_GRIMM")

    gl2GRIMM = files_address + '/glstring_for_GRIMM.txt'
    bin2GRIMM = files_address + '/binary_for_GRIMM.txt'

    move(gl2GRIMM, grimm_path + "/validation/simulation/data/input_test.txt")
    move(bin2GRIMM, grimm_path + "/validation/simulation/data/bin_input_test.txt")

    run_GRIMM(res_100)

    input2post = grimm_path + "/validation/output/don.hap.freqs"

    results, run_again = run_Post_GRIMM(input2post, input2GRAMM, alleles_names, files_address,
                                        aux_tools, errors_in_families)

    visualize(results, aux_tools['parent_has_empty_hap'], files_address)

    return results, files_address + '/errors.txt', run_again

