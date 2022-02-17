import os
from shutil import move
from GR_code.GG_GRAMM.code.runfile_update import run_GRAMM
from GR_code.GG_GRIMM.validation.runfile import run_GRIMM
from GR_code.GG_Post.code.create_results_updated import run_Post_GRIMM
from GR_code.GG_Post.code.visualization_updated import visualize


def run_all(input2GRAMM, alleles_names, files_address, res_1000, is_serology, race_dict, open_ambiguity_sim):
    errors_in_families, aux_tools = run_GRAMM(input2GRAMM, files_address, alleles_names, race_dict,
                                              open_ambiguity_sim, is_serology)
    grimm_path = os.path.abspath("GR_code/GG_GRIMM")

    gl2GRIMM = files_address + '/glstring_for_GRIMM.txt'
    bin2GRIMM = files_address + '/binary_for_GRIMM.txt'

    move(gl2GRIMM, grimm_path + "/validation/simulation/data/input_test.txt")
    move(bin2GRIMM, grimm_path + "/validation/simulation/data/bin_input_test.txt")

    run_GRIMM(res_1000)

    input2post = grimm_path + "/validation/output/don.hap.freqs"

    results, run_again = run_Post_GRIMM(input2post, input2GRAMM, alleles_names, files_address,
                                        aux_tools, errors_in_families)

    visualize(results, aux_tools['parent_has_empty_hap'], files_address)

    return results, files_address + '/errors.txt', run_again

