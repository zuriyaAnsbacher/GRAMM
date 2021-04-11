import os
from shutil import move
from GR_code.GG_GRAMM.code.runfile import run_GRAMM
from GR_code.GG_GRIMM.validation.runfile import run_GRIMM
from GR_code.GG_Post.code.create_results import run_Post_GRIMM
from GR_code.GG_Post.code.visualization import visualize


def run_all(input2GRAMM, alleles, files_address, res_1000, is_serology, race_dict):
    gl2GRIMM, bin2GRIMM, errors, er_lst, double_als = run_GRAMM(input2GRAMM, alleles, files_address, is_serology, race_dict)
    grimm_path = os.path.abspath("GR_code/GG_GRIMM")

    move(gl2GRIMM, grimm_path + "/validation/simulation/data/input_test.txt")
    move(bin2GRIMM, grimm_path + "/validation/simulation/data/bin_input_test.txt")

    run_GRIMM(res_1000)

    input2post = grimm_path + "/validation/output/test.hap.freqs"

    results, run_again = run_Post_GRIMM(input2post, input2GRAMM, alleles, files_address, er_lst, double_als, is_serology)

    visualize(results, double_als, files_address)

    return results, errors, run_again

# path_res, path_er_pre = run_all("/home/dsi/zuriya/Projects/GR_Web/GR_code/GG_GRAMM/input/EASY.csv", ['A', 'B', 'C', 'DRB1', 'DQB1'], "/home/dsi/zuriya/Projects/GR_Web/output_train")


