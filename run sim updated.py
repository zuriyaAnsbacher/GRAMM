import os
import csv
from shutil import move
import pandas as pd
from GR_code.GG_GRAMM.code.runfile_update import run_GRAMM
from GR_code.GG_GRIMM.validation.runfile import run_GRIMM
from GR_code.GG_Post.code.create_results_updated import run_Post_GRIMM
from processing_data2input.processing_update import basic_csv2format


def clean_gl(hap):
    hap = hap.replace('g', '').split('~')
    hap.sort()  # check. need to be in this order: A,B,C,DQB1,DRB1
    hap = '~'.join(hap)
    return hap


def add_race_CAU(in_file, out_file):
    in_f = open(in_file, 'r')
    out_f = open(out_file, 'w')

    reader = csv.reader(in_f)
    writer = csv.writer(out_f)

    csv.field_size_limit(100000000)  # add it because limitation in csv

    for row in reader:
        writer.writerow([row[0], row[1], 'CAU', 'CAU'])

    in_f.close()
    out_f.close()


def add_children_haps(fam_dict, inheritance):
    # example for inheritance: C100001=F1~M1:C100002=F1~M2:C100003=F2~M2:C100004=F1~M1:C100005=F2~M2
    children = inheritance.split(";")  # [C100001=F1~M1, C100002=F1~M2 ...]
    for child in children:
        if 'X' not in child:  # todo: change it?
            id_child = child.split("=")[0].lstrip('C')  # if child is C100001=F1~M1, so his id is 100001
            f_hap_idx = child.split("=")[1].split("~")[0].lstrip('F')  # if child is C100001=F1~M1, so f_hap_idx is 1
            m_hap_idx = child.split("=")[1].split("~")[1].lstrip('M')  # if child is C100001=F1~M1, so m_hap_idx is 1

            hap_1 = fam_dict['F'][int(f_hap_idx) - 1]  # the "-1" because hap1 is in index 0, and hap2 is in index 1
            hap_2 = fam_dict['M'][int(m_hap_idx) - 1]

            fam_dict[id_child] = [hap_1, hap_2]


def check_if_exist_pred_equal_to_sim_GRAMM(sim_file, gramm_pred, skip_sim_headers):
    conclusion = {'right': 0, 'error': 0}

    sim, pred = open(sim_file, 'r'), open(gramm_pred, 'r')
    sim_reader, pred_reader = csv.reader(sim), csv.reader(pred)

    csv.field_size_limit(100000000)  # add it because limitation in csv

    next(pred_reader)  # skip headers
    if skip_sim_headers:  # in HARD and HARDER sim, there are no headers, so skip_sim_headers=False
        next(sim_reader)  # skip headers

    pred_dict = {}

    id_fam = None
    for line in pred_reader:
        if line[0] != id_fam:
            idx_option = 0
        id_fam = line[0]
        if id_fam not in pred_dict:
            pred_dict[id_fam] = []
        pred_dict[id_fam].append({})
        pred_dict[id_fam][idx_option]['F'] = [line[1], line[2]]
        pred_dict[id_fam][idx_option]['M'] = [line[3], line[4]]
        inheritance = line[5]
        add_children_haps(pred_dict[id_fam][idx_option], inheritance)
        idx_option += 1

    for line in sim_reader:
        id_fam = line[0].lstrip('0')  # lstrip because in the preds file, this is the foramt
        if id_fam in pred_dict:
            id_indv = line[1]

            if id_indv[0] == '0':
                id_indv = 'F' if int(id_indv.lstrip('0')) % 2 == 1 else 'M'

            correct_pred_exist = False
            hap1_sim, hap2_sim = clean_gl(line[6]), clean_gl(line[8])
            for option in pred_dict[id_fam]:
                hap1_pred, hap2_pred = clean_gl(option[id_indv][0]), clean_gl(option[id_indv][1])
                if pred_haps_VS_sim_haps(hap1_sim, hap2_sim, hap1_pred, hap2_pred):
                    conclusion['right'] += 1
                    correct_pred_exist = True
                    break  # we just need to check if there is 1 correct option at least
            if not correct_pred_exist:
                conclusion['error'] += 1

    return conclusion


def check_if_first_pred_equal_to_sim_GRAMM(sim_file, gramm_pred, skip_sim_headers):
    conclusion = {'right': 0, 'error': 0}

    sim, pred = open(sim_file, 'r'), open(gramm_pred, 'r')
    sim_reader, pred_reader = csv.reader(sim), csv.reader(pred)

    csv.field_size_limit(100000000)  # add it because limitation in csv

    next(pred_reader)  # skip headers
    if skip_sim_headers:  # in HARD and HARDER sim, there are no headers, so skip_sim_headers=False
        next(sim_reader)  # skip headers

    pred_dict = {}

    # build dict that contains the first result from preds
    for line in pred_reader:
        id_fam = line[0]
        if id_fam not in pred_dict:  # save first result only
            pred_dict[id_fam] = {}
            # TODO: check if always there are all parents' haps
            pred_dict[id_fam]['F'] = [line[1], line[2]]
            pred_dict[id_fam]['M'] = [line[3], line[4]]
            inheritance = line[5]
            add_children_haps(pred_dict[id_fam], inheritance)  # todo: maybe no need to check children?

    for line in sim_reader:
        id_fam = line[0].lstrip('0')  # lstrip because in the preds file, this is the foramt
        if id_fam in pred_dict:
            id_indv = line[1]
            """
            in father and mother, the id_indv not equal to their keys in the dict.
            so we check if the current person is father or mother by look in the first digit in their id_indv:
            in father anf mother it's 0, in children it's 1
            then, we check if it's father or mother by using modulo, cause fathers ID's are odd, and mother's are even 
            and then we convert id_indv to 'F' or 'M'
            """
            if id_indv[0] == '0':
                id_indv = 'F' if int(id_indv.lstrip('0')) % 2 == 1 else 'M'

            hap1_sim, hap2_sim = line[6], line[8]
            hap1_pred, hap2_pred = pred_dict[id_fam][id_indv][0], pred_dict[id_fam][id_indv][1]
            hap1_sim, hap2_sim, hap1_pred, hap2_pred = clean_gl(hap1_sim), clean_gl(hap2_sim), clean_gl(
                hap1_pred), clean_gl(hap2_pred)
            if pred_haps_VS_sim_haps(hap1_sim, hap2_sim, hap1_pred, hap2_pred):
                conclusion['right'] += 1
            else:
                conclusion['error'] += 1

    sim.close()
    pred.close()

    return conclusion


def pred_haps_VS_sim_haps(hap1_sim, hap2_sim, hap1_pred, hap2_pred):
    if hap1_sim == hap1_pred and hap2_sim == hap2_pred:
        return True
    elif hap1_sim == hap2_pred and hap2_sim == hap1_pred:
        return True
    elif hap1_sim == hap1_pred and hap2_pred == 'Unknown':
        return True
    elif hap1_sim == hap2_pred and hap1_pred == 'Unknown':
        return True
    elif hap2_sim == hap1_pred and hap2_pred == 'Unknown':
        return True
    elif hap2_sim == hap2_pred and hap1_pred == 'Unknown':
        return True
    else:
        return False


def check_if_first_pred_equal_to_sim_GRIMM(sim_file, grimm_pred, skip_sim_headers):
    conclusion = {'right': 0, 'error': 0}

    sim, pred = open(sim_file, 'r'), open(grimm_pred, 'r')
    sim_reader, pred_reader = csv.reader(sim), csv.reader(pred)

    csv.field_size_limit(100000000)  # add it because limitation in csv

    if skip_sim_headers:
        next(sim_reader, None)  # skip headers

    pred_dict = {}

    # build dict that contains the first result from preds
    for line in pred_reader:
        if line[0] not in pred_dict:  # save first result only
            pred_dict[line[0]] = [line[1], line[2]]  # hap1, hap2

    for line in sim_reader:
        id_fam_indv = ''.join([line[0].lstrip('0'), '~', line[1].lstrip('0')])

        hap1_sim, hap2_sim = clean_gl(line[6]), clean_gl(line[8])

        if id_fam_indv in pred_dict:
            hap1_pred, hap2_pred = pred_dict[id_fam_indv][0], pred_dict[id_fam_indv][1]

            if (hap1_sim == hap1_pred and hap2_sim == hap2_pred) or (hap1_sim == hap2_pred and hap2_sim == hap1_pred):
                conclusion['right'] += 1
            else:
                conclusion['error'] += 1


    sim.close()
    pred.close()

    return conclusion


def check_if_exist_pred_equal_to_sim_GRIMM(sim_file, grimm_pred, skip_sim_headers):
    conclusion = {'right': 0, 'error': 0}

    sim, pred = open(sim_file, 'r'), open(grimm_pred, 'r')
    sim_reader, pred_reader = csv.reader(sim), csv.reader(pred)

    csv.field_size_limit(100000000)  # add it because limitation in csv

    if skip_sim_headers:
        next(sim_reader, None)  # skip headers

    pred_dict = {}

    # build dict that contains the first result from preds
    for line in pred_reader:
        if line[0] not in pred_dict:  # save first result only
            pred_dict[line[0]] = [[line[1], line[2]]]  # hap1, hap2
        else:
            pred_dict[line[0]].append([line[1], line[2]])

    for line in sim_reader:
        id_fam_indv = ''.join([line[0].lstrip('0'), '~', line[1].lstrip('0')])

        # hap1_sim, hap2_sim = clean_gl(line[6]), clean_gl(line[8])
        # if id_fam_indv in pred_dict:
        #     hap1_pred, hap2_pred = pred_dict[id_fam_indv][0], pred_dict[id_fam_indv][1]
        #
        #     if (hap1_sim == hap1_pred and hap2_sim == hap2_pred) or (hap1_sim == hap2_pred and hap2_sim == hap1_pred):
        #         conclusion['right'] += 1
        #     else:
        #         conclusion['error'] += 1

        if id_fam_indv in pred_dict:
            correct_pred_exist = False
            hap1_sim, hap2_sim = clean_gl(line[6]), clean_gl(line[8])
            for option in pred_dict[id_fam_indv]:
                hap1_pred, hap2_pred = clean_gl(option[0]), clean_gl(option[1])
                if (hap1_sim == hap1_pred and hap2_sim == hap2_pred) or (hap1_sim == hap2_pred and hap2_sim == hap1_pred):
                    conclusion['right'] += 1
                    correct_pred_exist = True
                    break  # we just need to check if there is 1 correct option at least
            if not correct_pred_exist:
                conclusion['error'] += 1

    sim.close()
    pred.close()

    return conclusion


def compare_to_true_what_freq_of_right_res(sim_file, gramm_pred, skip_sim_headers):
    pass


def run_gramm_code(sim_path):
    reformat_file = f'running for GRAMM paper/files in runtime/reformat_{os.path.basename(sim_path)}'  # check path
    alleles_names = ['A', 'B', 'C', 'DRB1', 'DQB1']
    files_address = 'running for GRAMM paper/files in runtime'
    race_dict = {}

    basic_csv2format(sim_path, reformat_file, race_dict)  # if update_preprocessing, add race_dict

    errors_in_families, aux_tools = run_GRAMM(reformat_file, files_address, alleles_names, race_dict,
                                              is_serology=False)

    grimm_path = "GR_code/GG_GRIMM"

    gl2GRIMM = files_address + '/glstring_for_GRIMM.txt'
    bin2GRIMM = files_address + '/binary_for_GRIMM.txt'

    move(gl2GRIMM, grimm_path + "/validation/simulation/data/input_test.txt")
    move(bin2GRIMM, grimm_path + "/validation/simulation/data/bin_input_test.txt")

    run_GRIMM(res_1000=False, sim=False)

    input2post = grimm_path + "/validation/output/test.hap.freqs"

    results, run_again = run_Post_GRIMM(input2post, reformat_file, alleles_names, files_address, aux_tools,
                                        errors_in_families)

    return results


def get_sim_file_with_persons_passed_successfully_gramm_running(orig_sim, out_from_run, new_for_grimm_run,
                                                                skip_headers):
    sim = open(orig_sim, 'r')
    sim_reader = csv.reader(sim)
    new_file = open(new_for_grimm_run, 'w')
    writer = csv.writer(new_file)
    df = pd.read_csv(out_from_run)
    success_idx = list(set(list(df['FAMCODE'])))

    if skip_headers:
        next(sim_reader)  # skip headers

    csv.field_size_limit(100000000)  # add it because limitation in csv

    for row in sim_reader:
        fam_id = row[0].lstrip('0')
        if int(fam_id) in success_idx:
            indv_id = row[1].lstrip('0')
            id_fam_indv = ''.join([fam_id, '~', indv_id])
            writer.writerow([id_fam_indv, row[5], 'CAU', 'CAU'])

    sim.close()
    new_file.close()


def main():
    cur_sim = "MM"
    skip_headers = True
    # if cur_sim in ["HARD", "HARDER", 'HARD_MINI']:
    #     skip_headers = False  # no header in HARD and HARDER simulations files
    # deal the error - flag about not check them
    # results = run_gramm_code(f'running for GRAMM paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv')
    # results = run_gramm_code(f'running for GRAMM paper/cordmom/cordmom_mini/cordmom_{cur_sim}.csv')

    # get_sim_file_with_persons_passed_successfully_gramm_running(f'running for GRAMM paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv',
    #                                                                        f'running for GRAMM paper/output_from_web_running/new_version/results_{cur_sim}.csv',
    #                                                                        f'running for GRAMM paper/files_to_sapir_for_grimm_run/new_version/{cur_sim}.csv', skip_headers)

    # add_race_CAU(f'./files_to_sapir_for_grimm_run/{cur_sim}.csv',  f'./files_to_sapir_for_grimm_run/{cur_sim}_new.csv')

    # move('running for GRAMM paper/files in runtime/results.csv',
    #      f'running for GRAMM paper/output_from_web_running/new_version/results_{cur_sim}.csv')
    # move('running for GRAMM paper/files in runtime/errors.txt',
    #      f'running for GRAMM paper/output_from_web_running/new_version/errors_{cur_sim}.txt')

    # move('running for GRAMM paper/files in runtime/results.csv',
    #      f'running for GRAMM paper/output_from_web_running/cordmom_mini/results_cordmom_{cur_sim}.csv')
    # move('running for GRAMM paper/files in runtime/errors.txt',
    #      f'running for GRAMM paper/output_from_web_running/cordmom_mini/errors_cordmom_{cur_sim}.txt')


    """
    >>> GRAMM - first and exist
    """
    # right_and_error_gramm = check_if_first_pred_equal_to_sim_GRAMM(
    #     f'running for GRAMM paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv',
    #     f'running for GRAMM paper/output_from_web_running/new_version/results_{cur_sim}.csv', skip_headers)
    # #
    # right_and_error_gramm_exist = check_if_exist_pred_equal_to_sim_GRAMM(
    #     f'running for GRAMM paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv',
    #     f'running for GRAMM paper/output_from_web_running/new_version/results_{cur_sim}.csv', skip_headers)
    #
    # print(right_and_error_gramm)
    # print(right_and_error_gramm_exist)


    """ 
    >> GRIMM - first and exist
    """
    right_and_error_grimm = check_if_first_pred_equal_to_sim_GRIMM(f'running for GRAMM paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv',
                                                             f'running for GRAMM paper/grimm_output_from_sapir/{cur_sim}.hap.freqs',
                                                             skip_headers)
    right_and_error_grimm_exist = check_if_exist_pred_equal_to_sim_GRIMM(f'running for GRAMM paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv',
                                                             f'running for GRAMM paper/grimm_output_from_sapir/{cur_sim}.hap.freqs',
                                                             skip_headers)
    print(right_and_error_grimm)
    print(right_and_error_grimm_exist)

main()


