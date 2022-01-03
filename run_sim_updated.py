import os
import csv
import time
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
            id_child = child.split("=")[0].lstrip('C')  # if _child_ is C100001=F1~M1, so his id is 100001
            f_hap_idx = child.split("=")[1].split("~")[0].lstrip('F')  # if _child_ is C100001=F1~M1, so f_hap_idx is 1
            m_hap_idx = child.split("=")[1].split("~")[1].lstrip('M')  # if _child_ is C100001=F1~M1, so m_hap_idx is 1

            hap_1 = fam_dict['F'][int(f_hap_idx) - 1]  # the "-1" because hap1 is in index 0, and hap2 is in index 1
            hap_2 = fam_dict['M'][int(m_hap_idx) - 1]

            fam_dict[id_child] = [hap_1, hap_2]


def check_if_exist_pred_equal_to_sim_GRAMM(sim_file, gramm_pred, skip_sim_headers, genotype=False):
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

                if genotype:  # compare as genotypes
                    if pred_VS_sim_genotypes(hap1_sim, hap2_sim, hap1_pred, hap2_pred):
                        conclusion['right'] += 1
                        correct_pred_exist = True
                        break  # we just need to check if there is 1 correct option at least

                else:  # compare as haplotypes
                    if pred_haps_VS_sim_haps(hap1_sim, hap2_sim, hap1_pred, hap2_pred):
                        conclusion['right'] += 1
                        correct_pred_exist = True
                        break  # we just need to check if there is 1 correct option at least
            if not correct_pred_exist:
                conclusion['error'] += 1

    return conclusion


def check_if_first_pred_equal_to_sim_GRAMM(sim_file, gramm_pred, skip_sim_headers, genotype=False):
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

            if genotype:  # compare as genotypes
                if pred_VS_sim_genotypes(hap1_sim, hap2_sim, hap1_pred, hap2_pred):
                    conclusion['right'] += 1
                else:
                    conclusion['error'] += 1

            else:  # compare as haplotypes
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


# def pred_VS_sim_genotypes(hap1_sim, hap2_sim, hap1_pred, hap2_pred):
#     sim_genotype = []
#     for alleles1, alleles2 in zip(hap1_sim.split('~'), hap2_sim.split('~')):  # [A*01, B*02..] [A*03, B*04..]
#         sim_genotype.append([alleles1.split('*')[1], alleles2.split('*')[1]])  # [[01, 03], [02, 04] ..]
#     for idx_alleles, (alleles1, alleles2) in enumerate(zip(hap1_pred.split('~'), hap2_pred.split('~'))):
#         if alleles1.split('*')[1] not in sim_genotype[idx_alleles] or alleles2.split('*')[1] not in sim_genotype[idx_alleles]:
#             return False
#     return True


def pred_VS_sim_genotypes(hap1_sim, hap2_sim, hap1_pred, hap2_pred):

    # create genotype of sim data
    sim_genotype = []
    for alleles1, alleles2 in zip(hap1_sim.split('~'), hap2_sim.split('~')):  # [A*01, B*02..] [A*03, B*04..]
        sim_genotype.append([alleles1.split('*')[1], alleles2.split('*')[1]])  # [[01, 03], [02, 04] ..]

    # check only hap2_pred if hap1_pred is Unknown
    if hap1_pred == 'Unknown':
        for idx_alleles, alleles in enumerate(hap2_pred.split('~')):
            if alleles.split('*')[1] not in sim_genotype[idx_alleles]:
                return False

    # check only hap1_pred if hap2_pred is Unknown
    elif hap2_pred == 'Unknown':
        for idx_alleles, alleles in enumerate(hap1_pred.split('~')):
            if alleles.split('*')[1] not in sim_genotype[idx_alleles]:
                return False

    # check two haps preds
    else:
        for idx_alleles, (alleles1, alleles2) in enumerate(zip(hap1_pred.split('~'), hap2_pred.split('~'))):
            if alleles1.split('*')[1] not in sim_genotype[idx_alleles] or alleles2.split('*')[1] not in sim_genotype[idx_alleles]:
                return False

    return True


def check_if_first_pred_equal_to_sim_GRIMM(sim_file, grimm_pred, skip_sim_headers, genotype=False):
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

            if genotype:  # compare as genotypes
                if pred_VS_sim_genotypes(hap1_sim, hap2_sim, hap1_pred, hap2_pred):
                    conclusion['right'] += 1
                else:
                    conclusion['error'] += 1

            else:  # compare as haplotypes
                if (hap1_sim == hap1_pred and hap2_sim == hap2_pred) or (hap1_sim == hap2_pred and hap2_sim == hap1_pred):
                    conclusion['right'] += 1
                else:
                    conclusion['error'] += 1


    sim.close()
    pred.close()

    return conclusion


def check_if_exist_pred_equal_to_sim_GRIMM(sim_file, grimm_pred, skip_sim_headers, genotype=False):
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

                if genotype:  # compare as genotypes
                    if pred_VS_sim_genotypes(hap1_sim, hap2_sim, hap1_pred, hap2_pred):
                        conclusion['right'] += 1
                        correct_pred_exist = True
                        break  # we just need to check if there is 1 correct option at least

                else:  # compare as haplotypes
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


def run_gramm_code(sim_path, is_ser=False):
    reformat_file = f'running_for_GRAMM_paper/files_in_runtime/reformat_{os.path.basename(sim_path)}'  # check path
    alleles_names = ['A', 'B', 'C', 'DRB1', 'DQB1']
    files_address = 'running_for_GRAMM_paper/files_in_runtime'
    race_dict = {}

    open_ambiguity_sim = basic_csv2format(sim_path, reformat_file, race_dict)  # if update_preprocessing, add race_dict

    errors_in_families, aux_tools = run_GRAMM(reformat_file, files_address, alleles_names, race_dict, open_ambiguity_sim,
                                              is_serology=is_ser)

    grimm_path = "GR_code/GG_GRIMM"

    gl2GRIMM = files_address + '/glstring_for_GRIMM.txt'
    bin2GRIMM = files_address + '/binary_for_GRIMM.txt'

    move(gl2GRIMM, grimm_path + "/validation/simulation/data/input_test.txt")
    move(bin2GRIMM, grimm_path + "/validation/simulation/data/bin_input_test.txt")

    run_GRIMM(res_1000=False, sim=False)

    input2post = grimm_path + "/validation/output/don.hap.freqs"

    results, run_again = run_Post_GRIMM(input2post, reformat_file, alleles_names, files_address, aux_tools,
                                        errors_in_families)

    return results


def get_sim_file_with_persons_passed_successfully_gramm_running(orig_sim, out_from_run, new_for_grimm_run, skip_headers):
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
    cur_sim = "HARD"
    version = 'new_simulations'
    skip_headers = True

    SIMULATION = False
    NEW_SIMULATION = True
    Israel = False
    DD = False
    CORD_MOM = False
    part = ''

    Gen_or_Ser = 'serology'
    RUN_GRAMM = True
    RUN_GRIMM = True

    is_genotype = False

    if not os.path.exists(f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}'):
        os.mkdir(f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}')

    if not os.path.exists(f'running_for_GRAMM_paper/input_for_grimm_running/{version}'):
        os.mkdir(f'running_for_GRAMM_paper/input_for_grimm_running/{version}')

    if not os.path.exists(f'running_for_GRAMM_paper/output_from_grimm_after_update/{version}'):
        os.mkdir(f'running_for_GRAMM_paper/output_from_grimm_after_update/{version}')

    if SIMULATION:
        if RUN_GRAMM:

            _ = run_gramm_code(f'running_for_GRAMM_paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv')

            move('running_for_GRAMM_paper/files_in_runtime/results.csv',
                 f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}/results_{cur_sim}.csv')
            move('running_for_GRAMM_paper/files_in_runtime/errors.txt',
                 f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}/errors_{cur_sim}.txt')

            """
            >>> GRAMM - first and exist
            """
            right_and_error_gramm = check_if_first_pred_equal_to_sim_GRAMM(
                f'running_for_GRAMM_paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv',
                f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}/results_{cur_sim}.csv', skip_headers, genotype=is_genotype)

            right_and_error_gramm_exist = check_if_exist_pred_equal_to_sim_GRAMM(
                f'running_for_GRAMM_paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv',
                f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}/results_{cur_sim}.csv', skip_headers, genotype=is_genotype)

        if RUN_GRIMM:

            """ 
            >> GRIMM - first and exist
            """

            get_sim_file_with_persons_passed_successfully_gramm_running(
                f'running_for_GRAMM_paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv',
                f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}/results_{cur_sim}.csv',
                f'running_for_GRAMM_paper/input_for_grimm_running/{version}/{cur_sim}.csv', skip_headers)

            move(f'running_for_GRAMM_paper/input_for_grimm_running/{version}/{cur_sim}.csv',
                 'GR_code/GG_GRIMM/validation/simulation/data/input_test.txt')
            if os.path.exists("GR_code/GG_GRIMM/validation/simulation/data/bin_input_test.txt"):  # in grimm we dont use bin file
                os.remove("GR_code/GG_GRIMM/validation/simulation/data/bin_input_test.txt")

            run_GRIMM(res_1000=False, sim=False)

            right_and_error_grimm = check_if_first_pred_equal_to_sim_GRIMM(
                f'running_for_GRAMM_paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv',
                'GR_code/GG_GRIMM/validation/output/don.hap.freqs',
                skip_headers, genotype=is_genotype)

            right_and_error_grimm_exist = check_if_exist_pred_equal_to_sim_GRIMM(
                f'running_for_GRAMM_paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv',
                'GR_code/GG_GRIMM/validation/output/don.hap.freqs',
                skip_headers, genotype=is_genotype)

        # print results in the end of running
        if RUN_GRAMM:
            print('GRAMM: first result and exists result')
            print(right_and_error_gramm)
            print(right_and_error_gramm_exist)

        if RUN_GRIMM:
            print('GRIMM: first result and exists result')
            print(right_and_error_grimm)
            print(right_and_error_grimm_exist)

    if NEW_SIMULATION:
        if RUN_GRAMM:
            _ = run_gramm_code(f'running_for_GRAMM_paper/simulations updated/PEDIGREE_HAPLO_MASTER_CAU_{cur_sim}.csv')

            move('running_for_GRAMM_paper/files_in_runtime/results.csv',
                 f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}/results_{cur_sim}.csv')
            move('running_for_GRAMM_paper/files_in_runtime/errors.txt',
                 f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}/errors_{cur_sim}.txt')

            """
            >>> GRAMM - first and exist
            """
            # ************* haplotypes *************
            right_and_error_gramm_hap = check_if_first_pred_equal_to_sim_GRAMM(
                f'running_for_GRAMM_paper/simulations updated/PEDIGREE_HAPLO_MASTER_CAU_{cur_sim}.csv',
                f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}/results_{cur_sim}.csv',
                skip_headers, genotype=False)

            right_and_error_gramm_exist_hap = check_if_exist_pred_equal_to_sim_GRAMM(
                f'running_for_GRAMM_paper/simulations updated/PEDIGREE_HAPLO_MASTER_CAU_{cur_sim}.csv',
                f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}/results_{cur_sim}.csv',
                skip_headers, genotype=False)

            # ************* genotype *************
            right_and_error_gramm_geno = check_if_first_pred_equal_to_sim_GRAMM(
                f'running_for_GRAMM_paper/simulations updated/PEDIGREE_HAPLO_MASTER_CAU_{cur_sim}.csv',
                f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}/results_{cur_sim}.csv',
                skip_headers, genotype=True)

            right_and_error_gramm_exist_geno = check_if_exist_pred_equal_to_sim_GRAMM(
                f'running_for_GRAMM_paper/simulations updated/PEDIGREE_HAPLO_MASTER_CAU_{cur_sim}.csv',
                f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}/results_{cur_sim}.csv',
                skip_headers, genotype=True)

        if RUN_GRIMM:

            """ 
            >> GRIMM - first and exist
            """

            get_sim_file_with_persons_passed_successfully_gramm_running(
                f'running_for_GRAMM_paper/simulations updated/PEDIGREE_HAPLO_MASTER_CAU_{cur_sim}.csv',
                f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}/results_{cur_sim}.csv',
                f'running_for_GRAMM_paper/input_for_grimm_running/{version}/{cur_sim}.csv', skip_headers)

            move(f'running_for_GRAMM_paper/input_for_grimm_running/{version}/{cur_sim}.csv',
                 'GR_code/GG_GRIMM/validation/simulation/data/input_test.txt')
            if os.path.exists(
                    "GR_code/GG_GRIMM/validation/simulation/data/bin_input_test.txt"):  # in grimm run we dont use bin file
                os.remove("GR_code/GG_GRIMM/validation/simulation/data/bin_input_test.txt")

            run_GRIMM(res_1000=False, sim=False)

            # ************* haplotypes *************
            right_and_error_grimm_hap = check_if_first_pred_equal_to_sim_GRIMM(
                f'running_for_GRAMM_paper/simulations updated/PEDIGREE_HAPLO_MASTER_CAU_{cur_sim}.csv',
                'GR_code/GG_GRIMM/validation/output/don.hap.freqs',
                skip_headers, genotype=False)

            right_and_error_grimm_exist_hap = check_if_exist_pred_equal_to_sim_GRIMM(
                f'running_for_GRAMM_paper/simulations updated/PEDIGREE_HAPLO_MASTER_CAU_{cur_sim}.csv',
                'GR_code/GG_GRIMM/validation/output/don.hap.freqs', skip_headers, genotype=False)

            # ************* genotype *************
            right_and_error_grimm_geno = check_if_first_pred_equal_to_sim_GRIMM(
                f'running_for_GRAMM_paper/simulations updated/PEDIGREE_HAPLO_MASTER_CAU_{cur_sim}.csv',
                'GR_code/GG_GRIMM/validation/output/don.hap.freqs',
                skip_headers, genotype=True)

            right_and_error_grimm_exist_geno = check_if_exist_pred_equal_to_sim_GRIMM(
                f'running_for_GRAMM_paper/simulations updated/PEDIGREE_HAPLO_MASTER_CAU_{cur_sim}.csv',
                'GR_code/GG_GRIMM/validation/output/don.hap.freqs', skip_headers, genotype=True)

            # save output from gramm (for calculate pi^2, freqs..)
            move('GR_code/GG_GRIMM/validation/output/don.hap.freqs',
                 f'running_for_GRAMM_paper/output_from_grimm_after_update/{version}/{cur_sim}.txt')

        # print results in the end of running                                                                         s
        if RUN_GRAMM:
            print('GRAMM:\nHaplotypes:\nFirst result:')
            print(right_and_error_gramm_hap)
            print('Exists result')
            print(right_and_error_gramm_exist_hap)

            print('Genotypes:\nFirst result:')
            print(right_and_error_gramm_geno)
            print('Exists result')
            print(right_and_error_gramm_exist_geno)

        if RUN_GRIMM:
            print('GRIMM:\nHaplotypes:\nFirst result:')
            print(right_and_error_grimm_hap)
            print('Exists result')
            print(right_and_error_grimm_exist_hap)

            print('Genotypes:\nFirst result:')
            print(right_and_error_grimm_geno)
            print('Exists result')
            print(right_and_error_grimm_exist_geno)


    if Israel:
        if RUN_GRAMM:
            _ = run_gramm_code(f'running_for_GRAMM_paper/Australian_Israeli_data/processed_files_by_pyard/Israel_serology.csv', False)

            move('running_for_GRAMM_paper/files_in_runtime/results.csv',
                 f'running_for_GRAMM_paper/Australian_Israeli_data/results_and_errors/results_Israel.csv')
            move('running_for_GRAMM_paper/files_in_runtime/errors.txt',
                 f'running_for_GRAMM_paper/Australian_Israeli_data/results_and_errors/errors_Israel.txt')

    if DD:
        if RUN_GRAMM:

            if Gen_or_Ser == 'genetic':
                _ = run_gramm_code(
                    f'running_for_GRAMM_paper/Australian_Israeli_data/DD{part}_{Gen_or_Ser}_full.csv')
            elif Gen_or_Ser == 'serology':
                _ = run_gramm_code(
                    f'running_for_GRAMM_paper/Australian_Israeli_data/processed_files_by_pyard/DD_serology_full.csv', False)
            else:
                print('error: Gen_or_Ser must be "genetic" or "serology"')
                return


            move('running_for_GRAMM_paper/files_in_runtime/results.csv',
                 f'running_for_GRAMM_paper/Australian_Israeli_data/results_and_errors/results_DD{part}_{Gen_or_Ser}.csv')
            move('running_for_GRAMM_paper/files_in_runtime/errors.txt',
                 f'running_for_GRAMM_paper/Australian_Israeli_data/results_and_errors/errors_DD{part}_{Gen_or_Ser}.txt')

    if CORD_MOM:
        _ = run_gramm_code(f'running_for_GRAMM_paper/cordmom/cordmom_to_gramm_corrected.csv')

        move('running_for_GRAMM_paper/files_in_runtime/results.csv',
             f'running_for_GRAMM_paper/cordmom/results_and_errors/results_cord_mom.csv')
        move('running_for_GRAMM_paper/files_in_runtime/errors.txt',
             f'running_for_GRAMM_paper/cordmom/results_and_errors/errors_cord_mom.txt')


def main1():  # run grimm for save the haplotypes anf freqs it return, for calculate pi^2
    cur_sim = 'HARDEST'
    version = 'version_after_update_grimm'
    skip_headers = True

    get_sim_file_with_persons_passed_successfully_gramm_running(
        f'running_for_GRAMM_paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv',
        f'running_for_GRAMM_paper/results_and_errors_simulations_gramm/{version}/results_{cur_sim}.csv',
        f'running_for_GRAMM_paper/input_for_grimm_running/{version}/{cur_sim}.csv', skip_headers)

    move(f'running_for_GRAMM_paper/input_for_grimm_running/{version}/{cur_sim}.csv',
         'GR_code/GG_GRIMM/validation/simulation/data/input_test.txt')
    if os.path.exists(
            "GR_code/GG_GRIMM/validation/simulation/data/bin_input_test.txt"):  # in grimm we dont use bin file
        os.remove("GR_code/GG_GRIMM/validation/simulation/data/bin_input_test.txt")

    run_GRIMM(res_1000=False, sim=False)

    right_and_error_grimm = check_if_first_pred_equal_to_sim_GRIMM(
        f'running_for_GRAMM_paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv',
        'GR_code/GG_GRIMM/validation/output/don.hap.freqs',
        skip_headers)

    right_and_error_grimm_exist = check_if_exist_pred_equal_to_sim_GRIMM(
        f'running_for_GRAMM_paper/simulations/PEDIGREE_HAPLO_MASTER_{cur_sim}.csv',
        'GR_code/GG_GRIMM/validation/output/don.hap.freqs', skip_headers)

    print(right_and_error_grimm)
    print(right_and_error_grimm_exist)

    move('GR_code/GG_GRIMM/validation/output/don.hap.freqs', f'running_for_GRAMM_paper/output_from_grimm_after_update/{cur_sim}.txt')


main()
