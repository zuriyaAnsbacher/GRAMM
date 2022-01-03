import csv
import os
import json
import itertools
from os import remove

from GR_code.GG_GRAMM.code.utils_oldVersion import sum_errors, convert_to_serology
from GR_code.GG_Post.code.readers.sim_reader_oldVersion import simReader


def is_equal(al_1, al_2):
    if al_1 == al_2 or al_1.startswith(al_2) or al_2.startswith(al_1) or al_1 == '' or al_2 == '':
        return True
    return False


# def validate(hap_1, hap_2, member, alleles_names):
#     hap_1 = string_list(hap_1, alleles_names)
#     hap_2 = string_list(hap_2, alleles_names)
#     for j in range(len(hap_1)):
#         if not (is_equal(member[j][0], hap_1[j]) and is_equal(member[j][1], hap_2[j])) \
#                 and not (is_equal(member[j][0], hap_2[j]) and is_equal(member[j][1], hap_1[j])):
#             return False
#     return True


def validate(hap_1, hap_2, member, alleles):
    hap_1 = string_list(hap_1, alleles)
    hap_2 = string_list(hap_2, alleles)
    # pairs = list(itertools.combinations(member, 2))
    # pairs_consistent = [True] * len(pairs)
    # for k, pair in enumerate(pairs):
    for j in range(len(hap_1)):
        pairs = list(itertools.combinations(member[j], 2))  # todo: is it for serology case?
        pairs_consistent = [True] * len(pairs)
        for k, pair in enumerate(pairs):
            if not (is_equal(pair[0], hap_1[j]) and is_equal(pair[1], hap_2[j])) \
                    and not (is_equal(pair[0], hap_2[j]) and is_equal(pair[1], hap_1[j])):
                pairs_consistent[k] = False
        pairs_consistent_opposite = [not element for element in pairs_consistent]
        if any(pairs_consistent_opposite):
            return False

    return True


        # is_consistency_per_allele = []
        # for als_member in member[j]:
        #     is_consistency_per_allele.append(not (is_equal(als_member, hap_1[j]) and is_equal(als_member, hap_2[j])))
        # if all(is_consistency_per_allele):
    #         return False
    # return True


def string_list(gl_str, alleles):
    als = gl_str.split('~')
    mydict = {}
    for al in als:
        key, val = al.split('*')
        mydict[key] = val
    return [mydict[x] for x in alleles]


def valid_family(key, family, p_1, p_2, alleles, double_als, writer):
    hap_1 = p_1[0]
    hap_2 = p_1[1]
    hap_3 = p_2[0]
    hap_4 = p_2[1]
    freq_1 = float(p_1[2])
    freq_2 = float(p_2[2])
    children = ''
    valid = []
    not_valid = []
    for j in range(len(family[0])):
        fmals = family[1][j]
        if family[0][j] == 'F' or family[0][j] == 'M':
            if validate(hap_1, hap_2, fmals, alleles) or validate(hap_3, hap_4, fmals, alleles):
                valid.append(family[0][j])
            else:
                not_valid.append(family[0][j])
        else:
            if validate(hap_1, hap_3, fmals, alleles):
                valid.append("C" + family[0][j] + '=F1~M1')
            elif validate(hap_1, hap_4, fmals, alleles):
                valid.append("C" + family[0][j] + '=F1~M2')
            elif validate(hap_2, hap_3, fmals, alleles):
                valid.append("C" + family[0][j] + '=F2~M1')
            elif validate(hap_2, hap_4, fmals, alleles):
                valid.append("C" + family[0][j] + '=F2~M2')
            else:
                not_valid.append(family[0][j])
    if len(not_valid) == 0:
        match_child2haps = [member for member in valid if '~' in member]  # example: [1=F2~M2, 2=F1~M2 ...]
        match_child2haps = ":".join(match_child2haps)

        # when the value of double_als[key] is 'F' or 'M', it means that 1 haplotype of this parent is unknown
        if key in double_als:
            if double_als[key] == 'F':
                hap_2 = 'unknown'
            elif double_als[key] == 'M':
                hap_4 = 'unknown'

        row = [key, hap_1, hap_2, hap_3, hap_4, match_child2haps, freq_1 * freq_2]
        writer.writerows([row])
        return True
    return False


# in 2 list: [01:02, 03:05] ; [01:03, 04:05] -- > return False because there is no common allele
# here: [01:02, 03:05] ; [01:02, 04:05] -- > return True because 01:02 is common
def common_allele(als_list1, als_list2):
    for al1 in als_list1:
        for al2 in als_list2:
            if is_equal(al1, al2):
                return True
    return False


def consistent_haps_par_and_geno_child(pair1, family, alleles):
    hap1 = string_list(pair1[0], alleles)
    hap2 = string_list(pair1[1], alleles)
    for ind_member, member in enumerate(family[0]):
        if member != 'F' and member != 'M':
            als_member = family[1][ind_member]
            for ind_als, als in enumerate(als_member):
                if not common_allele([hap1[ind_als], hap2[ind_als]], als):
                    return False
    return True


def open_serology_options_children(child, types, s_dict):
    group_of_antigen_to_remove = []
    for i, alleles in enumerate(child):
        for al in alleles:
            group = types[i] + '*' + al
            # check if the value (group) is in dict (low res only)
            # if it exists, open to all the serology options and remove the original value
            # (it's only group name, not an antigen)
            # example: A*10 appears, so add A*25, A*25, A*34, A*66 (only the digits) and remove A*10
            if ":" not in al and group in s_dict:
                for ser_value in s_dict[group]:
                    child[i].append(ser_value.split('*')[1])
                if group != "B*14" and group != "B*15" and group != "DRB1*03":
                    group_of_antigen_to_remove.append(group.split('*')[1])

        for gr_to_remove in group_of_antigen_to_remove:
            child[i].remove(gr_to_remove)
        group_of_antigen_to_remove = []

    # for chromosome in [c1f, c2f, c1m, c2m]:
    #     for type_ in alleles_names:
    #         for allele in chromosome[type_]:
    #             group = type_ + '*' + allele
    #             # check if the value (group) is in dict (low res only)
    #             # if it exists, open to all the serology options and remove the original value
    #             # (it's only group name, not an antigen)
    #             # example: A*10 appears, so add A*25, A*25, A*34, A*66 (only the digits) and remove A*10
    #             if ":" not in allele and group in s_dict:
    #                 for ser_value in s_dict[group]:
    #                     chromosome[type_].append(ser_value.split('*')[1])
    #                 chromosome[type_].remove(group.split('*')[1])


def run_Post_GRIMM(input_file, orig_input, alleles, output_path, er_lst, double_als, is_ser):
    min_freq = 1e-15  # change?

    wf = open(output_path + '/results.csv', 'w', newline='')
    writer = csv.writer(wf)
    row = ['FAMCODE', 'HAP_F1', 'HAP_F2', 'HAP_M1', 'HAP_M2', 'HAPS_INHERITANCE', 'PROBABILITY']
    writer.writerows([row])

    grimm_out = input_file
    grimm_out = open(grimm_out)
    grimm_out = grimm_out.readlines()
    data_dict = {}

    # create dict with the data from GG_GRIMM
    for i in range(len(grimm_out)):
        grimm_out[i] = grimm_out[i].strip().split(",")
        if grimm_out[i][0] in data_dict:
            data_dict[grimm_out[i][0]].append(grimm_out[i][1:])
        else:
            data_dict[grimm_out[i][0]] = [grimm_out[i][1:]]

    ddr = simReader(orig_input)

    fam_dict = {}
    family_ids, family_als, orig, _ = ddr.get_family()
    # create dict of original families (from input of user)
    while family_ids is not None:
        id = str(family_ids[0])
        s_id = id.split('.')
        fam_dict[s_id[0]] = [family_ids[1:], family_als]
        family_ids, family_als, orig, _ = ddr.get_family()

    if is_ser:
        gram_path = os.path.abspath("GR_code/GG_GRAMM")  # check that

        # with open(gram_path + '/data/ser_dict_group2antigen.json') as ser_dict_path_group2anti:
        #     ser_dict_group2anti = json.load(ser_dict_path_group2anti)
        # with open(gram_path + '/data/ser_dict_antigen2group.json') as ser_dict_path_anti2group:
        #     ser_dict_anti2group = json.load(ser_dict_path_anti2group)

        with open(
                "/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_GRAMM/data/ser_dict_group2antigen.json") as ser_dict_path_group2anti:
            ser_dict_group2anti = json.load(ser_dict_path_group2anti)
        with open(
                "/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_GRAMM/data/ser_dict_antigen2group.json") as ser_dict_path_anti2group:
            ser_dict_anti2group = json.load(ser_dict_path_anti2group)

        for key_ser, fam_val in fam_dict.items():
            for ind_c, child in enumerate(fam_val[1]):
                if fam_val[0][ind_c] != 'F' and fam_val[0][ind_c] != 'M':
                    for i, als in enumerate(child):
                        for j, al in enumerate(als):
                            fam_dict[key_ser][1][ind_c][i][j] = convert_to_serology(ser_dict_anti2group, alleles[i], al)
                    open_serology_options_children(child, alleles, ser_dict_group2anti)

    valid = 0
    invalid = 0
    valid_list = []
    invalid_list = []

    # compare between original data about families and the data from GG_GRIMM
    # by checking if children (from input) are consistent with haps of parents (from GG_GRIMM).
    for key in fam_dict:
        if key + '.0' in data_dict.keys() and key + '.1' in data_dict.keys():
            p1_lst = data_dict[key + '.0']
            p2_lst = data_dict[key + '.1']
            # convert freq to float, then sort the list (begin in the max freq)
            for p in p1_lst:
                p[2] = float(p[2])
            for p in p2_lst:
                p[2] = float(p[2])
            p1_lst = sorted(p1_lst, key=lambda x: (x[2]))
            p2_lst = sorted(p2_lst, key=lambda x: (x[2]))
            p1_lst = p1_lst[::-1]
            p2_lst = p2_lst[::-1]
            max_res2fam = 0
            incons_haps = 0
            for pair1 in p1_lst:
                if max_res2fam == 100:
                    break

                if not consistent_haps_par_and_geno_child(pair1, fam_dict[key], alleles):
                    incons_haps += 1
                    continue

                for pair2 in p2_lst:
                    if max_res2fam == 100:
                        break  # max results of family is 100

                    if not consistent_haps_par_and_geno_child(pair1, fam_dict[key], alleles):
                        incons_haps += 1
                        continue

                    # if fam key already appears in results and the freqs is too low, it won't be added to results:
                    if key not in valid_list or float(pair1[2]) * float(pair2[2]) > min_freq:
                        if valid_family(key, fam_dict[key], pair1, pair2, alleles, double_als, writer):
                            max_res2fam += 1
                            if key not in valid_list:
                                valid += 1
                                valid_list.append(key)
            if key not in valid_list:
                invalid += 1
                invalid_list.append(key)
            # print("num of incons haps in parents:" + str(incons_haps))

    # summarize all the errors (in pre processing and in post) to one file: "errors.txt"
    num_err = int(er_lst[0] + er_lst[1] + er_lst[2] + er_lst[3] / 2)
    all_err = int(num_err) + len(invalid_list)

    with open(output_path + '/errors.txt', 'a', newline='') as errors_d:
        for inv in invalid_list:
            errors_d.write(
                'family "' + str(inv) + '": ' + 'Inconsistency between the received parents\' haplotypes and the'
                                                ' children\'s HLA (The common reason: recombination)\n')

        sum_errors(errors_d, er_lst, len(fam_dict), len(invalid_list), all_err)

    wf.close()

    # run_again is a flag that give sign if there is 1 family and it failed_count because inconsistent.
    # in this case, we run again the process, with 1000 results in grimm (instead of 10).
    # it can lead to results (when the alleles_names are rare)
    run_again = True if len(fam_dict) == 1 and len(invalid_list) == 1 else False
    return output_path + '/results.csv', run_again

# run_Post_GRIMM("/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_Post/input_from_GRIMM/EASY_1.freqs", "/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_Post/original_input_from_user/EASY_1.csv", ['A', 'B', 'C', 'DRB1', 'DQB1'], "/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_Post/output")
