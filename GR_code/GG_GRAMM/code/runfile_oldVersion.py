import os
import copy
import json
from GR_code.GG_GRAMM.code.readers.sim_reader_oldVersion import simReader
from GR_code.GG_GRAMM.code.utils_oldVersion import list_to_dict, is_valid, options_num_fm, errors_report, sum_errors
from GR_code.GG_GRAMM.code.DualChromosome_oldVersion import DualChromosome
from GR_code.GG_GRAMM.code.emb_with_par_oldVersion import emb_wp
from GR_code.GG_GRAMM.code.emb_no_par_oldVersion import emb_np, emb_4als_in_pars
from GR_code.GG_GRAMM.code.add_data_by_children_oldVersion import add_data_child
from GR_code.GG_GRAMM.code.chroms_to_GL_oldVersion import create_gl


# old version

def run_GRAMM(input_path, alleles, output_path, is_serology, race_dict, sim=False):
    cur_path = os.path.abspath("GR_code/GG_GRAMM")  # check that

    if sim:
        cur_path = '../GR_code/GG_GRAMM'  # remove. only for simulations

    gl_file = open(output_path + "/input_test.txt", 'w+')
    bin_path = output_path + "/bin_input_test.txt"
    errors_d = open(output_path + '/errors.txt', 'w')

    # 'double_als' is a dict that map between *families* whose F or M have 1 empty haplotype
    # and the *parent* with the empty haplotype. (happen in case of insertion 1 parent and 1 _child_)
    # i.g: {"2": "empty1f", "15": "empty1m"}.
    double_als = {}

    with open(cur_path + '/data/low2high.txt') as json_file:
        d_low2high = json.load(json_file)

    with open(cur_path + '/data/ser_dict_antigen2group.json') as ser_dict_path_anti2group:
        ser_dict_anti2group = json.load(ser_dict_path_anti2group)

    with open(cur_path + '/data/ser_dict_group2antigen.json') as ser_dict_path_group2anti:
        ser_dict_group2anti = json.load(ser_dict_path_group2anti)

    # with open("/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_GRAMM/data/low2high.txt") as json_file:
    #     d_low2high = json.load(json_file)
    # with open("/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_GRAMM/data/ser_dict_antigen2group.json") as ser_dict_path_anti2group:
    #     ser_dict_anti2group = json.load(ser_dict_path_anti2group)
    # with open("/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_GRAMM/data/ser_dict_group2antigen.json") as ser_dict_path_group2anti:
    #     ser_dict_group2anti = json.load(ser_dict_path_group2anti)

    error1 = ''  # ?
    error2 = 'Less than 4 different alleles_names. Algorithm can not be executed'
    error3 = 'Failed in adding data. Contradiction between parents and _child_'
    error4 = 'Failed in creating GL string'
    errors_count = [0, 0, 0, 0]

    # reading data from xlsx file in format described for DDReader class
    ddr = simReader(input_path)

    # get the data of family from data file
    family_ids, family_als, par_num, amb = ddr.get_family(is_serology, ser_dict_anti2group)
    fam_count = 0
    bin_d = {}

    while family_ids:
        al_types = copy.deepcopy(alleles)
        fam_count += 1
        if fam_count > int(family_ids[0]):
            break
        # convert data from list to dict
        fam_dict = list_to_dict(family_ids[1:], family_als, al_types)
        # validation test
        validation = is_valid(fam_dict, al_types, par_num)
        if validation != 0:
            errors_report(errors_d, errors_count, 0, validation, family_ids[0])
            family_ids, family_als, par_num, amb = ddr.get_family(is_serology, ser_dict_anti2group)
            continue
        else:
            chF = DualChromosome(al_types)  # create dual chromosome of father (without data, yet)
            chM = DualChromosome(al_types)  # create dual chromosome of mother (without data, yet)
            if par_num > 0:  # one or two parents. embed data in parents chromosomes
                if 'F' in family_ids:
                    chF.emb_data(fam_dict, 'F', par_num)
                if 'M' in family_ids:
                    chM.emb_data(fam_dict, 'M', par_num)
            # create children list
            child_lst = []
            for member in family_ids[1:]:
                if member != 'F' and member != 'M':
                    child_lst.append(member)
            children_d = {}
            for child in child_lst:
                children_d[child] = copy.deepcopy(fam_dict[child])  # like regular d, without pars
            children_d_c = copy.deepcopy(children_d)
            if par_num == 0:
                # emb 4 als in pars chromosomes
                type_emb, is_deter = emb_4als_in_pars(chF, chM, child_lst, children_d, al_types)
                if not type_emb:
                    errors_report(errors_d, errors_count, 1, error2, family_ids[0])
                    family_ids, family_als, par_num, amb = ddr.get_family(is_serology, ser_dict_anti2group)
                    continue
            adding_child_failed = False
            for child in child_lst:
                if par_num > 0:
                    child_d = copy.deepcopy(fam_dict[child])  # d for specific _child_: {'A':[..], 'B':[..], ...}  # TODO: should be deepcopy?
                    emb_FM = emb_wp(chF, chM, child_d, al_types)
                elif type_emb:  # no parents
                    # als_of_ty = children_d[_child_][type_emb].copy()
                    als_of_ty = copy.deepcopy(children_d[child][type_emb])
                    f1, f2, m1, m2 = chF.ch1, chF.ch2, chM.ch1, chM.ch2
                    emb_FM = emb_np(type_emb, als_of_ty, is_deter, f1, f2, m1, m2)
                if not emb_FM:
                    break
                father, mother = add_data_child(child, chF, chM, emb_FM, al_types, children_d_c)
                if not father:  # failed_count in adding data -> father = None
                    errors_report(errors_d, errors_count, 2, error3, family_ids[0])
                    family_ids, family_als, par_num, amb = ddr.get_family(is_serology, ser_dict_anti2group)
                    adding_child_failed = True
                    break
        if adding_child_failed:
            continue
        if validation == 0 and father:
            options = options_num_fm(father, mother)

            temp, empty1f, empty1m = create_gl(father, mother, al_types, family_ids[0], amb, gl_file, bin_d, d_low2high,
                                               race_dict, is_serology, ser_dict_group2anti)
            if temp == -1:  # failed_count in creating gl string to father
                errors_report(errors_d, errors_count, 3, error4, family_ids[0])
                family_ids, family_als, par_num, amb = ddr.get_family(is_serology, ser_dict_anti2group)
                continue
            elif temp == -2:  # failed_count in creating gl string to mother
                errors_report(errors_d, errors_count, 3, error4, str(family_ids[0]).replace('.0', '.1'))
                family_ids, family_als, par_num, amb = ddr.get_family(is_serology, ser_dict_anti2group)
                continue
            else:
                if empty1f:
                    double_als[family_ids[0]] = "F"
                elif empty1m:
                    double_als[family_ids[0]] = "M"
                al_types = copy.deepcopy(alleles)
        family_ids, family_als, par_num, amb = ddr.get_family(is_serology, ser_dict_anti2group)

    with open(bin_path, 'w') as bin_gl_f:
        bin_j = json.dump(bin_d, bin_gl_f)

    gl_file.close()
    errors_d.close()

    return output_path + "/input_test.txt", bin_path, output_path + "/errors.txt", errors_count, double_als


# from processing_data2input.processing import basic_csv2format
# basic_csv2format("/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_GRAMM/input/EASY_SIM.csv", "/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_GRAMM/input/EASY_SIM_old_format.csv")
# run_GRAMM("/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_GRAMM/input/EASY_SIM_old_format.csv", ['A', 'B', 'C', 'DRB1', 'DQB1'],
#           "/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_GRAMM/output", False, {})
