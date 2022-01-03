import os
import json

from GR_code.GG_GRAMM.code.als import Als


def is_valid(fam_d, al_types, par_num):
    """
    validation tests of alleles_names data
    return: 0 if valid. string with error description otherwise.
    """
    for types in al_types:
        for key in fam_d:
            for single_al in fam_d[key][types]:
                if single_al == "" or single_al == " ":
                    continue
                elif ":" in single_al:
                    parts = single_al.split(":")
                    for part in parts:
                        if not (part.isnumeric() or part.isupper() or part == ''):
                            return "There is an allele which contains an invalid character. " \
                                   "(It must contain numbers or uppercase only)."
                else:
                    if not (single_al.isnumeric() or single_al.isupper()):
                        return "There is an allele which contains an invalid character. " \
                               "(It must contain numbers or uppercase only)."
    if len(fam_d) < 2:
        return "Missing data about the family (less than two people)."
    if len(fam_d) == 2 and par_num == 2:
        return "Missing data about the family - there are no children."
    for types in al_types:
        lst = Als()
        for key in fam_d:  # F, M, 1, 2...
            if any(fam_d[key][types]):  # not empty
                # lst = fam_d[key][alleles_names].merge(lst)
                lst = lst.merge(fam_d[key][types])
        if len(lst) > 4:
            return "Too many alleles_names in the family."
    if par_num == 2:
        for types in al_types:
            fm_als = fam_d['F'][types] + fam_d['M'][types]
            for key in fam_d:
                if key != 'F' and key != 'M' and len(fm_als) == 4 and all(fm_als):
                    in_fm = fam_d[key][types].sub_lst(fm_als)
                    if not in_fm:
                        return "There is an allele in a _child_ that does not exist in the parents."

    return 0


def list_to_dict(family_ids, family_als, als_names):
    """
    convert family data from list to dict
    :param family_ids: alleles_names of family members (F, M, 1, 2...)
    :param family_als: alleles_names data of family
    :param als_names: alleles_names alleles_names (A, B, C, DR, DQ)
    :return: dict with data about the family
    """
    fam_dict = {}
    for i in range(len(family_ids)):
        fam_member = family_ids[i]
        fam_dict[fam_member] = {}
        for j in range(len(family_als[0])):
            als_data = Als()
            als_data.append(family_als[i][j][0])
            als_data.append(family_als[i][j][1])
            fam_dict[fam_member][als_names[j]] = als_data
    return fam_dict


def child_ls_merge(child_ls, child_d, types):
    """
    merge alleles_names of children to list (ex. [01:03, 08, 01, 08] -> [01:03, 08])
    @param child_ls: _child_ list
    @param child_d: _child_ dictionary
    @param types: A/B/C/DR/DB
    @return: merged list
    """
    lst = Als()
    for child in child_ls:
        if any(child_d[child][types]):
            # lst = child_d[_child_][alleles_names].merge(lst)
            lst = lst.merge(child_d[child][types])
    return lst


# convert parents chromosomes to list (for writing in file)
def options_num_fm(father, mother):
    op_lst = []
    i = 0
    for chro in [father.ch1, father.ch2, mother.ch1, mother.ch2]:
        cur_op = 1
        for value in chro.values():
            if len(value) > 1:
                cur_op = cur_op * len(value)
        op_lst.append(cur_op)
        i += 1
    return op_lst


# write the current error to file
def errors_report(f_er, count, ind, error_name, fam_id):
    count[ind] += 1
    f_er.write('family "' + str(fam_id) + '": ' + error_name + '\n')  # write error to errors file


def sum_errors(f, er_count, fam_count, inconcis_er, all_err):
    f.write('\n\n***** Errors Summary *****' + '\n')
    f.write('Errors in pre-processing. Invalid families: ' + str(er_count[0]) +
            ' families. ' + "{:.1f}".format(er_count[0] / fam_count * 100) + '%' + ' from data\n')
    f.write('Errors in matching. Can not match the children\'s data to the parents\' data: ' + str(er_count[1]) +
            ' families. ' + "{:.1f}".format(er_count[1] / fam_count * 100) + '%' + ' from data\n')
    f.write('Errors in adding data. Contradiction between the children\'s data and the parents\' data: ' + str(er_count[2]) +
            ' families. ' + "{:.1f}".format(er_count[2] / fam_count * 100) + '%' + ' from data\n')
    f.write('Errors in creating GL string: ' + str(er_count[3]) +
            ' parents. ' + "{:.1f}".format(er_count[3] / fam_count * 50) + '%' + ' from data\n')
    f.write('Errors in the results. Inconsistency between the children\'s haplotypes and the parents\' haplotypes: ' + str(inconcis_er) + ' families. ' +
            "{:.1f}".format(inconcis_er / int(fam_count) * 100) + '% from data\n')
    f.write('\nTotal errors:\n' + str(all_err) + ' from ' + str(fam_count) + ' families\n')


# ser_dict = {"A*23": "A*09", "A*24": "A*09"..} so if we get allele = A, single_al = 23, change to 09
# that only if is_serology (in "extract_als") is True
def convert_to_serology(ser_dict, allele_type, single_al):
    # case of low res (A*23 -> A*09)
    if allele_type + '*' + single_al in ser_dict:
        return ser_dict[allele_type + '*' + single_al].split('*')[1]
    # case of high res (A*23:01 -> A*09)
    if ':' in single_al and allele_type + '*' + single_al.split(':')[0] in ser_dict:
        return ser_dict[allele_type + '*' + single_al.split(':')[0]].split('*')[1]
    return single_al

