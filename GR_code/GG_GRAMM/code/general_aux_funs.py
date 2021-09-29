from GR_code.GG_GRAMM.code.als_update import Als


def equal_als(al_1, al_2):
    """
    :param al_1: first allele
    :param al_2: second allele
    :return: true if the allele are identical or including, false otherwise
    """
    if al_1 == al_2 or al_1 is al_2 or al_1.startswith(al_2) or al_2.startswith(al_1):
        if bool(al_1) == bool(al_2):
            return True
    return False


def empty_hap(hap):
    """
    check if haplotype is empty (no data, just: {'A':[], 'B':[], ..}
    @param hap: dict represents haplotype
    @return: True if empty, else False
    """
    for value in hap.values():
        if value:
            return False
    return True


def convert_to_serology(ser_dict, allele_name, single_al):
    """
    when data is serology, the alleles need to be converted, for example:
    ser_dict = {"A*23": "A*09", "A*24": "A*09"..} so if we get allele = A, single_al = 23, change to 09
    :param ser_dict: serology dict
    :param allele_name: allele name
    :param single_al: single allele
    :return:
    """
    # case of low res (A*23 -> A*09)
    if allele_name + '*' + single_al in ser_dict:
        return ser_dict[allele_name + '*' + single_al].split('*')[1]
    # case of high res (A*23:01 -> A*09)
    if ':' in single_al and allele_name + '*' + single_al.split(':')[0] in ser_dict:
        return ser_dict[allele_name + '*' + single_al.split(':')[0]].split('*')[1]
    return single_al

# # unused functions (?)
# def child_ls_merge(child_ls, child_d, types):
#     """
#     merge alleles of children to list (ex. [01:03, 08, 01, 08] -> [01:03, 08])
#     @param child_ls: child list
#     @param child_d: child dictionary
#     @param types: A/B/C/DR/DB
#     @return: merged list
#     """
#     lst = Als()
#     for child in child_ls:
#         if any(child_d[child][types]):
#             # lst = child_d[child][types].merge(lst)
#             lst = lst.merge(child_d[child][types])
#     return lst
#
#
# # convert parents chromosomes to list (for writing in file)
# def options_num_fm(father, mother):
#     op_lst = []
#     i = 0
#     for chro in [father.ch1, father.ch2, mother.ch1, mother.ch2]:
#         cur_op = 1
#         for value in chro.values():
#             if len(value) > 1:
#                 cur_op = cur_op * len(value)
#         op_lst.append(cur_op)
#         i += 1
#     return op_lst
#
#
# # write the current error to file
# def errors_report(f_er, count, ind, error_name, fam_id):
#     count[ind] += 1
#     f_er.write('family "' + str(fam_id) + '": ' + error_name + '\n')  # write error to errors file
#
#
# def sum_errors(f, er_count, fam_count, inconcis_er, all_err):
#     f.write('\n\n***** Errors Summary *****' + '\n')
#     f.write('Errors in pre-processing. Invalid families: ' + str(er_count[0]) +
#             ' families. ' + "{:.1f}".format(er_count[0] / fam_count * 100) + '%' + ' from data\n')
#     f.write('Errors in matching. Can not match the children\'s data to the parents\' data: ' + str(er_count[1]) +
#             ' families. ' + "{:.1f}".format(er_count[1] / fam_count * 100) + '%' + ' from data\n')
#     f.write('Errors in adding data. Contradiction between the children\'s data and the parents\' data: ' + str(er_count[2]) +
#             ' families. ' + "{:.1f}".format(er_count[2] / fam_count * 100) + '%' + ' from data\n')
#     f.write('Errors in creating GL string: ' + str(er_count[3]) +
#             ' parents. ' + "{:.1f}".format(er_count[3] / fam_count * 50) + '%' + ' from data\n')
#     f.write('Errors in the results. Inconsistency between the children\'s haplotypes and the parents\' haplotypes: ' + str(inconcis_er) + ' families. ' +
#             "{:.1f}".format(inconcis_er / int(fam_count) * 100) + '% from data\n')
#     f.write('\nTotal errors:\n' + str(all_err) + ' from ' + str(fam_count) + ' families\n')

