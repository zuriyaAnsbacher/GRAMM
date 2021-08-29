import os
import json


def get_files_path(out_files_path):

    out_GLstr = os.path.join(out_files_path, 'glstring_for_GRIMM.txt')
    out_binary = os.path.join(out_files_path, 'binary_for_GRIMM.txt')

    return out_GLstr, out_binary


def load_jsons():
    # cur_path = os.path.abspath("GR_code/GG_GRAMM")
    data_path = "../data"

    with open(os.path.join(data_path, 'low2high.txt')) as json1:
        low2high = json.load(json1)

    with open(os.path.join(data_path, 'ser_dict_antigen2group.json')) as json2:
        antigen2group = json.load(json2)

    with open(os.path.join(data_path, 'ser_dict_group2antigen.json')) as json3:
        group2antigen = json.load(json3)

    return low2high, antigen2group, group2antigen


def remove_duplicate_children(family, alleles_names):
    """
    if there are two children with identical alleles data, remove one
    :param family: family dict
    :param alleles_names: alleles_names
    """

    def duplicate(child1, child2):
        for al_name in alleles_names:
            data1 = family[child1][al_name]
            data2 = family[child2][al_name]
            if data1 != data2:
                return False
        return True

    tested_children = []  # insert to this list the children we checked
    to_remove = []  # children to remove. (don't del them into the loop because it crash,
    # cause can't change dict while we run on it)
    for child in family:
        if child not in ['F', 'M']:
            if not tested_children:  # empty list
                tested_children.append(child)
            else:
                for other_child in tested_children:
                    if other_child != child:  # check it's not the same child
                        are_identical = duplicate(child, other_child)
                        if are_identical:
                            to_remove.append(child)
                            break
                        else:
                            tested_children.append(child)

    for dup_child in to_remove:
        del family[dup_child]


def insert_parents_data_to_their_haps(family, hapF, hapM, par_num):
    """
    when parents exist (al least one), we insert their data to their haplotypes
    (explanation about the insertion manner in documentation of 'insert_parents_data')
    :param family: family dict
    :param hapF: father haplotypes
    :param hapM: mother haplotypes
    :param par_num: how many parents (0/1/2)
    """
    if 'F' in family:
        hapF.insert_parents_data(family, 'F', par_num)
    if 'M' in family:
        hapM.insert_parents_data(family, 'M', par_num)


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



# def main():
#     fam = {
#         '1': {
#             'A': ['02'], 'B': ['03']
#         },
#         '2': {
#             'A': ['02'], 'B': ['04']
#         },
#         '3': {
#             'A': ['02'], 'B': ['03']
#         },
#         '4': {
#             'A': ['02'], 'B': ['03']
#         }
#     }
#     remove_duplicate_children(fam, ['A', 'B'])



