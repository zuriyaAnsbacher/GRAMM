from GR_code.GG_GRAMM.code.read_families_from_file import get_families
from GR_code.GG_GRAMM.code.aux_functions_update import get_files_path, load_jsons, remove_duplicate_children, insert_parents_data_to_their_haps
from GR_code.GG_GRAMM.code.validation import is_valid_family
from GR_code.GG_GRAMM.code.utils_update_part2 import check_num_parents, convert_data_to_Als
from GR_code.GG_GRAMM.code.dualHaplotype import DualHaplotype
from GR_code.GG_GRAMM.code.no_parents import insert_data_to_single_allele_in_par_haps, \
    associate_children_to_par_haps_while_parents_dont_exist


def run_GRAMM(input_path, out_files_path, alleles_names, races_dict, is_serology):
    out_GLstr, out_binary = get_files_path(out_files_path)
    low2high, antigen2group, group2antigen = load_jsons()

    aux_tools = {'problematic_child': None, 'is_serology': is_serology, 'antigen2group': antigen2group, 'amb': {}}
    one_parent_empty = {}
    errors_in_families = {}

    families_dict = get_families(input_path, aux_tools)

    # sorted keys of families_dict, for process the families in the order (dict is not ordered)
    sorted_keys = sorted(families_dict.keys(), key=lambda my_key: int(my_key))  # check that it sorted

    # TODO: check if 'alleles_names' is changed in the loop and need to create a copy
    for count_fam in range(len(families_dict)):
        family = families_dict[sorted_keys[count_fam]]  # get one family
        par_num = check_num_parents(family)
        children_keys = [k for k in family.keys() if k not in ['F', 'M']]
        associating_children = {}

        # convert the data structure of each pair alleles from list to Als
        # 'Als' is a class which inherited from 'list', with adjusted attributes for alleles
        convert_data_to_Als(family)

        remove_duplicate_children(family, alleles_names)

        valid_family = is_valid_family(family, alleles_names, par_num, count_fam, aux_tools, errors_in_families)  # TODO: add checking format? (mistake in writing)
        if not valid_family:
            continue  # we do not want to analyze this family, so continue to the next family
        # TODO: in post-grimm, need to add this child to result file and visualization file

        # create dual haplotype for father and mother (without data, yet)
        hapF, hapM = DualHaplotype(alleles_names), DualHaplotype(alleles_names)

        if par_num > 0:
            insert_parents_data_to_their_haps(family, hapF, hapM, par_num)
        else:  # num_par == 0
            success_insert_single_allele, allele_inserted_to_parents = \
                insert_data_to_single_allele_in_par_haps(family, alleles_names, hapF, hapM, count_fam, errors_in_families)
            if not success_insert_single_allele:
                continue

        planA, planB = [], []  # TODO
        for child in children_keys:
            if par_num > 0:
                associate_children_to_par_haps_while_parents_exist()  # should separate to cases of 1 par and 2 pars?
            else:  # num_par == 0
                success_associate, associating = \
                    associate_children_to_par_haps_while_parents_dont_exist(hapF, hapM, family[child], allele_inserted_to_parents)
                if success_associate:
                    planA.append(child)
                else:
                    planB.append(child)


        # TODO:
        #  1.   if parents exist:
        #           a. insert their data to their haplotypes.
        #           b. associate each child with the 2 haplotypes he inherited
        #       else:
        #           a. add data to one allele from the children (from validation test must be at least one allele
        #           with 4 different values, if there is not parents)
        #           b. associate each child with the 2 haplotypes he inherited
        #  2.   add data from children to parents haplotypes
        #  3.   create gl strings from parents haplotypes and write to file for GRIMM








# def main():
#     a = {'1': {'F': ['3'], 'M': ['2'], '1': ['5']}, '2': {'F': ['4'], 'M': ['2']}}
#     b = a['1']
#     del b['F']
#     del a['1']['M']
#     print('')
#
# main()