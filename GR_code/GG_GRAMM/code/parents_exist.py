from aux_functions_update import empty_hap


def associate_children_to_par_haps_while_parents_exist(hapF, hapM, child, alleles_names):
    """
    associate child to the parents haplotypes he inherited. this function is for scenario that parents exist in data.
    for each parent, we go over the alleles and try to associate. the different cases explained in the function.
    :param hapF: father haplotype
    :param hapM: mother haplotype
    :param child: child dict
    :param alleles_names: alleles names
    :return: True if the associating succeed completely (two parents), False otherwise ; and 'associating' list
    """
    associating = [None, None]  # first index: father, second: mother

    for idx_parent, parent_haps in enumerate([hapF, hapM]):
        if empty_hap(parent_haps.hap1) and empty_hap(parent_haps.hap2):  # 2 haplotypes are empty
            associating[idx_parent] = 1  # associate (randomly) in first haplotype
            continue

        for allele_name in alleles_names:
            current_alleles_hap1, current_alleles_hap2 = parent_haps.hap1[allele_name], parent_haps.hap2[allele_name]
            current_allele_child = child[allele_name]

            if associating[idx_parent] is not None:
                break  # associating of current parent succeed

            elif not any(current_allele_child):  # child has no data in this allele (["", ""])
                continue  # to next allele

            # 1 haplotype empty, 1 with data (in parent)
            # (it enough to check that at least 1 hap is empty, because case of 2 empty was treated before)
            elif empty_hap(parent_haps.hap1) or empty_hap(parent_haps.hap2):
                full_hap, full_hap_idx = (parent_haps.hap1, 1) if not empty_hap(parent_haps.hap1) else (parent_haps.hap2, 2)
                # in this case, we can associate only if full haplotype contradicts child data
                # (so associate with the empty).
                if not current_allele_child.sub_lst(full_hap[allele_name]):
                    # if 'full_hap_idx' = 1, the empty hap is 2, and if 'full_hap_idx' = 2, the empty hap is 1
                    associating[idx_parent] = 3 - full_hap_idx

            # two haplotypes of parent are full
            # parent has 1 value in each haplotype, and no homozygous in this haplotype. (e.g. hap1: A:[01], hap2: [02])
            elif len(current_alleles_hap1) == len(current_alleles_hap2) == 1 and \
                    current_alleles_hap1 != current_alleles_hap2:
                # check intersection between alleles of child and parent
                alleles_2_haps_parent = current_alleles_hap1 + current_alleles_hap2
                inter, idx_inter = current_allele_child.intersection(alleles_2_haps_parent)

                if inter == 1:  # single match between child and specific parent haplotype
                    associating[idx_parent] = idx_inter + 1

                # when inter = 2, it means that parents have common allele, like: f:[01, 02], m:[02, 03] child:[01, 02]
                # so f has 2 intersections with child. (could be also if parent homomz', but condition before deny it)
                elif inter == 2:
                    # we check inter' with *other* parent. if there is 1 inter', we find the common allele
                    # (02, in the example above), then we associate child with the hap that doesnt contain the common
                    # because the common allele is the one that connect the child to the *other* parent
                    # (like case '4b' in function 'associate_children_to_par_haps_while_parents_dont_exist')
                    other_parent = hapM if idx_parent == 0 else hapF
                    other_parent_alleles = other_parent.hap1[allele_name] + other_parent.hap2[allele_name]
                    inter_other, idx_inter_other = current_allele_child.intersection(other_parent_alleles)
                    if inter_other == 1:
                        common_allele_parents = other_parent_alleles[idx_inter_other]
                        if common_allele_parents in current_alleles_hap1:
                            associating[idx_parent] = 2
                        elif common_allele_parents in current_alleles_hap2:
                            associating[idx_parent] = 1

                # any other value of inter (0 or >2) doesn't give us information for the associating

    # because we go over each parent separately, could be partial associating (only for father or mother).
    # success is defined as associating with 2 parents
    complete_success = True if all(associating) else False

    return complete_success, associating
