import copy
from general_aux_funs import equal_als


def get_inherited_haps(hapF, hapM, associating):
    """
    interpret the list 'associating' to pair of the haplotypes the child inherited, and pair of not inherited.
    for example, if associating = [2, 1], so inherited haps are hapF.hap2, hapM.hap1. (and not inherited are hapF.hap1, hapM.hap2)
    if the associating is unknown, the value in 'associating' is 'None', so if associating = [None, 2] so inherited are
    None and hapM.hap2
    :param hapF: father haplotype
    :param hapM: mother haplotype
    :param associating: associating list
    :return: hapF_inherited, hapF_no_inherited, hapM_inherited, hapM_no_inherited
    """
    inheritance = [None, None, None, None]  # order:hapF_inherited, hapF_no_inherited, hapM_inherited, hapM_no_inherited
    for idx_par, hap, associate in zip([0, 2], [hapF, hapM], associating):
        if associate == 1:  # hap1 is inherited
            inheritance[idx_par] = hap.hap1
            inheritance[idx_par + 1] = hap.hap2
        elif associate == 2:  # hap2 is inherited
            inheritance[idx_par] = hap.hap2
            inheritance[idx_par + 1] = hap.hap1
    return inheritance[0], inheritance[1], inheritance[2], inheritance[3]


def get_high_res_allele(str1, str2):
    """
    return allele with high res. if str1 = '01', str2='01:02', return str2
    :param str1: string of first allele
    :param str2: atring of second allele
    :return: allele with high resolution
    """
    if len(str1) >= len(str2):
        return str1
    return str2


def compare_with_connected_hap_and_remove_if_needed(p_als, no_inher_hap, inher_hap_second_par, idx_common, allele_name):
    """
    there are 2 cases that because of them, this function is needed:
    1. when parents exist, and we insert their data to their haplotypes. except of the first allele, all the insertion
       are non deterministic (2 values in each allele).
       for example: hap1: A:[01], B:[03, 04] ; hap2: A:[02], B:[03, 04]. so if we find, in this stage, by comparison to
       children, that in hap1, B is 03, so automatically we want that this value will be removed from hap2.
    2. when parents don't exist, and we try to add data from children to them, can be a scenario that child inherited
       hap1 from father and hap1 from mother, and his values in A are [01, 02], so we insert these values to hap1 of
       father and hap1 of mother (in allele A). and again, if now we determine that this father has only 01 in hap1 (allele A)
       (by comparison to another child, for example), we want to remove '01' from hap2 of mother.
    we call this function if we want to determine 1 value from the 2, so we check to cases above:
    if the allele in the inherited haplotype (that we want to determine 1 value) equal to the allele in the haplotype of
    the current parent (case 1) or to the haplotype that the child inherited from the second parent (case 2) - we remove
    from them the value that we want to determine in the inherited haplotype (p_als).

    :param p_als: the allele from current parent, that child inherited
    :param no_inher_hap: the allele from the current parent, that the child didn't inherit
    :param inher_hap_second_par: the allele that child inherited from second parent
    :param idx_common: idx of the common value between parent and child (that we want determine in p_als and remove
    from the connected haplotype
    :param allele_name: allele name
    """
    no_inherited_hap_in_par = no_inher_hap[allele_name] if no_inher_hap else None
    inherited_second_par = inher_hap_second_par[allele_name] if inher_hap_second_par else None

    # first condition (if no_inherited.., elif inherited_second..) for check it's not None
    # second (p_als == ..) for check that they have same values (were inserted together, probably)

    if no_inherited_hap_in_par and p_als == no_inherited_hap_in_par:
        no_inherited_hap_in_par.remove_a(p_als[idx_common])
    elif inherited_second_par and p_als == inherited_second_par and p_als:
        inherited_second_par.remove_a(p_als[idx_common])


def add_child_data(hapF, hapM, child, idx_child, associating, alleles_names, aux_tools, idx_fam, errors_in_families):
    """
    after we associated a child with the haplotypes he inherited from parents, we go over on the alleles of child and
    parents in these haplotypes and try to add data to the parents, based on the comparison between them.
    we handle the options :
        empty child (no data in allele): do nothing
        empty parent (no data in allele): add child data to parent
        child and parent both have 1 value in allele: if contradictory - error, else: insert the value with high resolution
        child has 2 parent has 1: if contradictory - error, else: insert the value with high res, and remove from child
        the common allele
        child has 1 parent has 2: if contradictory - error, else: insert the value with high res, and remove the common
        from parent (+ maybe remove from another haplotype, explained in 'compare_with..._remove_if_needed')
        child has 2 parent has 2: like above, with remove the common allele from child

    about errors: if only one child is problematic - continue, if more than one - stop the running.
    :param hapF: father haplotypes
    :param hapM: mother haplotypes
    :param child: child dict
    :param idx_child: child idx (1, 2..)
    :param associating: associating list
    :param alleles_names: alleles names
    :param aux_tools: dict with auxiliary tools
    :param idx_fam: family count (1, 2..)
    :param errors_in_families: dict contains errors in families
    :return: if adding was succeed or not
    """
    _child_ = copy.deepcopy(child)  # copy because we change '_child_'
    hapF_inherited, hapF_no_inherited, hapM_inherited, hapM_no_inherited = get_inherited_haps(hapF, hapM, associating)
    inherited_haps, no_inherited_haps = [hapF_inherited, hapM_inherited], [hapF_no_inherited, hapM_no_inherited]
    continue_running = True

    def child_error():
        if aux_tools['problematic_child']:  # more than one problematic child exist (the value in the dict != None)
            return False  # stop the running
        else:
            aux_tools['problematic_child'] = idx_child  # this is the first problematic child, so we allow it
            return True  # continue the running

    # we run also on the haplotypes from the second parent for calling to function 'compare_with..remove_if_needed'
    for inher_hap, no_inher_hap, inher_hap_second_par in zip(inherited_haps, no_inherited_haps, inherited_haps[::-1]):

        if not inher_hap:  # hap = None, we have no data about the inheritance of this parent
            continue
        for allele_name in alleles_names:
            c_als = _child_[allele_name]
            p_als = inher_hap[allele_name]

            empty_child = True if c_als.empty_Als() else False  # child: ['', '']
            empty_parent = True if p_als.empty_Als() else False  # parent: ['', '']
            child_has_1_parent_has_1 = True if len(c_als) == len(p_als) == 1 else False  # c:[01] p:[01:05]
            child_has_2_parent_has_1 = True if (len(c_als) == 2 and len(p_als) == 1) else False   # c:[01:05, 02] p:[01]
            child_has_1_parent_has_2 = True if (len(c_als) == 1 and len(p_als) == 2) else False  # c:[01] p:[01:05, 02]
            child_has_2_parent_has_2 = True if len(c_als) == len(p_als) == 2 else False  # c:[01:05, 02] p:[01, 03]

            if empty_child:
                continue  # to next allele

            elif empty_parent:  # tested
                # if child homozygous, remove one value and insert only one to the parent
                if len(c_als) == 2 and c_als[0] == c_als[1]:
                    c_als.remove_a(c_als[1])
                p_als.clear()  # p_als: ['', ''] - > []. for not contain empty value when we add new values
                p_als.extend(c_als)  # insert child data to parent

            elif child_has_1_parent_has_1:
                if c_als != p_als:  # the values are contradictory, so it's an error
                    continue_running = child_error()
                else:
                    p_als[0] = get_high_res_allele(p_als[0], c_als[0])  # insert to parent the value with high res

            elif child_has_2_parent_has_1:
                if p_als[0] not in c_als:  # the values are contradictory, so it's an error
                    continue_running = child_error()
                else:
                    idx_common = c_als.index_a(p_als[0])
                    p_als[0] = get_high_res_allele(p_als[0], c_als[idx_common])
                    c_als.remove_a(c_als[idx_common]) # remove common allele from child- because it's inserted to parent

            elif child_has_1_parent_has_2:  # tested
                if c_als[0] not in p_als:
                    continue_running = child_error()  # the values are contradictory, so it's an error
                else:
                    idx_common = p_als.index_a(c_als[0])
                    # insert the value with the high res
                    p_als[idx_common] = get_high_res_allele(p_als[idx_common], c_als[0])
                    compare_with_connected_hap_and_remove_if_needed(p_als, no_inher_hap, inher_hap_second_par,
                                                                    idx_common, allele_name)
                    p_als.remove_a(p_als[1 - idx_common])  # remove the non-common allele between parent and child

            elif child_has_2_parent_has_2:  # tested
                par_inter, idx_par_inter = c_als.intersection(p_als)
                if par_inter == 1:  # one intersection between values of child and parent
                    idx_common_child = c_als.index_a(p_als[idx_par_inter])
                    p_als[idx_par_inter] = get_high_res_allele(p_als[idx_par_inter], c_als[idx_common_child])

                    compare_with_connected_hap_and_remove_if_needed(p_als, no_inher_hap, inher_hap_second_par,
                                                                    idx_par_inter, allele_name)
                    p_als.remove_a(p_als[1 - idx_par_inter])  # remove the non-common allele between parent and child
                    c_als.remove_a(c_als[idx_common_child])  # remove the common allele between parent and child

                elif par_inter == 0:
                    continue_running = child_error()

            if not continue_running:  # more than 1 problematic child exist
                break

        if not continue_running:  # more than 1 problematic child exist
            break

    success = True if continue_running else False  # unnecessary variable, but I added for the readability

    if not success:
        errors_in_families[idx_fam] = ['All', '7']

    return success







