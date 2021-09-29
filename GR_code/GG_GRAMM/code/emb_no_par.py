from GR_code.GG_GRAMM.code.aux_functions import sep_child_chro, d_per_type, exist_homoz_d, exist_homoz_l, \
    child_from_spec_d, append_to_fm12, fix_indexes
from GR_code.GG_GRAMM.code.als import Als
from GR_code.GG_GRAMM.code.utils import child_ls_merge


def div_2als(gr0, gr1, als, after):
    """
    divide 2 als to 2 existing groups (with data)
    @param gr0: first group
    @param gr1: second group
    @param als: allele to emb
    @param after: list for case need more running to this als
    @return: 2 groups
    """
    gr = gr0 + gr1
    match, ind = gr.match_par(als)
    if match == 0:
        after.append(als)
    elif match == 1:
        if gr0.index_a(als[ind]) != -1 and gr1.index_a(als[ind]) == -1:
            gr1.append(als[1 - ind])
        else:
            gr0.append(als[1 - ind])
    return gr0, gr1


def div_4als(spec_d):
    """
    divide 4 alleles to 2 groups, with keep the logic of inheritance from parents
    for example : 01+02, 02+03, 03+04 --> [01, 03] [02, 04]
    @param spec_d: dict with data of children alleles
    @return: 2 groups
    """
    gr0 = Als()
    gr1 = Als()
    is_ex, homoz = exist_homoz_d(spec_d)
    if is_ex:  # there is homozygous
        gr0.append(homoz)
        gr1.append(homoz)
        for als in spec_d.values():
            for al in als:
                if al not in gr0 and al not in gr1:  # al not in groups. add to group that small than 2 als
                    if len(gr0) < 2:
                        gr0.append(al)
                    elif len(gr1) < 2:
                        gr1.append(al)
            if len(gr0) == len(gr1) == 2:  # the groups are full
                break
    else:  # no homozygous
        """
        "after" is a list help to make more running if at first the als didnt embed
        for exa: first run: gr0[01], gr1[02]. als = [03, 04], so cant emb
                 second run: g0[01, 03], gr1[02] so can emb als
        """
        after = []
        i = 0
        for child, als in spec_d.items():
            if i == 0:
                gr0.append(als[0])
                gr1.append(als[1])
            else:
                gr0, gr1 = div_2als(gr0, gr1, als, after)
            i += 1
        if len(after) > 0:
            gr0, gr1 = div_2als(gr0, gr1, after[0], [])
    return gr0, gr1


def emb_4als_in_pars(chF, chM, child_ls, child_d, al_types):
    """
    in case that in some allele type there are 4 different numbers - embed in parents chromosome.
    must be coherent with data about children
    e.g: c1:[01, 02], c2:[02, 03], c3:[01, 04], so -> par1: 01~03, par2: 02~04
    @param chF: chromosome of father
    @param chM: chromosome of mother
    @param child_ls: children list
    @param child_d: children dict
    @param al_types: all types: A, B, C, DR, DB
    @return: the division of 4 alleles to parents chromosomes
    """
    f1, f2, m1, m2 = chF.ch1, chF.ch2, chM.ch1, chM.ch2
    find_4 = False
    is_deter = True
    type_emb = None
    for types in al_types:
        if find_4:
            break
        lst = child_ls_merge(child_ls, child_d, types)
        if len(lst) == 4:  # 4 different als
            spec_d = d_per_type(child_d, types)  # dict without repetitions in specific allele
            cur_ls = child_from_spec_d(spec_d)  # list of keys (children) in spec_d
            if len(cur_ls) == 2:  # 2 children
                c_00, c_01, c_10, c_11 = sep_child_chro(spec_d, cur_ls)  # the 4 different alleles
                append_to_fm12(f1, m1, f2, m2, types, c_00, c_01, c_10, c_11)
                if len(set(lst)) == 4:  # the second emb is probability, so add more option
                    f2[types].append(c_11)
                    m2[types].append(c_10)
                    is_deter = False
                find_4 = True
                type_emb = types
            else:  # num children > 2
                gr1, gr2 = div_4als(spec_d)  # divide 4 als to 2 groups
                if len(gr1) < 2 or len(gr2) < 2:
                    return None, None
                append_to_fm12(f1, f2, m1, m2, types, gr1[0], gr1[1], gr2[0], gr2[1])
                find_4 = True
                type_emb = types
    return type_emb, is_deter


def emb_np(types, als_of_ty, is_deter, f1, f2, m1, m2):
    """
    embed child to chromosomes from parents that he inherited. case of no data about parents
    @param types: one type from A/B/C/DR/DB
    @param als_of_ty: dict of allles of specific type
    @param is_deter: False in case of only 2 children, who inherited different chromosomes. True otherwise.
    @param f1: first chromosome of father
    @param f2: second chromosome of father
    @param m1: first chromosome of mother
    @param m2: second chromosome of mother
    @return: list (len 2). first index - which chromosome child inherited from f (1/2), second - from m.
    """
    # TODO: what about children without data in this allele?
    emb_FM = [0, 0]  # emb_FM signs which chromo of parents the child emb. if chF2, chM1: [2,1]
    f_als = f1[types] + f2[types]
    m_als = m1[types] + m2[types]
    als_par = f_als + m_als
    is_ex, homoz = exist_homoz_l(als_par)  # TODO: but in als_par there are 4 alleles, no?
    if is_ex and homoz in als_of_ty:  # if was embed homozygous and this al exist in child
        del als_of_ty[als_of_ty.index_a(homoz)]  # remove homoz' and emb the second al
        if als_of_ty[0] in f_als:
            emb_FM = [f_als.index_a(als_of_ty[0]) + 1, m_als.index_a(homoz) + 1]
        else:
            emb_FM = [f_als.index_a(homoz) + 1, m_als.index_a(als_of_ty[0]) + 1]
    else:
        ind1 = als_par.index_a(als_of_ty[0])
        ind2 = als_par.index_a(als_of_ty[1])
        emb_FM = [min(ind1, ind2) + 1, max(ind1, ind2) - 1]
        if not is_deter:  # e.g: f[01][02,03] m[08][02,03]
            ind1, ind2 = fix_indexes(ind1, ind2)  # TODO: it doesnt work well in some cases
            emb_FM = [ind1, ind2]
    was_emd = True
    return emb_FM
