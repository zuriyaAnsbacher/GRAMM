import copy
from GR_code.GG_GRAMM.code.aux_functions import common_chr, who_fuller, high_res, equal_al
from GR_code.GG_GRAMM.code.als import Als


def add_or_remove(al_p, al_p2, c2_inpar, c2_outpar, checker):
    new = al_p + al_p2
    if new == c2_outpar and new != checker:
        c2_outpar.remove_a(al_p[0])
    elif new == c2_inpar and new != checker:
        c2_inpar.remove_a(al_p[0])


def data_pc_22(al_p, al_c, j, fm_in, fm_out, types):
    match, ind_c = al_p.match_par(al_c)
    if match == 1:  # 1 pair same, 1 pair different
        ind_p = al_p.index_a(al_c[ind_c])
        al_p[ind_p] = high_res(al_p[ind_p], al_c[ind_c])  # emb chromo with high res between par and child
        al_p2 = Als()
        al_p2.append(al_p[1 - ind_p])
        al_p.remove_a(al_p2[0])  # remove from par different allele
        al_c.remove_a(al_p[0])  # remove from child same allele
        c2_inpar = fm_out[j][types] if fm_out[j] else None  # the second chromo in cur par
        c2_outpar = fm_in[1 - j][types] if fm_in[1 - j] else None  # the second chromo in other par (that child inheritance)
        # try
        checker = fm_out[1 - j][types] if fm_out[1 - j] else None  # avoid case of incorrect deletion
        if c2_inpar and c2_outpar:
            add_or_remove(al_p, al_p2, c2_inpar, c2_outpar, checker)  # add or remove allele to the connected chromo


def add_data_child(child, chF, chM, emb_FM, al_types, d_chi):
    f_in, f_out, m_in, m_out = common_chr(chF, chM, emb_FM)  # f_in/m_in = chromos the child inheritance
    fm_in = [f_in, m_in]
    fm_out = [f_out, m_out]
    # if who_fuller(f_in, m_in) == 2:
    #     fm_in, fm_out = [m_in, f_in], [m_out, f_out]
    j = 0  # j = which iteration in pars loop
    for chr_par in fm_in:
        if not chr_par:  # TODO: Im not sure it could be True (maybe yes..)
            j += 1
            continue
        for types in al_types:
            al_c = d_chi[child][types]
            al_p = fm_in[j][types]
            if al_c.is_empty_a():
                continue
            if al_p.is_empty_a():
                if len(al_c) == 2 and al_c[0] == al_c[1]:
                    al_c.remove_a(al_c[1])
                chr_par[types] = copy.deepcopy(al_c)  # emb data from child to par
            elif len(al_c) == len(al_p) == 1:  # par and child has 1 chromo
                if al_c != al_p:
                    return None, None
                al_p[0] = high_res(al_p[0], al_c[0])  # emb chromo with high res between par and child
            elif len(al_c) == 2 and len(al_p) == 1:  # child has 2, par has 1
                # e.g. par: 01 , child: 01:08+02 ; so- par: 01:08, child: 02
                if al_p[0] in al_c:
                    al_p[0] = high_res(al_p[0], al_c[al_c.index_a(al_p[0])])
                    al_c.remove_a(al_p[0])
                else:
                    return None, None
            elif len(al_c) == len(al_p) == 2:  # par and child has 2 chromos
                data_pc_22(al_p, al_c, j, fm_in, fm_out, types)
            elif len(al_c) == 1 and len(al_p) == 2:
                if al_c[0] not in al_p:
                    return None, None
                if equal_al(al_p[0], al_p[1]) and al_p[0] != al_p[1]:
                    if len(al_c[0]) == 2 and equal_al(al_p[0], al_c[0]):  # e.g 01:02, 01:03  # TODO: it's not a good example! because they aren't equal
                        continue
                ind_p = al_p.index_a(al_c[0])
                al_p[ind_p] = high_res(al_p[ind_p], al_c[0])  # high res  # TODO: why not chr_par[types] = ... ?
                al_p2 = Als()
                al_p2.append(al_p[1 - ind_p])
                al_p.remove_a(al_p2[0])
                c2_inpar = fm_out[j][types] if fm_out[j] else None  # the second chromo in cur par
                c2_outpar = fm_in[1 - j][types] if fm_in[1 - j] else None # the second chromo in other par (that child inheritance)
                # try
                checker = fm_out[1 - j][types] if fm_out[1 - j] else None # avoid case of incorrect deletion
                if c2_inpar and c2_outpar:
                    add_or_remove(al_p, al_p2, c2_inpar, c2_outpar, checker)  # add or remove allele to the connected chromo

        j += 1
    return chF, chM



