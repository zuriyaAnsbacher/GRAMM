from GR_code.GG_GRAMM.code.aux_functions import equal_al, empty_dict
from GR_code.GG_GRAMM.code.als import Als


def emb_wp(chF, chM, child_d, al_types):
    """
    embed children to parents chromosomes, in case that parents exist
    @param chF: father chromosome
    @param chM: mother chromosome
    @param child_d: dict of child
    @param al_types: alleles types (A, B ...)
    @return: list with the embedding of the child
    """
    i = 0
    m_before = False
    flag = 0
    first = True
    emb_FM = [0, 0]  # emb_FM signs which chromo of parents the child emb. if chF2, chM1: [2,1]
    pars = [chF, chM]
    if not (empty_dict(chM.ch1) or empty_dict(chM.ch2)):  # if M full and F empty, begin with M
        pars = [chM, chF]
        m_before = True  # sign the opposite list. in the end of code - flip it
    for par in pars:
        was_emb = False
        if empty_dict(par.ch1) and empty_dict(par.ch2):
            emb_FM[i] = 1  # emb (random) in first chromo
            was_emb = True
        elif not (empty_dict(par.ch1) or empty_dict(par.ch2)):  # 2 chrom of par are full
            d1 = par.ch1
            d2 = par.ch2
            for types in al_types:
                if not was_emb and len(d1[types]) == len(d2[types]) == 1 and not \
                        equal_al(d1[types][0], d2[types][0]):  # not embedded, par has singles and not homozygous
                    par_pair = Als()
                    par_pair.extend([d1[types][0], d2[types][0]])
                    match, m_par = child_d[types].match_par(par_pair)
                    if len(child_d[types]) == 2 and match == 1:
                        emb_FM[i] = m_par + 1
                        was_emb = True
                        del_al = child_d[types][0] if child_d[types][0] in par_pair else child_d[types][1]
                        del child_d[types][child_d[types].index_a(del_al)]  # remove the emb allele
                    elif len(child_d[types]) == 1:  # child has 1 option to this allele (the first allele was embedded)
                        ind = par_pair.index_a(child_d[types][0])
                        if ind != -1:
                            emb_FM[i] = ind + 1
                            was_emb = True
            if not was_emb and flag != 1 and first == 0:
                flag = 1
                i = i - 1
                pars.append(par)  # add to more one running on this parents, to success the emb

        else:  # 1 chrom full, one empty
            full_d = par.ch1 if not empty_dict(par.ch1) else par.ch2
            full_ind = 0 if not empty_dict(par.ch1) else 1
            for types in al_types:
                if not was_emb and not child_d[types].sub_lst(full_d[types]):
                    emb_FM[i] = 2 - full_ind
                    was_emb = True
        first = 1
        i += 1
    if m_before:  # indexes opposite
        emb_FM.reverse()
    if flag == 1:  # and emb_FM[0] != 0:
        emb_FM.reverse()

    return emb_FM
