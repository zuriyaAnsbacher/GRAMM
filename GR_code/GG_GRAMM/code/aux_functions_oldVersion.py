
# old version

def equal_al(al_1, al_2):
    """
    :param al_1: first allele
    :param al_2: second allele
    :return: true if the allele are identical or including, false otherwise
    """
    if al_1 == al_2 or al_1 is al_2 or al_1.startswith(al_2) or al_2.startswith(al_1):
        if bool(al_1) == bool(al_2):
            return True
    return False


def empty_dict(d):
    """
    check if dict is empty
    @param d: dict
    @return: True if empty, else False
    """
    for value in d.values():
        if value:
            return False
    return True


def who_fuller(d1, d2):
    if not d1:
        return 2
    if not d2:
        return 1
    bin_d1 = []
    bin_d2 = []
    for d in [d1, d2]:
        for value in d.values():
            if d is d1:
                bin_d1.append(bool(value))
            else:
                bin_d2.append(bool(value))
    if bin_d1.count(True) >= bin_d2.count(True):
        return 1
    else:
        return 2


def sep_child_chro(spec_d, cur_ls):
    """
    organize alleles_names to 4 parameters, instead of list
    @param spec_d: dict children
    @param cur_ls: list children
    @return: the 4 parameters
    """
    c_00 = spec_d[cur_ls[0]][0]
    c_01 = spec_d[cur_ls[0]][1]
    c_10 = spec_d[cur_ls[1]][0]
    c_11 = spec_d[cur_ls[1]][1]
    return c_00, c_01, c_10, c_11


def d_per_type(child_d, types):
    exist = []
    d_type = {}
    for child in child_d.keys():
        if child_d[child][types] not in exist:
            d_type[child] = child_d[child][types]
        exist.append(child_d[child][types])
    return d_type


def exist_homoz_d(spec_d):
    for als in spec_d.values():
        if als[0] == als[1]:
            return True, als[0]
    return False, None


def exist_homoz_l(par_alls):
    exis = []
    for al in par_alls:
        if al in exis:
            return True, al
        exis.append(al)
    return False, None


def child_from_spec_d(spec_d):
    cur_ls = []
    for key in spec_d:
        cur_ls.append(key)
    return cur_ls


def append_to_fm12(f1, f2, m1, m2, types, a1, a2, a3, a4):
    f1[types].append(a1)
    f2[types].append(a2)
    m1[types].append(a3)
    m2[types].append(a4)


def fix_indexes(ind1, ind2):
    inds =[]
    for ind in [ind1, ind2]:
        if ind1 == 3:
            ind = 2
        elif ind == 0 or ind == 3:
            ind = 1
        inds.append(ind)
    return inds[0], inds[1]


def high_res(str1, str2):
    if len(str1) >= len(str2):
        return str1
    return str2


def common_chr(chF, chM, emb_FM):
    if emb_FM[0] == 1:
        f_child = chF.ch1
        f_no_child = chF.ch2
    elif emb_FM[0] == 2:
        f_child = chF.ch2
        f_no_child = chF.ch1
    else:
        f_child = None
        f_no_child = None
    if emb_FM[1] == 1:
        m_child = chM.ch1
        m_no_child = chM.ch2
    elif emb_FM[1] == 2:
        m_child = chM.ch2
        m_no_child = chM.ch1
    else:
        m_child = None
        m_no_child = None
    return f_child, f_no_child, m_child, m_no_child


