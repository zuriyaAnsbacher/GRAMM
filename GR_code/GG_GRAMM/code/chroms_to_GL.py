import os
import copy
import itertools
import json
from GR_code.GG_GRAMM.code.aux_functions import empty_dict


# if par has 1 full 1 empty (chromo), duplicate the full (for GG_GRIMM) and sign it (by f_dupl/m_dupl)
def copy_full2empty(c1f, c2f, c1m, c2m):
    empty1f = False
    empty1m = False
    if empty_dict(c1f):
        c1f = copy.deepcopy(c2f)
        empty1f = True
    if empty_dict(c2f):
        c2f = copy.deepcopy(c1f)
        empty1f = True
    if empty_dict(c1m):
        c1m = copy.deepcopy(c2m)
        empty1m = True
    if empty_dict(c2m):
        c2m = copy.deepcopy(c1m)
        empty1m = True
    return c1f, c2f, c1m, c2m, empty1f, empty1m


# ex. A*02: + [01, 05] -> A*02:01/A*02:05
def options2gl(als, options):
    gl = als + options[0]
    if len(options) == 1:
        return gl
    for op in options[1:]:
        gl += '/' + als + str(op)
    return gl


def open_serology_options(c1f, c2f, c1m, c2m, types, s_dict):
    for chromosome in [c1f, c2f, c1m, c2m]:
        for type_ in types:
            for allele in chromosome[type_]:
                group = type_ + '*' + allele
                # check if the value (group) is in dict (low res only)
                # if it exists, open to all the serology options and remove the original value
                # (it's only group name, not an antigen)
                # example: A*10 appears, so add A*25, A*25, A*34, A*66 (only the digits) and remove A*10
                if ":" not in allele and group in s_dict:
                    for ser_value in s_dict[group]:
                        chromosome[type_].append(ser_value.split('*')[1])
                    if group != "B*14" and group != "B*15" and group != "DRB1*03":
                        chromosome[type_].remove(group.split('*')[1])


# create gl string with high res of specific allele
def single2high(als, types, amb_d, d_low2high):
    in_low2high = True
    # empty data
    if als == 'UUUU' or als == '':
        return ''.join([str(types), '*UUUU'])

    # in high res. i.g: 01:02, return it without processing
    if als.count(':') == 1 and not str(als).endswith(':'):
        return ''.join([str(types), '*', str(als)])

    # low res, i.g: 02. check options in d_low2high (dict with the options)
    if ':' not in als:
        try:
            options = d_low2high[''.join([str(types), "*", str(als)])]
        except KeyError:
            return ''.join([str(types), '*', str(als)])  # if not in the dict, return allele with low res

    # ambiguity, i.g: 02:APC (in the code its 02:)
    elif str(als).endswith(':'):
        # amb_d contain the ambiguity of alleles in cur family that removed before
        # i.g:  in data: A*02:APC , after remove: 02: , in amb_d: {"A*02": APC}
        ambiguity = amb_d[types + '*' + als.rstrip(':')]
        try:
            # there are values in d_low2high in type list and there are in type dict
            if isinstance(d_low2high[ambiguity], list):
                options = d_low2high[ambiguity]
            else:  # its dict
                options = list(d_low2high[ambiguity].keys())
        # this ambiguity (seq of letters) not in d_low2high.
        # so we open all the options by the low res
        # (i.g: we have A*02:APC but APC not in d_low2high, so we look for all options of 02--> 02:01/03....)
        except KeyError:
            try:
                options = d_low2high[''.join([str(types), "*", str(als).rstrip(':')])]
            except KeyError:
                return ''.join(([str(types), '*', str(als).rstrip(':')]))  # return allele with low res

    elif ':' in als and len(als) > 4:
        return ''.join([str(types), '*', str(als[:5])])

    else:
        return -1

    # als = 'A*02:' for example
    als = ''.join([str(types), "*", str(als)]) if str(als).endswith(':') else ''.join([str(types), "*", str(als), ':'])
    if in_low2high:
        high_gl = options2gl(als, options)
        return high_gl
    return -1


def gl_cartesian_product(c1, c2, al_types):
    lst_options = []
    gl_options = []
    for types in al_types:
        lst_types = []
        for al1 in c1[types]:
            for al2 in c2[types]:
                lst_types.append(str(al1) + '+' + str(al2))
        lst_options.append(lst_types)
    for op in itertools.product(*lst_options):
        gl_options.append(op)
    return gl_options


def gl_cartesian_product_1empty(c, al_types):
    lst_options = []
    gl_options = []
    for types in al_types:
        al_lst = []
        for al in c[types]:
            al_lst.append(al)
        lst_options.append(al_lst)
    for op in itertools.product(*lst_options):
        op_lst = []
        for al in op:
            al1 = str(al) + "+" + str(al)
            op_lst.append(al1)
        gl_options.append(op_lst)
    return gl_options


# create gl string to one person
def gl_for_person(c1, c2, al_types, amb_d, is1empty, d_low2high):
    ignore_in_grimm = {}
    for types in al_types:
        if (c1[types] == [] or c1[types] == ['']) and (c2[types] == [] or c2[types] == ['']):
            # no data about allele in 2 chroms
            del c1[types]
            del c2[types]
            al_types.remove(types)
            continue
        if c1[types] == [] or c1[types] == ['']:  # no data in this chrom but exist data in other chrom
            c1[types] = ['UUUU']  # sign to GG_GRIMM
            ignore_in_grimm[types] = -1
        elif c2[types] == [] or c2[types] == ['']:
            c2[types] = ['UUUU']
            ignore_in_grimm[types] = -2
        if len(c1[types]) == 2 and c1[types] == c2[types] and not is1empty:
            # case of 2 same als in 2 chromo. divide 1 for each chromo
            # ex: c1[05, 07] c2[05, 07] --> c1[05] c2[07]. in GG_GRIMM it sign with open phase
            second = c1[types][1]
            c1[types].remove(second)
            c2[types].clear()
            c2[types].append(second)
    if not is1empty:
        options = gl_cartesian_product(c1, c2, al_types)
    else:
        options = gl_cartesian_product_1empty(c1, al_types)
    new_options = []

    for op in options:
        j = 0  # index to types
        new_op = []
        for pair in op:
            al1 = single2high(pair.split('+')[0], al_types[j], amb_d, d_low2high)  # each allele: from low to high
            al2 = single2high(pair.split('+')[1], al_types[j], amb_d, d_low2high)  # each allele: from low to high
            if al1 == -1 or al2 == -1:  # error in convert low res to high res
                return -1
            new_pair = ''.join([al1, "+", al2])
            new_op.append(new_pair)
            j += 1
        new_op = '^'.join(new_op)
        new_options.append(new_op)
    return new_options


def create_bin(al_types, c1, c2):
    # binary1 is 0/1 for each allele (al1, al2..)
    # binary2 is 0/1 for pair alleles (al1->al2, al2->al3...)
    binary1 = [0] * len(al_types)
    i = 0
    for types in al_types:
        if len(c1[types]) == 2 and c1[types] == c2[types]:
            binary1[i] = 1
        i += 1
    a = binary1[4]
    binary1[4] = binary1[3]
    binary1[3] = a   # replace 2 last because in GG_GRIMM, DQ before DR
    return binary1


# CHECK THIS FUNCTION
# if in some allele, there is data only in one chromosome, delete the data from the full allele
# that because GRIMM can not get data in this way (A*02:01+A*UUUU)
# UUUU sign empty (no data)
def delete_data_from_one_chromosome(c1, c2):
    for key, value in c1.items():
        if value == [''] and c2[key] != []:
            c2[key] = ['']
        elif c2[key] == [] and value != ['']:
            c1[key] = ['']


# def create_bin(al_types, c1, c2):
#     # binary1 is 0/1 for each allele (al1, al2..)
#     # binary2 is 0/1 for pair alleles (al1->al2, al2->al3...)
#     binary1 = [0] * len(al_types)
#     binary2 = [0] * (len(al_types) - 1)
#     i = 0
#     for types in al_types:
#         if len(c1[types]) == 2 and c1[types] == c2[types]:
#             binary1[i] = 1
#         i += 1
#     for j in range(len(binary2)):
#         if any([binary1[j], binary1[j+1]]):
#             binary2[j] = 1
#     return binary2


def write2gl_file(gl_file, bin_d, id_, gl, binary, race_dict):
    if gl != -1:
        bin_d[id_] = binary

        for single_gl in gl:
            if 'all' in race_dict:  # Loren simulations. all races are CAU
                gl_file.write(id_ + ',' + single_gl + ',' + 'CAU' + ',\n')
            elif str(id_[0]) in race_dict:  # add the races of the family to the file for GRIMM
                races = ';'.join(race_dict[str(id_[0])])
                gl_file.write(id_ + ',' + single_gl + ',' + races + ',\n')
            else:  # no races to this family
                gl_file.write(id_ + ',' + single_gl + '\n')


# create gl string for insert to GG_GRIMM
def create_gl(father, mother, al_types, id_fam, amb_d, gl_file, bin_d, d_low2high, race_dict, is_ser, ser_dict_g2a):
    c1f, c2f, c1m, c2m = father.ch1, father.ch2, mother.ch1, mother.ch2

    for c1, c2 in [(c1f, c2f), (c1m, c2m)]:
        delete_data_from_one_chromosome(c1, c2)

    id_f = str(id_fam) if '.' in str(id_fam) else str(id_fam) + '.0'
    id_m = str(id_fam).replace('.0', '.1') if '.' in str(id_fam) else str(id_fam) + '.1'
    bin_f = create_bin(al_types, c1f, c2f)
    bin_m = create_bin(al_types, c1m, c2m)

    if is_ser:
        open_serology_options(c1f, c2f, c1m, c2m, al_types, ser_dict_g2a)

    c1f, c2f, c1m, c2m, empty1f, empty1m = copy_full2empty(c1f, c2f, c1m, c2m)
    gl_f = gl_for_person(c1f, c2f, al_types, amb_d, empty1f, d_low2high)
    gl_m = gl_for_person(c1m, c2m, al_types, amb_d, empty1m, d_low2high)
    write2gl_file(gl_file, bin_d, id_f, gl_f, bin_f, race_dict)
    write2gl_file(gl_file, bin_d, id_m, gl_m, bin_m, race_dict)
    if gl_f == -1:
        return -1, 0, 0
    if gl_m == -1:
        return -2, 0, 0
    return 0, empty1f, empty1m




