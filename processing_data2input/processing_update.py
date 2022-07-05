import csv
import os.path
from os import path
import sys
import numpy as np


def add2file(files_folder, id_person, A1_data, B1_data, C1_data, DRB11_data, DQB11_data, A2_data, B2_data, C2_data,
             DRB12_data, DQB12_data, temp_f):
    """
    create a file from the data that the user inserted manually to web (each person at time)
    :param files_folder: files folder in server
    :param id_person: Father, Mother, or digit (if _child_)
    :param temp_f: name of the file for write the data from the user
    :param A1_data, .., DQB12_data : data of alleles_names
    """
    if id_person == "Father" or id_person == "Mother":
        id_person = id_person[0]  # Father -> F, Mother -> M

    # in first iteration, need to create the file and add headers
    if not path.exists(os.path.join(files_folder, temp_f)):
        with open(os.path.join(files_folder, temp_f), 'w') as f1:
            writer1 = csv.writer(f1)
            writer1.writerow(["FAMCODE", "BIRTHSEQ", "A1", "A2", "B1", "B2", "C1", "C2", "DRB11", "DRB12", "DQB11", "DQB12"])

    # from the second iteration onwards, add a line at time
    with open(os.path.join(files_folder, temp_f), 'a+') as f2:
        writer2 = csv.writer(f2)
        row = ["1", id_person, A1_data, A2_data, B1_data, B2_data, C1_data, C2_data, DRB11_data, DRB12_data, DQB11_data, DQB12_data]
        writer2.writerow(row)


def convert_gl_file_to_columns_file(reader, writer, match_als2col, headers, race_dict):
    """
    convert the inserted file, with format of gl string (A*02:01+A*30:04^B*40+B*51...)
    to columns format (02:01,30:04,40,51...)
    :param reader: reader of input file
    :param writer: csv writer, for the output file
    :param match_als2col: dict {"A": 2, "B": 4, "C": 6, "DRB1": 8, "DQB1": 10} for the indexes in output file
    :param headers: input file headers. for checking if races exist.
    :param race_dict: dict of races, for case that there are races in input file (so keep in this dict)
    """
    for line in reader:
        output_line = [''] * 12
        output_line[0] = line[0]  # id family.
        output_line[1] = line[1]  # id person

        all_als = line[2].split("^")  # A*02:01+A*30:04, B*40:01+B*51, C*03:04+C*14:02...
        for pair_als in all_als:  # A*02:01+A*30:04, for example
            al1, al2 = pair_als.split("+")[0], pair_als.split("+")[1]  # A*02:01, A*30:04

            if al1.split('*')[0] in match_als2col:
                col_index = match_als2col[al1.split('*')[0]]
            else:
                raise KeyError('The name of the inserted allele is not a name from: A, B, C, DRB1, DQB1')

            for ind, al in enumerate([al1, al2]):
                if "/" in al:  # ambiguity: A*03:01/A*03:02/...
                    output_line[col_index + ind] = al.split("*")[1].split(":")[0]
                else:  # high or low res: A*02:01, A*02
                    output_line[col_index + ind] = al.split("*")[1]

        writer.writerow(output_line)

        # in case of insertion races to input file, we put them on a dict
        # the key is the id of the family and the value is a list with the races.
        # example: {'1': ['CAU', 'AFA'], '2': ['CARB']}
        # and when we write the data to file (for GRIMM), we insert the races
        if 'RACE' in headers:
            # using try-except so the run does not crash if there is a mistake in races (the races aren't not critical)
            try:
                if line[0].isdigit() and line[0].lstrip('0') not in race_dict and len(line) == 4:
                    races = list(set(line[3].split('+')))
                    race_dict[line[0].lstrip('0')] = races
            except:  # races not in valid format/location.
                race_dict.clear()


def save_ambiguity_to_dict(pair_als, open_ambiguity, id_fam, id_person):  # A*03:01/A*03:02/A*03:08...
    """
    in some simulations files, the genotype is given with ambiguity, e.g: A*03:01/A*03:02/A*03:08...
    we remove all the ambiguities options (because in GRAMM we cant handle them), save them, and in the end of
    GRAMM (before GRIMM), we write all the options again.
    :param pair_als: pairs of allelles
    :param open_ambiguity: the dict that keep all the options of the ambiguity
    :param id_fam: family index
    :param id_person: person index (F, M, 1, 2..)
    """
    if id_person not in ['F', 'M']:
        id_person = 'children'

    if id_fam not in open_ambiguity:
        open_ambiguity[id_fam] = {}
    if id_person not in open_ambiguity[id_fam]:
        open_ambiguity[id_fam][id_person] = {'A': {}, 'B': {}, 'C': {}, 'DRB1': {}, 'DQB1': {}}

    # assumes the pair_als is 'A*03:01/A*03:02... + A*01:01/A*01:02...'
    for allele in pair_als.split('+'):  # A*03:01/A*03:02 ... , A*01:01/A*01:02...
        if allele == '' or 'UUUU' in allele:
            continue
        allele_name = allele.split('*')[0]  # A
        first2_digits = allele.split(':')[0].split('*')[1]  # 03
        ambiguities = [single.split(':')[1].rstrip() for single in allele.split('/')] # A*03:01/A*03:02/A*03:08 -> [01, 02, 08..]

        if first2_digits in open_ambiguity[id_fam][id_person][allele_name]:
            open_ambiguity[id_fam][id_person][allele_name][first2_digits].extend(ambiguities)
        else:
            open_ambiguity[id_fam][id_person][allele_name][first2_digits] = ambiguities

        # remove duplicates, while keep the order
        open_ambiguity[id_fam][id_person][allele_name][first2_digits] = \
            list(dict.fromkeys(open_ambiguity[id_fam][id_person][allele_name][first2_digits]))


def convert_PED_file_to_columns_file(reader, writer, match_als2col, race_dict):
    """
    convert the inserted file, with format of PED (Lorem simulations)
    to columns format (02:01,30:04,40,51...)
    :param reader: reader of input file
    :param writer: csv writer, for the output file
    :param match_als2col: dict {"A": 2, "B": 4, "C": 6, "DRB1": 8, "DQB1": 10} for the indexes in output file
    :param race_dict: dict of races. in this case' of PED, it always CAU
    """
    csv.field_size_limit(sys.maxsize)  # avoid from error when reading file with huge lines (in some of the sim files)
    output_lines_lst = []
    open_ambiguity = {'is_amb': False}  # save the ambiguity in this dict for re-write before the insertion to GRIMM

    for line in reader:
        output_line = [''] * 12
        id_fam = line[0].lstrip('0')
        id_person = line[1]

        if line[2] == 'NULL':  # NULL exist in parents rows only
            id_person = 'F' if line[4] == '1' else 'M'

        output_line[0] = id_fam
        output_line[1] = id_person

        all_als = line[5].split("^")
        for pair_als in all_als:
            if "|" in line[5]:  # ambiguity (in parents): A*03:01+A*68:02|A*03:01+A*68:18N|...
                # get al1 = A*03, al2 = A*68
                al1, al2 = pair_als.split('|')[0].split("+")[0].split(":")[0], \
                           pair_als.split('|')[0].split("+")[1].split(":")[0]
            elif "/" in line[5]:  # ambiguity (in children): A*03:01/A*03:02/A*03:08 ...
                open_ambiguity['is_amb'] = True
                # if id_person in ['F', 'M']:  # save amb only for parents
                save_ambiguity_to_dict(pair_als, open_ambiguity, id_fam=id_fam, id_person=id_person)
                al1, al2 = pair_als.split("+")[0].split("/")[0].split(":")[0], \
                           pair_als.split("+")[1].split("/")[0].split(":")[0]
            else:  # no ambiguity: A*03:02+A*01:01 (or A*03+A*01)
                al1, al2 = pair_als.split("+")[0], pair_als.split("+")[1]

            al1, al2 = al1.replace('g', '').replace('p', ''), al2.replace('g', '').replace('p', '')

            if al1.split('*')[0] in match_als2col:
                col_index = match_als2col[al1.split('*')[0]]
            else:
                raise KeyError('The name of the inserted allele is not a name from: A, B, C, DRB1, DQB1')

            for ind, al in enumerate([al1, al2]):
                output_line[col_index + ind] = al.split("*")[1]

        output_lines_lst.append(output_line)

    race_dict['all'] = 'CAU'  # in simulations, all the races are CAU. GRAMM know to handle that (by the key "all")

    # simulation file (PED) need to sort by family id (because parents in start, children in the end)
    # so: sort and then write to file
    output_lines_lst = sorted(output_lines_lst, key=lambda x: int(x[0]))  # int for numeric sort (1, 2.. , not 1, 10..)
    writer.writerows(output_lines_lst)

    if open_ambiguity['is_amb']:
        return open_ambiguity
    else:
        return None


def convert_old_columns_file_to_new_columns_file(reader, writer, headers, race_dict):
    """
    remove from file the races, if exist. if not it the same format of columns
    :param reader: reader of input file
    :param writer: csv writer, for the output file
    :param headers: input file headers. for checking if races exist.
    :param race_dict: dict of races, for case that there are races in input file (so keep in this dict)
    """
    for line in reader:
        output_line = line[:12]  # output_line is identical to line, only without races (if exist, in index 12)
        writer.writerow(output_line)

        # explanation to race_dict in "convert_gl_file_to_columns_file" function
        if 'RACE' in headers:
            # using try-except so the run does not crash if there is a mistake in races (the races are not critical)
            try:
                if line[0].isdigit() and line[0].lstrip('0') not in race_dict and len(line) == 13:
                    races = list(set(line[12].split('+')))
                    race_dict[line[0].lstrip('0')] = races
            except:
                race_dict.clear()


def convert_DD_or_Israel_file(reader, writer, match_als2col):
    open_ambiguity = {'is_amb': True}
    for line in reader:
        output_line = [''] * 12
        id_fam = line[0].lstrip('0')
        id_person = line[1]
        output_line[0] = id_fam
        output_line[1] = id_person

        for i in range(2, 12, 2):
            pair_als = line[i] + '+' + line[i + 1]
            # if id_person in ['F', 'M']:
            save_ambiguity_to_dict(pair_als, open_ambiguity, id_fam=id_fam, id_person=id_person)

            al1, al2 = pair_als.split("+")[0].split("/")[0].split(":")[0], \
                       pair_als.split("+")[1].split("/")[0].split(":")[0]

            if al1 != '' and al2 == '':
                al2 = al1

            if '*' in al1:  # else: al1 = ''
                if al1.split('*')[0] in match_als2col:
                    col_index = match_als2col[al1.split('*')[0]]
                else:
                    raise KeyError('The name of the inserted allele is not a name from: A, B, C, DRB1, DQB1')

                for ind, al in enumerate([al1, al2]):
                    output_line[col_index + ind] = al.split("*")[1]

        writer.writerow(output_line)

    return open_ambiguity


def basic_csv2format(input_path, output_path, race_dict):
    """
    Process the basic file to a valid format (that the code know how to work with it).
    Basic file could arrives in 3 option:
    1. data in columns. It doesnt need change
    2. data in gl string (ex: 1,F, A*02+A*03^B*01:01+B*03:01^...)
    3. Loren simulations (PED format)
    :param input_path: input_path
    :param output_path: output_path
    :param race_dict: empty dict for races.
    """
    # keys: alleles_names, values: column index in file (each have the index i, i+1)
    match_als2col = {"A": 2, "B": 4, "C": 6, "DRB1": 8, "DQB1": 10}
    open_ambiguity_exist = False

    wf = open(output_path, 'w')
    writer = csv.writer(wf)
    writer.writerow(
        ["FAMCODE", "BIRTHSEQ", "A1", "A2", "B1", "B2", "C1", "C2", "DRB11", "DRB12", "DQB11", "DQB12"])

    csv.field_size_limit(100000000)

    with open(input_path, 'r') as input_p:
        temp_reader = csv.reader(input_p)  # reader just for split by second_line (to files alleles_names: gl, PED, columns)
        next(temp_reader, None)
        second_line = next(temp_reader)

    with open(input_path, 'r') as input_p:
        reader = csv.reader(input_p)
        headers = next(reader)

        if "^" in second_line[2]:  # the input file is in format gl string
            convert_gl_file_to_columns_file(reader, writer, match_als2col, headers, race_dict)
        elif "^" in second_line[5]:  # the input is in format PED (Loren simulations)
            open_ambiguity = convert_PED_file_to_columns_file(reader, writer, match_als2col, race_dict)
            open_ambiguity_exist = True
        # Austrilian or Israel which were converted from serology to genetic by py-ard
        elif "/" in second_line[2] or "/" in second_line[8]:
            open_ambiguity = convert_DD_or_Israel_file(reader, writer, match_als2col)
            open_ambiguity_exist = True
        else:  # the input is in format of columns
            convert_old_columns_file_to_new_columns_file(reader, writer, headers, race_dict)
    wf.close()

    if open_ambiguity_exist:
        return open_ambiguity
    else:
        return None


def split_gl(gl_str):
    """
    split gl string to separate data : A*02:01+A*30:04^B*40:01+B*51:22 ... ---> 02:01, 30:04, 40:01, 51:22 ...
    :param gl_str: gl string
    :return: separate data
    """
    gl_str = gl_str.strip()
    pairs = gl_str.split("^")
    d = {"A": [], "B": [], "C": [], "DRB1": [], "DQB1": []}
    for pair in pairs:
        # example: pair = A*02:01+A*03:01; allele = A; val1 = 02:01; val2 = 03:01
        allele = pair.split("*")[0]
        val1 = pair.split("+")[0].split("*")[1]
        val2 = pair.split("+")[1].split("*")[1]
        d[allele].extend([val1, val2])
    for key in d.keys():
        if not d[key]:  # [] change to ["", ""]
            d[key] = ["", ""]
    return d["A"][0], d["B"][0], d["C"][0], d["DRB1"][0], d["DQB1"][0], d["A"][1], d["B"][1], d["C"][1], \
           d["DRB1"][1], d["DQB1"][1]


def create_glstring(A1, A2, B1, B2, C1, C2, DRB11, DRB12, DQB11, DQB12):
    """
    create gl string from split data (to show to user on the screen)
    e.g: 02:01, 30:04, 40:01, 51:22 ... ---> A*02:01+A*30:04^B*40:01+B*51:22 ...
    :return: gl_string
    """
    als_names = {A1: "A", A2: "A", B1: "B", B2: "B", C1: "C", C2: "C", DRB11: "DRB1", DRB12: "DRB1",
                 DQB11: "DQB1", DQB12: "DQB1"}
    gl_lst = []
    for pair in [[A1, A2], [B1, B2], [C1, C2], [DRB11, DRB12], [DQB11, DQB12]]:
        pair_lst = []
        for single in pair:
            if single == "":
                continue
            single = str(single)
            full_single = als_names[single] + '*' + single  # e.g: A --> A*02:08
            pair_lst.append(full_single)
        pair_gl = '+'.join(pair_lst)    # e.g: [A*02:08, A*03:01] --> A*02:08+A*03:01
        if pair_gl:  # if not empty
            gl_lst.append(pair_gl)
    gl_string = '^'.join(gl_lst)
    return gl_string


def add_races_from_manual_insertion(race_str, race_dict):
    """
    in case of insertion races when the insertion is manually, we put them on a dict
    because the manual insertion is possible by one family each time,
    there is one key only ('1') and the races are the value, in list (e.g: CAU;NAM;AAFA)
    :param race_str: the races, separate by ';'
    :param race_dict: race dict
    """
    try:
        if race_str == '':  # no races
            return
        race_str = race_str.split(';')
        race_dict['1'] = race_str
    except:
        race_dict.clear()
