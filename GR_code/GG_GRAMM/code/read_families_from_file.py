import re
import csv
import json

from utils import convert_to_serology


def get_families(user_file, aux_tools):
    """
    get families data from user file
    :param user_file: user file
    :param aux_tools: tools for special scenarios. it contains:
        'amb' dict, that keep the ambiguity data, if exists.
        'is_serology' flag, to convert the allele to serology interpretation, if the user choose that
        'ser_dict' that contains the serology conversion rules
    :return: dict with data about the families, in this structure:
        keys: families indexes (1, 2..)
        values: dict for each family
            in these sub-dicts (for each family):
                keys: indexes of family members: F, M, 1, 2 ...
                values: alleles data
                    in these sub-sub dict (of alleles data for each family member):
                    keys: alleles names (A, B..)
                    values: list of numbers of each alleles (e.g: [02:01, 30:04])
    """
    file = open(user_file, 'r')
    reader = csv.reader(file)
    next(reader, None)  # skip headers
    families_dict = {}

    for line in reader:
        id_family = line[0]
        id_person = line[1]
        indv_alleles, data_exist = get_individual(line, aux_tools)
        if data_exist:
            if id_family not in families_dict:
                families_dict[id_family] = {}  # create a sub dict to the family in families_dict
            families_dict[id_family][id_person] = indv_alleles

    file.close()

    return families_dict


def get_individual(line, aux_tools):
    """
    get individual data
    :param line: line in file
    :param aux_tools: explained in "get_families" function
    :return: individual data, and flag that equal to False if there is no data about this individual
    """
    indv_alleles = {}
    data_exist = True

    alleles_map = {2: 'A', 4: 'B', 6: 'C', 8: 'DRB1', 10: 'DQB1'}
    for i in range(2, 11, 2):
        alleles_pair = line[i:i+2]
        new_alleles_pair = process_alleles(alleles_pair, alleles_map[i], aux_tools)
        indv_alleles[alleles_map[i]] = new_alleles_pair
    if not any(item for sublist in indv_alleles.values() for item in sublist):  # completely empty: {('',''),('','')..}
        data_exist = False

    return indv_alleles, data_exist


def process_alleles(pair, allele, aux_tools):
    """
    process alleles: replace irrelevant characters, remove ambiguity (if exists) and save it to a dict,
    convert to serology if needed
    :param pair: alleles pair
    :param allele: the specific allele of this pair (A/B..)
    :param aux_tools: explained in "get_families" function
    :return: alleles pair, after process
    """
    # TODO: add try-except and error message or error may not occurs?
    amb = aux_tools['amb']  # amb is a dict of ambiguity of alleles
    is_serology = aux_tools['is_serology']  # is_serology is boolean (serology samples or genetic)
    ser_dict = aux_tools['ser_dict']  # ser_dict is a dict for the conversion to serology interpretation

    al1, al2 = pair[0], pair[1]

    # in simulation files, al1='02:01', al2='', mean that it homozygous (they both 02:01)
    if al2 == '':
        al2 = al1

    new_als = []
    for al in [al1, al2]:
        al = str(al).replace('p', '')

        """
            Save the ambiguity in dict and remove it from al1, al2.
            That because in the stage of comparing between parents and children, it could disturb.
            For example: parent: A*02:BJFV, child:A*02:02, we will not recognize it could match.
            Before the insertion to GRIMM, we restore it to data.
        """
        # TODO: right now, the saving of the amb in dict doesn't have identification about the person
        #  (tiny probability that will be 2 relatives with same ambiguity).
        #  Maybe need - if I read the whole file at once?
        if ":" in al and al.split(':')[1].isupper():
            amb[allele + '*' + al.split(':')[0]] = al.split(':')[1]

        al = re.sub('[a-zA-Z]', '', al)  # remove ambiguity

        if len(al) == 1 and al.isdigit():  # one digit to two digits (2 -> 02)
            al = '0' + al

        if is_serology:
            al = convert_to_serology(ser_dict, allele, al)

        new_als.append(al)

    return new_als


def run_script():
    aux_tools_dict = {}
    aux_tools_dict['amb'] = {}
    aux_tools_dict['is_serology'] = True
    with open('/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_GRAMM/data/ser_dict_antigen2group.json') as ser_dict_path_anti2group:
        ser_dict_anti2group = json.load(ser_dict_path_anti2group)
    aux_tools_dict['ser_dict'] = ser_dict_anti2group

    get_families('/home/zuriya/PycharmProjects/GR_Web/static/example_file1.csv', aux_tools_dict)
