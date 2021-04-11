import json


# create dict with serology map, according the file "rel_ser_ser.txt"
def antigen_with_same_interpretation(file_ser_ser):
    dict_allele_antigen = {}
    with open(file_ser_ser) as info_file:
        for line in info_file:
            if not line.startswith('#'):
                line = line.strip().split(';')
                if len(line[1]) == 1:
                    line[1] = '0' + line[1]
                for i in range(2, len(line)):
                    if line[i] != '':
                        alleles = line[i].split('/')
                        for allele in alleles:
                            key, value = line[0] + '*' + line[1], line[0] + '*' + allele
                            if key not in dict_allele_antigen:
                                dict_allele_antigen[key] = [value]
                            else:
                                dict_allele_antigen[key].append(value)
    return dict_allele_antigen


ser_file = "../GR_code/GG_GRAMM/data/rel_ser_ser.txt"
dict = antigen_with_same_interpretation(ser_file)
with open('serology_dict_old.json', 'w') as output_file:
    json.dump(dict, output_file)