import json

input_file = "/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_GRAMM/data/serology_dict_new.json"
output_file = "../GR_code/GG_GRAMM/data/ser_dict_antigen2group.json"

new_dict = {}

with open(input_file) as input_f:
    orig_dict = json.load(input_f)

for key, value in orig_dict.items():
    for val in value:
        if len(val.split('*')[1]) == 1:
            val = val.split('*')[0] + '*' + '0' + val.split('*')[1]
        new_dict[val] = key

with open(output_file, 'w') as output_f:
    json.dump(new_dict, output_f)




