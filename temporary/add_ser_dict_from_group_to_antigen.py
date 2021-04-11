import json

input_file = "/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_GRAMM/data/ser_dict_antigen2group.json"
output_file = "/home/zuriya/PycharmProjects/GR_Web/GR_code/GG_GRAMM/data/ser_dict_group2antigen.json"

new_dict = {}

with open(input_file) as input_f:
    orig_dict = json.load(input_f)

for key, value in orig_dict.items():
    if value not in new_dict:
        new_dict[value] = [key]
    else:
        new_dict[value].append(key)

with open(output_file, 'w') as output_f:
    json.dump(new_dict, output_f)