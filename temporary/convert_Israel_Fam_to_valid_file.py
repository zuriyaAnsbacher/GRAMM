import csv


def convert_row(row_):
    new_row_ = [0]*10
    for ind, cell in enumerate(row_):
        if 'Fam' in cell:
            continue
        elif cell == '':
            new_row_[ind - 1] = cell
        else:
            new_row_[ind - 1] = cell.split('*')[1]
    return new_row_


input_file = "Israel_Fam.csv"
output_file = "Israel_Fam_valid_format.csv"

with open(input_file, 'r') as input_f:
    with open(output_file, 'w') as output_f:
        reader = csv.reader(input_f, delimiter=',')
        writer = csv.writer(output_f, delimiter=',')

        writer.writerow(['FAMCODE', 'BIRTHSEQ', 'A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'DRB11', 'DRB12', 'DQB11', 'DQB12'])

        fam_ind = 1
        child_ind = 1

        for row in reader:
            famcode = int(row[0].split('m')[1])
            if famcode > fam_ind:
                fam_ind += 1
                child_ind = 1
            alleles_row = convert_row(row)
            indexes_row = [str(famcode), str(child_ind)]
            indexes_row.extend(alleles_row)
            writer.writerow(indexes_row)
            child_ind += 1


