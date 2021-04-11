import re


class simReader:
    """
    this class reads file to family list of IDs and list of alleles
    file columns: family ID, seq in family (F/M/child num), 4 cols for each gene: allele-1*, allele-2*,allele-1,allele-2
    """
    def __init__(self, data_path):
        self._f = open(data_path, 'rt')
        dummy = next(self._f)
        self._curr_line = (next(self._f)).split(',')
        self._currIndex = 1
        self._end = False

    @staticmethod
    def extract_als(al_list):
        """
        extract als in format to work with, ignoring letters that represent ambiguity.
        first 2 are updated and will be taken if exist.
        :param al_list: from file
        :return: the 2 alleles that were extracted as the best of 4 in file.
        """
        al1 = al_list[0]
        al2 = al_list[1]
        if al_list[0] == '':  # no high res
            al1 = al_list[2]
            al2 = al_list[3]
        if al2 == '':  # homozygous
            al2 = al1
        al1 = str(al1).replace('p', '')
        al2 = str(al2).replace('p', '')
        new_al1 = re.sub('[a-zA-Z]', '', al1)  # remove ambiguity
        new_al2 = re.sub('[a-zA-Z]', '', al2)  # remove ambiguity
        if new_al1 == '\n':
            new_al1 = ''
        if new_al2 == '\n':
            new_al2 = ''

        if len(new_al1) == 1 and new_al1.isdigit():  # one digit to two digits (2 -> 02)
            new_al1 = '0' + new_al1
        if len(new_al2) == 1 and new_al2.isdigit():
            new_al2 = '0' + new_al2

        return new_al1, new_al2

    @staticmethod
    def parents_num(family_ids):
        """
        :param family_ids: ids list of family
        :return: num of parents
        """
        if 'M' in family_ids and 'F' in family_ids:
            return 2
        elif 'M' in family_ids or 'F' in family_ids:
            return 1
        return 0

    def get_family(self):
        """
        read from current place in file all rows with the family id
        :return: list of ids and all lists of alleles for each
        """
        if self._end or self._curr_line[0] == '\n':
            print('end of file')
            return None, None, None, None
        fam_id = self._curr_line[0]
        i = self._currIndex
        family_ids = [fam_id]
        family_als = []

        seq, ind = self.get_individual(self._curr_line)
        family_als.append(ind)
        family_ids.append(seq)
        self._end = True
        for line in self._f:
            self._end = False
            line_list = line.split(',')
            if line_list[0] == fam_id:
                seq, ind = self.get_individual(line_list)
                if seq != "empty_member":  # if its not an empty list
                    # identical child - ignore from one of them
                    index_ind = family_als.index(ind) if ind in family_als else -1
                    if index_ind == -1 or not (isinstance(seq, (int, float))) \
                            or not (isinstance(family_ids[index_ind + 1], (int, float))):
                        family_als.append(ind)
                        family_ids.append(seq)
            else:
                self._curr_line = line_list
                par_num = self.parents_num(family_ids)
                return family_ids, family_als, par_num, None
        par_num = self.parents_num(family_ids)
        if self._end:
            return None, None, None, None
        return family_ids, family_als, par_num, None

    def get_individual(self, line):
        """
        create list of alleles for individual from row in sheet.
        :param row: of individual's data.
        :return: individuals id and list of alleles.
        """
        seq = line[1]  # individual seq number ID
        ind = []
        if len(line[-1]) > 1 and line[-1][-1] == '\n':
            line[-1] = line[-1][:-1]
        for i in range(2, 21, 4):
            val = line[i:i+4]
            ind.append(list(self.extract_als(val)))  # adding individual's alleles
        if not any(item for sublist in ind for item in sublist):  # completely empty: [("",""),("","")..]
            seq = None  # no data about this member
        return seq, ind
