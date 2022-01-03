from GR_code.GG_GRAMM.code.als import Als


# old version

class DualChromosome:
    """
    this class represent 2 chromosome of father or mother
    each chromosome is a dictionary, which contains data of alleles_names
    """
    def __init__(self, al_types):
        self.ch1 = {}
        self.ch2 = {}
        for al in al_types:
            self.ch1[al] = Als()
            self.ch2[al] = Als()

    def emb_data(self, fam_d, cur_par, par_num):
        """
        embed alleles_names data in the chromosomes of parents
        one- permanent, other- two options (if there is homozygous allele- it permanent, too)
        :param fam_d: dict with allele data
        """
        perm = 0
        for types in fam_d[cur_par]:
            al1 = fam_d[cur_par][types][0]
            al2 = fam_d[cur_par][types][1]
            if al1 == al2:  # homozygous
                self.ch1[types].append(al1)
                self.ch2[types].append(al2)
            else:
                if perm == 0 and not(par_num == 2 and fam_d['F'][types] == fam_d['M'][types]):
                    self.ch1[types].append(al1)
                    self.ch2[types].append(al2)
                    perm = 1
                else:
                    self.ch1[types].extend([al1, al2])
                    self.ch2[types].extend([al1, al2])
