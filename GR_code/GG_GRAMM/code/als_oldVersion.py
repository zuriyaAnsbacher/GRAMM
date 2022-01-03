from GR_code.GG_GRAMM.code.aux_functions_oldVersion import equal_al

# old version

class Als(list):
    def __init__(self):
        super().__init__()

    def __contains__(self, value):
        if value == "":
            return True
        lst1 = [equal_al(ele, value) for ele in self]
        return True if any(lst1) else False

    def __eq__(self, other):
        if len(self) == len(self):  # TODO: should be len(self) == len(other)
            lst1 = [ele in other for ele in self]
            lst2 = [ele in self for ele in other]
            return True if all(lst1 + lst2) else False
        return False

    def __ne__(self, other):
        return not (self == other)

    def __add__(self, other):  # TODO: why use this implementation and not just 'append' from 'other' to 'self'? or - why it need to be overrided?
        new = Als()
        for ele in self:
            new.append(ele)
        for ele in other:
            new.append(ele)
        return new

    def is_empty_a(self):
        for ele in self:
            if ele:
                return False
        return True

    def copy_a(self):
        copy_l = Als()
        for ele in self:
            copy_l.append(ele)
        return copy_l

    def sub_lst(self, other):
        lst1 = [ele in other for ele in self]
        return True if all(lst1) else False

    def index_a(self, al):
        for i in range(len(self)):
            if equal_al(self[i], al):
                return i
        return -1

    def remove_a(self, value):
        if self.index_a(value) != -1:
            del self[self.index_a(value)]

    # def merge(self, other):
    #     if not any(self):  # TODO: what is this condition?
    #         return self
    #     if len(self) == 0 or len(other) == 0:
    #         lst = other + self
    #         return lst
    #     if self[0] not in other:
    #         other.append(self[0])
    #     lst_c = other.copy()
    #     del lst_c[lst_c.index_a(self[0])]
    #     if self[1] not in lst_c:
    #         lst_c.append(self[1])
    #     lst_c.append(self[0])
    #     return lst_c

    def merge(self, other):
        # if not any(self):  # TODO: what is this condition?
        #     return self
        if len(self) == 0 or len(other) == 0:
            lst = other + self
            return lst
        if other[0] not in self:
            self.append(other[0])
        lst_c = self.copy_a()
        del lst_c[lst_c.index_a(other[0])]
        if other[1] not in lst_c:
            lst_c.append(other[1])
        lst_c.append(other[0])
        return lst_c

    def match_par(self, par):  # TODO: Im not sure what this function does. can match be only 0/1 or also 2? if it could be 2, it strange that m_par contain 1 index
        match = 0
        m_par = None
        for ele in par:
            if ele in self:
                match += 1
                m_par = par.index_a(ele)
        return match, m_par

