from GR_code.GG_GRAMM.code.aux_functions import equal_al


class Als(list):
    """
    This class illustrates a list of alleles.
    It's a different class because the alleles are required difference treatment:
    Two alleles considered as equal if they are not contradictory
    For example, 02 and 02:01 are equal. 03:01 and 03:02 are not.
    It causes to different implementations of some of the functions
    For example:
        "contain": checking if allele exists in alleles list (e.g. 02:01 exists in [03, 02])
        "sublist": checking if one list is a sublist of another
        and more
    In additional, there are more methods we want to add.
    Therefore, we build a class that inherited from "List", but adds and overrides some of the function
    """
    def __init__(self):
        super().__init__()

    def __contains__(self, value):
        """
        check if self contains value (override "in")
        :param value: value
        :return: true if contain, false otherwise
        """
        if value == "":  # empty value considers as contained in each Als list
            return True
        lst1 = [equal_al(ele, value) for ele in self]  # check equality with every item in Als
        return True if any(lst1) else False  # return True if equal to one at least

    def __eq__(self, other):
        """
        check equality (override "==")
        :param other: other Als list
        :return: true if equal, false otherwise
        """
        if len(self) == len(other):  # a necessary condition for equality
            lst1 = [ele in other for ele in self]
            lst2 = [ele in self for ele in other]
            return True if all(lst1 + lst2) else False  # all item in self are in other, and vice versa
        return False

    def __ne__(self, other):
        """
        override "!="
        :param other: other Als
        :return: true if not equal, false else
        """
        return not (self == other)

    def __add__(self, other):
        """
        override "+". need to be override, cause otherwise created a list, not an Als
        :param other: other Als
        :return: new Als, with elements from self and other
        """
        new = Als()
        for ele in self:
            new.append(ele)
        for ele in other:
            new.append(ele)
        return new

    def is_empty_a(self):  # (exists in "add_data_by_children" in lines 45, 47)
        """
        check if the Als is empty.
        this function is actually unnecessary, cause we can use: if len(obj) == 0, or: if not obj, to check if empty.
        meanwhile I leave it for emphasize that this is an Als and not a list
        # TODO: maybe it is necessary, for consider also this: ["", ""] as empty
        :return: true if empty, false else
        """
        for ele in self:
            if ele:
                return False
        return True

    def copy_a(self):
        """
        create a copy to an Als.
        use this function (and not use regular copy()) because we want that an Als() will be created, not a list
        :return: copied Als list
        """
        copy_l = Als()
        for ele in self:
            copy_l.append(ele)
        return copy_l

    def sub_lst(self, other):
        """
        check if self is a sub list of other
        :param other: other Als list
        :return: true if self is sub list of other
        """
        lst1 = [ele in other for ele in self]
        return True if all(lst1) else False  # return True if all the items in self are in other

    def index_a(self, value):
        """
        return the item index, if exists. if not exist: return -1
        :param value: value we want its index
        :return: the value index if exist. otherwise, -1
        """
        for i in range(len(self)):
            if equal_al(self[i], value):
                return i
        return -1

    def remove_a(self, value):
        """
        check if a value exists in als list, and if exist: remove from the list.
        pay attention that it does not return error if the value is not exist, unlike remove() in list
        :param value: value to remove
        """
        if self.index_a(value) != -1:
            self.remove(value)

    def merge(self, other):
        """
        merge 'other' to 'self', with no repetitions.
        for example, self = [01:02, 04], other = [01, 03:01], so return [01:02, 04, 03:01]
        >> Note: we do not want repetitions **between** the two Als, but within each one - it could be.
        for example, self = [02, 03], other = [01, 01] so return [02, 03, 01, 01]
        (because the first 01 and the second 01 are different: they were inherited from different parents)
        :param other: Als
        :return: merged Als
        """
        if len(self) == 0 or len(other) == 0:  # one empty, at least
            lst = self + other
            return lst

        # 'other' has two identical values, so we want to add they both (see explanation in the function) comment
        if other[0] == other[1] and other[0]:  # second condition for check that the values in 'other' are not empty
            lst = self + other
            return lst

        lst = self.copy_a()  # create a copy for not changing 'self'
        for item in other:
            if item not in lst and item:  # second condition for not add empty value ("") to the merged lst
                lst.append(item)
        return lst

    # TODO: Im not sure what this function does. can match be only 0/1 or also 2? if it could be 2, it strange that m_par contain 1 index.
    #  check it when I work on the functions that uses this.
    def match_par(self, par):
        match = 0
        m_par = None
        for ele in par:
            if ele in self:
                match += 1
                m_par = par.index_a(ele)
        return match, m_par

    def intersection(self, other):
        intersection_count = 0
        idx_intersection_in_other = None
        for allele in other:
            if allele in self:
                intersection_count += 1
                idx_intersection_in_other = other.index_a(allele)  # if 2 intersections, idx will be the second
        return intersection_count, idx_intersection_in_other


