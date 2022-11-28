class CypherQueryPlanB(object):
    """
    classdocs
    """

    def __init__(self, loci_map):  # , indexes):
        """
        Constructor
        """

        self.loc_map = loci_map
        # self.index_dict =  indexes

    def indexToLocus(self, indexes):
        loci = []

        # Sort if the indexes are not sorted
        indexes = sorted(indexes)

        for index in indexes:
            # Index to locus
            loci.append(str(index))

        return loci

    def findTypeFromIndexes(self, indexes):
        full_form_loci = self.indexToLocus(indexes)
        return self.findFullFormType(full_form_loci)

    def findFullFormType(self, full_form_loci):
        type = "".join([allele for allele in full_form_loci])
        return type

    def getLocus(self, allele):
        """Get the locus back from any given allele"""
        loc_allele = allele.split("*")
        return self.loc_map[loc_allele[0]]

    def findLinkType(self, typing):
        loci = sorted([self.getLocus(a) for a in typing[0].split("~")])
        loci_lc = "".join([allele for allele in loci])
        return loci_lc
