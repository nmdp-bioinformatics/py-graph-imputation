from itertools import combinations
import logging

from typing import Dict


class CypherQuery(object):
    """
    classdocs
    """

    def __init__(self, loc_map, verbose: bool = False):
        """
        Constructor
        """
        self.verbose = verbose
        self.logger = logging.getLogger("Logger." + __name__)
        # haplo = "A*02:01~C*03:03~B*15:01~DRB1*11:01~DQB1*03:01"
        self.loc_map = loc_map
        for locus, val in self.loc_map.items():
            self.loc_map[locus] = str(val)

    def getLocus(self, allele):
        """Get the locus back from any given allele"""
        loc_allele = allele.split("*")
        return self.loc_map[loc_allele[0]]

    """def getMatch(self, typing):
        loci = sorted([self.getLocus(a) for a in typing[0].split("~")])
        loci_lc = "".join([allele.lower() for allele in loci])
        loci_uc = loci_lc.upper()
        match = "MATCH ("+loci_lc+":"+loci_uc+")"
        return match

    def makeWhereClause(self, haplo):
        # ** This is only used by plan B ** #
        alleles = haplo.split("~")
        where_clause = "WHERE"
        numcp = 5 - len(alleles)
        for i in range(0, numcp):
            c_ind = i + 1
            where_clause = where_clause + " c" + str(c_ind) + ".CP<20000.0 AND"
        where_clause = where_clause + "max(abcqr.frequency) > ep"
        return where_clause

    def haplosIn(self, haplos):
        # Takes in ambigous list of haplotypes and creates WHERE clause#
        loci = sorted([self.getLocus(a) for a in haplos[0].split("~")])
        haplo_quoted = ",".join(["\"" + haplo + "\"" for haplo in haplos])
        loci_lc = "".join([allele.lower() for allele in loci])
        where_clause = "WHERE " + loci_lc +  ".name IN [ " + \
            haplo_quoted + " ]"
        return where_clause

    def getPath(self, typing):

        loci_t = sorted([self.getLocus(x) for x in typing[0].split("~")])
        missing = list()
        for locus in self.loci_a:
            if locus not in loci_t:
                missing.append(locus)

        nodes = list()
        for missing_locus in missing:
            new_loci = loci_t
            new_loci.append(missing_locus)
            loci_lc = "".join(sorted([allele.lower() for allele in new_loci]))
            loci_up = "".join(sorted(new_loci))
            nodeName = "("+loci_lc+":"+loci_up+")"
            #
            if len(loci_up) == len(self.loci_a):
                nodes.append(nodeName)
        return "-->"+nodes[0]

    def allHaplos(self, typing):
        uc = "".join(sorted(typing[0].split("~")))
        lc = uc.lower()
        return ("MATCH ("+lc+":"+uc+") RETURN "+lc +".name", lc)

    def buildQuery(self, typing):

        alleles = typing[0].split("~")
        LenAlleles = len(alleles)
        # returns all of the options
        if(len(alleles[0]) == 1):
            return self.allHaplos(typing)
        else:
        # If there are missing loci...
            if LenAlleles < len(self.loci_a):
                match = self.getMatch(typing)
                path = self.getPath(typing)
                hapsin = self.haplosIn(typing)
                return(match + path + " " + hapsin + " " + self.return_freq)
            else:
                match = self.getMatch(typing)
                hapsin = self.haplosIn(typing)
                return(match + " " + hapsin + " " + self.return_freq)"""
