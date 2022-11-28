import cython


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef open_ambiguities(list hap, unsigned char loc, tuple split_loc):
    cdef unsigned int k, i, p, j #hap_len, haps_len, splits_len
    cdef Py_ssize_t hap_len, haps_len, splits_len
    cdef list hap_new, hap1
    # cdef np.ndarray[STR, ndim=1] hap_new, hap1
    p = 0
    if len(split_loc) > 1:
        # This opens all allele ambiguities
        hap_len = len(hap[0])
        haps_len = len(hap)
        splits_len = len(split_loc)
        hap_new = [None] * (haps_len * splits_len)
        # hap_new = np.empty(haps_len * splits_len, dtype=np.object)  # produces an empty list of haplotypes
        hap1 = [None] * hap_len
        # hap1 = np.empty(haps_len, dtype=np.object)
        for k in range(haps_len):  # split a given locus in all haps.

            for j in range(hap_len):
                hap1[j] = hap[k][j]

            for i in range(splits_len):
                hap1[loc] = split_loc[i]
                hap_new[p] = hap1[:]
                p += 1
        return hap_new
    return hap

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef create_hap_list(list all_haps, dict optionDict, unsigned int N_Loc):
    cdef unsigned int i, j, count
    cdef list hap_list = []
    cdef list all_hap_split

    for i in range(len(all_haps)):
        all_hap_split = all_haps[i].split('~')
        count = 0
        for j in range(len(all_hap_split)):
            if all_hap_split[j] not in optionDict:
                break
            else:
                count += 1

        if count == N_Loc:
            hap_list.append(all_hap_split)
    return hap_list

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef deepcopy_list(list l):
    cdef list copy_l
    cdef unsigned int i, length
    length = len(l)
    copy_l = [None] * length
    for i in range(length):
        if isinstance(l[i], list):
            copy_l[i] = deepcopy_list(l[i])
        else:
            copy_l[i] = l[i]
    return copy_l
