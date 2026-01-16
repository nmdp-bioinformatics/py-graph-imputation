def extract_locuses(haplotype):
    """Get one haplotype and return its locuses"""
    return [
        element[:4] if element[:3].isalpha() else element[0]
        for element in haplotype.split("^")
    ]


def get_3_dominant_locuses(locuses, threshold=3):
    """
    prioritize locuses by the list Sapir gave us:
    1. A
    2. B
    3. DRB1
    4. C
    5. DQB1
    6. DPB1
    7. DRB3/4/5
    8. DPA1 (but if you already have DPB1 but not DQB1 then switch between 8 and 9)
    9. DQA1
    """
    locuses = list(set(locuses))
    dominant_locuses = []
    if "A" in locuses:
        dominant_locuses.append("A")
    if "B" in locuses:
        dominant_locuses.append("B")
    if "DRB1" in locuses:
        dominant_locuses.append("DRB1")
    if len(dominant_locuses) == threshold:
        return dominant_locuses
    if "C" in locuses:
        dominant_locuses.append("C")
    if len(dominant_locuses) == threshold:
        return dominant_locuses
    if "DQB1" in locuses:
        dominant_locuses.append("DQB1")
    if len(dominant_locuses) == threshold:
        return dominant_locuses
    if "DPB1" in locuses:
        dominant_locuses.append("DPB1")
    if len(dominant_locuses) == threshold:
        return dominant_locuses
    if "DRB3" in locuses or "DRB4" in locuses or "DRB5" in locuses:
        dominant_locuses.append("DRB3/4/5")
    if len(dominant_locuses) == threshold:
        return dominant_locuses
    if "DPB1" in locuses and "DQB1" not in locuses:
        if "DQA1" in locuses:
            dominant_locuses.append("DQA1")
        if len(dominant_locuses) == threshold:
            return dominant_locuses
        if "DPA1" in locuses:
            dominant_locuses.append("DPA1")
        if len(dominant_locuses) == threshold:
            return dominant_locuses
    else:
        if "DPA1" in locuses:
            dominant_locuses.append("DPA1")
        if len(dominant_locuses) == threshold:
            return dominant_locuses
        if "DQA1" in locuses:
            dominant_locuses.append("DQA1")
        if len(dominant_locuses) == threshold:
            return dominant_locuses
    return dominant_locuses


def filter_haplotype(haplotype, dominant_locuses):
    """Filter haplotype by the dominant locuses"""
    haplotype = haplotype.split("^")
    return "^".join(
        [
            element
            for element in haplotype
            if (element[:4] in dominant_locuses) or (element[0] in dominant_locuses)
        ]
    )


def split_gl(subject_gl):
    """split the gl string into the top-3 dominant part of the haplotype and the less dominant."""
    subject_locuses = extract_locuses(subject_gl)
    dominant_loc = get_3_dominant_locuses(subject_locuses)
    short_gl = filter_haplotype(subject_gl, dominant_loc)

    extra_locuses = [locus for locus in subject_locuses if locus not in dominant_loc]
    extra_gl = filter_haplotype(subject_gl, extra_locuses)

    return short_gl, extra_gl


def change_donor_file(path):
    # Read the original file and store its contents
    with open(path, "r") as file:
        lines = file.readlines()

    modified_lines = []
    gls = {"subject_id": [], "short_gl": [], "extra_gl": []}
    for name_gl in lines:
        name_gl = name_gl.rstrip()  # remove trailing whitespace
        if "," in name_gl:
            list_gl = name_gl.split(",")
        else:
            list_gl = name_gl.split("%")

        if len(list_gl) > 1:
            id_gl = list_gl[0]
            subject_gl = list_gl[1]

            short_gl, extra_gl = split_gl(subject_gl)
            extra_gl = extra_gl.replace("g", "")
            # gl = gl.replace('N', '')
            extra_gl = extra_gl.replace("L", "")

            gls["subject_id"].append(id_gl)
            gls["short_gl"].append(short_gl)
            gls["extra_gl"].append(extra_gl)
            list_gl[1] = short_gl

        modified_lines.append(",".join(list_gl))

    # Write the modified contents back to the file
    with open(path, "w") as file:
        for line in modified_lines:
            file.write(line + "\n")
    file.close()
    return gls, lines
