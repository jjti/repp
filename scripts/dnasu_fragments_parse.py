import json
import os
import re

os.chdir(os.path.join("assets", "dnasu"))

# split up into separate parts. each a tuple with (id: string, year: int, seq: string)
def get_parts():
    id_to_ref = {}
    id_to_year = {}
    ref_to_seq = {}

    # col1 is id, col4 is year
    with open("DNASU-clonedates.csv") as clone_years:
        for p in clone_years.readlines()[1:]:
            cols = p.split(",")

            p_id = cols[0]
            year = cols[4]  # year
            id_to_year[p_id] = int(year)

    # col1 is id, col33 is seqrefid
    clone_data = open("DNASU-CloneData.csv", "r")
    for p in clone_data.readlines()[1:]:
        cols = p.split(",")

        p_id = cols[0]
        # bb_ref = cols[17]  # bb ref
        in_ref = cols[18]  # insert ref
        # seq_ref = cols[33]  # seq ref

        id_to_ref[p_id] = [in_ref]
        ref_to_seq[in_ref] = ["" for _ in range(0, 20)]
        # ref_to_seq[bb_ref] = ""
    clone_data.close()

    seq_data = open("DNASU-InsertSeq.csv", "r")
    for p in seq_data.readlines()[1:]:
        _, ref, _, order, seq = p.split(",")
        order = int(order)
        if ref in ref_to_seq and not ref_to_seq[ref][order]:
            ref_to_seq[ref][order] = seq
    seq_data.close()

    parts = []
    for p_id, refs in id_to_ref.items():
        in_ref = refs[0]

        if in_ref in ref_to_seq and ref_to_seq[in_ref] and p_id in id_to_year:
            parts.append((p_id, "".join(ref_to_seq[in_ref]), id_to_year[p_id]))
    return parts


def write_parts(parts):
    os.chdir(os.path.join("db"))

    year_to_count = {year: 0 for year in range(1990, 2020)}

    # write to the local filesystem
    with open("dnasu", "w") as out_file:
        for (name, seq, year) in parts:
            if seq:
                year_to_count[year] += 1
                out_file.write(">{} circular-{}\n{}".format(name, year, seq))

    for year, count in year_to_count.items():
        print(year, count)


write_parts(get_parts())
