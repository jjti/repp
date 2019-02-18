import json
import os
import re

os.chdir(os.path.join("assets", "dnasu"))

# split up into separate parts. each a tuple with (id: string, seq: string)
def get_parts():
    id_to_ref = {}
    ref_to_seq = {}

    # col1 is id, col33 is seqrefid
    clone_data = open("DNASU-CloneData.csv", "r")
    for p in clone_data.readlines()[1:]:
        cols = p.split(",")

        p_id = cols[0]
        bb_ref = cols[17]  # bb ref
        in_ref = cols[18]  # insert ref

        id_to_ref[p_id] = [bb_ref, in_ref]
        ref_to_seq[in_ref] = ""
    clone_data.close()

    seq_data = open("DNASU-InsertSeq.csv", "r")
    for p in seq_data.readlines()[1:]:
        ref, _, seq = p.split(",")
        if ref in ref_to_seq:
            ref_to_seq[ref] += seq
    seq_data.close()

    parts = []
    for p_id, refs in id_to_ref.items():
        bb_ref, in_ref = refs

        if bb_ref in ref_to_seq and ref_to_seq[bb_ref]:
            parts.append((p_id + ".bb", ref_to_seq[bb_ref]))

        if in_ref in ref_to_seq and ref_to_seq[in_ref]:
            parts.append((p_id + ".in", ref_to_seq[in_ref]))
    return parts


def write_parts(parts):
    os.chdir(os.path.join("db"))

    # write to the local filesystem
    with open("dnasu", "w") as out_file:
        for (name, seq) in parts:
            out_file.write(">gnl|dnasu|{} circular\n{}".format(name, seq))


write_parts(get_parts())
