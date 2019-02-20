import json
import os
import re

os.chdir(os.path.join("assets", "addgene"))


def id_to_year(p_id):
    # 2005	0
    # 2006	11000
    # 2007	12000
    # 2008	17000
    # 2009	20000
    # 2010	22000
    # 2011	25000
    # 2012	31000
    # 2013	41000
    # 2014	50000
    # 2015	56000
    # 2016	72000
    # 2017	85000
    # 2018	102000
    # 2019	114000
    return {
        p_id < 11000: 2005,
        p_id >= 11000 and p_id < 12000: 2006,
        p_id >= 12000 and p_id < 17000: 2007,
        p_id >= 17000 and p_id < 20000: 2008,
        p_id >= 20000 and p_id < 22000: 2009,
        p_id >= 22000 and p_id < 25000: 2010,
        p_id >= 25000 and p_id < 31000: 2011,
        p_id >= 31000 and p_id < 41000: 2012,
        p_id >= 41000 and p_id < 50000: 2013,
        p_id >= 50000 and p_id < 56000: 2014,
        p_id >= 56000 and p_id < 72000: 2015,
        p_id >= 72000 and p_id < 85000: 2016,
        p_id >= 85000 and p_id < 102000: 2017,
        p_id >= 102000 and p_id < 114000: 2018,
        p_id >= 114000: 2019,
    }[1]


# split up into separate parts. each a tuple with (addgene_id, seq, circular: bool, year: int, name)
def get_parts(log_composites=True):
    # create element tree object
    part_json = json.load(open("addgene.json", "r"))
    plasmids = part_json["plasmids"]

    parts = []  # list of name, seq, circular tuples
    full_seq_count = 0
    for p in plasmids:
        part_id = p["id"]
        # name = p["name"]

        # vectors
        seqs = p["sequences"]
        for s in (
            seqs["public_addgene_full_sequences"] + seqs["public_user_full_sequences"]
        ):
            full_seq_count += 1
            parts.append((part_id, s["sequence"], True, id_to_year(part_id)))
            break

        # linear fragments
        sindex = 1
        for s in (
            seqs["public_addgene_partial_sequences"]
            + seqs["public_user_partial_sequences"]
        ):
            parts.append(
                (
                    "{}.{}".format(part_id, sindex),
                    s["sequence"],
                    False,
                    id_to_year(part_id),
                )
            )
            sindex += 1
    print(full_seq_count)
    return parts


def write_parts(parts):
    os.chdir(os.path.join("db"))

    # write to the local filesystem
    with open("addgene", "w") as out_file:
        seen_ids = set()
        for (p_id, seq, circular, year) in parts:
            if p_id in seen_ids or len(seq) < 30 or ">" in seq:
                continue

            topo = "circular" if circular else "linear"
            out_seq = seq + seq if circular else seq

            out_file.write(">{} {}-{}\n{}\n".format(p_id, topo, year, out_seq))

            seen_ids.add(p_id)


write_parts(get_parts())
