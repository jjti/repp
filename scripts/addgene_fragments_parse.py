import json
import os
import re
import statistics

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(os.path.join(FILE_DIR, "..", "assets", "addgene"))


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
def get_parts(log_composites=True, from_year=-1, log_year_count=False):
    # create element tree object
    part_json = json.load(open("addgene.json", "r"))
    plasmids = part_json["plasmids"]

    def fix(s):
        """some seqs have multilines, split and join here"""

        return s.replace("\n", "").replace("\r", "")

    parts = []  # list of name, seq, circular tuples
    year_count = {year: 0 for year in range(2004, 2020)}
    year_seq_length = {year: [] for year in range(2003, 2020)}
    full_seq_count = 0
    for p in plasmids:
        part_id = p["id"]
        year = id_to_year(part_id)
        # name = p["name"]

        if log_year_count:
            year_count[year] += 1

        # vectors
        seqs = p["sequences"]
        for s in (
            seqs["public_addgene_full_sequences"] + seqs["public_user_full_sequences"]
        ):
            full_seq_count += 1
            parts.append((part_id, fix(s["sequence"]), True, year))
            year_seq_length[year].append(len(fix(s["sequence"])))
            break

        # linear fragments
        sindex = 1
        for s in (
            seqs["public_addgene_partial_sequences"]
            + seqs["public_user_partial_sequences"]
        ):
            parts.append(
                (f"{part_id}.{sindex}", fix(s["sequence"]), False, year)
            )
            sindex += 1
    
    if log_year_count:
        for year, count in year_count.items():
            print(year, count)
        for year, lengths in year_seq_length.items():
            if lengths:
                print(year, statistics.median(lengths))

    if from_year > 0:
        parts = [p for p in parts if p[3] == from_year]

    # remove ambiguous bp
    parts = [p for p in parts if "N" not in p[1]]
    parts = [p for p in parts if "\t" not in p[1]]

    return parts


def write_parts(parts, to="", only_topo=""):
    output = os.path.join("db", "addgene")
    if to:
        output = to

    # write to the local filesystem
    with open(output, "w") as out_file:
        seen_ids = set()
        for (p_id, seq, circular, year) in parts:
            if p_id in seen_ids or len(seq) < 30 or ">" in seq:
                continue

            topo = "circular" if circular else "linear"
            if only_topo and topo != only_topo:
                continue

            out_seq = seq + seq if circular else seq

            out_file.write(f">{p_id} {topo}-{year}\n{out_seq}\n")

            seen_ids.add(p_id)


if __name__ == "__main__":
    # write_parts(get_parts(from_year=2018), to="2018.addgene.fasta", only_topo="circular", )

    # write_parts(get_parts())

    get_parts(log_year_count=True)
