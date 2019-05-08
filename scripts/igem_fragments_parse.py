import os
import re
import statistics

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(os.path.join(FILE_DIR, "..", "assets", "igem"))


def clean_file():
    with open("xml_parts.xml", "rb") as part_file:
        parsed = part_file.read().decode("utf-8", errors="ignore")

        with open("xml_parts.parsed.xml", "w") as part_file:
            part_file.write(parsed)


# clean_file()

# split up into separate parts. each a tuple with (BBa_name, seq, circular: bool, year: int)
def get_parts(log_composites=True, log_year_count=False):
    # create element tree object
    parts_string = open("xml_parts.parsed.xml", "r").read()
    seen_names = set()
    composites = []
    parts = []  # list of name, seq, circular tuples
    year_count = {year: 0 for year in range(2003, 2020)}
    year_seq_length = {year: [] for year in range(2003, 2020)}
    for row in parts_string.split("<row>"):
        if '"discontinued">1<' in row or '"sample_status">Not in stock' in row:
            continue

        name = ""
        type = ""
        seq = ""
        year = ""

        # get the name match
        name_match = re.search(r"<field name=\"part_name\">(\w*)</field>", row)
        if name_match:
            name = name_match[1]

        # get the part type
        type_match = re.search(r"<field name=\"part_type\">(\w*)</field>", row)
        if type_match:
            type = type_match[1].lower()

        # get the seq
        seq_match = re.search(r"<field name=\"sequence\">(\w*)</field>", row)
        if seq_match:
            seq = seq_match[1]

        # get the creation year
        year_match = re.search(r"<field name=\"creation_date\">(\d*)-.*</field>", row)
        if year_match:
            year = year_match[1]

        if log_year_count and year and seq:
            year_count[int(year)] += 1
            year_seq_length[int(year)].append(len(seq))

        if name and type and seq and len(seq) > 10 and name not in seen_names and year:
            is_circular = "backbone" in type or "plasmid" in type
            parts.append((name, seq, is_circular, year, type))
            seen_names.add(name)

        if log_composites and "composite" in type:
            composites.append(name)
    
    if log_year_count:
        for year, count in year_count.items():
            if count:
                print(year, count)
        for year, lengths in year_seq_length.items():
            if lengths:
                print(year, statistics.mean(lengths))

    # print(" ".join(composites))
    return parts


def write_parts(parts):
    os.chdir(os.path.join("db"))

    # write to the local filesystem
    with open("igem", "w") as out_file:
        for (name, seq, circular, year, type) in parts:
            topo = "circular" if circular else "linear"
            out_seq = seq + seq if circular else seq

            out_file.write(">{} {}-{}-{}\n{}\n".format(name, topo, year, type, out_seq))


# write_parts(get_parts())
get_parts(log_year_count=True)
