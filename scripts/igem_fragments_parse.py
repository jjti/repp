import os
import re

# output directory for igem
os.chdir(os.path.join("assets", "igem"))


def clean_file():
    with open("xml_parts.xml", "rb") as part_file:
        parsed = part_file.read().decode("utf-8", errors="ignore")

        with open("xml_parts.parsed.xml", "w") as part_file:
            part_file.write(parsed)


# clean_file()

# split up into separate parts. each a tuple with (BBa_name, seq, circular: bool, year: int)
def get_parts(log_composites=True):
    # create element tree object
    parts_string = open("xml_parts.parsed.xml", "r").read()
    seen_names = set()
    composites = []
    parts = []  # list of name, seq, circular tuples
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

        if name and type and seq and len(seq) > 10 and name not in seen_names and year:
            is_circular = "backbone" in type or "plasmid" in type
            parts.append((name, seq, is_circular, year, type))
            seen_names.add(name)

        if log_composites and "composite" in type:
            composites.append(name)

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


write_parts(get_parts())
