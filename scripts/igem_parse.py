import os
import re

# output directory for igem
os.chdir(os.path.join("assets", "igem"))


def clean_file():
    with open("xml_parts.xml", "rb") as part_file:
        parsed = part_file.read().decode("utf-8", errors="ignore")

        with open("xml_parts.parsed.xml", "w") as part_file:
            part_file.write(parsed)


clean_file()

# split up into separate parts. each a tuple with (BBa_name, seq, circular: bool)
def get_parts():
    # create element tree object
    parts_string = open("xml_parts.parsed.xml", "r").read()

    seen_names = set()

    parts = []  # list of name, seq, circular tuples
    for row in parts_string.split("<row>"):
        if '"discontinued">1<' in row or '"sample_status">Not in stock' in row:
            continue

        name = ""
        type = ""
        seq = ""

        # get the name match
        name_match = re.search(r"<field name=\"part_name\">(\w*)</field>", row)
        if name_match:
            name = name_match[1]

        # get the part type
        type_match = re.search(r"<field name=\"part_type\">(\w*)</field>", row)
        if type_match:
            type = type_match[1]

        # get the seq
        seq_match = re.search(r"<field name=\"sequence\">(\w*)</field>", row)
        if seq_match:
            seq = seq_match[1]

        if name and type and seq and len(seq) > 10 and name not in seen_names:
            is_circular = "backbone" in type.lower() or "plasmid" in type.lower()
            parts.append((name, seq, is_circular))
            seen_names.add(name)
    return parts


def write_parts(parts):
    os.chdir(os.path.join("db"))

    # write to the local filesystem
    with open("igem", "w") as output_parts:
        for (name, seq, circular) in parts:
            if circular:
                output_parts.write(">gnl|igem|" + name + " circular fwd\n" + seq + "\n")
            else:
                output_parts.write(">gnl|igem|" + name + " linear fwd\n" + seq + "\n")


write_parts(get_parts())
