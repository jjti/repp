import io
import os
import re
from xml.etree import ElementTree as ET

import requests
from lxml import etree, html

# output directory for igem
os.chdir(os.path.join("assets", "igem", "files"))


def clean_file():
    parsed = ""
    with open("xml_parts.xml", "r") as part_file:
        parsed = unicode(part_file.read(-1), errors="ignore")

    parsed = re.sub(r"<field name=\"sequence_sha1\">.*</field>\n", "", parsed)

    with open("xml_parts.parsed.xml", "w") as part_file:
        part_file.write(parsed)


# clean_file()

# split up into separate parts. each a tuple with (BBa_name, seq, circular: bool)
def get_parts():
    tree = ET.parse('xml_parts.parsed.xml')
    root = tree.getroot()

    print root

    def invalid_row(row):
        return row.find("ok") != 1 or row.find("discontinued") == 1

    parts = []  # list of name, seq, circular tuples
    for row in root.iterfind("row"):
        if invalid_row(row):
            continue

        # relevant rows: deprecated, part_type, part_name, status == "Deleted"?, discontinued == 1, ok == 1, sequence
        name = row.find("part_name").text
        type = row.find("part_type").text
        seq = row.find("sequence").text
        parts.append((name, seq, "backbone" in type.lower()))
    return parts


print get_parts()[:5]


def write_parts(parts):
    # write to the local filesystem
    with open("igem.fa", "w") as output_parts:
        for (name, seq, circular) in parts:
            if circular:
                output_parts.write(">gnl|igem|" + name + " circular fwd\n" + seq + "\n")
            else:
                output_parts.write(">gnl|igem|" + name + " linear fwd\n" + seq + "\n")


def parse_to_db(parts_file):
    combined_output = ""

    for part in parts_file.split(">")[1:]:
        lines = part.split("\n")

        header = lines[0]

        part_id = header.split(" ")[0]

        circular = "backbone" in header.lower()

        seq_lines = [s.strip() for s in lines[1:]]
        seq = "".join(seq_lines)

        if len(seq) < 20:
            continue

        if circular:
            combined_output += ">gnl|igem|" + part_id + " circular fwd\n" + seq + "\n"
        else:
            combined_output += ">gnl|igem|" + part_id + " linear fwd\n" + seq + "\n"

    os.chdir(os.path.join("..", "db"))

    with open("db", "w") as igem_db:
        igem_db.write(combined_output)


# write_parts(get_parts()) # write to the fs, copying the text from the webpage

# with open("igem.fa", "r") as igem_text:
#     parse_to_db(igem_text.read(-1))
