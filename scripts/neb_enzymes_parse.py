import os
import re

input_file = os.path.join("assets", "neb", "enzymes.txt")
output_file = os.path.join("assets", "neb", "enzymes.tsv")

enzyme_recogs = []  # (name, seq)

non_char = re.compile("[^a-zA-Z0-9]")

with open(input_file, "r") as enzymes_txt:
    enzymes_string = enzymes_txt.read()

    for line in enzymes_string.split("\n"):
        cols = [c.strip() for c in line.split(":")]

        name = non_char.sub("", cols[0])

        seq = cols[2].split(",")[0]
        seq = non_char.sub("", seq)

        cut_ind = int(non_char.sub("", cols[4]))
        hang_ind = int(cols[3].split(",")[0].strip())

        # ^ for cut index, _ for hang index
        cut_placed = False
        hang_placed = False
        new_seq = ""
        i = 0
        while not cut_placed or not hang_placed or i < len(seq):
            c = "N"
            if i < len(seq):
                c = seq[i]
            if i == cut_ind:
                cut_placed = True
                new_seq += "^"
            if i == hang_ind:
                hang_placed = True
                new_seq += "_"
            new_seq += c
            i += 1

        if name == "SfaNI":
            assert new_seq == "GCATCNNNNN^NNNN_N"
        if name == "MslI":
            assert new_seq == "CAYNN^_NNRTG"
        if name == "BsmAI":
            assert new_seq == "GTCTCN^NNNN_N"

        enzyme_recogs.append((name, new_seq))

with open(output_file, "w") as output:
    for name, seq in enzyme_recogs:
        output.write("{}\t{}\n".format(name, seq))
