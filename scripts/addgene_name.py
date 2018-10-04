import os

os.chdir(os.path.join("assets", "addgene", "repo"))


"""
using name scheme here: https://ncbi.github.io/cxx-toolkit/pages/ch_demo#ch_demo.T5
"""
seen_ids = set()
count = 0

for f in sorted(os.listdir("."))[:100]:
    if ".DS" in f:
        continue

    # get lines
    file_lines = []
    with open(f, "r") as orig_file:
        file_lines = orig_file.readlines()

    duplicate_entry = False

    # rename each
    for name_line_index in range(0, len(file_lines), 2):
        line = file_lines[name_line_index]
        line = line.replace(">", "")

        # is it circular? if yes, we're going to double the sequence
        circular = "circular" in line

        # addgene id
        id = f.split("_")[0]

        if not circular:
            id += "." + str(name_line_index / 2 + 1)
        else:
            # brittle, but am storing whether it was circular in the FASTA ID
            id += "(circular)"

        id = id.replace("\n", "")

        id_line = ">gnl|addgene|" + id + "\n"

        print id_line

        if id_line in seen_ids:
            duplicate_entry = True
            break
        else:
            seen_ids.add(id_line)

        # update its name
        file_lines[name_line_index] = id_line

        # double seq if circular to match sequence across zero-index
        if circular:
            file_lines[name_line_index + 1] = (
                file_lines[name_line_index + 1] + file_lines[name_line_index + 1]
            ).replace("\n", "")

    # write it back
    if not duplicate_entry:
        with open(f, "w") as new_file:
            new_file.writelines(file_lines)
