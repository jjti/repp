import os

os.chdir(os.path.join("assets", "addgene", "repo"))


"""
using name scheme here: https://ncbi.github.io/cxx-toolkit/pages/ch_demo#ch_demo.T5
"""
seen_names = set()
count = 0

for f in sorted(os.listdir("."))[:10]:
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

        # addgene index
        index = f.split("_")[0]

        new_line = ">gnl|addgene|" + index.replace("\n", "")

        # is it circular? if yes, we're going to double the sequence
        circular = "circular" in line

        # add a secondary index
        if not circular:
            new_line += "." + str(name_line_index / 2 + 1)

        new_line += "\n"

        if new_line in seen_names:
            duplicate_entry = True
            break
        else:
            seen_names.add(new_line)

        # update its name
        file_lines[name_line_index] = new_line

        # double seq if circular to match sequence across zero-index
        if circular:
            file_lines[name_line_index + 1] = (
                file_lines[name_line_index + 1] + file_lines[name_line_index + 1]
            ).replace("\n", "")

    # write it back
    # if not duplicate_entry:
    # with open(f, "w") as new_file:
    #     new_file.writelines(file_lines)

    print "".join(file_lines)
