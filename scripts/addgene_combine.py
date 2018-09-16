import os

os.chdir(os.path.join("assets", "addgene"))

# open a single FASTA for a combined db
with open("addgene.fa", "w") as combined_fasta:
    # move into repo full of files
    os.chdir("repo")

    # read in all the addgene fasta files and combine into one that makeblastdb can act on
    for f in os.listdir("."):
        combined_fasta.write(open(f, "r").read().strip() + "\n")
