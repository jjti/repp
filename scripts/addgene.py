
import concurrent.futures
import os
import string
from concurrent.futures import ThreadPoolExecutor, wait

import requests
from lxml import html

INDEX_LIMIT = 150000

# output test directory for addgene
os.chdir(os.path.join("assets", "addgene"))

# valid file characters
valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)

# find the saved files and Addgene pages already seen
max_page = max([int(f.split("_")[0]) for f in os.listdir(".") if "_" in f])

page_indicies = [i for i in range(1, INDEX_LIMIT) if i > max_page]


def clean_name(addgene_name):
    """
    convert an addgene name to one that's file name compliant
    """
    return "".join([c for c in addgene_name if c in valid_chars])


"""
Scrape the AddGene database into the local test/addgene directory for development/testing

AddGene pages are all indexed by the vector's id like: https://www.addgene.org/8389/sequences/

create FASTA files for each in test/data/addgene

if Addgene has full sequence information (prefered):
    save it as <vector_name>_circular in the fasta file
else:
    save it as <vector_name>_linear in the fasta file

save the file as <vector_id>_<vector_name>.fa
"""


def parse(page_index):
    """
    get a page from addgene and save its fragments/vectors if it's a defined id
    """
    page = requests.get("https://www.addgene.org/" + str(page_index) + "/sequences/")
    tree = html.fromstring(page.content)

    # XPath to the title
    name = tree.xpath('//*[@id="page-body"]/main/div/h1/span/text()')
    if not name:
        # many ids are unused/not on the site
        return
    name = clean_name(name[0])

    # XPath to a seq block for full sequence information
    full_seq_blocks = tree.xpath(
        '//*[@id="addgene-full"]/ul/li/div/div/div[2]/textarea/text()'
    )

    # create the file, we'll delete if afterwards if there was only
    # partial DNA available (had N's in it)
    file_path = str(page_index) + "_" + name + ".fa"
    file_written = False
    with open(file_path, "w") as f:
        if full_seq_blocks:
            # full seq information available, write that
            seq = full_seq_blocks[0].split("\n")[1:]
            seq = "".join(seq)
            seq = "".join(seq.split())
            f.write(">" + name + "_circular" + "\n" + seq)
            file_written = True
        else:
            # partial seq information available, write that
            partial_seq_blocks = tree.xpath(
                '//*[@id="addgene-partial"]/ul/li/div/div/div[2]/textarea/text()'
            )

            count = 1
            for partial_seq in partial_seq_blocks:
                seq = partial_seq.split("\n")[1:]
                seq = "".join(seq)
                seq = "".join(seq.split())

                # don't save files with unknown basepairs
                if "N" not in seq:
                    f.write(">" + name + "_linear_" + str(count) + "\n" + seq + "\n")
                    count += 1
                    file_written = True

    # was all partial data
    if not file_written:
        os.remove(file_path)


futures = []
with ThreadPoolExecutor(max_workers=7) as executor:
    """
    try and scrape with 7-threads
    """
    for page in page_indicies:
        futures.append(executor.submit(parse, page))
wait(futures)

# recalculate the seen indicies and save to the fs so we don't scan again in future
seen_indices = set([int(f.split("_")[0]) for f in os.listdir(".") if "_" in f])
empty_indices = [i for i in range(1, INDEX_LIMIT) if i not in seen_indices]
with open(".addgeneignore", "w") as f:
    for ind in empty_indices:
        f.write(str(ind) + "\n")

print "max page index: ", max(seen_indices)
