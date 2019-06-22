"""Parse the iGEM database into a series of DBs with all future years' parts removed."""

if __name__ == "__main__":
    for year in range(2005, 2019):
        with open("igem", "r") as db:
            avoid = set([str(y) for y in range(year, 2019)])

            keep_names = []
            keep_seqs = []

            header = True
            keep = False
            for line in db.readlines():
                if header:
                    h = line.index("-")
                    part_year = line[h + 1 : h + 5]
                    keep = part_year not in avoid
                    if keep:
                        keep_names.append(line.strip())
                    header = False
                elif keep:
                    keep_seqs.append(line.strip())
                    header = True
                else:
                    header = True

            print(year, len(keep_names))

            with open(f"igem-{year}", "w") as new_db:
                for name, seq in zip(keep_names, keep_seqs):
                    new_db.write(f"{name}\n{seq}\n")

