"""Parse out the 2018 plasmids from the addgene database, make a new one."""


if __name__ == "__main__":

    with open("addgene", "r") as db:
        kept_names = []
        kept_seqs = []
        header = True
        keep = True
        for line in db.readlines():
            if header:
                keep = (
                    "-2018" not in line and "-2019" not in line and "linear" not in line
                )
                if keep:
                    kept_names.append(line.strip())
                header = False
            elif keep:
                kept_seqs.append(line.strip())
                header = True
            else:
                header = True

        print(len(kept_names))

        with open("addgene.2", "w") as new_db:
            for (name, seq) in zip(kept_names, kept_seqs):
                new_db.write(f"{name}\n{seq}\n")
