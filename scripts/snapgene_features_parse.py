import io
import os

# install dgparse from https://github.com/DeskGen/dgparse
import dgparse


def read_features():
    """read all SnapGene features from all feature collections
    
    Returns:
        dict -- key = feature name, value = feature seq
    """

    os.chdir(os.path.join("..", "assets", "snapgene"))
    feat_files = [os.path.join(d, f) for d, _, fs in os.walk(".") for f in fs]
    feat_files = [
        f
        for f in feat_files
        if not (f.startswith("./.DS") or f.startswith("./features"))
    ]

    # a dict with key = name and value = seq (1 for each feature)
    failures = []
    features = {}
    for file in feat_files:
        # with io.open(file, "r+", encoding="utf8", errors="replace") as f:
        #     data = f.read()
        #     f.seek(0)
        #     f.write(data)
        #     f.truncate()

        with open(file, "rb") as snapgene_file:
            try:
                output = dgparse.snapgene.parse(snapgene_file)
                for feature in output["dnafeatures"]:
                    dna_feature = feature["dnafeature"]

                    name = dna_feature["name"].encode("utf-8").strip()
                    seq = dna_feature["pattern"]["bases"]
                    if len(seq) > 1 and "™" not in name:
                        features[name] = seq
            except:
                # raise
                failures.append(file)
                continue
    print(failures)
    return features


def write_features():
    """get features from snapgene files and write them to a TSV file
    """

    features = read_features()

    with open("features.tsv", "w") as output:
        for name, seq in features.iteritems():
            output.write(name + "\t" + seq + "\n")


write_features()
