<img src="https://user-images.githubusercontent.com/13923102/59981196-28186600-95ce-11e9-81c8-1ebd8e239499.png" width="650" margin-bottom="10px" />

`Repp` is a software application for automating plasmid design. It takes a target plasmid and finds the least expensive combination of fragments from user and public repositories to create it via Gibson Assembly.

Biologists profit when they can re-use DNA during plasmid design: it enables cheaper designs and faster builds. But parsing through all re-usable DNA is completely infeasible. For example, there are over 75,000 plasmids in Addgene -- the likelihood of knowing the best combination and ordering of sub-sequences from Addgene for a given plasmid design is low.

`Repp` enables such plasmid design. It turns plasmid specifications into designs using the least expensive design with both existing DNA fragments (PCR) and newly synthesized DNA fragments. Plasmids are specifiable using their target sequence, features, or sub-fragments.

## Installation

Download links are available at SourceForge: [https://sourceforge.net/projects/repplasmid/files/](https://sourceforge.net/projects/repplasmid/files/)

### MacOS/Linux

```bash
wget -O repp_src_0.1.0.tar.gz 'https://sourceforge.net/projects/repplasmid/files/repp_src_0.1.0.tar.gz/download'
tar xzf repp_src_0.1.0.tar.gz
cd repp_src_0.1.0
make install
```

### Windows

1. Download the most recent `repp_windows.*.zip` from [SourceForge](https://sourceforge.net/projects/repplasmid/files/)
2. Unzip
3. Run `repp_install.exe`

## Documentation

See [the docs](https://jjtimmons.github.io/repp/) or use `--help` on any command.

## Plasmid Design

### Sequence

To design a plasmid based on its expected sequence save it to a FASTA or Genbank file. For example:

```
>2ndVal_mScarlet-I
CAACCTTACCAGAGGGCGCCCCAGCTGGCAATTCCGACGTCTAAGAAACCATTATTATCA...
```

Then call `repp make sequence` to design it. The following example uses Addgene and a local BLAST database `parts_library.fa` as fragment sources:

```bash
repp make sequence --in "./2ndVal_mScarlet-I.fa" --addgene --dbs "parts_library.fa"
```

### Databases

`Repp` includes three embedded databases from large public repositories: [Addgene](https://www.addgene.org/), [iGEM](http://parts.igem.org/Main_Page), and [DNASU](https://dnasu.org/DNASU/Home.do). Each embedded database and its file path after installation are as follows:

- Addgene, `--addgene`, `~/.repp/addgene`
- DNASU, `--dnasu`, `~/.repp/dnasu`
- iGEM, `--igem`, `~/.repp/igem`

Users can also use their or their lab's fragment databases through the `--dbs` as a list of comma-separated fragment [BLAST databases](https://www.ncbi.nlm.nih.gov/books/NBK279688/). An example of a plasmid design using Addgene, DNASU, and multiple user-defined BLAST repositories is below:

```bash
repp make sequence --in "./2ndVal_mScarlet-I.fa" --addgene --dnasu --dbs "proteins.fa,backbones.fa"
```

### Configuration

The default settings file used by `Repp` is in `~/.repp/config.yaml`. The maximum number of fragments in an assembly, the minimum overlap between adjacent fragments, and cost curves for synthesis are all defined there. Editing this file directly will change the default values used during plasmid designs. For more details, see [configuration](https://jjtimmons.github.io/repp/configuration).

To overwrite some `Repp` settings on a per-design basis, create another YAML file:

```yaml
# custom_settings.yaml
fragments-min-junction-length: 25
synthetic-fragment-cost:
  1800: # max length in bp
    fixed: false
    cost: 0.07 # per bp cost
```

And reference it during plasmid design:

```bash
repp make sequence --in "./2ndVal_mScarlet-I.fa" --addgene --settings "./custom_settings.yaml"
```

### Backbones and Enzymes

The plasmid sequence in the input file is designed as a circular plasmid by default. In other words, Repp assumes that the sequence includes an insert sequence as well as a backbone. To use the sequence in the input file as an insert sequence but another fragment as a backbone, use the `--backbone` and `--enzyme` command in combination. This will lookup `--backbone` in the fragment databases and digest it with the enzyme selected through the `--enzyme` flag. The linearized backbone will be concatenated to the insert sequence. For example, to insert a `GFP_CDS` sequence into iGEM's `pSB1A3` backbone after linearizing it with `PstI`:

```bash
repp make sequence --in "./GFP_CDS.fa" --addgene --igem --backbone pSB1A3 --enzyme PstI
```

### Output

`Repp` saves plasmid designs to JSON files at the path specified through the `--out` flag. Below is an abbreviated example of plasmid design output:

```json
{
  "target": "2ndVal_mScarlet-I",
  "seq": "CAACCTTACCAGAGGGCGCCCCAG...",
  "time": "2019/06/24 11:51:39",
  "solutions": [
    {
      "count": 2,
      "cost": 236.65,
      "fragments": [
        {
          "type": "pcr",
          "cost": 94.67,
          "url": "https://www.addgene.org/103998/",
          "seq": "ACAAATAAATGTCCAGACCTGCAG...",
          "pcrSeq": "ACAAATAAATGTCCAGACCTGCAG...",
          "primers": [
            {
              "seq": "ACAAATAAATGTCCAGACCTGCA",
              "strand": true,
              "penalty": 4.406456,
              "pairPenalty": 18.046225,
              "tm": 58.594,
              "gc": 39.13
            },
            {
              "seq": "CATATGTATATCTCCTTCTTAAATCT",
              "strand": false,
              "penalty": 13.639769,
              "pairPenalty": 18.046225,
              "tm": 52.36,
              "gc": 26.923
            }
          ]
        },
        {
          "id": "103998-103998-synthesis-1",
          "type": "synthetic",
          "cost": 129,
          "seq": "AGGAGATATACATATGGTGAGTAA..."
        }
      ]
    }
  ]
}
```
