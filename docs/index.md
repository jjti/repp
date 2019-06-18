# Plasmid Defragger

Biologists benefit when they re-use DNA during plasmid design: it is makes plasmid designs cheaper and faster to build. But parsing through all combinations of re-usable DNA is completely infeasible. For example, there are over 75,000 plasmids in Addgene -- the likelihood of knowing exactly which combination and ordering of sub-sequences to order from Addgene for a new plasmid design is low.

Plasmid Defragger is made to enable such plasmid design. It turns arbitrary plasmid sequences into designs using the least expensive design with both existing DNA fragments (PCR) and newly synthesized DNA fragments. Plasmid Defragger was written to minimize plasmid design costs and simply the design process.

## Installation

Download links are available at SourceForge: [https://sourceforge.net/projects/plasmid-defragger/files/](https://sourceforge.net/projects/plasmid-defragger/files/)

### Windows

1. Download the most recent `plasmid-defragger_windows.*.zip` file
2. Unzip
3. Click and run the `plade_install.exe` installer

### MacOS/Linux

```bash
wget -O plasmid_defragger_src.tar.gz 'https://sourceforge.net/projects/plasmid-defragger/files/plasmid_defragger_src_0.1.0.tar.gz/download'
tar xzf plasmid_defragger_src.tar.gz
cd plasmid-defragger-code
make install
```

## Design

See `plade --help` and the [application documentation](https://jjtimmons.github.io/plade/) for more depth and examples.

### Sequence-based design

To design a plasmid sequence in a local FASTA file called `desired_plasmid.fa` using Addgene and a local parts database `parts_library.fa`:

```bash
plade make seq -in ./desired_plasmid.fa --addgene --dbs parts_library.fa
```

To construct a plasmid in `desired_plasmid.fa`, with less consideration for sequence exactness (the default is 98% percentage identity) and while using all available deposits in Addgene, DNASU, and iGEM, use:

```bash
plade make seq -in ./desired_plasmid.fa --addgene --dnasu --igem --identity 94
```

### SEE ALSO

- [plade annotate](plade_annotate.md) - Annotate a circular sequence using features (comma-separated)
- [plade delete](plade_delete.md) - Delete a feature or enzyme
- [plade find](plade_find.md) - Find features or enzymes
- [plade make](plade_make.md) - make a plasmid from its fragments, features or sequence
- [plade set](plade_set.md) - Set a feature or enzyme
