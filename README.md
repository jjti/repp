# Repp

Biologists profit when they can re-use DNA during plasmid design: it enables cheaper designs and faster builds. But parsing through all re-usable DNA is completely infeasible. For example, there are over 75,000 plasmids in Addgene -- the likelihood of knowing the best combination and ordering of sub-sequences from Addgene for a given plasmid design is low.

Repp is made to enable such plasmid design. It turns arbitrary plasmid sequences into designs using the least expensive design with both existing DNA fragments (PCR) and newly synthesized DNA fragments. Repp was written to minimize plasmid design costs and simply the design process.

## Installation

Download links are available at SourceForge: [https://sourceforge.net/projects/repp/files/](https://sourceforge.net/projects/repp/files/)

### MacOS/Linux

```bash
wget -O repp_src_0.1.0.tar.gz 'https://sourceforge.net/projects/repp/files/repp_src_0.1.0.tar.gz/download'
tar xzf repp_src_0.1.0.tar.gz
cd repp_src_0.1.0
make install
```

### Windows

1. Download the most recent `repp_windows.*.zip` file
2. Unzip
3. Run `repp_install.exe`

## Design

See `repp --help` and the [application documentation](https://jjtimmons.github.io/repp/) for more depth and examples.

### Sequence-based design

To design a plasmid sequence in a local FASTA file called `desired_plasmid.fa` using Addgene and a local parts database `parts_library.fa`:

```bash
repp make seq -in ./desired_plasmid.fa --addgene --dbs parts_library.fa
```

To construct a plasmid in `desired_plasmid.fa`, with less consideration for sequence exactness and using all sequences in Addgene, DNASU, and iGEM:

```bash
repp make seq -in ./desired_plasmid.fa --addgene --dnasu --igem --identity 94
```
