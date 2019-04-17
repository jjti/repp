# rvec

> _rvec_: (of software) ~reduce the fragmentation of (a file)~ _create vectors_ by concatenating parts stored in separate locations on a disk.

Finds pareto optimal solutions for vector construction using existing DNA fragments in local and/or remote repositories.

## Commands

Assemble a list of fragments via Gibson Assembly. Generate primers and synthetic sequences (to avoid complexities).

```bash
rvec fragments -i building_fragments.fa -o vector.json
```

Find fragments to assemble a vector sequence. Use local and/or remote repositories. Rank assemblies by cost and number of fragments.

```bash
rvec vector -i seq.fa -o vectors.json --dbs "local-repo-1.fa local-repo-2.fa" --addgene --igem
```

Find fragments to assemble a vector from a list of features by name. iGEM part names or accession numbers.

```bash
rvec features pSB1C3 FJ172221 -o vector.json
```

## Features

Checks for and avoids primer off-targets (BLAST) with an annealing temp of >40tm

Lets primer3 choose the best primers available within a range of sequence if allowable given neighboring fragments:

```txt
In example below, A and C are too far away from B for annealing.
Need to synthesize from A to B and from B to C.
We let primer3_core pick the best primers available on B (within a 50bp range at the start and the end of B).

------ :A         ------ :C
          ------ :B
```
