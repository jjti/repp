---
layout: default
title: Configuration
nav_order: 1
permalink: /configuration
---

## Configuration

The default settings file used by `Repp` is in `~/.repp/config.yaml`. The maximum number of fragments in an assembly, the minimum overlap between adjacent fragments, and cost curves for synthesis are all defined there. Editing this file directly will change the default values used during plasmid designs.

## Per-Design Configuration

To overwrite some `Repp` settings on a per-design basis, create another YAML file:

```yaml
# custom_settings.yaml
fragments-min-junction-length: 25
synthetic-fragment-cost:
  1800: # max length in bp
    fixed: false
    cost: 0.07 # per bp cost
```

And reference it during plasmid design via `--settings`:

```bash
repp make sequence --in "./2ndVal_mScarlet-I.fa" --addgene --settings "./custom_settings.yaml"
```

## Parameters

| Name                           |  Default | Description                                                                                                                                                                                                                                                                                                                        |
| ------------------------------ | -------: | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| fragments-max-count            |        6 | Maximum number of fragments allowed in a plasmid design. Larger numbers of fragments limit assembly efficiency.                                                                                                                                                                                                                    |
| fragments-min-junction-length  |       15 | Minimum length of overlap between adjacent fragments in bp.                                                                                                                                                                                                                                                                        |
| fragments-max-junction-length  |      120 | Maximum length of overlap between adjacent fragments in bp.                                                                                                                                                                                                                                                                        |
| fragments-max-junction-hairpin |       47 | Maximum annealing temperature allowed in primers and at the ends of synthetic fragments.                                                                                                                                                                                                                                           |
| gibson-assembly-cost­          |    12.98 | The per reaction dollar cost of each Gibon Assembly reaction. Based upon the per reaction cost of NEB’s Gibson Assembly Master Mix.                                                                                                                                                                                                |
| gibson-assembly-time-cost      |        0 | The per reaction cost of human hours for the assembly. Depends on researcher’s value of time and the length required per assembly.                                                                                                                                                                                                 |
| pcr-bp-cost                    |      0.6 | The per bp cost of each primer bp. Used in estimating the final assembly cost of each assembly. Cost is based upon IDT’s primer bp cost for 100nmol of single-stranded DNA as of February 2019.                                                                                                                                    |
| pcr-rxn-cost                   |     0.27 | The per reaction cost of PCR. Estimated using the per reaction cost of ThermoFisher’s Taq DNA Polymerase PCR Buffer (10X).                                                                                                                                                                                                         |
| pcr-time-cost                  |        0 | The per reaction of human time for each PCR reaction. This cost is applied across each assembly. So an \$85 human cost for a PCR assembly include all PCRs necessary for that assembly.                                                                                                                                            |
| pcr-min-length                 |       60 | The minimum number of bp necessary for a fragment to be PCR’ed. Fragment matches less than this length are not considered.                                                                                                                                                                                                         |
| pcr-primer-max-pair-penalty    |       30 | The maximum pair penalty for primers generated via Primer3. The configuration penalty is related to Primer3’s PRIMER*PAIR*\*\_PENALTY score and is used to filter out poor primer combinations with large mismatches in annealing temperature or heterodimers.                                                                     |
| pcr-primer-max-embed-length    |       20 | The maximum length of embedded sequence at the end of a fragment via mutation in a primer.                                                                                                                                                                                                                                         |
| pcr-primer-max-ectopic-tm      |       55 | The maximum tolerable primer annealing temperature against an ectopic binding site. Calculated via the “ntthal” binary in Primer3. 2 PCR products with primers whose ectopic binding tm exceed this value are ignored.                                                                                                             |
| pcr-buffer-length              |       20 | The allowable range in which Plasmid Defragger lets Primer3 optimize primer pairs. Used when a PCR fragments neighbor is synthetic. The synthetic fragment can be expanded to overlap whatever range the PCR fragment winds up spanning, so Primer3 is given a range in which to generate primer pairs, rather than a fixed start. |
| synthetic-min-length           |      125 | The minimum length of a fragment to be considered or synthesized.                                                                                                                                                                                                                                                                  |
| synthetic-max-length           |     3000 | The maximum length of a fragment to be considered for synthesis. Synthetic spans of DNA larger than this are fragmented into smaller synthetic fragments with overlap for one another.                                                                                                                                             |
| synthetic-fragment-cost        | cost-map | A synthesis cost map. Default costs correspond to IDT’s “gBlocks” product as of February 2019.                                                                                                                                                                                                                                     |
| synthetic-plasmid-cost         | cost-map | A synthesis cost map. Default costs correspond to IDT’s “Custom gene synthesis” service as of February 2019.                                                                                                                                                                                                                       |
| addgene-cost                   |       65 | The cost of procuring a plasmid from Addgene.                                                                                                                                                                                                                                                                                      |
| igem-cost                      |        0 | The cost of procuring an iGEM part from iGEM.                                                                                                                                                                                                                                                                                      |
| dnasu-cost                     |       55 | The cost of procuring a plasmid from DNASU.                                                                                                                                                                                                                                                                                        |

### Synthesis Cost Maps

`Repp` supports cost maps for synthetic fragments (`synthetic-fragment-cost`) and synthetic pre-cloned plasmids (`synthetic-plasmid-cost`). The default costs in `Repp` are 1:1 with IDT synthesis prices for synthetic fragments and synthetic genes. To customize the cost curve of synthesis, define a new cost curve as a YAML map. Each key is an integer which is the maximum length fragment/gene to be synthesized at that cost, and the value contains two keys: `fixed` and `cost`. If `fixed` is true, the `cost` of a fragment of that length is used as is. If `fixed` is false, the `cost` of a fragment of that length is per basepair.

For example, in the cost curve below, the cost of the fragment's synthesis is 55.0 up to an including fragments that are 500 bp. For fragments greater than 500 bp and no more than 1800 bp, the cost of synthesis is 0.07 \* _fragment bp length_. No fragments can be produced that are greater than 1800 bp.

```yaml
synthetic-fragment-cost:
  500:
    fixed: true
    cost: 35.0
  1800:
    fixed: false
    cost: 0.07
```

### SEE ALSO

- [repp](/) - Repp
