package defrag

import (
	"math"
	"sort"
)

// build builds up circular assemblies with less fragments than the build limit
//
// maxNodes is the maximum number of frags in a single assembly
// seq is the target sequence we're trying to build up
//
// It is created by traversing a DAG in forward order:
// 	foreach this.Frag (sorted in forward order):
// 	  foreach that.Frag that this Frag overlaps with + synthCount:
//	 	foreach assembly on that.Frag:
//    	    add this.Frag to the assembly to create a new assembly, store on this.Frag
func build(frags []*Frag, maxNodes int, seq string) (assemblies []assembly) {
	// number of additional frags try synthesizing to, in addition to those that
	// already have enough homology for overlap without any modifications for each Frag
	synthCount := int(math.Max(5, 0.05*float64(len(frags)))) // 5 of 5%, whichever is greater

	// sort with increasing start index
	sort.Slice(frags, func(i, j int) bool {
		return frags[i].start < frags[j].start
	})

	// create a starting assembly on each Frag including just itself
	for i, f := range frags {
		// edge case where the Frag spans the entire target vector... 100% match
		// it is the target vector. just return that as the assembly
		if len(f.Seq) >= len(seq) {
			assemblies = append(assemblies, assembly{
				frags:  []*Frag{f}, // just self
				cost:   f.Cost,     // just cost of procurement
				synths: 0,
			})
			break // can't top that
		}

		frags[i].assemblies = []assembly{
			assembly{
				frags:  []*Frag{f},  // just self
				cost:   f.costTo(f), // just PCR,
				synths: 0,           // no synthetic frags at start
			},
		}
	}

	// for every Frag in the list of increasing start index frags
	for i, f := range frags {
		// for every overlapping fragment + synthCount more
		for _, j := range f.reach(frags, i, synthCount) {
			// for every assembly on the reaching fragment
			for _, a := range frags[i].assemblies {
				// see if we can create a new assembly with this Frag included
				if newAssembly, created, complete := a.add(frags[j], maxNodes); created {
					if complete {
						// we've completed a circlular plasmid assembly
						// it has wrapped back onto itself
						assemblies = append(assemblies, newAssembly)
						// TODO: check if we can break here
					} else {
						// add to this Frag's list of assemblies
						frags[j].assemblies = append(frags[j].assemblies, newAssembly)
					}
				}
			}
		}
	}

	return
}
