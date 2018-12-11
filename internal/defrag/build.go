package defrag

import (
	"sort"
)

// build builds up circular assemblies with less fragments than the build limit
//
// maxNodes is the maximum number of nodes in a single assembly
//
// It is created by traversing a DAG in reverse order:
// 	foreach this.node (sorted in reverse order):
// 	  foreach that.node that this node overlaps with + synth_count:
//	 	foreach assembly on that.node:
//    	    add this.node to the assembly to create a new assembly, store on this.node
func build(nodes []node, maxNodes int) (assemblies []assembly) {
	// number of additional nodes try synthesizing to, in addition to those that
	// already have enough homology for overlap without any modifications for each node
	// synth_count = math.max(5, 0.05 * len(nodes)); 5 of 5%, whichever is greater
	synthCount := 5

	// sort with increasing start index
	sort.Slice(nodes, func(i, j int) bool {
		return nodes[i].start < nodes[j].start
	})

	// create a starting assembly on each node including just itself
	for i, n := range nodes {
		nodes[i].assemblies = []assembly{
			assembly{
				nodes:  []node{n},   // just self
				cost:   n.costTo(n), // just PCR,
				synths: 0,           // no synthetic nodes at start
			},
		}
	}

	// for every node in the list of increasing start index nodes
	for i, n := range nodes {
		// for every overlapping fragment + synthCount more
		for _, j := range n.reach(nodes, i, synthCount) {
			// for every assembly on the reaching fragment
			for _, a := range nodes[i].assemblies {
				// see if we can create a new assembly with this node included
				if newAss, created, complete := a.add(nodes[j], maxNodes); created {
					if complete {
						// we've completed a circlular plasmid assembly
						// it has wrapped back onto itself
						assemblies = append(assemblies, newAss)
						// TODO: check if we can break here
					} else {
						// add to this node's list of assemblies
						nodes[j].assemblies = append(nodes[j].assemblies, newAss)
					}
				}
			}
		}
	}

	return
}
