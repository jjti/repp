package assemble

import (
	"math"
	"sort"
)

// build is for building up countMap
//
// countMap is a map from the number of fragments in an assembly to
// a list of assemblies that have that number of fragments within them
//
// It is created by traversing a DAG in reverse order:
//
// number of additional nodes to consider for each node
// (ie additional number of fragments to try synthing to)
// synth_count = math.max(5, 0.05 * len(nodes))
//
// traverse the nodes
// 	foreach this.node (sorted in reverse order):
// 	  foreach that.node that this node overlaps with + synth_count:
//	 	foreach assembly on that.node:
//    	    add this.node to the assembly to create a new assembly, store on this.node
//
// create a map from the number of fragments in each assembly to a list with assemblies
// 		containing that many assemblies
func build(nodes []node) (assemblies []assembly) {
	// number of nodes to try to synthesize to from each node (plus the natural overlap)
	synthCount := int(math.Max(5.0, 0.05*float64(len(nodes))))

	// sort with increasing start index
	sort.Slice(nodes, func(i, j int) bool {
		return nodes[i].start < nodes[j].start
	})

	// for every node in the list of reverse sorted nodes
	for i, n := range nodes {
		// for every overlapping fragment + synthCount more
		for _, r := range n.reach(nodes, i, synthCount) {
			// for every assembly on the reachable fragment
			for _, a := range r.assemblies {
				// see if we can create a new assembly with this node included
				if newAss, created, complete := a.add(n); created {
					if complete {
						// we've completed a circlular plasmid assembly
						// it has wrapped back onto itself
						assemblies = append(assemblies, newAss)
						// TODO: check if we can break here
					} else {
						// and add to this node's list of assemblies
						n.assemblies = append(n.assemblies, newAss)
					}
				}
			}
		}
	}

	return
}
