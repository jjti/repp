package assemble

import (
	"math"
	"sort"
)

// build is for building up an ASSEMBLIES_LIST
//
// ASSEMBLIES_LIST is a map from the number of fragments in an assembly to
// a list of assemblies that have that number of fragments within them
//
// It is created by traversing a DAG in reverse order:
//
// number of additional nodes to consider for each node
// (ie additional number of fragments to try synthing to)
// synth_count = math.max(5, 0.05 * len(nodes))
//
// traverse the nodes
//
// 	foreach this.node (sorted in reverse order):
// 	  foreach that.node that this node overlaps with + synth_count:
//	 	foreach assembly on that.node:
//    	    add this.node to the assembly to create a new assembly, store on this.node
// after the last node is traversed, add an empty node and synth to synth_count other nodes
//
// // store in a map
// create a map from the number of fragments in each assembly to a list with assemblies
// 		containing that many assemblies
func build(nodes []node) (countMap map[int][]assembly) {
	// newly created and complete assemblies (circular)
	var assemblies []assembly
	// number of nodes to try to synthesize to from each node (plus the natural overlap)
	synthCount := int(math.Max(5.0, 0.05*float64(len(nodes))))

	// sort in descending order of end index
	sort.Slice(nodes, func(i, j int) bool {
		return nodes[i].end > nodes[j].end
	})

	// for every node in the list of reverse sorted nodes
	for i, n := range nodes {
		// for every overlapping fragment + synthCount more
		for _, r := range n.reach(nodes, i, synthCount) {
			// for every assembly on the reachable fragment
			for _, a := range r.assemblies {
				// see if we can create a new one with this node included
				if newa, made, complete := a.add(n); made {
					if complete {
						// we've completed a circlular plasmid, has mapped
						// back onto itself
						assemblies = append(assemblies, newa)

					} else {
						// and add to this node's list of assemblies
						n.assemblies = append(n.assemblies, newa)
					}
				}
			}
		}
	}

	// build up the fragment count to assemblies list map
	for _, a := range assemblies {
		if as, ok := countMap[a.len()]; ok {
			countMap[a.len()] = append(as, a)
		}
	}

	return countMap
}
