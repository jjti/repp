package defrag

import (
	"github.com/jjtimmons/defrag/config"
)

// Assemble the BLAST matches into assemblies that span the target sequence
//
// First build up assemblies, creating all possible assemblies that are
// beneath the upper-bound limit on the number of fragments
//
// Then find the pareto optimal solutions that minimize either the cost
// of the assembly or the number of fragments (relative to the other
// assembly plans)
//
// Then, for each pareto optimal solution, traverse the assembly and
// "fill-in" the nodes. Create primers on the node if it's a PCR Fragment
// or create a sequence to be synthesized if it's a synthetic fragment.
// Error out and repeat the build stage if a node fails to be filled
func Assemble(matches []Match, seq string, conf *config.Config) [][]Fragment {
	// map fragment Matches to nodes
	var nodes []node
	for _, m := range matches {
		nodes = append(nodes, new(m, len(seq), conf))
	}

	// build up slice of assemblies that could, within the upper-limit on
	// fragment count, be assembled to make the target vecto
	assemblies := build(nodes, conf.Fragments.MaxCount)

	// build up a map from fragment count to a sorted list of assemblies with that number
	paretos := pareto(assemblies)

	// convert and fill the fragments
	var found [][]Fragment
	for _, assemblies := range paretos {
		// get the first assembly that fills properly (cheapest workable solution)
		for _, singleAssembly := range assemblies {
			filledFrags := singleAssembly.fill(seq, conf)

			// if a node in the assembly fails to be prepared,
			// remove all assemblies with the node and try again
			if filledFrags != nil {
				// add this list of fragments to the list of such
				found = append(found, filledFrags)
				break
			}
		}
	}
	return found
}
