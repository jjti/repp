// Package assemble turns blast Matches into building Fragments
package assemble

import (
	"github.com/jjtimmons/decvec/config"
	"github.com/jjtimmons/decvec/internal/dvec"
)

var (
	conf = config.New()
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
func Assemble(matches []dvec.Match, seq string) [][]dvec.Fragment {
	// map fragment Matches to nodes
	var nodes []node
	for _, m := range matches {
		nodes = append(nodes, new(m, len(seq)))
	}

	// build up slice of assemblies that could, within the upper-limit on
	// fragment count, be assembled to make the target vecto
	assemblies := build(nodes)

	// assemble into a list of fragment lists (building fragments)
	return assemble(assemblies, seq, [][]dvec.Fragment{})
}

// assemble the nodes up against the vector's sequence,
// create lists of building fragment lists (ways to get to the target seq)
func assemble(assemblies []assembly, seq string, found [][]dvec.Fragment) [][]dvec.Fragment {
	// find the pareto optimal set (best of either cost/fragment# or both)
	paretos := pareto(assemblies)

	// convert and fill the fragments
	for _, p := range paretos {
		failedNode, filledFrags := p.fill(seq)

		// if a node in the assembly fails to be prepared,
		// remove all assemblies with the node and try again
		if filledFrags == nil {
			filteredAssemblies := filterNode(failedNode, assemblies)
			return assemble(filteredAssemblies, seq, found)
		}

		// add this list of fragments to the list of such
		found = append(found, filledFrags)
	}

	return found
}

// filterNode returns a list of assemblies, sans the node that's failed
// for some reason, because we know it will fail again in the future
func filterNode(black node, old []assembly) (new []assembly) {
	for _, a := range old {
		if !a.contains(black) {
			new = append(new, a)
		}
	}
	return
}
