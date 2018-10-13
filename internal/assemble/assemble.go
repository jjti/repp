// Package assemble turns blast Matches into building Fragments
package assemble

import (
	"fmt"

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
func Assemble(matches []dvec.Match, seq string) {
	// map fragment Matches to nodes
	// store whether they're entry or terminal nodes based on whether they
	// overlap the first or last bp of the search range, respectively
	var nodes []node
	for _, m := range matches {
		nodes = append(nodes, node{
			id:       m.Entry,
			uniqueID: string(m.Start%len(seq)) + m.Entry,
			start:    m.Start,
			end:      m.End,
		})
	}

	// build up a the map from fragment to assemblies with that fragment count
	assemblies := build(nodes)

	// assemble into a list of fragment lists (building fragments)
	assemble(assemblies, seq, [][]dvec.Fragment{})
}

// assemble the nodes up against the vector's sequence,
// create lists of building fragment lists (ways to get to the target seq)
func assemble(assemblies []assembly, seq string, found [][]dvec.Fragment) [][]dvec.Fragment {
	// find the pareto optimal set (best of either cost/fragment# or both)
	paretos := pareto(assemblies)

	// convert and fill the fragments
	for _, p := range paretos {
		n, fs := p.fill(seq)

		// if a node in the assembly fails to be prepared,
		// remove all assemblies with the node and try again
		if fs == nil {
			filteredAssemblies := filterNode(n, assemblies)
			return assemble(filteredAssemblies, seq, found)
		}

		// add this list of fragments to the list of such
		found = append(found, fs)

		fmt.Printf("%+v", p)
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
