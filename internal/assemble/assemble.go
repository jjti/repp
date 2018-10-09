// Package assemble turns blast Matches into building Fragments
package assemble

import (
	"fmt"

	"github.com/jjtimmons/decvec/config"
	"github.com/jjtimmons/decvec/internal/dvec"
)

var (
	conf = config.NewConfig()
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
func Assemble(matches []dvec.Match, seqL int) {
	// map fragment Matches to nodes
	// store whether they're entry or terminal nodes based on whether they
	// overlap the first or last bp of the search range, respectively
	var nodes []node
	for _, m := range matches {
		nodes = append(nodes, node{
			id:       m.ID,
			uniqueID: string(m.Start%seqL) + m.ID,
			start:    m.Start,
			end:      m.End,
		})
	}

	assemble(nodes)

}

//
func assemble(nodes []node) {
	// build up a the map from fragment to assemblies with that fragment count
	assemblies := build(nodes)

	// find the pareto optimal set (best of either cost/fragment# or both)
	paretos := pareto(assemblies)

	//
	fmt.Println(paretos)

}
