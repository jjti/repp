// Package assemble is for traversing the target matches and
// creating a list of possible assemblies for making the target
// sequence
package assemble

import (
	"github.com/jjtimmons/decvec/config"
	"github.com/jjtimmons/decvec/internal/frag"
)

var (
	conf = config.NewConfig()
)

// Assemble the BLAST matches, as nodes, into an assembly
// that spans the template/target sequence and has a minimum
// total cost
func Assemble(matches []frag.Match, seqL int, from int, to int) {
	// map fragment Matches to nodes
	// store whether they're entry or terminal nodes,
	// based on whether they overlap the first or last bp of
	// the search range, respectively
	var nodes []node
	for _, m := range matches {
		nodes = append(nodes, node{
			id:       m.ID,
			start:    m.Start,
			end:      m.End,
			entry:    m.Start <= from,
			terminal: m.End >= to,
		})
	}

	// remove nodes that will never be in an assembly
	// beneath the upper-limit in settings and from the CLI
	nodes = limit(nodes, seqL)

	// add a single mock-starting node that's zero-bp long
	// this is used to force a synthesis from this node to the
	// rest and is used in-case there isn't already a node that's
	// reachable from this one
	nodes = append(nodes, node{
		id:    "mock-entry-node",
		start: from,
		end:   from,
		entry: true,
	})

	// traverse the nodes by finding their costs (minimum cost to be
	// included in an assembly that spans the whole vector) and
	// "filling in" each node by creating primers on it and synthesizing
	// and synthetic fragments if need be
	traverse(nodes)
}
