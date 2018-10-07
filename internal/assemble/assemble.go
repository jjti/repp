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
			terminus: m.End >= to,
		})
	}
}
