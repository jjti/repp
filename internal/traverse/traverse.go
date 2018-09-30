// Package traverse is for traversing the target matches and
// creating a list of possible assemblies for making the target
// sequence
package traverse

import (
	"github.com/jjtimmons/decvec/config"
	"github.com/jjtimmons/decvec/internal/frag"
)

var (
	conf = config.NewConfig()
)

// Traverse the matches on the target fragment
// and build up a list of possible assemblies
func Traverse(f *frag.Fragment) {
	// we search from 1x-2x the vector's sequence length
	entryBP := len(f.Seq)
	terminalBP := 2 * len(f.Seq)

	// map fragment Matches to nodes
	// store whether they're entry or terminal nodes,
	// based on whether they overlap the first or last bp of
	// the search range, respectively
	var nodes []node
	for _, m := range f.Matches {
		nodes = append(nodes, node{
			id:       m.ID,
			start:    m.Start,
			end:      m.End,
			entry:    m.Start < entryBP,
			terminal: m.End > terminalBP,
		})
	}

	// remove nodes that will never be in an assembly
	// beneath the upper-limit in settings and from the CLI
	// shortNodes := upperLimit(nodes, len(f.Seq))
}
