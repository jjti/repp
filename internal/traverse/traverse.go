// Package traverse is for traversing the target matches and
// creating a list of possible assemblies for making the target
// sequence
package traverse

import (
	"github.com/jjtimmons/decvec/internal/frag"
)

// node is a single node within the DP tree for building up
// the vector from smaller building fragments
type node struct {
	// id of the node's source in the database (will be used to avoid off-target primers in it)
	id string

	// start of this node on the target vector (which has been 3x'ed for BLAST)
	start int

	// end of this node on the target vector
	end int
}

// Traverse the matches on the target fragment
// and build up a list of possible assemblies
func Traverse(f *frag.Fragment) {
	// map fragment Matches to nodes
	var nodes []node
	for _, m := range f.Matches {
		nodes = append(nodes, node{
			id:    m.ID,
			start: m.Start,
			end:   m.End,
		})
	}

	// remove nodes that will never be in an assembly
	// beneath the upper-limit in settings and from the CLI
	// shortNodes := upperLimit(nodes, len(f.Seq))
}
