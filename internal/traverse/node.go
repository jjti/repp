package traverse

import (
	"math"
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

	// cost to prepare this node for assembly
	cost int

	// entry signifies whether this node is a valid entry node -- ie one that could be the first
	// node in a vector assembly
	entry bool

	// terminal signifies whether this node is a valid terminal node -- ie whether it overlaps
	// with the last bp in the assembly
	terminal bool
}

// synthDist returns the number of synthesized fragments
// that would need to be created between one node and another
// if the two were to be joined, with no existing fragments in between, in an assembly
func (n *node) synthDist(other node) int {
	distToNext := 0

	// this fragment is before the other
	if n.start < other.start {
		distToNext = other.start - n.end
	} else {
		distToNext = n.start - other.end
	}

	if distToNext > 0 {
		// split up the distance between them by the maximum synthesized
		// fragment size
		return int(math.Ceil(float64(distToNext) / float64(conf.Synthesis.MaxLength)))
	}

	// they already overlap
	return 0
}
