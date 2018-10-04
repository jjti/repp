package assemble

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

	// entry signifies whether this node is a valid entry node -- ie one that could be the first
	// node in a vector assembly
	entry bool

	// terminal signifies whether this node is a valid terminal node -- ie whether it overlaps
	// with the last bp in the assembly
	terminal bool

	// next is the index of the next fragment after this one in the minimum-cost-assembly
	// involving this node
	next int
}

// distToNext returns the number of bps between the end of the earlier of this
// and the other fragment and the start of the other
func (n *node) distTo(other node) int {
	dist := 0

	// this fragment is before the other
	if n.start < other.start {
		dist = other.start - n.end
	} else {
		dist = n.start - other.end
	}

	return dist
}

// synthDist returns the number of synthesized fragments
// that would need to be created between one node and another
// if the two were to be joined, with no existing fragments in between, in an assembly
func (n *node) synthDist(other node) int {
	dist := n.distTo(other)

	if dist > 0 {
		// split up the distance between them by the maximum synthesized
		// fragment size
		return int(math.Ceil(float64(dist) / float64(conf.Synthesis.MaxLength)))
	}

	// they already overlap
	return 0
}

// costTo calculates the $ amount needed to get from this fragment
// to the other node passed, either by PCR or synthesis
//
// PCR is of course preferred if we can add homology between this and
// the other fragment with PCR homology arms alone
//
// Otherwise we find the total synthesis distance between this and
// the other fragmetn and divide that by the cost per bp of synthesized DNA
func (n *node) costTo(other node) float32 {
	dist := n.distTo(other)

	if dist < 20 {
		// there's already overlap between this node and the one being tested
		// or little enough of a distance for it to be added via PCR
		return 60 * conf.PCR.BPCost
	}

	return float32(dist) * conf.Synthesis.BPCost
}
