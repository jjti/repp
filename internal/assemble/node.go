package assemble

import (
	"math"
)

// node is a single node within the DP tree for building up
// the vector from smaller building fragments
type node struct {
	// id of the node's source in the database (will be used to avoid off-targets within it)
	id string

	// unique-id of a match, the  start index % seq-length + id
	// (unique identified that catches nodes that cross the zero-index)
	uniqueID string

	// start of this node on the target vector (which has been 3x'ed for BLAST)
	start int

	// end of this node on the target vector
	end int

	// assemblies that span from this node to the end of the vector
	assemblies []assembly
}

// distTo returns the distance between the start of this node and the end of the other.
// assumes that this node starts before the other
// will return a negative number if this node overlaps with the other.
// will return a positive number if this node doesn't overlap with the other.
func (n *node) distTo(other node) (bpDist int) {
	return other.start - n.end
}

// synthDist returns the number of synthesized fragments that would need to be created
// between one node and another if the two were to be joined, with no existing
// fragments/nodes in-between, in an assembly
func (n *node) synthDist(other node) (synthCount int) {
	dist := n.distTo(other)

	if dist > -(conf.Fragments.MinHomology) {
		// can't be negative before the ceil
		floatDist := math.Max(1.0, float64(dist))
		// split up the distance between them by the max synthesized fragment size
		return int(math.Ceil(floatDist / float64(conf.Synthesis.MaxLength)))
	}

	// if the dist is <-20, there's enough overlap already for PCR
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
func (n *node) costTo(other node) (cost float32) {
	dist := n.distTo(other)

	if dist <= -(conf.Fragments.MinHomology) {
		// there's already overlap between this node and the one being tested
		// and enough for exiting homology to be enough.
		// estimating two primers, 20bp each.
		return 40 * conf.PCR.BPCost
	}

	// we need to create a new synthetic fragment to get from this fragment to the next
	// account for both the bps between them as well as the additional bps we need to add
	// for homology between the two
	return (float32(conf.Fragments.MinHomology) + float32(dist)) * conf.Synthesis.BPCost
}

// reach returns a slice of nodes that overlap with, or are the first synth_count nodes
// away from this one within a slice of ordered nodes
//
// nodes are nodes sorted in ascending start index order
// i is this node's index in the slice of nodes
// synthCount is the number of nodes to try to synthesize to, in addition to the
// 	nodes that are reachable with just existing homology
func (n *node) reach(nodes []node, i, synthCount int) (reachable []node) {
	reachable = []node{}

	// accumulate the nodes that overlap with this one
	for true {
		i++

		// we've run out of nodes
		if i >= len(nodes) {
			return reachable
		}

		// these nodes overlap by enough for assembly without PCR
		if n.distTo(nodes[i]) <= -(conf.Fragments.MinHomology) {
			reachable = append(reachable, nodes[i])
		} else if synthCount > 0 {
			// there's not enough existing overlap, but we can synthesize to it
			synthCount--
			reachable = append(reachable, nodes[i])
		} else {
			break
		}
	}

	return
}
