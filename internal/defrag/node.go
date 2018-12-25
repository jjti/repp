package defrag

import (
	"fmt"
	"math"
	"regexp"
	"strconv"
	"strings"

	"github.com/jjtimmons/defrag/config"
)

// node is a single node within the DP tree for building up
// the vector from smaller building fragments
type node struct {
	// id of the node's source in the database (will be used to avoid off-targets within it)
	id string

	// uniqueID of a match, the start index % seq-length + id
	// used identified that catches nodes that cross the zero-index
	uniqueID string

	// seq of the node against the range on the target sequence
	seq string

	// start of this node on the target vector (which has been 3x'ed for BLAST)
	start int

	// end of this node on the target vector
	end int

	// db that the node was derived from
	db string

	// assemblies that span from this node to the end of the vector
	assemblies []assembly

	// cost of the node. eg: fragments from Addgene cost $65
	cost float64

	// url to get the node, eg: link for addgene page
	url string

	// primers for amplifying a node, set with setPrimers
	primers []Primer

	// assembly configuration
	conf *config.Config
}

// newNode creates a node from a match
func newNode(m match, seqL int, conf *config.Config) *node {
	cost := 0.0
	url := ""
	if strings.Contains(m.entry, "addgene") {
		re := regexp.MustCompile("^.*addgene\\|(\\d*)")
		match := re.FindStringSubmatch(m.entry)

		cost = conf.AddGeneVectorCost
		url = fmt.Sprintf("https://www.addgene.org/%s/", match[1])
	}

	return &node{
		id:       m.entry,
		uniqueID: strconv.Itoa(m.start%seqL) + m.entry,
		seq:      strings.ToUpper(m.seq),
		start:    m.start,
		end:      m.end,
		db:       m.db,
		cost:     cost,
		url:      url,
		conf:     conf,
	}
}

// copy returns a deep dopy of a node. used because nodes are mutated
// during assembly filling, and we don't want primers being shared between
// nodes in different assemblies
func (n *node) copy() *node {
	return &node{
		id:       n.id,
		seq:      n.seq,
		uniqueID: n.uniqueID,
		start:    n.start,
		end:      n.end,
		db:       n.db,
		conf:     n.conf,
		cost:     n.cost,
		url:      n.url,
	}
}

// fragment converts a node into a fragment
func (n *node) fragment() Fragment {
	// should have primers by this point, add up their expected cost
	fragType := vector
	cost := n.cost
	for _, p := range n.primers {
		fragType = pcr // has primers, is a PCR fragment
		cost += conf.PCR.BPCost * float64(len(p.Seq))
	}

	return Fragment{
		ID:      n.id,
		Seq:     strings.ToUpper(n.seq),
		Entry:   n.id,
		Type:    fragType,
		URL:     n.url,
		Primers: n.primers,
		Cost:    cost,
	}
}

// distTo returns the distance between the start of this node and the end of the other.
// assumes that this node starts before the other
// will return a negative number if this node overlaps with the other and positive otherwise
func (n *node) distTo(other *node) (bpDist int) {
	return other.start - n.end
}

// synthDist returns the number of synthesized fragments that would need to be created
// between one node and another if the two were to be joined, with no existing
// fragments/nodes in-between, in an assembly
func (n *node) synthDist(other *node) (synthCount int) {
	dist := n.distTo(other)

	if dist <= n.conf.PCR.MaxEmbedLength {
		// if the dist is <MaxEmbedLength, we can PCR our way there
		// and add the mutated bp between the nodes with PCR
		return 0
	}

	// if it's > 5, we have to synthesize to get there (arbitatry, TODO: move to settings)
	// can't be negative before the ceil
	floatDist := math.Max(1.0, float64(dist))

	// split up the distance between them by the max synthesized fragment size
	return int(math.Ceil(floatDist / float64(n.conf.Synthesis.MaxLength)))
}

// costTo calculates the $ amount needed to get from this fragment
// to the other node passed, either by PCR or synthesis
//
// PCR is, of course, preferred if we can add homology between this and
// the other fragment with PCR homology arms alone
//
// Otherwise we find the total synthesis distance between this and
// the other fragment and divide that by the cost per bp of synthesized DNA
//
// This does not add in the cost of procurement, which is added to the assembly cost
// in assembly.add()
func (n *node) costTo(other *node) (cost float64) {
	dist := n.distTo(other)

	if dist <= n.conf.PCR.MaxEmbedLength {
		if dist < -(n.conf.Fragments.MinHomology) {
			// there's already enough overlap between this node and the one being tested
			// estimating two primers, 20bp each
			return 46 * n.conf.PCR.BPCost
		}

		// we have to create some additional primer sequence to reach the next fragment
		// guessing 46bp (2x primer3 target primer length) plus half MinHomology on each primer
		return float64(46+n.conf.Fragments.MinHomology) * n.conf.PCR.BPCost
	}

	// we need to create a new synthetic fragment to get from this fragment to the next
	// to account for both the bps between them as well as the additional bps we need to add
	// for homology between the two
	fragLength := n.conf.Fragments.MinHomology + dist
	return n.conf.SynthCost(fragLength)
}

// reach returns a slice of node indexes that overlap with, or are the first synth_count nodes
// away from this one within a slice of ordered nodes
func (n *node) reach(nodes []*node, i, synthCount int) (reachable []int) {
	reachable = []int{}

	// accumulate the nodes that overlap with this one
	for true {
		i++

		// we've run out of nodes
		if i >= len(nodes) {
			return reachable
		}

		// these nodes overlap by enough for assembly without PCR
		if n.distTo(nodes[i]) <= -(n.conf.Fragments.MinHomology) {
			reachable = append(reachable, i)
		} else if synthCount > 0 {
			// there's not enough existing overlap, but we can synthesize to it
			synthCount--
			reachable = append(reachable, i)
		} else {
			break
		}
	}
	return
}

// junction checks for and returns any 100% identical homology between the end of this
// node and the start of the other. returns an empty string if there's no junction between them
func (n *node) junction(other *node, minHomology, maxHomology int) (junction string) {
	thisSeq := strings.ToUpper(n.seq)
	otherSeq := strings.ToUpper(other.seq)

	//      v-maxHomology from end    v-minHomology from end
	// ------------------------------------
	//                    -----------------------------
	start := len(thisSeq) - maxHomology
	end := len(thisSeq) - minHomology

	if start < 0 {
		return ""
	}

	// for every possible start index
	for i := start; i <= end; i++ {
		// traverse from that index to the end of the seq
		j := 0
		for k := i; k < len(thisSeq); k++ {
			if thisSeq[k] != otherSeq[j] {
				break
			} else {
				j++
			}

			// we made it to the end of the sequence, there's a junction
			if k == len(thisSeq)-1 {
				return thisSeq[i:]
			}
		}
	}
	return
}

// synthTo returns synthetic fragments to get this node to the next.
// It creates a slice of building fragments that have homology against
// one another and are within the upper and lower synthesis bounds.
// seq is the vector's sequence. We need it to build up the target
// vector's sequence
func (n *node) synthTo(next *node, seq string) (synthedFrags []Fragment) {
	// check whether we need to make synthetic fragments to get
	// to the next fragment in the assembly
	fragC := n.synthDist(next)
	if fragC == 0 {
		return nil
	}

	// make the slice
	synthedFrags = []Fragment{}

	// length of each synthesized fragment
	fragL := n.distTo(next) / fragC
	if n.conf.Synthesis.MinLength > fragL {
		// need to synthesize at least Synthesis.MinLength bps
		fragL = n.conf.Synthesis.MinLength
	}

	// account for homology on either end of each synthetic fragment
	fragL += n.conf.Fragments.MinHomology * 2
	seq += seq // double to account for sequence across the zero-index

	// slide along the range of sequence to create synthetic fragments
	// and create one at each point, each w/ MinHomology for the fragment
	// before and after it
	for fragIndex := 0; fragIndex < int(fragC); fragIndex++ {
		start := n.end - n.conf.Fragments.MinHomology // start w/ homology
		start += fragIndex * fragL                    // slide along the range to cover
		end := start + fragL + n.conf.Fragments.MinHomology

		sFrag := Fragment{
			ID:   fmt.Sprintf("%s-synthetic-%d", n.id, fragIndex+1),
			Seq:  strings.ToUpper(seq[start:end]),
			Type: synthetic,
			Cost: n.conf.SynthCost(len(n.seq)) + n.cost,
		}
		synthedFrags = append(synthedFrags, sFrag)
	}
	return
}
