package defrag

import (
	"fmt"
	"strings"

	"github.com/jjtimmons/defrag/config"
)

// assembly is a slice of nodes ordered by the nodes
// distance from the end of the target vector
type assembly struct {
	// nodes, ordered by distance from the "end" of the vector
	nodes []*node

	// estimated cost of making this assembly
	cost float64

	// total number of synthetic nodes that will be needed to make this
	synths int
}

// add a node to the start of this assembly.
//
// Update the cost of the assembly to include the link between the new first node and the one after it.
// Store the node's id in the list of node ids.
// Complete  an assembly if a node has matched up onto itself across the zero-index.
func (a *assembly) add(n *node, maxCount int) (newAssembly assembly, created, complete bool) {
	// check if we could complete an assembly with this new node
	complete = n.uniqueID == a.nodes[0].uniqueID
	// calc the number of synthesis fragments needed to get to this next node
	synths := a.nodes[len(a.nodes)-1].synthDist(n)

	newCount := a.len() + synths
	if !complete {
		// only adding a new node if not annealing to starting node
		newCount++
	}

	// stay beneath upper limit
	if newCount > maxCount {
		return assembly{}, false, false
	}

	// we will create a new assembly
	created = true

	// calc the estimated dollar cost of getting to the next node
	annealCost := a.nodes[len(a.nodes)-1].costTo(n)

	// check whether the node is already contained in the assembly
	nodeContained := false
	for _, included := range a.nodes {
		if included.id == n.id {
			nodeContained = true
			break
		}
	}

	if !nodeContained {
		// add the cost of procuring this node to the total assembly cost
		annealCost += n.cost
	}

	if complete {
		if synths < 1 {
			// costs nothing to anneal node to self, already been PCR'ed
			annealCost = 0
		}

		return assembly{
			nodes:  append(a.nodes, n.copy()),
			cost:   a.cost + annealCost,
			synths: a.synths + synths,
		}, created, complete
	}

	return assembly{
		nodes:  append(a.nodes, n.copy()),
		cost:   a.cost + annealCost,
		synths: a.synths + synths,
	}, created, false
}

// contains returns if the id of the node has already been seen in this assembly
func (a *assembly) contains(n node) (hasNode bool) {
	for _, otherN := range a.nodes {
		// they're the same if they have the same id and start index
		// id isn't enough by itself because there may be multiple with the same
		// entry id in the BLAST db
		if otherN.uniqueID == n.uniqueID {
			return true
		}
	}
	return false
}

// len returns len(assembly.nodes) + the synthesis fragment count
func (a *assembly) len() int {
	return len(a.nodes) + a.synths
}

// fill traverses the nodes in an assembly and converts them into fragments
//
// - Vector fragments if a node spans the whole target range (entirely contains it)
//
// - PCR fragments if the node subselects a region (BLAST match)
//
// - Synthetic fragments if there wasn't a match for a region
//
// It can fail. For example, a PCR Fragment may have off-targets in
// the parent vector
func (a *assembly) fill(seq string, conf *config.Config) (frags []Fragment, err error) {
	// edge case where a single node fills the whole target vector. Return just a single
	// "fragment" (of Vector type... misnomer) that matches the target sequence 100%
	if a.len() == 1 && len(a.nodes[0].seq) >= len(seq) {
		n := a.nodes[0]
		return []Fragment{
			Fragment{
				ID:    n.id,
				Seq:   strings.ToUpper(n.seq)[0:len(seq)],
				Entry: n.id,
				Type:  vector,
				URL:   n.url,
				Cost:  n.cost, // only the ordering cost, no PCR/Synth etc
			},
		}, nil
	}

	// do two loops. the first is to fill in primers. let each node create primers for
	// itself that will span it to the last and next fragments (if reachable). this has to be done
	// in two loops because synthesis depends on nodes' ranges, and nodes' ranges
	// may change during setPrimers (we let primer3_core pick from a range of primer options)
	for i, n := range a.nodes {
		// last node, do nothing
		// here only to allow for vector "circularization" if we need to synthesize
		// from a.nodes[len(a.nodes)-2] to a.nodes[len(a.nodes)-1]
		if i > 0 && n.uniqueID == a.nodes[0].uniqueID {
			break
		}

		// try and make primers for the fragment (need last and next nodes)
		var last *node
		if i == 0 {
			// mock up a last fragment that's to the left of this starting node
			final := a.nodes[len(a.nodes)-2]
			last = &node{
				start: final.start - len(seq),
				end:   final.end - len(seq),
				conf:  conf,
			}
		} else {
			last = a.nodes[i-1]
		}

		var next *node
		if i < len(a.nodes)-1 {
			next = a.nodes[i+1]
		} else {
			// mock up a next fragment that's to the right of this terminal node
			first := a.nodes[0]
			next = &node{
				start: first.start + len(seq),
				end:   first.end + len(seq),
				conf:  conf,
			}
		}

		if err := n.setPrimers(last, next, seq, conf); err != nil || len(n.primers) < 2 {
			return nil, fmt.Errorf("failed to fill %s: %v", n.id, err)
		}
	}

	// do another loop to covert the nodes with primers to fragments and
	// synthesize the sequence between nodes
	for i, n := range a.nodes {
		// last node, do nothing
		// here only to allow for vector "circularization" if we need to synthesize
		// from a.nodes[len(a.nodes)-2] to a.nodes[len(a.nodes)-1]
		if i > 0 && n.uniqueID == a.nodes[0].uniqueID {
			break
		}

		// convert to a fragment from a node, store this to the list of building fragments
		// cost is calculated here as the summed cost of both primers (based on length)
		frags = append(frags, n.fragment())

		// add synthesized fragments between this node and the next (if necessary)
		if synthedFrags := n.synthTo(a.nodes[i+1], seq); synthedFrags != nil {
			frags = append(frags, synthedFrags...)
		}
	}
	return
}
