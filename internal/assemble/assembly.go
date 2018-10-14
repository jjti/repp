package assemble

import (
	"fmt"

	"github.com/jjtimmons/decvec/internal/dvec"
)

// assembly is a slice of nodes ordered by the nodes
// distance from the end of the target vector
type assembly struct {
	// nodes, ordered by distance from the "end" of the vector
	nodes []node

	// estimated cost of making this assembly
	cost float32

	// total number of synthetic nodes that will be needed to make this
	synths int
}

// add a node to the start of this assembly.
//
// Update the cost of the assembly to include the link between the new first node and the one after it.
// Store the node's id in the list of node ids.
// Complete  an assembly if a node has matched up onto itself across the zero-index.
func (a *assembly) add(n node) (newAssembly assembly, created, complete bool) {
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
	if newCount > conf.Fragments.MaxCount {
		return newAssembly, false, false
	}

	// calc the estimated dollar cost of getting to the next node
	annealCost := a.nodes[len(a.nodes)-1].costTo(n)

	if complete {
		if synths < 1 {
			// costs nothing to anneal node to self, already been PCR'ed
			annealCost = 0
		}

		return assembly{
			nodes:  append(a.nodes, n),
			cost:   a.cost + annealCost,
			synths: a.synths + synths,
		}, true, true
	}

	return assembly{
		nodes:  append(a.nodes, n),
		cost:   a.cost + annealCost,
		synths: a.synths + synths,
	}, true, false
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

// fill traverses the nodes in an assembly and converts them into
// fragments -- either pcr fragments or synthetic fragments -- that
// will match the sequences that come together in vitro during assembly
//
// it can fail out. For example, a PCR Fragment may have off-targets in
// the parent vector. If that happens, we return the problem node and nil
// building fragments
func (a *assembly) fill(seq string) (blacklist node, frags []dvec.Fragment) {
	for i, n := range a.nodes {
		// last node, do nothing
		// here only to allow for vector "circularization" if we need to synthesize
		// from a.nodes[len(a.nodes)-2] to a.nodes[len(a.nodes)-1]
		if i == len(a.nodes)-1 {
			break
		}

		// convert
		frag := n.fragment()

		// try and make primers for the fragment
		fragPrimers, err := primers(frag)
		if err != nil {
			fmt.Printf("node %s failed: %v", n.id, err)
			// return the node as a blackmailed node if pcr fails
			return n, nil
		}

		// set primers and store this to the list of building fragments
		frag.Primers = fragPrimers
		frags = append(frags, frag)

		// add synthesized fragments between the two if necessary
		if synthedFrags := n.synthTo(a.nodes[i+1], seq); synthedFrags != nil {
			frags = append(frags, synthedFrags...)
		}
	}

	return
}
