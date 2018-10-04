package assemble

import (
	"sort"
)

// cost is for building up the minimum assembly cost of each node
//
// This finds each node's $ amount from the "end of the vector" to
// get the minumum cost of an assembly including the node
//
// For example, if a node overlaps the end of a vector target sequence
// but needs to be PCR'ed to create it, the cost of that node is
// the cost per bp of DNA times the summed length of both primers
//
// But if the node is two "nodes" away from the end of the target sequence,
// it is its own cost to prepare plus that of the other nodes in between
// it and the end of the vector
func cost(nodes []node) map[node]float32 {
	// sort by start index
	sort.Slice(nodes, func(i, j int) bool {
		return nodes[i].start < nodes[j].start
	})

	costs := make(map[node]float32)

	// costOf is for finding the minimum cost of a single node.
	// its index is used to quickly figure out which fragments
	// could potentially come after this one
	var costOf func(int, node) float32
	costOf = func(i int, n node) float32 {
		// already memoized node's cost, return
		if nCost, cached := costs[n]; cached {
			return nCost
		}

		// if this overlaps with the end of the vector, return
		if n.terminal {
			// guestimate a fixed cost to PCR this fragment and cache
			nCost := conf.PCR.BPCost * 60
			costs[n] = nCost
			return nCost
		}

		// for each node that's further along in the list of nodes,
		// estimate the cost of getting from this to those
		minCost := n.costTo(nodes[i+1])
		for j, nn := range nodes[i+1:] {
			if costTo := n.costTo(nn) + costOf(i+j+1, nn); costTo < minCost {
				minCost = costTo
				n.next = i + j + 1
			}
		}
		costs[n] = minCost
		return minCost
	}

	// fill costs
	costOf(0, nodes[0])

	return costs
}
