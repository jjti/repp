package traverse

import (
	"sort"

	"github.com/jjtimmons/decvec/config"
)

// upperLimit is for removing all fragments that are greater than the target distance
// from the end of the vector, only keep those that are beneath a
// hard limit in the number of fragments in a vector
func upperLimit(nodes []node, seqL int) []node {
	// get the config to figure out the maximum length of DNA synthesis
	c := config.NewConfig()

	// get each fragments minimum number of fragments in an assembly
	dists := distanceToEnd(nodes, seqL, c.Synthesis.MaxLength)

	var shortNodes []node
	for n, dist := range dists {
		// if it overlaps with the start of the search range
		// 	compare it against the max-distance from the end
		// else
		// 	check if it's max-distance - 1 from end and remove
		// 	it if it's not
		if n.entry {
			if dist < c.Fragments.MaxCount {
				shortNodes = append(shortNodes, n)
			}
		} else {
			if dist+1 < c.Fragments.MaxCount {
				shortNodes = append(shortNodes, n)
			}
		}
	}

	return shortNodes
}

// distanceToEnd is for calculating each node's minimum distance to
// a vector "endpoint". used to weed out fragments that are too far
// away from the end of the target
//
// for each building node, find the minimum number of fragments
// needed to get from it to a "last-bp" (2x the length of the target
// node's sequence length)
//
// maxSynth is the maximum length of a synthesized stretch of DNA
func distanceToEnd(nodes []node, seqL int, maxSynth int) map[node]int {
	// last bp within range being scanned
	lastBP := 2 * seqL
	// map from a node to the minimum number of fragments
	// needed in an assembly including it
	dists := make(map[node]int)

	// sort by start
	sort.Slice(nodes, func(i, j int) bool {
		return nodes[i].start < nodes[j].start
	})

	// add a "sink" to ensure there's a building node just past
	// the end of the scanned range
	sink := node{
		start: lastBP + 1,
		end:   lastBP + 1,
	}
	nodes = append(nodes, sink)

	// number of fragments distance for this node
	// through to the target bp (2x the seqLength)
	//
	// if dist has been calculated already
	// 	return that
	// else if this node ends after the lastBP,
	// 	return 1 and store 1 to the cache for this node
	// else
	//	return the min-distance from this to all fragments
	// 	this could reasonably be assembled with
	var distFor func(int, node) int
	distFor = func(i int, n node) int {
		if dist, cached := dists[n]; cached {
			return dist
		}

		if n.terminal {
			dists[n] = 1
			return 1
		}

		// find the minimum distance among these options
		minDistNext := 1 + n.synthDist(nodes[i+1]) + distFor(i+1, nodes[i+1])
		for j, nn := range nodes[i+1:] {
			distToFrag := 1 + n.synthDist(nn) + distFor(i+j, nn)
			if distToFrag < minDistNext {
				minDistNext = distToFrag
			}
		}

		// store in cache in case this is referenced later
		dists[n] = minDistNext
		return minDistNext
	}

	// fill cache
	distFor(0, nodes[0])

	// remove the sink node
	delete(dists, sink)

	return dists
}
