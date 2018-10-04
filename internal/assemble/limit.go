package assemble

import (
	"sort"
)

// limit is for removing all fragments that cannot be included in any
// assembly with less than the upper limit on the number of fragments
func limit(nodes []node, seqL int) []node {
	// get each fragments' minimum number of fragments in an assembly
	dists := distance(nodes, seqL)

	var shortNodes []node
	for n, dist := range dists {
		// if it overlaps with the start of the search range
		// 	compare it against the max-distance from the end
		// else
		// 	check if it's max-distance - 1 from end and remove
		// 	it if it's not
		if n.entry {
			if dist < conf.Fragments.MaxCount {
				shortNodes = append(shortNodes, n)
			}
		} else {
			if dist+1 < conf.Fragments.MaxCount {
				shortNodes = append(shortNodes, n)
			}
		}
	}

	return shortNodes
}

// distance is for calculating each node's minimum distance to
// a vector "endpoint". used to weed out fragments that are too far
// away from the end of the target
//
// for each building node, find the minimum number of fragments
// needed to get from it to a "last-bp" (2x the length of the target
// node's sequence length)
func distance(nodes []node, seqL int) map[node]int {
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
		start:    lastBP + 1,
		end:      lastBP + 1,
		terminal: true,
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
		// we've already found distance from end of vector, return that
		if dist, cached := dists[n]; cached {
			return dist
		}

		// this node overlaps the end of the vector, return 1 (for it)
		if n.terminal {
			dists[n] = 1
			return 1
		}

		// find the minimum distance among these options
		minDist := 1000000000
		for j, nn := range nodes[i+1:] {
			if distTo := 1 + n.synthDist(nn) + distFor(i+j+1, nn); distTo < minDist {
				minDist = distTo
			}
		}

		// store in cache in case this is referenced later
		dists[n] = minDist
		return minDist
	}

	// fill cache
	distFor(0, nodes[0])

	// remove the sink node
	delete(dists, sink)

	return dists
}
