// Package traverse is for traversing the target matches and
// creating a list of possible assemblies for making the target
// sequence
package traverse

import (
	"fmt"
	"math"
	"sort"

	"github.com/jjtimmons/decvec/internal/frag"
)

type nodeType int

const (
	pcr   nodeType = 1
	synth nodeType = 2
)

// Traverse the matches on the target fragment
// and build up a list of possible assemblies
func Traverse(f *frag.Fragment) [][]frag.Fragment {
	var assemblies [][]frag.Fragment

	// get each fragments minimum number of fragments in an assembly
	// dists := calcFragDistance(f)

	return assemblies
}

// calcFragDistance is for calculating each fragment's distance from
// a vector "endpoint"
//
// for each building fragment, find the minimum number of fragments
// needed to get from it to a "last-bp" (2x the length of the target
// fragment's sequence length)
//
// maxSynth is the maximum length of a synthesized stretch of DNA
func calcFragDistance(target *frag.Fragment, maxSynth int) map[frag.Match]int {
	// last bp within range being scanned
	lastBP := 2 * len(target.Seq)
	// map from a fragment to the minimum number of fragments
	// needed in an assembly including it
	dists := make(map[frag.Match]int)

	// sort by start
	sort.Slice(target.Matches, func(i, j int) bool {
		return target.Matches[i].Start < target.Matches[j].Start
	})

	// add a "sink" to ensure there's a building fragment just past
	// the end of the scanned range
	matches := append(target.Matches, frag.Match{
		ID:    "sink",
		Start: lastBP + 1,
		End:   lastBP + 1,
	})

	// number of fragments distance for this fragment
	// through to the target bp (2x the seqLength)
	//
	// if dist has been calculated already
	// 	return that
	// else if this fragment ends after the lastBP,
	// 	return 1 and store 1 to the cache for this fragment
	// else
	//	return the min-distance from this to all fragments
	// 	this could reasonably be assembled with
	var distFor func(int, frag.Match) int
	distFor = func(i int, f frag.Match) int {
		if dist, cached := dists[f]; cached {
			return dist
		}

		fmt.Printf("%d %d %t", i, f.End, f.End >= lastBP)

		if f.End >= lastBP {
			dists[f] = 1
			return 1
		}

		// calculate how many synthesis fragments are necessary to get
		// from this fragment to each fragment after it
		var synthCount []int
		for _, m := range matches[i+1:] {
			synthsToNext := 0
			distToNext := m.Start - f.End
			if distToNext > 0 {
				synthsToNext = int(math.Floor(float64(distToNext) / float64(maxSynth)))
			}
			synthCount = append(synthCount, synthsToNext)
		}

		// find the minimum distance among these options, accounting
		// for the fact that we need to each synthesis is an additional (unseen) fragment
		minDistNext := distFor(i+1, matches[i+1]) + synthCount[i+1]
		for i, f := range matches[i+2:] {
			if synthCount[i]+distFor(i, f) < minDistNext {
				minDistNext = synthCount[i] + distFor(i, f)
			}
		}

		// store in cache in case this is referenced later
		dists[f] = minDistNext
		return minDistNext
	}

	// fill cache
	distFor(0, matches[0])

	return dists
}
