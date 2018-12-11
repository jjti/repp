package defrag

import (
	"sort"
)

// filter "proper-izes" the matches from BLAST
//
// TODO: filter on multiple BLAST databases. Keep ALL fragments from non-remote repositories (ie free
// fragments) but still properize out the remote fragments. So if a large fragment in the local repository
// completely emcompasses a smaller fragment in Addgene, remove the fragment from the Addgene database.
// Cannot do the same if there's a large fragment in the remote but not local db
//
// proper-izing fragment matches means removing those that are completely
// self-contained in other fragments: the larger of the available fragments
// will be the better one, since it covers a greater region and will almost
// always be preferable to the smaller one
//
// Circular-arc graph: https://en.wikipedia.org/wiki/Circular-arc_graph
//
// also remove small fragments here, that are too small to be useful during assembly
func filter(matches []Match, minSize int) (properized []Match) {
	properized = []Match{}

	// remove fragments that are shorter the minimum cut off size
	var largeEnough []Match
	for _, m := range matches {
		if m.Length() > minSize {
			largeEnough = append(largeEnough, m)
		}
	}

	// sort largeEnough by their start index
	// for fragments with equivelant starting indexes, put the larger one first
	sort.Slice(largeEnough, func(i, j int) bool {
		if largeEnough[i].Start != largeEnough[j].Start {
			return largeEnough[i].Start < largeEnough[j].Start
		}
		return largeEnough[i].Length() > largeEnough[j].Length()
	})

	// only include those that aren't encompassed by the one before it
	for _, m := range largeEnough {
		lastMatch := len(properized) - 1
		if lastMatch < 0 || m.End > properized[lastMatch].End {
			properized = append(properized, m)
		}
	}

	return properized
}
