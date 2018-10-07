package blast

import (
	"sort"

	"github.com/jjtimmons/decvec/internal/frag"
)

// filter is for "propertizing" the matches from BLAST and removing
// those that are outside the range of overlap we're testing
//
// propertizing fragment matches means removing those that are completely
// self-contained in other fragments: the larger of the available fragments
// will be the better one, since it covers a greater region and will almost
// always be preferable to the smaller one
//
// the search range is from 1x-2x the length of the target vector sequence
// and it used in an effort to find fragments that overlap with large regions
// of the target vector but do so in regions outside the central target range
//
// also remove small fragments here, that are too small to be useful during
// assembly
func filter(matches []frag.Match, from int, to int) []frag.Match {
	// sort matches by their start index
	// if they're same, put the larger one first
	sort.Slice(matches, func(i, j int) bool {
		if matches[i].Start != matches[j].Start {
			return matches[i].Start < matches[j].Start
		}
		return matches[i].Length() > matches[j].Length()
	})

	// only include those that aren't encompassed in the one before it
	var properMatches []frag.Match
	for _, m := range matches {
		lastMatch := len(properMatches) - 1
		if lastMatch < 0 || m.End > properMatches[lastMatch].End {
			properMatches = append(properMatches, m)
		}
	}

	// remove fragments that start end before 1x the target vector sequence's length
	var beforeEnd []frag.Match
	for _, m := range properMatches {
		if m.End > from {
			beforeEnd = append(beforeEnd, m)
		}
	}

	// remove fragments that start past 2x the target vector sequence's length
	var afterStart []frag.Match
	for _, m := range beforeEnd {
		if m.Start < to {
			afterStart = append(afterStart, m)
		}
	}

	// remove fragments that are larger the minimum cut off size
	var largeEnough []frag.Match
	for _, m := range afterStart {
		if m.Length() > 50 {
			largeEnough = append(largeEnough, m)
		}
	}

	return largeEnough
}
