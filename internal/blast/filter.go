package blast

import (
	"sort"

	"github.com/jjtimmons/decvec/internal/frag"
)

// filter is for "propertizing" the matches from BLAST and removing
// those that start past 2x the target sequence's length
//
// propertizing fragment matches means removing those that are completely
// self-contained in other fragments: the larger of the available fragments
// will be the better one, since it covers a greater region and will almost
// always be preferable to the smaller one
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
	return afterStart
}
