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
// will be the better one, since it covers a greater region
func filter(f *frag.Fragment) {
	// sort matches by their start index
	// if they're same, put the larger one first
	sort.Slice(f.Matches, func(i, j int) bool {
		if f.Matches[i].Start != f.Matches[j].Start {
			return f.Matches[i].Start < f.Matches[j].Start
		}
		return f.Matches[i].Length() > f.Matches[j].Length()
	})

	// only include those that aren't encompassed in the one before it
	var properMatches []frag.Match
	for _, m := range f.Matches {
		lastMatch := len(properMatches) - 1
		if lastMatch < 0 || m.End > properMatches[lastMatch].End {
			properMatches = append(properMatches, m)
		}
	}

	// remove fragments that start end before 1x the target vector sequence's length
	start := len(f.Seq)
	var beforeEnd []frag.Match
	for _, m := range properMatches {
		if m.End > start {
			beforeEnd = append(beforeEnd, m)
		}
	}

	// remove fragments that start past 2x the target vector sequence's length
	end := 2 * len(f.Seq)
	var afterStart []frag.Match
	for _, m := range beforeEnd {
		if m.Start < end {
			afterStart = append(afterStart, m)
		}
	}

	f.Matches = afterStart
}
