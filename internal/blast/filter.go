package blast

import (
	"sort"

	"github.com/jjtimmons/decvec/internal/frag"
)

// filter is for "propertizing" the matches from BLAST and removing
// those that are either 1) "truncated" or 2) overlap with the edges of the
// basted target sequence
//
// propertizing fragment matches means removing those that are completely
// self-contained in other fragments: the larger of the available fragments
// will be the better one, since it covers a greater region
func filter(f *frag.Fragment) {
	// we doubled target sequence before blast to cover the matches across
	// the "zero-index"
	lastbp := 2 * len(f.Seq)

	// remove matches that are truncated
	// ie, remove matches that start at 0 or end at the last bp
	var newMatches []frag.Match
	for _, m := range f.Matches {
		if m.Start > 0 && m.End < lastbp {
			newMatches = append(newMatches, m)
		}
	}
	f.Matches = newMatches

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

	f.Matches = properMatches
}
