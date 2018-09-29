// Package frag is for the models and io of Fragments.
// The target and the building pieces are all Fragments
package frag

// Match is a fragment that this fragment matches against
type Match struct {
	// ID of the matched fragment
	ID string

	// Start of the match
	Start int

	// End of the match
	End int
}

// Length returns the length of the match on the target fragment
func (m *Match) Length() int {
	return m.End - m.Start + 1 // it's inclusive
}

// Fragment is a single potential building block used to
// assemble the target sequence vector
type Fragment struct {
	// ID is a unique identifier for this fragment
	//
	// IDs are set before being sent to make FASTA files
	ID string

	// the fragment's sequence
	Seq string

	// slice of matches that this fragment has against others
	Matches []Match
}
