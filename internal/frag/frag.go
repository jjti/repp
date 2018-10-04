// Package frag is for the models and io of Fragments.
// The target and the building pieces are all Fragments
package frag

// Fragment is a single potential building block used to
// assemble the target sequence vector
type Fragment struct {
	// ID is a unique identifier for this fragment
	//
	// IDs are set before being sent to make FASTA files
	ID string

	// the fragment's sequence
	Seq string
}

// Match is a fragment that this fragment matches against
type Match struct {
	// ID of the matched fragment
	ID string

	// Start of the match
	Start int

	// End of the match
	End int

	// Whether it's a circular fragment
	Circular bool

	// The sequence of the template/parent fragment
	Template string
}

// Length returns the length of the match on the target fragment
func (m *Match) Length() int {
	return m.End - m.Start + 1 // it's inclusive
}
