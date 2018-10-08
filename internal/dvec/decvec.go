package dvec

// Match is a fragment that this fragment matches against
type Match struct {
	// ID of the matched fragment in the database
	ID string

	// Seq of the match on the target vector
	Seq string

	// Start of the fragment (0-indexed)
	Start int

	// End of the fragment (0-indexed)
	End int

	// Whether it's a circular fragment (vector, plasmid, etc)
	Circular bool

	// The number of mismatching bps in the match (for primer off-targets)
	Mismatch int
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
}
