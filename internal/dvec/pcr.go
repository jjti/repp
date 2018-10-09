package dvec

// Primer is a single Primer used to ampligy a parent fragment
type Primer struct {
	// Seq of the primer (In 5' to 3' direction)
	Seq string

	// Strand of the primer; true if template, false if complement
	Strand bool

	// Start of the fragment (0-index)
	Start int

	// End of the fragment (0-indexed)
	End int

	// Penalty score
	Penalty float32

	// PairPenalty score from primer3
	PairPenalty float32

	// Tm of the primer
	Tm float32

	// GC % max
	GC float32
}

// PCR is a PCR building fragment. It includes primers to create
// the fragment (using its template sequence/vector for addgene as a source)
// and the sequence of the resulting fragment
type PCR struct {
	// ID of this fragment
	ID string

	// Seq of this fragment
	Seq string

	// Entry of this fragment In the DB that it came from
	// Used to look for off-targets
	Entry string

	// Start of the fragment (0-indexed)
	Start int

	// End of the fragment (0-indexed)
	End int

	// Primers necessary to create this PCR Fragment from the template sequence
	Primers []Primer
}
