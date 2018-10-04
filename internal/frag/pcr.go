package frag

// primer is a single primer used to ampligy a parent fragment
type primer struct {
	// seq of the primer (in 5' to 3' direction)
	seq string

	// strand of the primer; true if template, false if complement
	strand bool

	// Start of the fragment (0-index)
	Start int

	// End of the fragment (0-indexed)
	End int

	// TODO: add more primer3 information here, like penalty score etc
}

// PCR is a PCR building fragment. It includes primers to create
// the fragment (using its template sequence/vector for addgene as a source)
// and the sequence of the resulting fragment
type PCR struct {
	Fragment

	// Start of the fragment (0-index)
	Start int

	// End of the fragment (0-indexed)
	End int

	// Primers necessary to create this PCR Fragment from the template sequence
	Primers []primer
}
