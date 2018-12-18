package defrag

// Fragment is a single building block stretch of DNA for assembly
type Fragment struct {
	// ID is a unique identifier for this fragment
	ID string `json:"-"`

	// URL, eg link to a vector's addgene page
	URL string `json:"url,omitempty"`

	// Cost to make the fragment
	Cost float64 `json:"costDollars"`

	// fragment's sequence (linear)
	Seq string `json:"seq,omitempty"`

	// primers necessary to create this (if pcr fragment)
	Primers []Primer `json:"primers,omitempty"`

	// Entry of this fragment In the DB that it came from
	// Used to look for off-targets
	Entry string `json:"-"`

	// Type of this fragment
	Type Type `json:"-"`
}

// Primer is a single Primer used to create a PCR fragment
type Primer struct {
	// Seq of the primer (In 5' to 3' direction)
	Seq string `json:"seq"`

	// Strand of the primer; true if template, false if complement
	Strand bool `json:"strand"`

	// Penalty score
	Penalty float64 `json:"penalty"`

	// PairPenalty score from primer3
	PairPenalty float64 `json:"pairPenalty"`

	// Tm of the primer
	Tm float64 `json:"tm"`

	// GC % max
	GC float64 `json:"gc"`
}
