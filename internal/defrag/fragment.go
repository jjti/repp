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
