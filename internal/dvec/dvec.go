package dvec

const (
	// PCR fragments are those prepared by PCR, often a subselection of their
	// parent vector
	PCR Type = 0

	// Synthetic fragments are those that will be fully synthesized (ex: gBlocks)
	Synthetic Type = 1
)

// Type is the Fragment building type to be used in the assembly
type Type int

// Fragment is a single building block stretch of DNA for assembly
type Fragment struct {
	// ID is a unique identifier for this fragment
	ID string `json:"id,omitempty"`

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

// Match is a blast "hit" in the blastdb
type Match struct {
	// Entry of the matched fragment in the database
	Entry string

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

// Primer is a single Primer used to create a PCR fragment
type Primer struct {
	// Seq of the primer (In 5' to 3' direction)
	Seq string `json:"seq"`

	// Strand of the primer; true if template, false if complement
	Strand bool `json:"strand"`

	// Start of the fragment (0-index)
	Start int `json:"-"`

	// End of the fragment (0-indexed)
	End int `json:"-"`

	// Penalty score
	Penalty float32 `json:"penalty"`

	// PairPenalty score from primer3
	PairPenalty float32 `json:"pairPenalty"`

	// Tm of the primer
	Tm float32 `json:"tm"`

	// GC % max
	GC float32 `json:"gc"`
}
