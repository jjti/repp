package frag

// Fragment is a single potential building block used to
// assemble the target sequence vector
type Fragment struct {
	// a unique identifier
	ID string

	// the fragment's sequence
	Seq string
}
