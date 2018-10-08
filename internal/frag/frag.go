// Package frag is for the models and io of Fragments.
// The target and the building sequences are all Fragments
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
