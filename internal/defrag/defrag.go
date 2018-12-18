package defrag

// Type is the Fragment building type to be used in the assembly
type Type int

const (
	// PCR fragments are those prepared by PCR, often a subselection of their
	// parent vector
	PCR Type = 0

	// Synthetic fragments are those that will be fully synthesized (ex: gBlocks)
	Synthetic Type = 1
)
