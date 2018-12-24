package defrag

// Type is the Fragment building type to be used in the assembly
type Type int

const (
	// Vector is a circular sequence of DNA, e.g.: many of Addgene's plasmids
	Vector Type = 0

	// PCR fragments are those prepared by PCR, often a subselection of their parent vector
	PCR Type = 1

	// Synthetic fragments are those that will be fully synthesized (ex: gBlocks)
	Synthetic Type = 2
)
