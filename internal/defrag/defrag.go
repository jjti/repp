package defrag

// Type is the Fragment building type to be used in the assembly
type Type int

const (
	// vector is a circular sequence of DNA, e.g.: many of Addgene's plasmids
	vector Type = 0

	// pcr fragments are those prepared by pcr, often a subselection of their parent vector
	pcr Type = 1

	// synthetic fragments are those that will be fully synthesized (ex: gBlocks)
	synthetic Type = 2

	// linear fragment, ie the type of a fragment as it was uploaded submitted and without PCR/synthesis
	existing Type = 3
)
