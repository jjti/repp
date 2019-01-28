package defrag

// Payload is the JSON shape of an input payload from the API endpoint
type Payload struct {
	// Target sequence that we want to build
	Target string `json:"target"`

	// Fragments that should be assembled
	Fragments []Frag `json:"fragments"`

	// Addgene is whether to use the Addgene repository
	Addgene bool `json:"addgene"`

	// IGEM is whether to use the IGEM repository
	IGEM bool `json:"igem"`

	// Backbone to insert the inserts into
	Backbone string `json:"backbone"`

	// Enzyme to cut the backbone with (name)
	Enzyme string `json:"enzyme"`
}

// Solution is a single solution to build up the target vector
type Solution struct {
	// Count is the number of fragments in this solution
	Count int `json:"count"`

	// Cost estimated from the primer and sequence lengths
	Cost float64 `json:"dollars"`

	// Fragments used to build this solution
	Fragments []*Frag `json:"fragments"`
}

// Output is a struct containing design results for the assembly
type Output struct {
	// Time, ex:
	// "2006-01-02 15:04:05.999999999 -0700 MST"
	// https://golang.org/pkg/time/#Time.String
	Time string `json:"time"`

	// Target sequence
	Target string `json:"target"`

	// Solutions builds
	Solutions []Solution `json:"solutions"`

	// SynthInsertCost cost of a full synthesis
	SynthInsertCost float64 `json:"synthInsertCost"`
}
