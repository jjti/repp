package defrag

import (
	"fmt"
	"strings"

	"github.com/jjtimmons/defrag/config"
	"github.com/spf13/cobra"
)

// FragmentFindCmd logs the building fragment with the name passed.
func FragmentFindCmd(cmd *cobra.Command, args []string) {
	if len(args) < 1 {
		cmd.Help()
		stderr.Fatalln("\nno fragment name passed.")
	}
	name := args[0]

	flags, _ := parseCmdFlags(cmd, args, false)
	frag, err := queryDatabases(name, flags.dbs)
	if err != nil {
		stderr.Fatalln(err)
	}
	fmt.Printf("%s\t%s\t%s\n", name, frag.db, frag.Seq)
}

// FragmentsCmd accepts a cobra commands and assembles a list of building fragments in order
func FragmentsCmd(cmd *cobra.Command, args []string) {
	if _, err := fragments(parseCmdFlags(cmd, args, true)); err != nil {
		stderr.Println(err)
	}
}

// fragments pieces together a list of fragments into a single vector
// with the fragments in the order and orientation specified
func fragments(input *Flags, conf *config.Config) (output []byte, err error) {
	// read in the constituent fragments
	inputFragments, err := read(input.in, false)
	if err != nil {
		return
	}

	// add in the backbone if it was provided
	if input.backbone.ID != "" {
		inputFragments = append([]Frag{input.backbone}, inputFragments...)
	}

	// piece together the adjacent fragments
	target, fragments, err := prepareFragments(inputFragments, conf)
	if err != nil {
		return
	}

	// write the single list of fragments as a possible solution to the output file
	return writeJSON(input.out, target.ID, target.Seq, [][]*Frag{fragments}, len(target.Seq), conf, 0)
}

// prepareFragments takes a list of Fragments and returns the Vector we assume the user is
// trying to build as well as the Fragments (possibly prepared via PCR)
func prepareFragments(targetFrags []Frag, conf *config.Config) (target Frag, solution []*Frag, err error) {
	if len(targetFrags) < 1 {
		return Frag{}, nil, fmt.Errorf("failed: no fragments to assemble")
	}

	// convert the fragments to frags (without a start and end and with the conf)
	frags := make([]*Frag, len(targetFrags))
	for i, f := range targetFrags {
		frags[i] = &Frag{
			ID:       f.ID,
			Seq:      f.Seq,
			fullSeq:  f.Seq,
			conf:     conf,
			start:    0,
			end:      0,
			fragType: existing,
		}
	}

	// find out how much overlap the *last* Frag has with its next one
	// set the start, end, and vector sequence based on that
	//
	// add all of each frags seq to the vector sequence, minus the region overlapping the next
	minHomology := conf.FragmentsMinHomology
	maxHomology := conf.FragmentsMaxHomology
	junction := frags[len(frags)-1].junction(frags[0], minHomology, maxHomology)
	var vectorSeq strings.Builder
	for i, n := range frags {
		// correct for this Frag's overlap with the last Frag
		n.start = vectorSeq.Len() - len(junction)
		n.end = n.start + len(n.Seq) - 1

		// find the junction between this Frag and the next (if there is one)
		junction = n.junction(frags[(i+1)%len(frags)], minHomology, maxHomology)

		// add this Frag's sequence onto the accumulated vector sequence
		vectorSeq.WriteString(n.Seq[0 : len(n.Seq)-len(junction)])
	}

	// create the assumed target vector object
	target = Frag{
		Seq:      vectorSeq.String(),
		fragType: circular,
	}

	// create an assembly out of the frags (to fill/convert to fragments with primers)
	a := assembly{frags: frags}
	solution, err = a.fill(target.Seq, conf)
	if err != nil {
		return Frag{}, nil, err
	}

	return target, solution, nil
}

// ValidateJunctions checks each fragment and confirms that it has sufficient homology
// with its adjacent fragments and that the match is exact. Largely for testing
func ValidateJunctions(frags []*Frag, conf *config.Config) {
	for i, f := range frags {
		next := frags[(i+1)%len(frags)]
		j := f.junction(next, conf.FragmentsMinHomology, conf.FragmentsMaxHomology)
		if j == "" {
			stderr.Fatalf("no junction found between %s and %s", f.ID, next.ID)
		}
	}
}
