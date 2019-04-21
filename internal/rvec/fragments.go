package rvec

import (
	"fmt"
	"os"
	"strings"
	"testing"

	"github.com/jjtimmons/rvec/config"
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
	if frag.fragType == circular {
		frag.Seq = frag.Seq[:len(frag.Seq)/2]
	}

	fmt.Printf("%s\t%s\n%s\n", name, frag.db, frag.Seq)
}

// FragmentsCmd accepts a cobra commands and assembles a list of building fragments in order
func FragmentsCmd(cmd *cobra.Command, args []string) {
	flags, conf := parseCmdFlags(cmd, args, true)

	// read in the constituent fragments
	var frags []*Frag
	frags, err := read(flags.in, false)
	if err != nil {
		stderr.Fatalln(err)
	}

	// add in the backbone if it was provided
	if flags.backbone.ID != "" {
		frags = append([]*Frag{flags.backbone}, frags...)
	}

	target, solution := fragments(frags, conf)

	// write the single list of fragments as a possible solution to the output file
	writeJSON(
		flags.out,
		target.ID,
		target.Seq,
		[][]*Frag{solution},
		len(target.Seq),
		0,
		flags.backboneMeta,
		conf,
	)

	os.Exit(0)
}

// fragments pieces together a list of fragments into a single vector
// with the fragments in the order and orientation specified
func fragments(frags []*Frag, conf *config.Config) (target *Frag, solution []*Frag) {
	// piece together the adjacent fragments
	if len(frags) < 1 {
		stderr.Fatalln("failed: no fragments to assemble")
	}

	// anneal the fragments together, shift their junctions and create the vector sequence
	vecSeq := annealFragments(conf.FragmentsMinHomology, conf.FragmentsMaxHomology, frags)

	// create the assumed target vector object
	target = &Frag{
		Seq:      vecSeq,
		fragType: circular,
	}

	// create an assembly out of the frags (to fill/convert to fragments with primers)
	a := assembly{frags: frags}
	solution, err := a.fill(target.Seq, conf)
	if err != nil {
		stderr.Fatalln(err)
	}

	return target, solution
}

// annealFragments shifts the start and end of junctions that overlap one another
func annealFragments(min, max int, frags []*Frag) (vec string) {
	// set the start, end, and vector sequence based on that
	//
	// add all of each frags seq to the vector sequence, minus the region overlapping the next
	var vecSeq strings.Builder
	for i, f := range frags {
		next := frags[(i+1)%len(frags)]
		// if we're on the last fragment, mock the first one further along the vector
		if i == len(frags)-1 {
			next = &Frag{
				Seq:   next.Seq,
				start: next.start + vecSeq.Len(),
				end:   next.end + vecSeq.Len(),
			}
		}

		jL := len(f.junction(next, min, max)) // junction length

		contrib := f.Seq[0 : len(f.Seq)-jL] // frag's contribution to vector

		// correct for this Frag's overlap with the next Frag
		f.start = vecSeq.Len()
		f.end = f.start + len(f.Seq) - 1

		// add this Frag's sequence onto the accumulated vector sequence
		vecSeq.WriteString(contrib)
	}

	return vecSeq.String()
}

// validateJunctions checks each fragment and confirms that it has sufficient homology
// with its adjacent fragments and that the match is exact. Largely for testing
func validateJunctions(name string, frags []*Frag, conf *config.Config, t *testing.T) {
	for i, f := range frags {
		next := frags[(i+1)%len(frags)]
		j := f.junction(next, conf.FragmentsMinHomology, conf.FragmentsMaxHomology)
		if j == "" {
			s1 := f.Seq
			if f.PCRSeq != "" {
				s1 = f.PCRSeq
			}

			s2 := next.Seq
			if next.PCRSeq != "" {
				s2 = next.PCRSeq
			}

			t.Errorf("%s -- no junction found between %s and %s\n%s\n\n%s", name, f.ID, next.ID, s1, s2)
		}
	}
}
