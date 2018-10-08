package pcr

import (
	"fmt"
	"log"
	"path"

	"github.com/jjtimmons/decvec/config"
	"github.com/jjtimmons/decvec/internal/blast"
)

// Primer is a single Primer used to ampligy a parent fragment
type Primer struct {
	// seq of the primer (In 5' to 3' direction)
	seq string

	// strand of the primer; true if template, false if complement
	strand bool

	// Start of the fragment (0-index)
	Start int

	// End of the fragment (0-indexed)
	End int

	// Penalty score
	Penalty float32

	// PairPenalty score from primer3
	PairPenalty float32

	// Tm of the primer
	Tm float32

	// GC % max
	GC float32
}

// PCR is a PCR building fragment. It includes primers to create
// the fragment (using its template sequence/vector for addgene as a source)
// and the sequence of the resulting fragment
type PCR struct {
	// ID of this fragment
	ID string

	// Seq of this fragment
	Seq string

	// Entry of this fragment In the DB that it came from
	// Used to look for off-targets
	Entry string

	// Start of the fragment (0-indexed)
	Start int

	// End of the fragment (0-indexed)
	End int

	// Primers necessary to create this PCR Fragment from the template sequence
	Primers []Primer
}

// SetPrimers creates primers on a PCR fragment and returns an error if
//	1. the primers have an unacceptably high primer3 penalty score
//	2. there are off-targets in the primers
func (p *PCR) SetPrimers() error {
	maxPairP := config.NewConfig().PCR.P3MaxPenalty

	handleP3 := func(err error) {
		// we should fail totally on any primer3 errors -- shouldn't happen
		if err != nil {
			log.Panic(err)
		}
	}

	exec := p3exec{
		Frag: p,
		In:   path.Join(p3Dir, p.ID+".in"),
		Out:  path.Join(p3Dir, p.ID+".out"),
	}

	// make input file
	err := exec.input()
	handleP3(err)

	// execute
	err = exec.run()
	handleP3(err)

	// parse the results into primers for storing on the fragment
	primers, err := exec.parse()
	handleP3(err)
	p.Primers = primers

	// 1. check for whether the primers have too have a pair penalty score
	if p.Primers[0].PairPenalty > maxPairP {
		return fmt.Errorf(
			"primers have pair primer3 penalty score of %f, should be less than %f:\n%+v\n%+v",
			p.Primers[0].PairPenalty,
			maxPairP,
			p.Primers[0],
			p.Primers[1],
		)
	}

	// 2. check for whether either of the primers have an off-target/mismatch
	for _, primer := range p.Primers {
		mismatchExists, mismatch, err := blast.Mismatch(primer.seq, p.Entry)
		handleP3(err) // shouldn't be erroring here either
		if mismatchExists {
			return fmt.Errorf(
				"found a mismatching sequence, %s, against the primer %s",
				mismatch.Seq,
				primer.seq,
			)
		}
	}

	return nil
}
