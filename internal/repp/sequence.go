package repp

import (
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"
	"text/tabwriter"
	"time"

	"github.com/jjtimmons/repp/config"
	"github.com/spf13/cobra"
)

// SequenceFindCmd is for BLAST'ing a sequence against the dbs and finding matches
func SequenceFindCmd(cmd *cobra.Command, args []string) {
	if len(args) < 1 {
		cmd.Help()
		stderr.Fatalln("\nno sequence passed.")
	}
	seq := args[0]

	flags, _ := parseCmdFlags(cmd, args, false)
	tw := blastWriter()
	matches, err := blast("find_cmd", seq, true, flags.dbs, flags.filters, flags.identity, tw)
	if err != nil {
		stderr.Fatalln(err)
	}

	if len(matches) == 0 {
		stderr.Fatalln("no matches found")
	}

	// sort so the largest matches are first
	sort.Slice(matches, func(i, j int) bool {
		return (matches[i].subjectEnd - matches[i].subjectStart) > (matches[j].queryEnd - matches[j].queryStart)
	})

	// to avoid logging the same matches multiple times
	key := func(m match) string {
		return m.entry + strconv.Itoa(m.subjectStart) + strconv.Itoa(m.subjectEnd)
	}

	seenIds := make(map[string]bool)
	writer := tabwriter.NewWriter(os.Stdout, 0, 4, 3, ' ', 0)
	fmt.Fprintf(writer, "entry\tqstart\tqend\tsstart\tsend\tdatabase\tURL\t\n")
	for _, m := range matches {
		if _, seen := seenIds[key(m)]; seen {
			continue
		}

		if m.subjectEnd-m.subjectStart < 20 {
			continue
		}

		fmt.Fprintf(writer, "%s\t%d\t%d\t%d\t%d\t%s\t%s\n", m.entry, m.queryStart, m.queryEnd, m.subjectStart, m.subjectEnd, m.db, parseURL(m.entry, m.db))
		seenIds[key(m)] = true
	}
	writer.Flush()
}

// SequenceCmd takes a cobra command (with its flags) and runs Vector.
func SequenceCmd(cmd *cobra.Command, args []string) {
	Sequence(parseCmdFlags(cmd, args, true))
}

// Sequence is for running an end to end vector design using a target sequence.
func Sequence(flags *Flags, conf *config.Config) [][]*Frag {
	start := time.Now()

	insert, target, solutions, err := sequence(flags, conf) // build up the assemblies that make the sequence
	if err != nil {
		stderr.Fatalln(err)
	}

	// write the results to a file
	elapsed := time.Since(start)
	_, err = writeJSON(
		flags.out,
		target.ID,
		target.Seq,
		solutions,
		len(insert.Seq),
		elapsed.Seconds(),
		flags.backboneMeta,
		conf,
	)
	if err != nil {
		stderr.Fatalln(err)
	}

	if conf.Verbose {
		fmt.Printf("%s\n\n", elapsed)
	}

	return solutions
}

// sequence builds a vector using a simple cost optimization scheme.
//
// The goal is to find an "optimal" assembly sequence with:
// 	1. the fewest fragments
// 	2. the lowest overall assembly cost ($)
// and, secondarily:
//	3. no duplicate end regions between Gibson fragments
// 	4. no hairpins in the junctions
// 	5. no off-target binding sites in the parent vectors
//	6. low primer3 penalty scores
//
// First build up assemblies, creating all possible assemblies that are
// beneath the upper-bound limit on the number of fragments fully covering
// the target sequence
//
// Then find the pareto optimal solutions that minimize either the cost
// of the assembly or the number of fragments (relative to the other
// assembly plans)
//
// Then, for each pareto optimal solution, traverse the assembly and
// "fill-in" the nodes. Create primers on the Frag if it's a PCR Frag
// or create a sequence to be synthesized if it's a synthetic fragment.
// Error out and repeat the build stage if a Frag fails to be filled
func sequence(input *Flags, conf *config.Config) (insert, target *Frag, solutions [][]*Frag, err error) {
	// read the target sequence (the first in the slice is used)
	fragments, err := read(input.in, false)
	if err != nil {
		return &Frag{}, &Frag{}, nil, fmt.Errorf("failed to read target sequence from %s: %v", input.in, err)
	}

	if len(fragments) > 1 {
		stderr.Printf(
			"warning: %d fragments were in %s. Only targeting the sequence of the first: %s\n",
			len(fragments),
			input.in,
			fragments[0].ID,
		)
	}

	target = fragments[0]
	if conf.Verbose {
		fmt.Printf("Building %s\n", target.ID)
	}

	// if a backbone was specified, add it to the sequence of the target frag
	insert = target.copy() // store a copy for logging later
	if input.backbone.ID != "" {
		target.Seq += input.backbone.Seq
	}

	// get all the matches against the target vector
	tw := blastWriter()
	matches, err := blast(target.ID, target.Seq, true, input.dbs, input.filters, input.identity, tw)
	if conf.Verbose {
		tw.Flush()
	}
	if err != nil {
		dbMessage := strings.Join(input.dbs, ", ")
		return &Frag{}, &Frag{}, nil, fmt.Errorf("failed to blast %s against the dbs %s: %v", target.ID, dbMessage, err)
	}

	// keep only "proper" arcs (non-self-contained)
	matches = cull(matches, len(target.Seq), conf.PCRMinLength)
	if conf.Verbose {
		fmt.Printf("%d matches after culling\n", len(matches)/2)
	}

	// map fragment Matches to nodes
	frags := newFrags(matches, conf)

	if input.backbone.ID != "" {
		// add the backbone in as fragments (copy twice)
		input.backbone.conf = conf
		input.backbone.start = len(insert.Seq)
		input.backbone.end = input.backbone.start + len(input.backbone.Seq) - 1
		input.backbone.uniqueID = "backbone" + strconv.Itoa(input.backbone.start)
		frags = append(frags, input.backbone)

		copiedBB := input.backbone.copy()
		copiedBB.start += len(target.Seq)
		copiedBB.end += len(target.Seq)
		copiedBB.uniqueID = input.backbone.uniqueID
		frags = append(frags, copiedBB)

		sort.Slice(frags, func(i, j int) bool {
			return frags[i].start < frags[j].start
		})
	}

	// build up a slice of assemblies that could, within the upper-limit on
	// fragment count, be assembled to make the target vector
	assemblies := createAssemblies(frags, target.Seq, len(target.Seq), false, conf)

	// build up a map from fragment count to a sorted list of assemblies with that number
	assemblyCounts, countToAssemblies := groupAssembliesByCount(assemblies)

	// fill in pareto optimal assembly solutions
	solutions = fillAssemblies(target.Seq, assemblyCounts, countToAssemblies, conf)

	return insert, target, solutions, nil
}
