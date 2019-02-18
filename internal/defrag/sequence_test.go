package defrag

import (
	"path"
	"reflect"
	"sort"
	"strings"
	"testing"

	"github.com/jjtimmons/defrag/config"
)

// if an input fragment being built is exactly the same as one in a DB, it should be used
// as is and without PCR or any preparation
func Test_vector_single_vector(t *testing.T) {
	fs, c := NewFlags(
		path.Join("..", "..", "test", "input", "109049.addgene.fa"),
		path.Join("..", "..", "test", "output", "109049.output.json"),
		"",
		"",
		"",
		[]string{},
		true,
		false,
	)

	_, _, assemblies, err := sequence(fs, c) // use addgene database
	if err != nil {
		t.Error(err)
	}

	if !strings.Contains(assemblies[0][0].ID, "109049") && !strings.Contains(assemblies[1][0].ID, "109049") {
		t.Fatal("failed to use 109049 to build the vector")
	}
}

func Test_build(t *testing.T) {
	seq := "GTCGACGGATCGGGAGATCTCCCGATCCCCTATGGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATAGTCGACGGATCGGGAGATCTCCCGATCCCCTATGGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATA"
	c := config.New()

	c.SynthesisMaxLength = 5
	c.FragmentsMinHomology = 1
	c.FragmentsMaxCount = 3
	c.PCRMaxEmbedLength = 0

	// ignore cost for now
	c.PCRBPCost = 0
	c.SynthesisFragmentCost = map[int]config.SynthCost{
		100000: {
			Fixed: true,
			Cost:  0.0,
		},
	}

	// input
	n21 := &Frag{
		ID:       "1",
		uniqueID: "1",
		start:    0,
		end:      10,
		conf:     c,
	}
	n22 := &Frag{
		ID:       "2",
		uniqueID: "2",
		start:    5,
		end:      15,
		conf:     c,
	}
	n23 := &Frag{
		ID:       "3",
		uniqueID: "1",
		start:    14,
		end:      30,
		conf:     c,
	}
	n24 := &Frag{
		ID:       "4",
		uniqueID: "2",
		start:    28,
		end:      35,
		conf:     c,
	}

	// output
	a1 := assembly{
		// wrap w/o synthesis
		frags:  []*Frag{n21, n22, n23},
		synths: 0,
	}
	a2 := assembly{
		// wrap w/ synthesis
		frags:  []*Frag{n21, n23},
		synths: 1,
	}
	a3 := assembly{
		// wrap w/o synthsis
		frags:  []*Frag{n22, n23, n24},
		synths: 0,
	}

	type args struct {
		nodes []*Frag
	}
	tests := []struct {
		name           string
		args           args
		wantAssemblies []assembly
	}{
		{
			"test building of assemblies",
			args{
				nodes: []*Frag{n21, n22, n23, n24},
			},
			[]assembly{a1, a2, a3},
		},
	}

	// turn list list of assemblies into a list of Frag ids
	getNodeSet := func(as []assembly) (assemblyIDs []string) {
		for _, a := range as {
			thisID := ""
			for _, n := range a.frags {
				thisID += n.ID
			}
			assemblyIDs = append(assemblyIDs, thisID)
		}

		sort.Strings(assemblyIDs)

		return
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			assemblies := createAssemblies(tt.args.nodes, len(seq), c)

			// concatenate Frag ids together
			actualIds := getNodeSet(assemblies)

			// wanted ID sets
			wantedIds := getNodeSet(tt.wantAssemblies)

			if !reflect.DeepEqual(actualIds, wantedIds) {
				t.Errorf("build() = %v, want %v", actualIds, wantedIds)
			}
		})
	}
}

// Test_sequence is for sequence from end to end
func Test_sequence(test *testing.T) {
	c := config.New()

	type testFlags struct {
		in       string
		out      string
		backbone string
		enzyme   string
		filters  string
		dbs      []string
		addgene  bool
		igem     bool
	}

	tests := []testFlags{
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K2602025.fa"),
			path.Join("..", "..", "test", "output", "BBa_K2602025.json"),
			"pSB1A3",
			"PstI",
			"",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K2779020.fa"),
			path.Join("..", "..", "test", "output", "BBa_K2779020.json"),
			"pSB1A3",
			"PstI",
			"BBa_K277",
			[]string{},
			false,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_E0610.fa"),
			path.Join("..", "..", "test", "output", "BBa_E0610.json"),
			"pSB1C3",
			"EcoRI",
			"2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,BBa_E061",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_I5310.fa"),
			path.Join("..", "..", "test", "output", "BBa_I5310.json"),
			"pSB1C3",
			"EcoRI",
			"",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K077557.fa"),
			path.Join("..", "..", "test", "output", "BBa_K077557.json"),
			"pSB1C3",
			"EcoRI",
			"2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,BBa_K0077",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K2651001.fa"),
			path.Join("..", "..", "test", "output", "BBa_K2651001.json"),
			"pSB1C3",
			"EcoRI",
			"BBa_K265",
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K2779020.fa"),
			path.Join("..", "..", "test", "output", "BBa_K2779020.json"),
			"pSB1A3",
			"PstI",
			"BBa_K277", // no year filters needed
			[]string{},
			true,
			true,
		},
		testFlags{
			path.Join("..", "..", "test", "input", "BBa_K1085023.fa"),
			path.Join("..", "..", "test", "output", "BBa_K1085023.json"),
			"pSB1C3",
			"EcoRI",
			"2009,2010,2011,2012,BBa_K108",
			[]string{},
			true,
			true,
		},
	}

	for _, t := range tests {
		sols := Sequence(NewFlags(t.in, t.out, t.backbone, t.enzyme, t.filters, t.dbs, t.addgene, t.igem))

		if len(sols) < 1 {
			test.Fail()
		}

		for _, s := range sols {
			ValidateJunctions(s, c)
		}
	}
}
