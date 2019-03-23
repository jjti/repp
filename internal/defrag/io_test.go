package defrag

import (
	"fmt"
	"path"
	"path/filepath"
	"reflect"
	"testing"
)

func Test_inputParser_dbPaths(t *testing.T) {
	parser := inputParser{}

	db1 := "../exampleDir/exampleDB.fa"
	dbAbs1, _ := filepath.Abs(db1)

	db2 := "otherBLASTDir/otherDB.fa"
	dbAbs2, _ := filepath.Abs(db2)

	type args struct {
		dbList string
	}
	tests := []struct {
		name      string
		args      args
		wantPaths []string
		wantError error
	}{
		{
			"single blast path",
			args{
				dbList: db1,
			},
			[]string{dbAbs1},
			nil,
		},
		{
			"multi fasta db paths",
			args{
				dbList: fmt.Sprintf("%s, %s", db1, db2),
			},
			[]string{dbAbs1, dbAbs2},
			nil,
		},
		{
			"empty",
			args{
				dbList: fmt.Sprintf(", %s", db1),
			},
			[]string{dbAbs1},
			nil,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotPaths, _ := parser.dbPaths(tt.args.dbList); !reflect.DeepEqual(gotPaths, tt.wantPaths) {
				t.Errorf("parseDBs() = %v, want %v", gotPaths, tt.wantPaths)
			}
		})
	}
}

func Test_inputParser_parseOut(t *testing.T) {
	parser := inputParser{}

	type args struct {
		in string
	}
	tests := []struct {
		name    string
		args    args
		wantOut string
	}{
		{
			"parse relative path to neighboring output path",
			args{
				in: "./test_file.fa",
			},
			"./test_file.output.json",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotOut := parser.guessOutput(tt.args.in); gotOut != tt.wantOut {
				t.Errorf("parseOut() = %v, want %v", gotOut, tt.wantOut)
			}
		})
	}
}

func Test_inputParser_getFilters(t *testing.T) {
	type args struct {
		filterFlag string
	}
	tests := []struct {
		name string
		p    *inputParser
		args args
		want []string
	}{
		{
			"biobrick separated from year by commas",
			&inputParser{},
			args{
				filterFlag: "tests,BBa_k222000,2004",
			},
			[]string{"TESTS", "BBA_K222000", "2004"},
		},
		{
			"single year",
			&inputParser{},
			args{
				filterFlag: "2004",
			},
			[]string{"2004"},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			p := &inputParser{}
			if got := p.getFilters(tt.args.filterFlag); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("inputParser.getFilters() = %v, want %v", got, tt.want)
			}
		})
	}
}

// Test reading of a FASTA file
func Test_read(t *testing.T) {
	type fileRead struct {
		name         string
		file         string
		fragCount    int
		readFeatures bool
	}

	files := []fileRead{
		{
			"113726(circular)",
			path.Join("..", "..", "test", "input", "113726(circular).parent"),
			1,
			false,
		},
		{
			"multi.fasta",
			path.Join("..", "..", "test", "input", "multi.fasta"),
			5,
			false,
		},
		{
			"genbank sequence",
			path.Join("..", "..", "test", "input", "genbank.gb"),
			1,
			false,
		},
		{
			"genbank features",
			path.Join("..", "..", "test", "input", "genbank.gb"),
			66,
			true,
		},
	}

	for _, f := range files {
		fragments, err := read(f.file, f.readFeatures)

		if err != nil {
			t.Error(err)
		}

		if len(fragments) != f.fragCount {
			t.Errorf("failed to load fragments, len=%d, expected=%d", len(fragments), f.fragCount)
		}

		for _, f := range fragments {
			// ensure we got an ID
			if len(f.ID) < 1 {
				t.Error("failed to load an ID for a Frag from FASTA")
			}

			// ensure we got a Seq
			if len(f.Seq) < 1 {
				t.Errorf("failed to parse a sequence for Frag %s", f.ID)
			}
		}
	}
}

func Test_writeGenbank(t *testing.T) {
	type args struct {
		filename string
		name     string
		seq      string
		frags    []*Frag
		feats    []match
	}
	tests := []struct {
		name string
		args args
	}{
		{
			"include forward and reverse fields",
			args{
				filepath.Join("..", "..", "test", "output", "writeGenbank.gb"),
				"mock part",
				"aattgtgagcggataacaattgacattgtgagcggataacaagatactgagcacatactagagaaagaggagaaatactagatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccttcggctacggcctgcaatgcttcgcccgctaccccgaccacatgaagctgcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagctaccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagaggcctgctgcaaacgacgaaaactacgctttagtagcttaataatactagagtcacactggctcaccttcgggtgggcctttctgcgtttatatactagagagagaatataaaaagccagattattaatccggcttttttattattt",
				[]*Frag{},
				[]match{
					match{
						entry:      "feature 1",
						queryStart: 0,
						queryEnd:   10,
						forward:    true,
					},
					match{
						entry:      "feature 2",
						queryStart: 15,
						queryEnd:   20,
						forward:    false,
					},
				},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			writeGenbank(tt.args.filename, tt.args.name, tt.args.seq, tt.args.frags, tt.args.feats)
		})
	}
}
