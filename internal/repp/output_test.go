package repp

import (
	"path/filepath"
	"testing"
)

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
