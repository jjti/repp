package defrag

import (
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/jjtimmons/defrag/config"
)

func TestNewFeatureDB(t *testing.T) {
	db := NewFeatureDB()

	if len(db.features) < 1 {
		t.Fail()
	}
}

func Test_features(t *testing.T) {
	test1, conf := NewFlags(
		"p10 promoter, mEGFP, T7 terminator",
		filepath.Join("..", "..", "test", "output", "features.json"),
		"pSB1A3",
		"EcoRI",
		"",
		[]string{},
		true,
		true,
	)

	type args struct {
		flags *Flags
		conf  *config.Config
	}
	tests := []struct {
		name string
		args args
	}{
		{
			"test end to end features creation",
			args{
				flags: test1,
				conf:  conf,
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			features(tt.args.flags, tt.args.conf)

			t.Fail()
		})
	}
}

// GACCTTTAATTCAACCCAACACAATATATTATAGTTAAATAAGAATTATTATCAAATCATTTGTATATTAATTAAAATACTATACTGTAAATTACATTTTATTTACAATC

func Test_queryFeatures(t *testing.T) {
	type args struct {
		flags *Flags
	}
	tests := []struct {
		name string
		args args
		want [][]string
	}{
		{
			"gather SV40 origin, p10 promoter, mEGFP",
			args{
				&Flags{
					in:  "SV40 origin,p10 promoter,mEGFP",
					dbs: []string{config.AddgeneDB, config.IGEMDB},
				},
			},
			[][]string{
				[]string{"SV40 origin", "ATCCCGCCCCTAACTCCGCCCAGTTCCGCCCATTCTCCGCCCCATGGCTGACTAATTTTTTTTATTTATGCAGAGGCCGAGGCCGCCTCGGCCTCTGAGCTATTCCAGAAGTAGTGAGGAGGCTTTTTTGGAGGCC"},
				[]string{"p10 promoter", "GACCTTTAATTCAACCCAACACAATATATTATAGTTAAATAAGAATTATTATCAAATCATTTGTATATTAATTAAAATACTATACTGTAAATTACATTTTATTTACAATC"},
				[]string{"mEGFP", "AGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGCGCGGCGAGGGCGAGGGCGATGCCACCAACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTCCTTCAAGGACGACGGCACCTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTTCAACAGCCACAACGTCTATATCACGGCCGACAAGCAGAAGAACGGCATCAAGGCGAACTTCAAGATCCGCCACAACGTCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCAAGCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAG"},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got, _ := queryFeatures(tt.args.flags); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("queryFeatures() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_blastFeatures(t *testing.T) {
	type args struct {
		flags          *Flags
		targetFeatures [][]string
	}
	tests := []struct {
		name string
		args args
	}{
		{
			"blast a feature against the part databases",
			args{
				flags: &Flags{
					dbs:      []string{config.AddgeneDB, config.IGEMDB},
					filters:  []string{},
					identity: 100.0,
				},
				targetFeatures: [][]string{
					[]string{"SV40 origin", "ATCCCGCCCCTAACTCCGCCCAGTTCCGCCCATTCTCCGCCCCATGGCTGACTAATTTTTTTTATTTATGCAGAGGCCGAGGCCGCCTCGGCCTCTGAGCTATTCCAGAAGTAGTGAGGAGGCTTTTTTGGAGGCC"},
				},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := blastFeatures(tt.args.flags, tt.args.targetFeatures)

			matches := []match{}
			for _, ms := range got {
				for _, m := range ms {
					matches = append(matches, m.match)
				}
			}

			// confirm that the returned fragments sequences contain at least the full queried sequence
			for _, m := range matches {
				containsTargetSeq := false
				for _, wantedSeq := range tt.args.targetFeatures {
					if strings.Contains(m.seq, wantedSeq[1]) {
						containsTargetSeq = true
					}
				}

				if !containsTargetSeq {
					t.Fatalf("match with seq %s doesn't contain any of the target features", m.seq)
				}
			}
		})
	}
}
