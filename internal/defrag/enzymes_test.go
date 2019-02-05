package defrag

import (
	"reflect"
	"testing"
)

func Test_recogRegex(t *testing.T) {
	type args struct {
		recog string
	}
	tests := []struct {
		name        string
		args        args
		wantDecoded string
	}{
		{
			"decode PpuMI: RGGWCCY",
			args{
				recog: "RGGWCCY",
			},
			"(A|G)GG(A|T)CC(C|T)",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotDecoded := recogRegex(tt.args.recog); gotDecoded != tt.wantDecoded {
				t.Errorf("recogRegex() = %v, want %v", gotDecoded, tt.wantDecoded)
			}
		})
	}

	// should be able to decode every recognition site without failing
	for _, enz := range NewEnzymeDB().enzymes {
		recogRegex(newEnzyme(enz).recog)
	}
}

func Test_digest(t *testing.T) {
	type args struct {
		frag Frag
		enz  enzyme
	}
	tests := []struct {
		name         string
		args         args
		wantDigested Frag
		wantErr      bool
	}{
		{
			"fail with no recognition sequence",
			args{
				Frag{
					Seq: "ATGAGGTTAGCCAAAAAAGCACGTGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				},
				enzyme{recog: "GAATTC", hangInd: 5, cutInd: 1},
			},
			Frag{},
			true,
		},
		{
			"digest in template sequence, no overhang",
			args{
				Frag{
					Seq: "ATGAGGTTAGCCAAAAAAGCACGTGAATTCGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				},
				enzyme{recog: "GAATTC", hangInd: 3, cutInd: 3},
			},
			Frag{
				Seq: "TTCGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCAATGAGGTTAGCCAAAAAAGCACGTGAA",
			},
			false,
		},
		{
			"digest in reverse complement sequence, no overhang",
			args{
				Frag{
					Seq: "ATGAGGTTAGCCAAAAAAGCACGCTGGGGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				},
				enzyme{recog: "CCCAGC", hangInd: 3, cutInd: 3}, // rev comp GCTGGG
			},
			Frag{
				Seq: "GGGGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCAATGAGGTTAGCCAAAAAAGCACGCT",
			},
			false,
		},
		{
			"digest in template sequence, positive overhang",
			args{
				Frag{
					Seq: "ATGAGGTTAGCCAAAAAAGCACGTGAATTCGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				},
				enzyme{recog: "GAATTC", hangInd: 1, cutInd: 5},
			},
			Frag{
				Seq: "CGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCAATGAGGTTAGCCAAAAAAGCACGTGAATT",
			},
			false,
		},
		{
			"digest in template sequence, negative overhang",
			args{
				Frag{
					Seq: "ATGAGGTTAGCCAAAAAAGCACGTGAATTCGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				},
				enzyme{recog: "GAATTC", hangInd: 5, cutInd: 1},
			},
			Frag{
				Seq: "CGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCAATGAGGTTAGCCAAAAAAGCACGTG",
			},
			false,
		},
		{
			"digest in reverse complement sequence, positive overhang",
			args{
				Frag{
					Seq: "ATGAGGTTAGCCAAAAAAGCACGTGCTGGGGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				},
				enzyme{recog: "CCCAGC", hangInd: 1, cutInd: 5}, // rev comp = GCTGGG
			},
			Frag{
				Seq: "GGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCAATGAGGTTAGCCAAAAAAGCACGTGCTGG",
			},
			false,
		},
		{
			"digest in reverse complement sequence, negative overhang",
			args{
				Frag{
					Seq: "ATGAGGTTAGCCAAAAAAGCACGTGCTGGGGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCA",
				},
				enzyme{recog: "CCCAGC", hangInd: 5, cutInd: 1}, // rev comp = GCTGGG
			},
			Frag{
				Seq: "GGGTGGCGCCCACCGACTGTTCCCAAACTGTAGCTCTTCGTTCCGTCAAGGCCCGACTTTCATCGCGGCCCATTCCAATGAGGTTAGCCAAAAAAGCACGTG",
			},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotDigested, err := digest(tt.args.frag, tt.args.enz)
			if (err != nil) != tt.wantErr {
				t.Errorf("digest() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(gotDigested, tt.wantDigested) {
				t.Errorf("digest() = %v, want %v", gotDigested, tt.wantDigested)
			}
		})
	}
}
