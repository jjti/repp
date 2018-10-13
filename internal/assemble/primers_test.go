package assemble

import (
	"path"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/jjtimmons/decvec/internal/dvec"
)

func Test_setPrimers(t *testing.T) {
	db, _ := filepath.Abs(path.Join("..", "..", "test", "blast", "db"))

	conf.PCR.P3MaxPenalty = 10.0
	conf.DB = db

	type args struct {
		p dvec.Fragment
	}
	tests := []struct {
		name        string
		args        args
		wantPrimers []dvec.Primer
		wantErr     bool
	}{
		{
			"create primers basic",
			args{
				p: dvec.Fragment{
					ID:    "gnl|addgene|39412.1",
					Seq:   "TGCTGACTGTGGCGGGTGAGCTTAGGGGGCCTCCGCTCCAGCTCGACACCGGGCAGCTGCTGAAGATCGCGAAGAGAGGGGGAGTAACAGCGGTAGAGGCAGTGCACGCCTGGCGCAATGCGCTCACCGGGGCCCCCTTGAACCTGACCCCAGACCAGGTAGTCGCAATCGCGAACAATAATGGGGGAAAGCAAGCCCTGGAAACCGTGCAAAGGTTGTTGCCGGTCCTTTGTCAAGACCACGGCCTTACACCGGAGCAAGTCGTGGCCATTGCAAGCAATGGGGGTGGCAAACAGGCTCTTGAGACGGTTCAGAGACTTCTCCCAGTTCTCTGTCAAGCCCACGGGCTGACTCCCGATCAAGTTGTAGCGATTGCGTCGCATGACGGAGGGAAACAAGCATTGGAGACTGTCCAACGGCTCCTTCCCGTGTTGTGTCAAGCCCACGGTTTGACGCCTGCACAAGTGGTCGCCATCGCCAGCCATGATGGCGGTAAGCAGGCGCTGGAAACAGTACAGCGCCTGCTGCCTGTACTGTGCCAGGATCATGGACTGACCCCAGACCAGGTAGTCGCAATCGCGAACAATAATGGGGGAAAGCAAGCCCTGGAAACCGTGCAAAGGTTGTTGCCGGTCCTTTGTCAAGACCACGGCCTTACACCGGAGCAAGTCGTGGCCATTGCAAATAATAACGGTGGCAAACAGGCTCTTGAGACGGTTCAGAGACTTCTCCCAGTTCTCTGTCAAGCCCACGGGCTGACTCCCGATCAAGTTGTAGCGATTGCGTCGCATGACGGAGGGAAACAAGCATTGGAGACTGTCCAACGGCTCCTTCCCGTGTTGTGTCAAGCCCACGGTTTGACGCCTGCACAAGTGGTCGCCATCGCCAACAACAACGGCGGTAAGCAGGCGCTGGAAACAGTACAGCGCCTGCTGCCTGTACTGTGCCAGGATCATGGACTGACCCCAGACCAGGTAGTCGCAATCGCGTCGAACATTGGGGGAAAGCAAGCCCTGGAAACCG",
					Entry: "gnl|addgene|39412.1",
				},
			},
			[]dvec.Primer{
				dvec.Primer{},
				dvec.Primer{},
			},
			false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			ps, err := primers(tt.args.p)
			if (err != nil) != tt.wantErr {
				t.Errorf("setPrimers() error = %v, wantErr %v", err, tt.wantErr)
			}
			if !reflect.DeepEqual(ps, tt.wantPrimers) {
				t.Errorf("setPrimers() primers = %v, wantPrimers %v", ps, tt.wantPrimers)
			}
		})
	}
}

func Test_p3exec_input(t *testing.T) {
	type fields struct {
		f   *dvec.Fragment
		in  string
		out string
	}
	tests := []struct {
		name    string
		fields  fields
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			p := &p3exec{
				f:   tt.fields.f,
				in:  tt.fields.in,
				out: tt.fields.out,
			}
			if err := p.input(); (err != nil) != tt.wantErr {
				t.Errorf("p3exec.input() error = %v, wantErr %v", err, tt.wantErr)
			}
		})
	}
}

func Test_p3exec_run(t *testing.T) {
	type fields struct {
		f   *dvec.Fragment
		in  string
		out string
	}
	tests := []struct {
		name    string
		fields  fields
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			p := &p3exec{
				f:   tt.fields.f,
				in:  tt.fields.in,
				out: tt.fields.out,
			}
			if err := p.run(); (err != nil) != tt.wantErr {
				t.Errorf("p3exec.run() error = %v, wantErr %v", err, tt.wantErr)
			}
		})
	}
}

func Test_p3exec_parse(t *testing.T) {
	type fields struct {
		f   *dvec.Fragment
		in  string
		out string
	}
	tests := []struct {
		name    string
		fields  fields
		want    []dvec.Primer
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			p := &p3exec{
				f:   tt.fields.f,
				in:  tt.fields.in,
				out: tt.fields.out,
			}
			got, err := p.parse()
			if (err != nil) != tt.wantErr {
				t.Errorf("p3exec.parse() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("p3exec.parse() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_primers(t *testing.T) {
	type args struct {
		p dvec.Fragment
	}
	tests := []struct {
		name        string
		args        args
		wantPrimers []dvec.Primer
		wantErr     bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotPrimers, err := primers(tt.args.p)
			if (err != nil) != tt.wantErr {
				t.Errorf("primers() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(gotPrimers, tt.wantPrimers) {
				t.Errorf("primers() = %v, want %v", gotPrimers, tt.wantPrimers)
			}
		})
	}
}