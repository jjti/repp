package defrag

import (
	"fmt"
	"os"
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
			"./test_file.defrag.json",
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

func Test_inputParser_getInput(t *testing.T) {
	// move into the test directory
	parser := inputParser{}
	os.Chdir(filepath.Join("..", "..", "test"))

	tests := []struct {
		name   string
		wantIn string
	}{
		{
			"get fasta file from directory alone",
			"109049.addgene.fa",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotIn := parser.guessInput(); gotIn != tt.wantIn {
				t.Errorf("getInput() = %v, want %v", gotIn, tt.wantIn)
			}
		})
	}
}
