package defrag

import (
	"fmt"
	"os"
	"path"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/spf13/cobra"
)

func Test_Vector(t *testing.T) {
	in, _ := filepath.Abs(path.Join("..", "..", "test", "109049.addgene.fa"))
	out, _ := filepath.Abs(path.Join("..", "..", "bin", "109049.addgene.json"))
	dbs, _ := filepath.Abs(path.Join("..", "..", "test", "db", "db"))

	// https://stackoverflow.com/a/50880663
	cmd := &cobra.Command{}
	cmd.Flags().String("in", in, "")
	cmd.Flags().String("out", out, "")
	cmd.Flags().String("dbs", dbs, "")
	cmd.Flags().Bool("addgene", true, "")

	type args struct {
		cmd  *cobra.Command
		args []string
	}
	tests := []struct {
		name string
		args args
	}{
		{
			"end to end design of a test vector",
			args{
				cmd: cmd,
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			Vector(tt.args.cmd, tt.args.args)
		})
	}
}

func Test_parseDBs(t *testing.T) {
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
			if gotPaths, _ := parseDBs(tt.args.dbList); !reflect.DeepEqual(gotPaths, tt.wantPaths) {
				t.Errorf("parseDBs() = %v, want %v", gotPaths, tt.wantPaths)
			}
		})
	}
}

// if an input fragment being built is exactly the same as one in a DB, it should be used
// as is and without PCR or any preparation
func Test_execute_single_vector(t *testing.T) {
	in := path.Join("..", "..", "test", "109049.addgene.fa")
	out := path.Join("..", "..", "bin", "109049.output.json")

	assemblies := execVector(in, out, "", true) // use addgene database

	if len(assemblies) != 1 {
		t.Fatal("failed to return the pareto optimal solution: 109049 alone")
	}

	if len(assemblies[0]) != 1 || !strings.Contains(assemblies[0][0].ID, "109049") {
		t.Fatal("failed to use 109049 to build the vector")
	}

	if assemblies[0][0].Type != vector {
		t.Fatalf("failed to recognize 109049 as a Type.Vector, was %d", assemblies[0][0].Type)
	}
}

func Test_parseOut(t *testing.T) {
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
			if gotOut := guessOutput(tt.args.in); gotOut != tt.wantOut {
				t.Errorf("parseOut() = %v, want %v", gotOut, tt.wantOut)
			}
		})
	}
}

func Test_getInput(t *testing.T) {
	// move into the test directory
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
			if gotIn := guessInput(); gotIn != tt.wantIn {
				t.Errorf("getInput() = %v, want %v", gotIn, tt.wantIn)
			}
		})
	}
}
