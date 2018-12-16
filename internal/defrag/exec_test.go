package defrag

import (
	"fmt"
	"path"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/spf13/cobra"
)

func TestExecute(t *testing.T) {
	in, _ := filepath.Abs(path.Join("..", "..", "test", "target.fa"))
	out, _ := filepath.Abs(path.Join("..", "..", "bin", "test_output.json"))
	// db, _ := filepath.Abs(path.Join("..", "..", "test", "mockDB", "mockDB"))

	// https://stackoverflow.com/a/50880663
	cmd := &cobra.Command{}
	cmd.Flags().String("in", in, "")
	cmd.Flags().String("out", out, "")
	cmd.Flags().String("dbs", "", "")
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
			Execute(tt.args.cmd, tt.args.args)
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
