package cmd

import (
	"path"
	"path/filepath"
	"testing"

	"github.com/spf13/cobra"
)

func Test_makeExec(t *testing.T) {
	target, _ := filepath.Abs(path.Join("..", "test", "target.fa"))
	db, _ := filepath.Abs(path.Join("..", "test", "mockDB", "mockDB"))
	output, _ := filepath.Abs(path.Join("..", "bin", "test.output.json"))

	MakeCmd.PersistentFlags().Set("target", target)
	MakeCmd.PersistentFlags().Set("dbs", db)
	MakeCmd.PersistentFlags().Set("out", output)
	MakeCmd.PersistentFlags().Set("addgene", "true")

	type args struct {
		cmd  *cobra.Command
		args []string
	}
	tests := []struct {
		name string
		args args
	}{
		{
			"end to end test",
			args{
				cmd:  MakeCmd,
				args: []string{},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if err := MakeCmd.Execute(); err != nil {
				t.Fatal(err)
			}
		})
	}

	// t.Fail()
}
