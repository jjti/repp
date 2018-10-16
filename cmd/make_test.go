package cmd

import (
	"path"
	"path/filepath"
	"testing"

	"github.com/spf13/cobra"
)

func Test_makeExec(t *testing.T) {
	target, _ := filepath.Abs(path.Join("..", "test", "target.fa"))
	db, _ := filepath.Abs(path.Join("..", "assets", "addgene", "db", "addgene"))

	makeCmd.PersistentFlags().Set("target", target)
	makeCmd.PersistentFlags().Set("db", db)

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
				cmd:  makeCmd,
				args: []string{},
			},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			makeExec(tt.args.cmd, tt.args.args)
		})
	}
}
