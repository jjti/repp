package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

// makeCmd represents the make command
var makeCmd = &cobra.Command{
	Use:   "make",
	Short: "Make a vector from its target sequence and a list of possible consituent fragments",
	Long:  ``,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("make called")
	},
}

func init() {
	rootCmd.AddCommand(makeCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// makeCmd.PersistentFlags().String("foo", "", "A help for foo")

	makeCmd.Flags().String("target", "t", "target sequence to make")
	makeCmd.Flags().String("input-fasta", "f", "multi-fasta input file with building sequences")
}
