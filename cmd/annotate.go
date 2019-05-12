package cmd

import (
	"github.com/jjtimmons/rvec/internal/rvec"
	"github.com/spf13/cobra"
)

// annotateCmd is for finding features or enzymes by their name.
var annotateCmd = &cobra.Command{
	Use:                        "annotate [seq]",
	Run:                        rvec.Annotate,
	Short:                      "Annotate a circular sequence using the feature db",
	SuggestionsMinimumDistance: 3,
}

// set flags
func init() {
	annotateCmd.Flags().StringP("in", "i", "", "input file name")
	annotateCmd.Flags().StringP("out", "o", "", "output file name")
	annotateCmd.Flags().StringP("exclude", "x", "", "keywords for excluding features")
	annotateCmd.Flags().IntP("identity", "p", 96, "match %-identity threshold (see 'blastn -help')")
	annotateCmd.Flags().BoolP("cull", "c", true, "remove features enclosed in others")
	annotateCmd.Flags().BoolP("names", "n", false, "log feature names to the console")

	rootCmd.AddCommand(annotateCmd)
}
