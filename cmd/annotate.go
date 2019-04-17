package cmd

import (
<<<<<<< HEAD
	"github.com/jjtimmons/rvec/internal/rvec"
=======
	"github.com/jjtimmons/defrag/internal/defrag"
>>>>>>> 6f3450c4125f17d3ff1136ff8c157a24a9b4d467
	"github.com/spf13/cobra"
)

// annotateCmd is for finding features or enzymes by their name.
var annotateCmd = &cobra.Command{
	Use:                        "annotate [seq]",
<<<<<<< HEAD
	Run:                        rvec.Annotate,
=======
	Run:                        defrag.Annotate,
>>>>>>> 6f3450c4125f17d3ff1136ff8c157a24a9b4d467
	Short:                      "Annotate a circular sequence using the feature db",
	SuggestionsMinimumDistance: 3,
}

// set flags
func init() {
	annotateCmd.Flags().StringP("in", "i", "", "input file name")
	annotateCmd.Flags().StringP("out", "o", "", "output file name")
	annotateCmd.Flags().StringP("exclude", "x", "", "keywords for excluding features")
	annotateCmd.Flags().IntP("identity", "t", 100, "match %-identity threshold (see 'blastn -help')")
	annotateCmd.Flags().BoolP("enclosed", "c", false, "annotate with features enclosed in others")

	rootCmd.AddCommand(annotateCmd)
}
