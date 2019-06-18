package cmd

import (
	"github.com/jjtimmons/plade/internal/plade"
	"github.com/spf13/cobra"
)

// annotateCmd is for finding features or enzymes by their name.
var annotateCmd = &cobra.Command{
	Use:                        "annotate [seq]",
	Run:                        plade.Annotate,
	Short:                      "Annotate a circular sequence using features (comma-separated)",
	SuggestionsMinimumDistance: 3,
}

// set flags
func init() {
	annotateCmd.Flags().StringP("in", "i", "", "input file name")
	annotateCmd.Flags().StringP("out", "o", "", "output file name")
	annotateCmd.Flags().StringP("exclude", "x", "", "keywords for excluding features")
	annotateCmd.Flags().BoolP("addgene", "a", false, "use the Addgene repository")
	annotateCmd.Flags().BoolP("igem", "g", false, "use the iGEM repository")
	annotateCmd.Flags().BoolP("dnasu", "u", false, "use the DNASU repository")
	annotateCmd.Flags().StringP("dbs", "d", "", "comma separated list of local fragment databases to consider as features")
	annotateCmd.Flags().IntP("identity", "p", 96, "match %-identity threshold (see 'blastn -help')")
	annotateCmd.Flags().BoolP("cull", "c", true, "remove features enclosed in others")
	annotateCmd.Flags().BoolP("names", "n", false, "log feature names to the console")

	RootCmd.AddCommand(annotateCmd)
}
