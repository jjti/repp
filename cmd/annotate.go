package cmd

import (
	"github.com/jjtimmons/repp/internal/repp"
	"github.com/spf13/cobra"
)

// annotateCmd is for finding features or enzymes by their name.
var annotateCmd = &cobra.Command{
	Use:                        "annotate [seq]",
	Run:                        repp.Annotate,
	Short:                      "Annotate a circular sequence using features (comma-separated)",
	SuggestionsMinimumDistance: 3,
	Long: `Accepts a sequence file as input and runs alignment against the
embedded feature database. Each alignment feature is included as
a feature in the output: a Genbank file. Individual databases
can be selected, in which case the entries in the database will
be used in the alignment _rather_ than the feature database.

The feature database and the default 96% identity are based on
information from [SnapGene](https://www.snapgene.com/resources/plasmid-files/)`,
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
