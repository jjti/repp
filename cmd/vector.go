package cmd

import (
	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
)

// vectorCmd is for assembling a vector (single circular sequence) after finding
// building fragments to piece it together efficiently
var vectorCmd = &cobra.Command{
	Use:                        "vector",
	Short:                      "Build a target vector from local and/or remote fragments",
	Run:                        defrag.VectorCmd,
	SuggestionsMinimumDistance: 4,
	Long: `
Build up a vector from its target sequence using a combination of existing and
synthesized fragments.

Solutions have either a minimum fragment count or assembly cost (or both).`,
}

// set flags
func init() {
	rootCmd.AddCommand(vectorCmd)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	vectorCmd.Flags().StringP("in", "i", "", "input FASTA with target sequence")
	vectorCmd.Flags().StringP("out", "o", "", "output file name")
	vectorCmd.Flags().StringP("dbs", "d", "", "delimited list of local fragment databases")
	vectorCmd.Flags().BoolP("addgene", "a", false, "use the Addgene repository")
	vectorCmd.Flags().BoolP("igem", "g", false, "use the iGEM repository")
	vectorCmd.Flags().StringP("backbone", "b", "", backboneHelp)
	vectorCmd.Flags().StringP("enzyme", "e", "", enzymeHelp)
	vectorCmd.Flags().StringP("filter", "f", "", "delimited keywords for removing fragments")
	vectorCmd.Flags().Float64P("identity", "t", 100.0, "match %-identity threshold (see 'blastn -help')")
}
