package cmd

import (
	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
)

// FragmentsCmd is for piecing together a list of input fragments into a vector
// and preparing the fragments to assemble into that vector
var FragmentsCmd = &cobra.Command{
	Use:   "fragments",
	Short: "Assemble a list of fragments via Gibson Assembly",
	Long: `Prepare a list of fragments into a vector via Gibson Assembly. Fragments are
	checked for existing homology with their neighbors (in the multi-FASTA input file).
	Fragments lacking the requisite homology are prepared for assembly with PCR.`,
	Run: defrag.Fragments,
}

// set flags
func init() {
	rootCmd.AddCommand(FragmentsCmd)

	// Flags for specifying the paths to the input file, input fragment files, and output file
	FragmentsCmd.Flags().StringP("in", "i", "", "Input file with a list of building fragments <FASTA>")
	FragmentsCmd.Flags().StringP("out", "o", "", "Output file name with fragments and primers <JSON>")
	FragmentsCmd.Flags().StringP("dbs", "d", "", "Comma separated list of databases with building fragments")
	VectorCmd.Flags().BoolP("addgene", "a", false, "Whether to use the Addgene repository as a source of building fragments")
	FragmentsCmd.Flags().StringP("backbone", "b", "",
		`Entry name of a backbone to insert the fragments into. Must match an entry 
	in one of the dbs (either those passed manually or AddGene, iGEM). If an
	enzyme chosen, the backbone will be linearized with that enzyme and
	the largest resulting fragment will be used as the backbone. If no enzyme
	is specified, defrag will chose one nearest the first bp of the backbone with a
	single recognition site`)
	FragmentsCmd.Flags().StringP("enzyme", "e", "",
		`Name of a enzyme to linearize the backbone with (backbone must be specified).
	The enzyme's name must be recognizable by defrag. Use 'defrag enzymes' for a list
	of recognized enzyme names and their recognition sequences`)

	// TODO: let the user run their assembly through `defrag vector` with the passed building fragments
	// FragmentsCmd.Flags().StringP("dbs", "d", "", "Comma separated list of BLAST databases (multi-FASTA)")
}
