package cmd

import (
	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
)

var featureDB = defrag.NewFeatureDB()

// featuresCmd is for building a vector from its list of contained features. Also
// includes sub-commands for updating the feature database
var featuresCmd = &cobra.Command{
	Use:                        "features",
	Short:                      "Find or build a vector from its features",
	Run:                        defrag.FragmentsCmd, // TODO
	SuggestionsMinimumDistance: 4,
	Long:                       "",
}

// featuresCreateCmd is for adding a new feature to the features db
var featuresCreateCmd = &cobra.Command{
	Use:                        "create [feature name] [feature sequence]",
	Short:                      "add a feature to the features database",
	Run:                        featureDB.Create, // TODO
	SuggestionsMinimumDistance: 2,
	Long:                       "",
	Aliases:                    []string{"add", "make", "new"},
	Example:                    "  defrag features create \"custom terminator 3\" CTAGCATAACAAGCTTGGGCACCTGTAAACGGGTCTTGAGGGGTTCCATTTTG",
}

// featuresReadCmd is for reading features (close to the one requested) from the db
var featuresReadCmd = &cobra.Command{
	Use:                        "read [feature name]",
	Short:                      "find a feature in the db by its name (returns it and similar features)",
	Run:                        featureDB.Read, // TODO
	SuggestionsMinimumDistance: 2,
	Long:                       "",
	Aliases:                    []string{"find"},
	Example:                    "  defrag features read terminator",
}

// featuresUpdateCmd is for updating features in the feature db
var featuresUpdateCmd = &cobra.Command{
	Use:                        "update [feature name] [feature sequence]",
	Short:                      "update a feature's sequence in the features database",
	Run:                        featureDB.Update, // TODO
	SuggestionsMinimumDistance: 2,
	Long:                       "",
	Aliases:                    []string{"change"},
	Example:                    "  defrag features update \"T7 terminator\" CTAGCATAACCCCTTGGGGCCTGTAAACGGGTCTTGAGGGGTTTTTTG",
}

// featuresDeleteCmd is for deleting features from the feature db
var featuresDeleteCmd = &cobra.Command{
	Use:                        "delete [feature name]",
	Short:                      "delete a feature from the features database",
	Run:                        featureDB.Delete, // TODO
	SuggestionsMinimumDistance: 2,
	Long:                       "",
	Aliases:                    []string{"remove"},
	Example:                    "  defrag features delete \"T7 terminator\"",
}

// set flags
func init() {
	featuresCmd.AddCommand(featuresCreateCmd)
	featuresCmd.AddCommand(featuresReadCmd)
	featuresCmd.AddCommand(featuresUpdateCmd)
	featuresCmd.AddCommand(featuresDeleteCmd)

	rootCmd.AddCommand(featuresCmd)

	// TODO: let the user run their assembly through `defrag vector` with the passed building fragments
	// FeaturesCmd.Flags().StringP("dbs", "d", "", "Comma separated list of BLAST databases (multi-FASTA)")
}
