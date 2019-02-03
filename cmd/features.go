package cmd

import (
	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
)

var featureDB = defrag.NewFeatureDB()

// featuresCmd is for building a vector from its list of contained features. Also
// includes sub-commands for updating the feature database
var featuresCmd = &cobra.Command{
	Use:                        "features [feature name ...]",
	Short:                      "Find or build a vector from its features",
	Run:                        defrag.FeaturesCmd, // TODO
	SuggestionsMinimumDistance: 4,
	Long:                       "",
	Example:                    "  defrag features create \"custom terminator 3\" CTAGCATAACAAGCTTGGGCACCTGTAAACGGGTCTTGAGGGGTTCCATTTTG",
}

// featuresCreateCmd is for adding a new feature to the features db
var featuresCreateCmd = &cobra.Command{
	Use:                        "create [feature name] [feature sequence]",
	Short:                      "Add a feature to the features database",
	Run:                        featureDB.CreateCmd,
	SuggestionsMinimumDistance: 2,
	Long:                       "\nAdd a new feature to the features database so it can be use used in 'defrag features'",
	Aliases:                    []string{"add"},
	Example:                    "  defrag features create \"custom terminator 3\" CTAGCATAACAAGCTTGGGCACCTGTAAACGGGTCTTGAGGGGTTCCATTTTG",
}

// featuresReadCmd is for reading features (close to the one requested) from the db
var featuresReadCmd = &cobra.Command{
	Use:                        "read [feature name]",
	Short:                      "Find features in the features database",
	Run:                        featureDB.ReadCmd,
	SuggestionsMinimumDistance: 2,
	Aliases:                    []string{"find"},
	Example:                    "  defrag features read terminator",
	Long: `
Find features in the features database that are similar to [feature name].
Writes each feature to the stdout with their name and sequence.
If multiple features contain the feature name sent, each are logged.
Otherwise, all features with names similar to the feature name are writen to stdout`,
}

// featuresUpdateCmd is for updating features in the feature db
var featuresUpdateCmd = &cobra.Command{
	Use:                        "update [feature name] [feature sequence]",
	Short:                      "Update a feature's sequence in the features database",
	Run:                        featureDB.UpdateCmd,
	SuggestionsMinimumDistance: 2,
	Aliases:                    []string{"change"},
	Example:                    "  defrag features update \"T7 terminator\" CTAGCATAACCCCTTGGGGCCTGTAAACGGGTCTTGAGGGGTTTTTTG",
	Long: `
Update a feature's sequence with the new one passed.
If no such feature name exists in the database, an error is logged to stdout.`,
}

// featuresDeleteCmd is for deleting features from the feature db
var featuresDeleteCmd = &cobra.Command{
	Use:                        "delete [feature name]",
	Short:                      "Delete a feature from the features database",
	Run:                        featureDB.DeleteCmd,
	SuggestionsMinimumDistance: 2,
	Aliases:                    []string{"remove"},
	Example:                    "  defrag features delete \"T7 terminator\"",
	Long: `
Delete a feature from the features database by its name.
If no such feature name exists in the database, an error is logged to stdout.`,
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
