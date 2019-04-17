package cmd

import (
	"github.com/spf13/cobra"
)

// deleteCmd is for finding features or enzymes by their name
var deleteCmd = &cobra.Command{
	Use:                        "delete [feature,enzyme]",
	Short:                      "Delete a feature or enzyme",
	SuggestionsMinimumDistance: 2,
	Long:                       `Delete a feature or enzyme by name.`,
	Aliases:                    []string{"rm", "remove"},
}

// featuresDeleteCmd is for deleting features from the feature db
var featuresDeleteCmd = &cobra.Command{
	Use:                        "feature [name]",
	Short:                      "Delete a feature from the features database",
	Run:                        featureDB.DeleteCmd,
	SuggestionsMinimumDistance: 2,
	Aliases:                    []string{"remove"},
<<<<<<< HEAD
	Example:                    "  rvec delete feature \"T7 terminator\"",
=======
	Example:                    "  defrag delete feature \"T7 terminator\"",
>>>>>>> 6f3450c4125f17d3ff1136ff8c157a24a9b4d467
	Long: `Delete a feature from the features database by its name.
If no such feature name exists in the database, an error is logged to stderr.`,
}

// set flags
func init() {
	deleteCmd.AddCommand(featuresDeleteCmd)

	rootCmd.AddCommand(deleteCmd)
}
