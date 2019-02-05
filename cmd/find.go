package cmd

import (
	"github.com/spf13/cobra"
)

// findCmd is for finding features or enzymes by their name
var findCmd = &cobra.Command{
	Use:                        "find [feature,enzyme]",
	Short:                      "Find features or enzymes",
	SuggestionsMinimumDistance: 2,
	Long: `
Find features or enzymes by name.
If there is no exact match, similar entries are returned`,
	Aliases: []string{"ls", "list"},
}

// featureFindCmd is for reading features (close to the one requested) from the db
var featureFindCmd = &cobra.Command{
	Use:                        "feature [name]",
	Short:                      "Find features in the features database",
	Run:                        featureDB.ReadCmd,
	SuggestionsMinimumDistance: 2,
	Example:                    "  defrag features find terminator",
	Long: `
Find features in the features database that are similar to [name].
Writes each feature to the stdout with their name and sequence.
If multiple features contain the feature name sent, each are logged.
Otherwise, all features with names similar to the feature name are writen to stdout`,
	Aliases: []string{"features"},
}

// enzymeFindCmd is for listing out all the available enzymes usable for digesting
// a backbone. Useful for if the user doesn't know which enzymes are available
var enzymeFindCmd = &cobra.Command{
	Use:                        "enzyme [name]",
	Short:                      "Find enzymes available for linearizing backbones by name",
	Run:                        enzymeDB.ReadCmd,
	SuggestionsMinimumDistance: 2,
	Long: `
Lists out all the enzymes with the same or a similar a similar name as the argument.

'defrag find enzyme' without any arguments logs all enzymes available.`,
	Aliases: []string{"enzymes"},
}

// set flags
func init() {
	findCmd.AddCommand(featureFindCmd)
	findCmd.AddCommand(enzymeFindCmd)

	rootCmd.AddCommand(findCmd)
}
