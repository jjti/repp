package cmd

import (
	"github.com/jjtimmons/defrag/internal/defrag"
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
	Aliases: []string{"ls"},
}

// featureFindCmd is for reading features (close to the one requested) from the db
var featureFindCmd = &cobra.Command{
	Use:                        "feature [name]",
	Short:                      "Find features in the features database",
	Run:                        featureDB.ReadCmd,
	SuggestionsMinimumDistance: 2,
	Example:                    "  defrag features read terminator",
	Long: `
Find features in the features database that are similar to [name].
Writes each feature to the stdout with their name and sequence.
If multiple features contain the feature name sent, each are logged.
Otherwise, all features with names similar to the feature name are writen to stdout`,
}

// enzymesCmd is for listing out all the available enzymes usable for digesting
// a backbone. Useful for if the user doesn't know which enzymes are available
var enzymesCmd = &cobra.Command{
	Use:                        "enzymes",
	Short:                      "List enzymes available for linearizing backbones",
	Run:                        defrag.Enzymes,
	SuggestionsMinimumDistance: 4,
	Long: `
Lists out all the enzymes in defrag by name along with their recognition sequence.
[name]	[recognition sequence]`,
}

// set flags
func init() {
	findCmd.AddCommand(featureFindCmd)
	findCmd.AddCommand(enzymesCmd)

	rootCmd.AddCommand(findCmd)
}
