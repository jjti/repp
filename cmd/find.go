package cmd

import (
<<<<<<< HEAD
	"github.com/jjtimmons/rvec/internal/rvec"
=======
	"github.com/jjtimmons/defrag/internal/defrag"
>>>>>>> 6f3450c4125f17d3ff1136ff8c157a24a9b4d467
	"github.com/spf13/cobra"
)

// findCmd is for finding features or enzymes by their name.
var findCmd = &cobra.Command{
	Use:                        "find",
	Short:                      "Find features or enzymes",
	SuggestionsMinimumDistance: 2,
	Long: `Find features or enzymes by name.
If there is no exact match, similar entries are returned`,
	Aliases: []string{"ls", "list"},
}

// featureFindCmd is for reading features (close to the one requested) from the db.
var featureFindCmd = &cobra.Command{
	Use:                        "feature [name]",
	Short:                      "Find features in the features database",
	Run:                        featureDB.ReadCmd,
	SuggestionsMinimumDistance: 2,
<<<<<<< HEAD
	Example:                    "  rvec features find terminator",
=======
	Example:                    "  defrag features find terminator",
>>>>>>> 6f3450c4125f17d3ff1136ff8c157a24a9b4d467
	Long: `Find features in the features database that are similar to [name].
Writes each feature to the stdout with their name and sequence.
If multiple features contain the feature name sent, each are logged.
Otherwise, all features with names similar to the feature name are writen to stdout`,
	Aliases: []string{"features"},
}

// enzymeFindCmd is for listing out all the available enzymes usable for digesting
// a backbone. Useful for if the user doesn't know which enzymes are available.
var enzymeFindCmd = &cobra.Command{
	Use:                        "enzyme [name]",
	Short:                      "Find enzymes available for linearizing backbones",
	Run:                        enzymeDB.ReadCmd,
	SuggestionsMinimumDistance: 2,
	Long: `List out all the enzymes with the same or a similar a similar name as the argument.

<<<<<<< HEAD
'rvec find enzyme' without any arguments logs all enzymes available.`,
=======
'defrag find enzyme' without any arguments logs all enzymes available.`,
>>>>>>> 6f3450c4125f17d3ff1136ff8c157a24a9b4d467
	Aliases: []string{"enzymes"},
}

// fragmentFindCmd is for finding a fragment by its name
var fragmentFindCmd = &cobra.Command{
	Use:                        "fragment [name]",
	Short:                      "Find a fragment in the databases",
<<<<<<< HEAD
	Run:                        rvec.FragmentFindCmd,
=======
	Run:                        defrag.FragmentFindCmd,
>>>>>>> 6f3450c4125f17d3ff1136ff8c157a24a9b4d467
	SuggestionsMinimumDistance: 2,
	Long:                       `Find a fragment with a given name in the databases requested.`,
}

// sequenceFindCmd is for finding a sequence in the dbs
var sequenceFindCmd = &cobra.Command{
	Use:                        "sequence [name]",
	Short:                      "Find a sequence in the databases",
<<<<<<< HEAD
	Run:                        rvec.SequenceFindCmd,
=======
	Run:                        defrag.SequenceFindCmd,
>>>>>>> 6f3450c4125f17d3ff1136ff8c157a24a9b4d467
	SuggestionsMinimumDistance: 2,
	Long:                       `Find a sequence's matches in the databases requested.`,
	Aliases:                    []string{"seq"},
}

// set flags
func init() {
	fragmentFindCmd.Flags().StringP("dbs", "d", "", "comma separated list of local fragment databases")
	fragmentFindCmd.Flags().BoolP("addgene", "a", false, "use the Addgene repository")
	fragmentFindCmd.Flags().BoolP("igem", "g", false, "use the iGEM repository")
	fragmentFindCmd.Flags().BoolP("dnasu", "u", false, "use the DNASU respository")

	sequenceFindCmd.Flags().StringP("dbs", "d", "", "comma separated list of local fragment databases")
	sequenceFindCmd.Flags().BoolP("addgene", "a", false, "use the Addgene repository")
	sequenceFindCmd.Flags().BoolP("igem", "g", false, "use the iGEM repository")
	sequenceFindCmd.Flags().BoolP("dnasu", "u", false, "use the DNASU respository")
	sequenceFindCmd.Flags().StringP("exclude", "x", "", "keywords for excluding fragments")
	sequenceFindCmd.Flags().IntP("identity", "t", 100, "match %-identity threshold (see 'blastn -help')")

	findCmd.AddCommand(featureFindCmd)
	findCmd.AddCommand(enzymeFindCmd)
	findCmd.AddCommand(fragmentFindCmd)
	findCmd.AddCommand(sequenceFindCmd)

	rootCmd.AddCommand(findCmd)
}
