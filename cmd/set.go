package cmd

import (
	"github.com/spf13/cobra"
)

// setCmd is for piecing together a list of input fragments into a plasmid
// and preparing the fragments to make into that plasmid
var setCmd = &cobra.Command{
	Use:                        "set [feature,enzyme]",
	Short:                      "Set a feature or enzyme",
	SuggestionsMinimumDistance: 1,
	Long: `Create/update a feature or enzyme with its name and sequence/recognition-site.
Set features can be passed to the 'repp build features' command and enzymes can
be passed to the --enzyme flag`,
	Aliases: []string{"add", "update"},
}

// featureCreateCmd is for adding a new feature to the features db
var featureCreateCmd = &cobra.Command{
	Use:                        "feature [name] [sequence]",
	Short:                      "Add a feature to the features database",
	Run:                        featureDB.SetCmd,
	SuggestionsMinimumDistance: 2,
	Long:                       "\nSet a feature in the features database so it can be use used in 'repp builde features'",
	Aliases:                    []string{"add", "update"},
	Example:                    "  repp set feature \"custom terminator 3\" CTAGCATAACAAGCTTGGGCACCTGTAAACGGGTCTTGAGGGGTTCCATTTTG",
}

// enzymeCreateCmd is for adding a new feature to the features db
var enzymeCreateCmd = &cobra.Command{
	Use:                        "enzyme [name] [sequence]",
	Short:                      "Add an enzyme to the enzymes database",
	Run:                        enzymeDB.SetCmd,
	SuggestionsMinimumDistance: 2,
	Long: `
Set an enzyme in the enzymes database so it can be used to linearize backbones.
Enzymes are passed to the build command, by name, with the --enzyme flag.

Valid recognition sequences have both a cut site in the template sequence: "^" and
a cut site in the complement sequence: "_". Use 'repp ls enzyme' for examples`,
	Aliases: []string{"add", "update"},
	Example: "  repp set enzyme BbvCI CC^TCA_GC",
}

func init() {
	setCmd.AddCommand(featureCreateCmd)
	setCmd.AddCommand(enzymeCreateCmd)

	RootCmd.AddCommand(setCmd)
}
