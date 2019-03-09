// Package cmd is for command line interactions with the defrag application
package cmd

import (
	"log"

	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
)

var (
	featureDB = defrag.NewFeatureDB()

	enzymeDB = defrag.NewEnzymeDB()
)

// rootCmd represents the base command when called without any subcommands.
var rootCmd = &cobra.Command{
	Use: "defrag",
	Short: `Build vectors using DNA sequences available in public repositories.
Specify vectors using their sequence, features, or fragments`,
	Long:    ``,
	Version: "0.1.0",
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		log.Fatalf("%v", err)
	}
}
