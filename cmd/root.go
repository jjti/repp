<<<<<<< HEAD
// Package cmd is for command line interactions with the rvec application
=======
// Package cmd is for command line interactions with the defrag application
>>>>>>> 6f3450c4125f17d3ff1136ff8c157a24a9b4d467
package cmd

import (
	"log"

<<<<<<< HEAD
	"github.com/jjtimmons/rvec/internal/rvec"
=======
	"github.com/jjtimmons/defrag/internal/defrag"
>>>>>>> 6f3450c4125f17d3ff1136ff8c157a24a9b4d467
	"github.com/spf13/cobra"
)

var (
<<<<<<< HEAD
	featureDB = rvec.NewFeatureDB()

	enzymeDB = rvec.NewEnzymeDB()
=======
	featureDB = defrag.NewFeatureDB()

	enzymeDB = defrag.NewEnzymeDB()
>>>>>>> 6f3450c4125f17d3ff1136ff8c157a24a9b4d467
)

// rootCmd represents the base command when called without any subcommands.
var rootCmd = &cobra.Command{
<<<<<<< HEAD
	Use: "rvec",
=======
	Use: "defrag",
>>>>>>> 6f3450c4125f17d3ff1136ff8c157a24a9b4d467
	Short: `Build vectors using DNA sequences available in public repositories.
Specify vectors using their sequence, features, or fragments`,
	Version: "0.1.0",
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		log.Fatalf("%v", err)
	}
}
