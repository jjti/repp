package cmd

import (
	"log"

	"github.com/jjtimmons/defrag/internal/defrag"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

var featureDB = defrag.NewFeatureDB()

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "defrag",
	Short: "Build vectors using their target sequence, constiuent features, or fragments",
	Long:  ``,
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		log.Fatalf("failed to setup defrag commands: %v", err)
	}
}

func init() {
	viper.BindPFlag("settings", rootCmd.PersistentFlags().Lookup("settings"))
}
