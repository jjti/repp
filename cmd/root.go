package cmd

import (
	"log"
	"path"
	"path/filepath"

	"github.com/spf13/viper"

	"github.com/spf13/cobra"
)

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "decvec",
	Short: "Build vectors from their sequence and a fragment database",
	Long:  ``,

	// Uncomment the following line if your bare application
	// has an action associated with it:
	//	Run: func(cmd *cobra.Command, args []string) { },
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		log.Fatalf("Failed to execute decvec: %v", err)
	}
}

func init() {
	settings, _ := filepath.Abs(path.Join("../config/settings"))

	// Here you will define your flags and configuration settings.
	// Cobra supports persistent flags, which, if defined here,
	// will be global for your application.
	rootCmd.Flags().StringP("config", "c", settings, "config file (default is /config/settings.yaml)")
	viper.BindPFlag("make.target", rootCmd.Flags().Lookup("config"))
}
