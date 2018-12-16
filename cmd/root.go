package cmd

import (
	"log"
	"path"
	"path/filepath"

	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "defrag",
	Short: "Build vectors from any of their sequence, features or fragments",
	Long:  ``,

	// Uncomment the following line if your bare application
	// has an action associated with it:
	//	Run: func(cmd *cobra.Command, args []string) { },
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		log.Fatalf("Failed to execute defrag: %v", err)
	}
}

func init() {
	settings, _ := filepath.Abs(path.Join("..", "config"))

	// Here you will define your flags and configuration settings.
	// Cobra supports persistent flags, which, if defined here,
	// will be global for your application.
	rootCmd.PersistentFlags().StringP("config", "c", settings, "config file")
	viper.BindPFlag("config", rootCmd.PersistentFlags().Lookup("config"))
}
