// Package config is for app wide settings that are unmarshalled
// from Viper (see: /cmd)
package config

import (
	"testing"
)

func TestConfig_SynthCost(t *testing.T) {
	type fields struct {
		SynthesisCost map[int]SynthCost
	}
	type args struct {
		fragLength int
	}

	configFields := fields{
		SynthesisCost: map[int]SynthCost{
			50: SynthCost{
				Fixed:   true,
				Dollars: 5,
			},
			200: SynthCost{
				Fixed:   false,
				Dollars: 0.1,
			},
			10000: SynthCost{
				Fixed:   false,
				Dollars: 0.50,
			},
		},
	}

	tests := []struct {
		name   string
		fields fields
		args   args
		want   float64
	}{
		{
			"fixed cost synthesis",
			configFields,
			args{
				fragLength: 20,
			},
			5.0,
		},
		{
			"variable cost synthesis (small)",
			configFields,
			args{
				fragLength: 100,
			},
			10.0,
		},
		{
			"variable cost synthesis (large)",
			configFields,
			args{
				fragLength: 1000,
			},
			500.0,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			c := Config{
				SynthesisCost: tt.fields.SynthesisCost,
			}
			if got := c.SynthCost(tt.args.fragLength); got != tt.want {
				t.Errorf("Config.SynthCost() = %v, want %v", got, tt.want)
			}
		})
	}
}
