// Package config is for app wide settings that are unmarshalled
// from Viper (see: /cmd)
package config

import (
	"reflect"
	"testing"
)

func TestNew(t *testing.T) {
	tests := []struct {
		name  string
		wantC Config
	}{
		{
			"create a map from fragment size upper limit to cost",
			Config{},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotC := New(); !reflect.DeepEqual(gotC, tt.wantC) {
				t.Errorf("New() = %+v, want %+v", gotC, tt.wantC)
			}
		})
	}
}

func TestConfig_SynthCost(t *testing.T) {
	type fields struct {
		Root      string
		DB        string
		Fragments FragmentConfig
		PCR       PCRConfig
		Synthesis SynthesisConfig
		filled    bool
	}
	type args struct {
		fragLength int
	}

	configFields := fields{
		Root:      "",
		DB:        "",
		Fragments: FragmentConfig{},
		PCR:       PCRConfig{},
		Synthesis: SynthesisConfig{
			Cost: map[int]SynthCost{
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
			MaxLength: 1000,
			MinLength: 20,
		},
		filled: true,
	}

	tests := []struct {
		name   string
		fields fields
		args   args
		want   float32
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
				Root:      tt.fields.Root,
				DB:        tt.fields.DB,
				Fragments: tt.fields.Fragments,
				PCR:       tt.fields.PCR,
				Synthesis: tt.fields.Synthesis,
				filled:    tt.fields.filled,
			}
			if got := c.SynthCost(tt.args.fragLength); got != tt.want {
				t.Errorf("Config.SynthCost() = %v, want %v", got, tt.want)
			}
		})
	}
}
