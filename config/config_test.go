package config

import (
	"testing"
)

func TestConfig_SynthCost(t *testing.T) {
	type fields struct {
		SynthesisCost      map[int]SynthCost
		SyntheticMaxLength int
	}
	type args struct {
		fragLength int
	}

	configFields := fields{
		SynthesisCost: map[int]SynthCost{
			50: SynthCost{
				Fixed: true,
				Cost:  5,
			},
			200: SynthCost{
				Fixed: false,
				Cost:  0.1,
			},
			10000: SynthCost{
				Fixed: false,
				Cost:  0.50,
			},
		},
		SyntheticMaxLength: 10000,
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
				CostSyntheticFragment: tt.fields.SynthesisCost,
				SyntheticMaxLength:    tt.fields.SyntheticMaxLength,
			}
			if got := c.SynthFragmentCost(tt.args.fragLength); got != tt.want {
				t.Errorf("Config.SynthCost() = %v, want %v", got, tt.want)
			}
		})
	}
}
