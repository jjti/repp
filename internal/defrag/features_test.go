package defrag

import (
	"testing"
)

func TestNewFeatureDB(t *testing.T) {
	db := NewFeatureDB()

	if len(db.features) < 1 {
		t.Fail()
	}
}
