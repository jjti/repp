package blast

import (
	"github.com/jjtimmons/decvec/internal/frag"
)

// Mismatch finds mismatching sequences between the query sequences (from primers
// right now) and the fragment's full parent sequence in the database
func Mismatch(query *frag.Fragment, subject string) (matches []frag.Match, err error) {
	matches, err = BLAST(query)

	// keep only those matches that are against the parent fragment
	return matches, err
}
