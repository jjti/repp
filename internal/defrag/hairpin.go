package defrag

import (
	"log"
	"os/exec"
	"regexp"
	"strconv"

	"github.com/jjtimmons/defrag/config"
)

// hairpin finds the melting temperature of a hairpin in a sequence
// returns 0 if there is none
func hairpin(seq string, conf *config.Config) (melt int) {
	// if it's longer than 60bp (max for ntthal) find the max between
	// the start and end of the sequence
	if len(seq) > 60 {
		startHairpin := hairpin(seq[:60], conf)
		endHairpin := hairpin(seq[len(seq)-60:], conf)

		if startHairpin > endHairpin {
			return startHairpin
		}
		return endHairpin
	}

	// see nnthal (no parameters) help. within primer3 distribution
	ntthalCmd := exec.Command(
		"ntthal",
		"-a", "HAIRPIN",
		"-t", "50",
		"-s1", seq,
		"-path", conf.Primer3Config,
	)

	ntthalOut, err := ntthalCmd.CombinedOutput()
	if err != nil {
		log.Fatal(err)
	}

	ntthalOutString := string(ntthalOut)
	tempRegex := regexp.MustCompile("t = (\\d*)")

	if tempRegex.MatchString(ntthalOutString) {
		tempString := tempRegex.FindAllStringSubmatch(ntthalOutString, 1)
		if len(tempString) < 1 {
			return 0
		}
		temp, err := strconv.Atoi(tempString[0][1])

		if err != nil {
			log.Fatal(err)
		}

		return temp
	}

	return 0
}
