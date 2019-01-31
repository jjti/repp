package defrag

import (
	"fmt"
	"log"
	"os"
	"regexp"
	"sort"
	"strings"

	"github.com/spf13/cobra"
)

// enzymes that can be used to linearize a backbone in preparatoin for an insert sequence
// map from the name of the enzyme to its recog, hangIndex, and cutIndex
var enzymes map[string]enzyme

// enzyme is a single enzyme that can be used to linearize a backbone before
// inserting a sequence
type enzyme struct {
	recog   string
	hangInd int
	cutInd  int
}

// Enzymes is a Cobra command that prints out all the enzymes to the stdout
// along with their recognition sequence
func Enzymes(cmd *cobra.Command, args []string) {
	var enzymeNameLog strings.Builder
	enzymeNameLog.WriteString("Enzymes in defrag:")

	// first get all the enzyme names and sort them alphabetically
	sortedNames := []string{}
	for name := range enzymes {
		sortedNames = append(sortedNames, name)
	}
	sort.Strings(sortedNames)

	// write all of them to the log
	for _, name := range sortedNames {
		enzymeNameLog.WriteString(
			fmt.Sprintf("\n%s: %s", name, enzymes[name].recog),
		)
	}
	log.Println(enzymeNameLog.String())

	os.Exit(0) // that's it, just log
}

// digest a Frag (backbone) with an enzyme's first recogition site
//
// remove the 5' end of the fragment post-cleaving. it will be degraded.
// keep exposed 3' ends. good visual explanation:
// https://warwick.ac.uk/study/csde/gsp/eportfolio/directory/pg/lsujcw/gibsonguide/
func digest(frag Frag, enz enzyme) (digested Frag, err error) {
	wrappedBp := 38 // largest current recognition site in the list of enzymes below
	if len(frag.Seq) < wrappedBp {
		return Frag{}, fmt.Errorf("%s is too short for digestion", frag.ID)
	}

	if frag.fragType == circular {
		// double check that it's circular, first 20 bp should show up at least twice
		firstBpsRegex := regexp.MustCompile(frag.Seq[:20])
		firstBpsCount := firstBpsRegex.FindAllStringIndex(frag.Seq, -1)

		if len(firstBpsCount) > 1 {
			frag.Seq = frag.Seq[:len(frag.Seq)/2] // undo the doubling of sequence for circular parts
		}
	}

	// turn recognition site (with ambigous bps) into a recognition seq
	reg := regexp.MustCompile(recogRegex(enz.recog))
	seq := frag.Seq + frag.Seq[0:wrappedBp]
	revCompSeq := revComp(frag.Seq) + revComp(frag.Seq[0:wrappedBp])

	// positive if template strand has overhang
	// negative if rev comp strand has overhang
	templateOverhangLength := enz.cutInd - enz.hangInd
	cutIndex := -1
	digestedSeq := ""
	if reg.MatchString(seq) {
		// template
		cutIndex = reg.FindStringIndex(seq)[0] // first int is the start of match
	}
	if reg.MatchString(revCompSeq) {
		// reverse complement
		revCutIndex := reg.FindStringIndex(revCompSeq)[0]
		revCutIndex = len(frag.Seq) - revCutIndex - len(enz.recog) // flip it to account for being on rev comp
		revCutIndex = (revCutIndex + len(frag.Seq)) % len(frag.Seq)
		if revCutIndex >= 0 && (revCutIndex < cutIndex || cutIndex < 0) {
			cutIndex = revCutIndex // take whichever occurs sooner in the sequence
		}
	}
	if cutIndex == -1 {
		// no valid cutsites in the sequence
		return Frag{}, fmt.Errorf("no %s cutsites found in %s", enz.recog, frag.ID)
	}

	if templateOverhangLength >= 0 {
		cutIndex = (cutIndex + enz.cutInd) % len(frag.Seq)
		digestedSeq = frag.Seq[cutIndex:] + frag.Seq[:cutIndex]
	} else {
		bottomIndex := (cutIndex + enz.cutInd) % len(frag.Seq)
		topIndex := (cutIndex + enz.hangInd) % len(frag.Seq)
		digestedSeq = frag.Seq[topIndex:] + frag.Seq[:bottomIndex]
	}

	return Frag{
		ID:  frag.ID,
		Seq: digestedSeq,
	}, nil
}

// recogRegex turns a recognition sequence into a regex sequence for searching
// sequence for searching the template sequence for digestion sites
func recogRegex(recog string) (decoded string) {
	regexDecode := map[rune]string{
		'A': "A",
		'C': "C",
		'G': "G",
		'T': "T",
		'M': "(A|C)",
		'R': "(A|G)",
		'W': "(A|T)",
		'Y': "(C|T)",
		'S': "(C|G)",
		'K': "(G|T)",
		'H': "(A|C|T)",
		'D': "(A|G|T)",
		'V': "(A|C|G)",
		'B': "(C|G|T)",
		'N': "(A|C|G|T)",
		'X': "(A|C|G|T)",
	}

	var regexDecoder strings.Builder
	for _, c := range recog {
		regexDecoder.WriteString(regexDecode[c])
	}

	return regexDecoder.String()
}

func init() {
	enzymes = map[string]enzyme{
		"PI-SceI":    enzyme{recog: "ATCTATGTCGGGTGCGGAGAAAGAGGTAATGAAATGG", hangInd: 11, cutInd: 15},
		"PI-PspI":    enzyme{recog: "TGGCAAACAGCTATTATGGGTATTATGGGT", hangInd: 13, cutInd: 17},
		"I-CeuI":     enzyme{recog: "TAACTATAACGGTCCTAAGGTAGCGAA", hangInd: 14, cutInd: 18},
		"I-SceI":     enzyme{recog: "TAGGGATAACAGGGTAAT", hangInd: 5, cutInd: 9},
		"AscI":       enzyme{recog: "GGCGCGCC", hangInd: 6, cutInd: 2},
		"AsiSI":      enzyme{recog: "GCGATCGC", hangInd: 3, cutInd: 5},
		"FseI":       enzyme{recog: "GGCCGGCC", hangInd: 2, cutInd: 6},
		"NotI":       enzyme{recog: "GCGGCCGC", hangInd: 6, cutInd: 2},
		"NotI-HF":    enzyme{recog: "GCGGCCGC", hangInd: 6, cutInd: 2},
		"PacI":       enzyme{recog: "TTAATTAA", hangInd: 3, cutInd: 5},
		"PmeI":       enzyme{recog: "GTTTAAAC", hangInd: 4, cutInd: 4},
		"PspXI":      enzyme{recog: "VCTCGAGB", hangInd: 6, cutInd: 2},
		"SbfI":       enzyme{recog: "CCTGCAGG", hangInd: 2, cutInd: 6},
		"SbfI-HF":    enzyme{recog: "CCTGCAGG", hangInd: 2, cutInd: 6},
		"SfiI":       enzyme{recog: "GGCCNNNNNGGCC", hangInd: 5, cutInd: 8},
		"SgrAI":      enzyme{recog: "CRCCGGYG", hangInd: 6, cutInd: 2},
		"SrfI":       enzyme{recog: "GCCCGGGC", hangInd: 4, cutInd: 4},
		"SwaI":       enzyme{recog: "ATTTAAAT", hangInd: 4, cutInd: 4},
		"BaeI":       enzyme{recog: "NNNNNNNNNNNNNNNACNNNNGTAYCNNNNNNNNNNNN", hangInd: 33, cutInd: 38},
		"BbvCI":      enzyme{recog: "CCTCAGC", hangInd: 5, cutInd: 2},
		"BspQI":      enzyme{recog: "GCTCTTCNNNN", hangInd: 11, cutInd: 8},
		"CspCI":      enzyme{recog: "NNNNNNNNNNNNNCAANNNNNGTGGNNNNNNNNNNNN", hangInd: 35, cutInd: 37},
		"Nb.BbvCI":   enzyme{recog: "CCTCAGC", hangInd: 5, cutInd: 0},
		"Nt.BbvCI":   enzyme{recog: "CCTCAGC", hangInd: 0, cutInd: 2},
		"Nt.BspQI":   enzyme{recog: "GCTCTTCN", hangInd: 0, cutInd: 8},
		"PpuMI":      enzyme{recog: "RGGWCCY", hangInd: 5, cutInd: 2},
		"RsrII":      enzyme{recog: "CGGWCCG", hangInd: 5, cutInd: 2},
		"SapI":       enzyme{recog: "GCTCTTCNNNN", hangInd: 11, cutInd: 8},
		"SexAI":      enzyme{recog: "ACCWGGT", hangInd: 6, cutInd: 1},
		"AatII":      enzyme{recog: "GACGTC", hangInd: 1, cutInd: 5},
		"Acc65I":     enzyme{recog: "GGTACC", hangInd: 5, cutInd: 1},
		"AccI":       enzyme{recog: "GTMKAC", hangInd: 4, cutInd: 2},
		"AclI":       enzyme{recog: "AACGTT", hangInd: 4, cutInd: 2},
		"AcuI":       enzyme{recog: "CTGAAGNNNNNNNNNNNNNNNN", hangInd: 20, cutInd: 22},
		"AfeI":       enzyme{recog: "AGCGCT", hangInd: 3, cutInd: 3},
		"AflII":      enzyme{recog: "CTTAAG", hangInd: 5, cutInd: 1},
		"AflIII":     enzyme{recog: "ACRYGT", hangInd: 5, cutInd: 1},
		"AgeI":       enzyme{recog: "ACCGGT", hangInd: 5, cutInd: 1},
		"AgeI-HF":    enzyme{recog: "ACCGGT", hangInd: 5, cutInd: 1},
		"AhdI":       enzyme{recog: "GACNNNNNGTC", hangInd: 5, cutInd: 6},
		"AleI":       enzyme{recog: "CACNNNNGTG", hangInd: 5, cutInd: 5},
		"AlwNI":      enzyme{recog: "CAGNNNCTG", hangInd: 3, cutInd: 6},
		"ApaI":       enzyme{recog: "GGGCCC", hangInd: 1, cutInd: 5},
		"ApaLI":      enzyme{recog: "GTGCAC", hangInd: 5, cutInd: 1},
		"ApoI":       enzyme{recog: "RAATTY", hangInd: 5, cutInd: 1},
		"ApoI-HF":    enzyme{recog: "RAATTY", hangInd: 5, cutInd: 1},
		"AseI":       enzyme{recog: "ATTAAT", hangInd: 4, cutInd: 2},
		"AvaI":       enzyme{recog: "CYCGRG", hangInd: 5, cutInd: 1},
		"AvrII":      enzyme{recog: "CCTAGG", hangInd: 5, cutInd: 1},
		"BaeGI":      enzyme{recog: "GKGCMC", hangInd: 1, cutInd: 5},
		"BamHI":      enzyme{recog: "GGATCC", hangInd: 5, cutInd: 1},
		"BamHI-HF":   enzyme{recog: "GGATCC", hangInd: 5, cutInd: 1},
		"BanI":       enzyme{recog: "GGYRCC", hangInd: 5, cutInd: 1},
		"BanII":      enzyme{recog: "GRGCYC", hangInd: 1, cutInd: 5},
		"BbsI":       enzyme{recog: "GAAGACNNNNNN", hangInd: 12, cutInd: 8},
		"BbsI-HF":    enzyme{recog: "GAAGACNNNNNN", hangInd: 12, cutInd: 8},
		"BcgI":       enzyme{recog: "NNNNNNNNNNNNCGANNNNNNTGCNNNNNNNNNNNN", hangInd: 34, cutInd: 36},
		"BciVI":      enzyme{recog: "GTATCCNNNNNN", hangInd: 11, cutInd: 12},
		"BclI":       enzyme{recog: "TGATCA", hangInd: 5, cutInd: 1},
		"BclI-HF":    enzyme{recog: "TGATCA", hangInd: 5, cutInd: 1},
		"BfuAI":      enzyme{recog: "ACCTGCNNNNNNNN", hangInd: 14, cutInd: 10},
		"BglI":       enzyme{recog: "GCCNNNNNGGC", hangInd: 4, cutInd: 7},
		"BglII":      enzyme{recog: "AGATCT", hangInd: 5, cutInd: 1},
		"BlpI":       enzyme{recog: "GCTNAGC", hangInd: 5, cutInd: 2},
		"BmgBI":      enzyme{recog: "CACGTC", hangInd: 3, cutInd: 3},
		"BmrI":       enzyme{recog: "ACTGGGNNNNN", hangInd: 10, cutInd: 11},
		"BmtI":       enzyme{recog: "GCTAGC", hangInd: 1, cutInd: 5},
		"BmtI-HF":    enzyme{recog: "GCTAGC", hangInd: 1, cutInd: 5},
		"BpmI":       enzyme{recog: "CTGGAGNNNNNNNNNNNNNNNN", hangInd: 20, cutInd: 22},
		"Bpu10I":     enzyme{recog: "CCTNAGC", hangInd: 5, cutInd: 2},
		"BpuEI":      enzyme{recog: "CTTGAGNNNNNNNNNNNNNNNN", hangInd: 20, cutInd: 22},
		"BsaAI":      enzyme{recog: "YACGTR", hangInd: 3, cutInd: 3},
		"BsaBI":      enzyme{recog: "GATNNNNATC", hangInd: 5, cutInd: 5},
		"BsaHI":      enzyme{recog: "GRCGYC", hangInd: 4, cutInd: 2},
		"BsaI":       enzyme{recog: "GGTCTCNNNNN", hangInd: 11, cutInd: 7},
		"BsaI-HF":    enzyme{recog: "GGTCTCNNNNN", hangInd: 11, cutInd: 7},
		"BsaI-HFv2":  enzyme{recog: "GGTCTCNNNNN", hangInd: 11, cutInd: 7},
		"BsaWI":      enzyme{recog: "WCCGGW", hangInd: 5, cutInd: 1},
		"BsaXI":      enzyme{recog: "NNNNNNNNNNNNACNNNNNCTCCNNNNNNNNNN", hangInd: 30, cutInd: 33},
		"BseRI":      enzyme{recog: "GAGGAGNNNNNNNNNN", hangInd: 14, cutInd: 16},
		"BseYI":      enzyme{recog: "CCCAGC", hangInd: 5, cutInd: 1},
		"BsgI":       enzyme{recog: "GTGCAGNNNNNNNNNNNNNNNN", hangInd: 20, cutInd: 22},
		"BsiEI":      enzyme{recog: "CGRYCG", hangInd: 2, cutInd: 4},
		"BsiHKAI":    enzyme{recog: "GWGCWC", hangInd: 1, cutInd: 5},
		"BsiWI":      enzyme{recog: "CGTACG", hangInd: 5, cutInd: 1},
		"BsiWI-HF":   enzyme{recog: "CGTACG", hangInd: 5, cutInd: 1},
		"BsmBI":      enzyme{recog: "CGTCTCNNNNN", hangInd: 11, cutInd: 7},
		"BsmI":       enzyme{recog: "GAATGCN", hangInd: 5, cutInd: 7},
		"BsoBI":      enzyme{recog: "CYCGRG", hangInd: 5, cutInd: 1},
		"Bsp1286I":   enzyme{recog: "GDGCHC", hangInd: 1, cutInd: 5},
		"BspDI":      enzyme{recog: "ATCGAT", hangInd: 4, cutInd: 2},
		"BspEI":      enzyme{recog: "TCCGGA", hangInd: 5, cutInd: 1},
		"BspHI":      enzyme{recog: "TCATGA", hangInd: 5, cutInd: 1},
		"BspMI":      enzyme{recog: "ACCTGCNNNNNNNN", hangInd: 14, cutInd: 10},
		"BsrBI":      enzyme{recog: "CCGCTC", hangInd: 3, cutInd: 3},
		"BsrDI":      enzyme{recog: "GCAATGNN", hangInd: 6, cutInd: 8},
		"BsrFI":      enzyme{recog: "RCCGGY", hangInd: 5, cutInd: 1},
		"BsrGI":      enzyme{recog: "TGTACA", hangInd: 5, cutInd: 1},
		"BsrGI-HF":   enzyme{recog: "TGTACA", hangInd: 5, cutInd: 1},
		"BssHII":     enzyme{recog: "GCGCGC", hangInd: 5, cutInd: 1},
		"BssSI":      enzyme{recog: "CACGAG", hangInd: 5, cutInd: 1},
		"BstAPI":     enzyme{recog: "GCANNNNNTGC", hangInd: 4, cutInd: 7},
		"BstBI":      enzyme{recog: "TTCGAA", hangInd: 4, cutInd: 2},
		"BstEII":     enzyme{recog: "GGTNACC", hangInd: 6, cutInd: 1},
		"BstEII-HF":  enzyme{recog: "GGTNACC", hangInd: 6, cutInd: 1},
		"BstXI":      enzyme{recog: "CCANNNNNNTGG", hangInd: 4, cutInd: 8},
		"BstYI":      enzyme{recog: "RGATCY", hangInd: 5, cutInd: 1},
		"BstZ17I-HF": enzyme{recog: "GTATAC", hangInd: 3, cutInd: 3},
		"Bsu36I":     enzyme{recog: "CCTNAGG", hangInd: 5, cutInd: 2},
		"BtgI":       enzyme{recog: "CCRYGG", hangInd: 5, cutInd: 1},
		"BtgZI":      enzyme{recog: "GCGATGNNNNNNNNNNNNNN", hangInd: 20, cutInd: 16},
		"BtsI":       enzyme{recog: "GCAGTGNN", hangInd: 6, cutInd: 8},
		"ClaI":       enzyme{recog: "ATCGAT", hangInd: 4, cutInd: 2},
		"DraI":       enzyme{recog: "TTTAAA", hangInd: 3, cutInd: 3},
		"DraIII-HF":  enzyme{recog: "CACNNNGTG", hangInd: 3, cutInd: 6},
		"DrdI":       enzyme{recog: "GACNNNNNNGTC", hangInd: 5, cutInd: 7},
		"EaeI":       enzyme{recog: "YGGCCR", hangInd: 5, cutInd: 1},
		"EagI":       enzyme{recog: "CGGCCG", hangInd: 5, cutInd: 1},
		"EagI-HF":    enzyme{recog: "CGGCCG", hangInd: 5, cutInd: 1},
		"EarI":       enzyme{recog: "CTCTTCNNNN", hangInd: 10, cutInd: 7},
		"EciI":       enzyme{recog: "GGCGGANNNNNNNNNNN", hangInd: 15, cutInd: 17},
		"Eco53kI":    enzyme{recog: "GAGCTC", hangInd: 3, cutInd: 3},
		"EcoNI":      enzyme{recog: "CCTNNNNNAGG", hangInd: 6, cutInd: 5},
		"EcoO109I":   enzyme{recog: "RGGNCCY", hangInd: 5, cutInd: 2},
		"EcoRI":      enzyme{recog: "GAATTC", hangInd: 5, cutInd: 1},
		"EcoRI-HF":   enzyme{recog: "GAATTC", hangInd: 5, cutInd: 1},
		"EcoRV":      enzyme{recog: "GATATC", hangInd: 3, cutInd: 3},
		"EcoRV-HF":   enzyme{recog: "GATATC", hangInd: 3, cutInd: 3},
		"Esp3I":      enzyme{recog: "CGTCTCNNNNN", hangInd: 11, cutInd: 7},
		"FspI":       enzyme{recog: "TGCGCA", hangInd: 3, cutInd: 3},
		"HaeII":      enzyme{recog: "RGCGCY", hangInd: 1, cutInd: 5},
		"HincII":     enzyme{recog: "GTYRAC", hangInd: 3, cutInd: 3},
		"HindIII":    enzyme{recog: "AAGCTT", hangInd: 5, cutInd: 1},
		"HindIII-HF": enzyme{recog: "AAGCTT", hangInd: 5, cutInd: 1},
		"HpaI":       enzyme{recog: "GTTAAC", hangInd: 3, cutInd: 3},
		"KasI":       enzyme{recog: "GGCGCC", hangInd: 5, cutInd: 1},
		"KpnI":       enzyme{recog: "GGTACC", hangInd: 1, cutInd: 5},
		"KpnI-HF":    enzyme{recog: "GGTACC", hangInd: 1, cutInd: 5},
		"MfeI":       enzyme{recog: "CAATTG", hangInd: 5, cutInd: 1},
		"MfeI-HF":    enzyme{recog: "CAATTG", hangInd: 5, cutInd: 1},
		"MluI":       enzyme{recog: "ACGCGT", hangInd: 5, cutInd: 1},
		"MluI-HF":    enzyme{recog: "ACGCGT", hangInd: 5, cutInd: 1},
		"MmeI":       enzyme{recog: "TCCRACNNNNNNNNNNNNNNNNNNNN", hangInd: 24, cutInd: 26},
		"MscI":       enzyme{recog: "TGGCCA", hangInd: 3, cutInd: 3},
		"MslI":       enzyme{recog: "CAYNNNNRTG", hangInd: 5, cutInd: 5},
		"MspA1I":     enzyme{recog: "CMGCKG", hangInd: 3, cutInd: 3},
		"NaeI":       enzyme{recog: "GCCGGC", hangInd: 3, cutInd: 3},
		"NarI":       enzyme{recog: "GGCGCC", hangInd: 4, cutInd: 2},
		"Nb.BsmI":    enzyme{recog: "GAATGC", hangInd: 5, cutInd: 0},
		"Nb.BsrDI":   enzyme{recog: "GCAATG", hangInd: 6, cutInd: 0},
		"Nb.BssSI":   enzyme{recog: "CACGAG", hangInd: 5, cutInd: 0},
		"Nb.BtsI":    enzyme{recog: "GCAGTG", hangInd: 6, cutInd: 0},
		"NcoI":       enzyme{recog: "CCATGG", hangInd: 5, cutInd: 1},
		"NcoI-HF":    enzyme{recog: "CCATGG", hangInd: 5, cutInd: 1},
		"NdeI":       enzyme{recog: "CATATG", hangInd: 4, cutInd: 2},
		"NgoMIV":     enzyme{recog: "GCCGGC", hangInd: 5, cutInd: 1},
		"NheI":       enzyme{recog: "GCTAGC", hangInd: 5, cutInd: 1},
		"NheI-HF":    enzyme{recog: "GCTAGC", hangInd: 5, cutInd: 1},
		"NmeAIII":    enzyme{recog: "GCCGAGNNNNNNNNNNNNNNNNNNNN", hangInd: 25, cutInd: 26},
		"NruI":       enzyme{recog: "TCGCGA", hangInd: 3, cutInd: 3},
		"NruI-HF":    enzyme{recog: "TCGCGA", hangInd: 3, cutInd: 3},
		"NsiI":       enzyme{recog: "ATGCAT", hangInd: 1, cutInd: 5},
		"NsiI-HF":    enzyme{recog: "ATGCAT", hangInd: 1, cutInd: 5},
		"NspI":       enzyme{recog: "RCATGY", hangInd: 1, cutInd: 5},
		"PaeR7I":     enzyme{recog: "CTCGAG", hangInd: 5, cutInd: 1},
		"PciI":       enzyme{recog: "ACATGT", hangInd: 5, cutInd: 1},
		"PflFI":      enzyme{recog: "GACNNNGTC", hangInd: 5, cutInd: 4},
		"PflMI":      enzyme{recog: "CCANNNNNTGG", hangInd: 4, cutInd: 7},
		"PluTI":      enzyme{recog: "GGCGCC", hangInd: 1, cutInd: 5},
		"PmlI":       enzyme{recog: "CACGTG", hangInd: 3, cutInd: 3},
		"PshAI":      enzyme{recog: "GACNNNNGTC", hangInd: 5, cutInd: 5},
		"PsiI":       enzyme{recog: "TTATAA", hangInd: 3, cutInd: 3},
		"PspOMI":     enzyme{recog: "GGGCCC", hangInd: 5, cutInd: 1},
		"PstI":       enzyme{recog: "CTGCAG", hangInd: 1, cutInd: 5},
		"PstI-HF":    enzyme{recog: "CTGCAG", hangInd: 1, cutInd: 5},
		"PvuI":       enzyme{recog: "CGATCG", hangInd: 2, cutInd: 4},
		"PvuI-HF":    enzyme{recog: "CGATCG", hangInd: 2, cutInd: 4},
		"PvuII":      enzyme{recog: "CAGCTG", hangInd: 3, cutInd: 3},
		"PvuII-HF":   enzyme{recog: "CAGCTG", hangInd: 3, cutInd: 3},
		"SacI":       enzyme{recog: "GAGCTC", hangInd: 1, cutInd: 5},
		"SacI-HF":    enzyme{recog: "GAGCTC", hangInd: 1, cutInd: 5},
		"SacII":      enzyme{recog: "CCGCGG", hangInd: 2, cutInd: 4},
		"SalI":       enzyme{recog: "GTCGAC", hangInd: 5, cutInd: 1},
		"SalI-HF":    enzyme{recog: "GTCGAC", hangInd: 5, cutInd: 1},
		"ScaI-HF":    enzyme{recog: "AGTACT", hangInd: 3, cutInd: 3},
		"SfcI":       enzyme{recog: "CTRYAG", hangInd: 5, cutInd: 1},
		"SfoI":       enzyme{recog: "GGCGCC", hangInd: 3, cutInd: 3},
		"SmaI":       enzyme{recog: "CCCGGG", hangInd: 3, cutInd: 3},
		"SmlI":       enzyme{recog: "CTYRAG", hangInd: 5, cutInd: 1},
		"SnaBI":      enzyme{recog: "TACGTA", hangInd: 3, cutInd: 3},
		"SpeI":       enzyme{recog: "ACTAGT", hangInd: 5, cutInd: 1},
		"SpeI-HF":    enzyme{recog: "ACTAGT", hangInd: 5, cutInd: 1},
		"SphI":       enzyme{recog: "GCATGC", hangInd: 1, cutInd: 5},
		"SphI-HF":    enzyme{recog: "GCATGC", hangInd: 1, cutInd: 5},
		"SspI":       enzyme{recog: "AATATT", hangInd: 3, cutInd: 3},
		"SspI-HF":    enzyme{recog: "AATATT", hangInd: 3, cutInd: 3},
		"StuI":       enzyme{recog: "AGGCCT", hangInd: 3, cutInd: 3},
		"StyI":       enzyme{recog: "CCWWGG", hangInd: 5, cutInd: 1},
		"StyI-HF":    enzyme{recog: "CCWWGG", hangInd: 5, cutInd: 1},
		"TspMI":      enzyme{recog: "CCCGGG", hangInd: 5, cutInd: 1},
		"Tth111I":    enzyme{recog: "GACNNNGTC", hangInd: 5, cutInd: 4},
		"XbaI":       enzyme{recog: "TCTAGA", hangInd: 5, cutInd: 1},
		"XcmI":       enzyme{recog: "CCANNNNNNNNNTGG", hangInd: 7, cutInd: 8},
		"XhoI":       enzyme{recog: "CTCGAG", hangInd: 5, cutInd: 1},
		"XmaI":       enzyme{recog: "CCCGGG", hangInd: 5, cutInd: 1},
		"XmnI":       enzyme{recog: "GAANNNNTTC", hangInd: 5, cutInd: 5},
		"ZraI":       enzyme{recog: "GACGTC", hangInd: 3, cutInd: 3},
		"AlwI":       enzyme{recog: "GGATCNNNNN", hangInd: 10, cutInd: 9},
		"ApeKI":      enzyme{recog: "GCWGC", hangInd: 4, cutInd: 1},
		"AvaII":      enzyme{recog: "GGWCC", hangInd: 4, cutInd: 1},
		"BbvI":       enzyme{recog: "GCAGCNNNNNNNNNNNN", hangInd: 17, cutInd: 13},
		"BccI":       enzyme{recog: "CCATCNNNNN", hangInd: 10, cutInd: 9},
		"BceAI":      enzyme{recog: "ACGGCNNNNNNNNNNNNNN", hangInd: 19, cutInd: 17},
		"BcoDI":      enzyme{recog: "GTCTCNNNNN", hangInd: 10, cutInd: 6},
		"BsmAI":      enzyme{recog: "GTCTCNNNNN", hangInd: 10, cutInd: 6},
		"BsmFI":      enzyme{recog: "GGGACNNNNNNNNNNNNNN", hangInd: 19, cutInd: 15},
		"BspCNI":     enzyme{recog: "CTCAGNNNNNNNNN", hangInd: 12, cutInd: 14},
		"BsrI":       enzyme{recog: "ACTGGN", hangInd: 4, cutInd: 6},
		"BstNI":      enzyme{recog: "CCWGG", hangInd: 3, cutInd: 2},
		"BtsCI":      enzyme{recog: "GGATGNN", hangInd: 5, cutInd: 7},
		"BtsIMutI":   enzyme{recog: "CAGTGNN", hangInd: 5, cutInd: 7},
		"FauI":       enzyme{recog: "CCCGCNNNNNN", hangInd: 11, cutInd: 9},
		"FokI":       enzyme{recog: "GGATGNNNNNNNNNNNNN", hangInd: 18, cutInd: 14},
		"HgaI":       enzyme{recog: "GACGCNNNNNNNNNN", hangInd: 15, cutInd: 10},
		"HphI":       enzyme{recog: "GGTGANNNNNNNN", hangInd: 12, cutInd: 13},
		"Hpy99I":     enzyme{recog: "CGWCG", hangInd: 0, cutInd: 5},
		"HpyAV":      enzyme{recog: "CCTTCNNNNNN", hangInd: 10, cutInd: 11},
		"MboII":      enzyme{recog: "GAAGANNNNNNNN", hangInd: 12, cutInd: 13},
		"MlyI":       enzyme{recog: "GAGTCNNNNN", hangInd: 10, cutInd: 10},
		"NciI":       enzyme{recog: "CCSGG", hangInd: 3, cutInd: 2},
		"Nt.AlwI":    enzyme{recog: "GGATCNNNN", hangInd: 0, cutInd: 9},
		"Nt.BsmAI":   enzyme{recog: "GTCTCN", hangInd: 0, cutInd: 6},
		"Nt.BstNBI":  enzyme{recog: "GAGTCNNNN", hangInd: 0, cutInd: 9},
		"PleI":       enzyme{recog: "GAGTCNNNNN", hangInd: 10, cutInd: 9},
		"PspGI":      enzyme{recog: "CCWGG", hangInd: 5, cutInd: 0},
		"SfaNI":      enzyme{recog: "GCATCNNNNNNNNN", hangInd: 14, cutInd: 10},
		"TfiI":       enzyme{recog: "GAWTC", hangInd: 4, cutInd: 1},
		"TseI":       enzyme{recog: "GCWGC", hangInd: 4, cutInd: 1},
		"Tsp45I":     enzyme{recog: "GTSAC", hangInd: 5, cutInd: 0},
		"TspRI":      enzyme{recog: "NNCASTGNN", hangInd: 0, cutInd: 9},
		"AciI":       enzyme{recog: "CCGC", hangInd: 3, cutInd: 1},
		"AluI":       enzyme{recog: "AGCT", hangInd: 2, cutInd: 2},
		"BfaI":       enzyme{recog: "CTAG", hangInd: 3, cutInd: 1},
		"BsaJI":      enzyme{recog: "CCNNGG", hangInd: 5, cutInd: 1},
		"BslI":       enzyme{recog: "CCNNNNNNNGG", hangInd: 4, cutInd: 7},
		"BstUI":      enzyme{recog: "CGCG", hangInd: 2, cutInd: 2},
		"Cac8I":      enzyme{recog: "GCNNGC", hangInd: 3, cutInd: 3},
		"CviAII":     enzyme{recog: "CATG", hangInd: 3, cutInd: 1},
		"CviKI-1":    enzyme{recog: "RGCY", hangInd: 2, cutInd: 2},
		"CviQI":      enzyme{recog: "GTAC", hangInd: 3, cutInd: 1},
		"DdeI":       enzyme{recog: "CTNAG", hangInd: 4, cutInd: 1},
		"DpnII":      enzyme{recog: "GATC", hangInd: 4, cutInd: 0},
		"FatI":       enzyme{recog: "CATG", hangInd: 4, cutInd: 0},
		"Fnu4HI":     enzyme{recog: "GCNGC", hangInd: 3, cutInd: 2},
		"HaeIII":     enzyme{recog: "GGCC", hangInd: 2, cutInd: 2},
		"HhaI":       enzyme{recog: "GCGC", hangInd: 1, cutInd: 3},
		"HinP1I":     enzyme{recog: "GCGC", hangInd: 3, cutInd: 1},
		"HinfI":      enzyme{recog: "GANTC", hangInd: 4, cutInd: 1},
		"HpaII":      enzyme{recog: "CCGG", hangInd: 3, cutInd: 1},
		"Hpy166II":   enzyme{recog: "GTNNAC", hangInd: 3, cutInd: 3},
		"Hpy188I":    enzyme{recog: "TCNGA", hangInd: 2, cutInd: 3},
		"Hpy188III":  enzyme{recog: "TCNNGA", hangInd: 4, cutInd: 2},
		"HpyCH4III":  enzyme{recog: "ACNGT", hangInd: 2, cutInd: 3},
		"HpyCH4IV":   enzyme{recog: "ACGT", hangInd: 3, cutInd: 1},
		"HpyCH4V":    enzyme{recog: "TGCA", hangInd: 2, cutInd: 2},
		"MboI":       enzyme{recog: "GATC", hangInd: 4, cutInd: 0},
		"MluCI":      enzyme{recog: "AATT", hangInd: 4, cutInd: 0},
		"MnlI":       enzyme{recog: "CCTCNNNNNNN", hangInd: 10, cutInd: 11},
		"MseI":       enzyme{recog: "TTAA", hangInd: 3, cutInd: 1},
		"MspI":       enzyme{recog: "CCGG", hangInd: 3, cutInd: 1},
		"MwoI":       enzyme{recog: "GCNNNNNNNGC", hangInd: 4, cutInd: 7},
		"NlaIII":     enzyme{recog: "CATG", hangInd: 0, cutInd: 4},
		"NlaIV":      enzyme{recog: "GGNNCC", hangInd: 3, cutInd: 3},
		"RsaI":       enzyme{recog: "GTAC", hangInd: 2, cutInd: 2},
		"Sau3AI":     enzyme{recog: "GATC", hangInd: 4, cutInd: 0},
		"Sau96I":     enzyme{recog: "GGNCC", hangInd: 4, cutInd: 1},
		"ScrFI":      enzyme{recog: "CCNGG", hangInd: 3, cutInd: 2},
		"StyD4I":     enzyme{recog: "CCNGG", hangInd: 5, cutInd: 0},
		"TaqI":       enzyme{recog: "TCGA", hangInd: 3, cutInd: 1},
		"Nt.CviPII":  enzyme{recog: "CCD", hangInd: 0, cutInd: 0},
	}
}
