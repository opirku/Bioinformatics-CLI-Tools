package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

type inputSequence struct {
	header    string
	sequences string
} // Struct to hold my input fasta file

func readFastaFile(fileName string) ([]inputSequence, error) {
	file, err := os.Open(fileName)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := bufio.NewReader(file)
	var records []inputSequence
	var header, sequences string

	for {
		line, err := reader.ReadString('\n')
		if err != nil && err.Error() != "EOF" {
			return nil, err
		}
		line = strings.TrimSpace(line) // Strip newlines and extra spaces

		// If we've reached the end of the file, break the loop
		if err != nil && line == "" {
			break
		}

		// If the line starts with '>', it's a header line
		if strings.HasPrefix(line, ">") {
			// If we already have a sequence, store the previous record
			if header != "" {
				records = append(records, inputSequence{header: header, sequences: sequences})
			}
			header = line
			sequences = ""
		} else {
			// Accumulate sequence lines
			sequences += line
		}
	}

	// Add the last sequence if any
	if header != "" {
		records = append(records, inputSequence{header: header, sequences: sequences})
	}

	return records, nil
}

type consensusProfile struct {
	profile [4][]int
}

func CalculateConsensusProfile(sequences []inputSequence) consensusProfile {
	if len(sequences) == 0 {
		return consensusProfile{}
	}
	sequenceLength := len(sequences[0].sequences)

	// Initialize the profile with zeros
	profile := consensusProfile{}
	for i := 0; i < 4; i++ {
		profile.profile[i] = make([]int, sequenceLength)
	}

	// Iterate through each sequence and count the occurrences of A, C, G, T at each position
	for _, sequence := range sequences {
		for i, nucleotide := range sequence.sequences {
			switch nucleotide {
			case 'A':
				profile.profile[0][i]++
			case 'C':
				profile.profile[1][i]++
			case 'G':
				profile.profile[2][i]++
			case 'T':
				profile.profile[3][i]++
			}
		}
	}

	return profile
}

func displayProfile(consensusProfile consensusProfile) {
	nucleotides := []string{"A", "C", "G", "T"}
	for i := 0; i < 4; i++ {
		fmt.Printf("%s: ", nucleotides[i])
		for _, count := range consensusProfile.profile[i] {
			fmt.Printf("%d ", count)
		}
		fmt.Println()
	}
}
func getConsensusSequence(consensusProfile consensusProfile) string {
	nucleotides := []string{"A", "C", "G", "T"}
	consensus := ""
	for i := 0; i < len(consensusProfile.profile[0]); i++ {
		maxCount := -1
		consensusNucleotide := ""

		for j, count := range consensusProfile.profile {
			if count[i] > maxCount {
				maxCount = count[i]
				consensusNucleotide = nucleotides[j]
			}
		}
		consensus += consensusNucleotide
	}
	return consensus
}

func main() {
	if len(os.Args) < 2 {
		fmt.Println("Usage: go run consensusProf.go <fasta-file>")
		return
	}
	fileName := os.Args[1]
	sequences, err := readFastaFile(fileName)
	if err != nil {
		fmt.Println("Error Reading FASTA file: ", err)
		return
	}
	profile := CalculateConsensusProfile(sequences)
	consensus := getConsensusSequence(profile)
	fmt.Println("The Consensus Profile is: %s\n", consensus)
	displayProfile(profile)
}
