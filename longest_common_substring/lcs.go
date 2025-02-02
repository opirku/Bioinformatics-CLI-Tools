package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"strings"
)

// readFastaFile reads the FASTA file and returns all sequences.
func readFastaFile(filename string) ([]string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var sequences []string
	var currentSequence strings.Builder
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			// If there is a current sequence, save it before starting a new one
			if currentSequence.Len() > 0 {
				sequences = append(sequences, currentSequence.String())
				currentSequence.Reset()
			}
			continue
		}
		currentSequence.WriteString(strings.TrimSpace(line))
	}

	if currentSequence.Len() > 0 {
		sequences = append(sequences, currentSequence.String())
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return sequences, nil
}

// checkIfExists checks if the substring exists in all strings.
func checkIfExists(substring string, sequences []string) bool {
	for _, seq := range sequences {
		if !strings.Contains(seq, substring) {
			return false
		}
	}
	return true
}

// findLCS finds the longest common substring of the collection.
func findLCS(sequences []string) string {
	// Binary search range
	minLen, maxLen := 0, len(sequences[0])

	var result string
	for length := minLen; length <= maxLen; length++ {
		var found bool
		// Check all substrings of current length from the first string
		for i := 0; i+length <= len(sequences[0]); i++ {
			substr := sequences[0][i : i+length]
			if checkIfExists(substr, sequences) {
				if len(substr) > len(result) {
					result = substr
				}
				found = true
			}
		}

		// If no substrings of the current length exist in all strings, stop
		if !found {
			break
		}
	}

	return result
}

func main() {
	// Define the FASTA file flag.
	fastaFile := flag.String("fasta", "", "Path to the FASTA file containing the DNA sequences")
	flag.Parse()

	// Ensure the FASTA file flag is provided.
	if *fastaFile == "" {
		fmt.Println("Please provide the path to a FASTA file using the -fasta flag.")
		return
	}

	// Read the DNA sequences from the FASTA file.
	sequences, err := readFastaFile(*fastaFile)
	if err != nil {
		fmt.Printf("Error reading FASTA file: %v\n", err)
		return
	}

	// Find the longest common substring
	lcs := findLCS(sequences)

	// Output the result
	fmt.Println("Longest Common Substring:", lcs)
}
