package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"sort"
	"strings"
)

const NUCLEOTIDES = 5

// readFastaFile reads the FASTA file and concatenates all non-header lines.
func readFastaFile(filename string) (string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return "", err
	}
	defer file.Close()

	var dnaSequence strings.Builder
	scanner := bufio.NewScanner(file)
	firstSequence := true

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			// If this is not the first sequence, insert a delimiter before starting the next sequence.
			if !firstSequence {
				dnaSequence.WriteString("$")
			}
			firstSequence = false
			continue
		}
		// Append the sequence part (remove newlines and spaces)
		dnaSequence.WriteString(strings.TrimSpace(line))
	}
	if err := scanner.Err(); err != nil {
		return "", err
	}
	return dnaSequence.String(), nil
}

// encodeDNASequence is kept for your original approach.
// (Not used in the naïve suffix array construction below.)
func encodeDNASequence(dnaSequence string, suffixes []int) {
	n := len(dnaSequence)
	if len(suffixes) < n {
		panic("suffixes slice is too small")
	}
	for i := 0; i < n; i++ {
		switch dnaSequence[i] {
		case 'A':
			suffixes[i] = 0
		case 'C':
			suffixes[i] = 1
		case 'G':
			suffixes[i] = 2
		case 'T':
			suffixes[i] = 3
		case '$':
			suffixes[i] = 4
		default:
			fmt.Printf("Invalid nucleotide: %c\n", dnaSequence[i])
			os.Exit(1)
		}
	}
}

// buildSuffixArrayNaive builds a suffix array by simply sorting all suffixes.
func buildSuffixArrayNaive(dnaSequence string) []int {
	n := len(dnaSequence)
	sa := make([]int, n)
	for i := 0; i < n; i++ {
		sa[i] = i
	}
	sort.Slice(sa, func(i, j int) bool {
		return dnaSequence[sa[i]:] < dnaSequence[sa[j]:]
	})
	return sa
}

// computeLCP computes the LCP array using Kasai's algorithm.
func computeLCP(suffixArray []int, dnaSequence string) []int {
	n := len(suffixArray)
	lcp := make([]int, n-1)
	rank := make([]int, n)
	for i, suffix := range suffixArray {
		rank[suffix] = i
	}

	lcpLength := 0
	// Loop through every index of the original string.
	for i := 0; i < n; i++ {
		if rank[i] > 0 {
			j := suffixArray[rank[i]-1]
			for i+lcpLength < n && j+lcpLength < n && dnaSequence[i+lcpLength] == dnaSequence[j+lcpLength] {
				lcpLength++
			}
			lcp[rank[i]-1] = lcpLength
			if lcpLength > 0 {
				lcpLength--
			}
		}
	}
	return lcp
}

// findLCS returns a slice containing all longest common substrings (LCS)
// found from the LCP array. If there is no common substring (max LCP 0), it returns an empty slice.
func findLCS(dnaSequence string, suffixArray []int) []string {
	lcp := computeLCP(suffixArray, dnaSequence)
	maxLCP := 0
	var lcsList []string
	for i, l := range lcp {
		if l > maxLCP {
			maxLCP = l
			// Start a new list of LCS.
			lcsList = []string{dnaSequence[suffixArray[i] : suffixArray[i]+maxLCP]}
		} else if l == maxLCP && maxLCP > 0 {
			// Append this LCS if it is of the maximum length.
			lcsList = append(lcsList, dnaSequence[suffixArray[i]:suffixArray[i]+maxLCP])
		}
	}
	// Optional: Deduplicate the list if the same substring appears multiple times.
	unique := make(map[string]bool)
	var deduped []string
	for _, s := range lcsList {
		if !unique[s] {
			unique[s] = true
			deduped = append(deduped, s)
		}
	}
	return deduped
}

func main() {
	// Define the FASTA file flag.
	fastaFile := flag.String("fasta", "", "Path to the FASTA file containing the DNA sequence")
	flag.Parse()

	// Ensure the FASTA file flag is provided.
	if *fastaFile == "" {
		fmt.Println("Please provide the path to a FASTA file using the -fasta flag.")
		return
	}

	// Read the DNA sequence from the FASTA file.
	dna, err := readFastaFile(*fastaFile)
	if err != nil {
		fmt.Printf("Error reading FASTA file: %v\n", err)
		return
	}

	// Build the suffix array using the naive method.
	sa := buildSuffixArrayNaive(dna)

	// Find all longest common substrings.
	lcs := findLCS(dna, sa)
	fmt.Println("Longest Common Substrings:")
	for _, seq := range lcs {
		fmt.Println(seq)
	}
}
