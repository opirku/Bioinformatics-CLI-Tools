package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"strings"
)

const NUCLEOTIDES = 5

func readFastaFile(filename string) (string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return "", err
	}
	defer file.Close()

	var dnaSequence string
	scanner := bufio.NewScanner(file)
	isFirstLine := true
	for scanner.Scan() {
		line := scanner.Text()
		if isFirstLine && strings.HasPrefix(line, ">") {
			isFirstLine = false
			continue // Skip header line
		}
		// Remove any spaces or newlines and add to dnaSequence
		dnaSequence += strings.ReplaceAll(line, "\n", "")
	}

	if err := scanner.Err(); err != nil {
		return "", err
	}

	return dnaSequence, nil
}

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

func computeBuckets(suffixes []int, n int, bucketStart, bucketEnd []int) {
	bucketSize := make([]int, NUCLEOTIDES)

	for i := 0; i < n; i++ {
		bucketSize[suffixes[i]]++
	}

	bucketStart[0] = 0
	bucketEnd[0] = bucketSize[0] - 1
	for i := 1; i < NUCLEOTIDES; i++ {
		bucketStart[i] = bucketStart[i-1] + bucketSize[i-1]
		bucketEnd[i] = bucketStart[i] + bucketSize[i] - 1
	}
}

func classifySuffixesLorS(suffixes []int, n int, suffixTypes []int) {
	suffixTypes[n-1] = 0

	for i := n - 2; i >= 0; i-- {
		if suffixes[i] > suffixes[i+1] || (suffixes[i] == suffixes[i+1] && suffixTypes[i+1] == 0) {
			suffixTypes[i] = 0 // S type
		} else {
			suffixTypes[i] = 1 // L type
		}
	}
}

func findLMSSuffixes(suffixes []int, n int, suffixTypes []int) []int {
	var lms []int
	for i := 1; i < n; i++ {
		if suffixTypes[i] == 0 && suffixTypes[i-1] == 1 {
			lms = append(lms, i)
		}
	}
	return lms
}

func induceSort(sa, suffixes, suffixTypes, bucketStartOrig, bucketEndOrig []int, n int) {
	bucketStart := append([]int(nil), bucketStartOrig...)
	bucketEnd := append([]int(nil), bucketEndOrig...)

	for i := 0; i < n; i++ {
		sa[i] = -1
	}

	for i := len(sa) - 1; i >= 0; i-- {
		if sa[i] >= 0 {
			ch := suffixes[sa[i]]
			sa[bucketEnd[ch]] = sa[i]
			bucketEnd[ch]--
		}
	}

	for i := 0; i < n; i++ {
		if sa[i] > 0 && suffixTypes[sa[i]-1] == 1 {
			ch := suffixes[sa[i]-1]
			sa[bucketStart[ch]] = sa[i] - 1
			bucketStart[ch]++
		}
	}

	for i := n - 1; i >= 0; i-- {
		if sa[i] > 0 && suffixTypes[sa[i]-1] == 0 {
			ch := suffixes[sa[i]-1]
			sa[bucketEnd[ch]] = sa[i] - 1
			bucketEnd[ch]--
		}
	}
}

func buildSuffixArray(dnaSequence string) []int {
	n := len(dnaSequence)
	suffixes := make([]int, n)
	encodeDNASequence(dnaSequence, suffixes)

	sa := make([]int, n)
	bucketStart := make([]int, NUCLEOTIDES)
	bucketEnd := make([]int, NUCLEOTIDES)
	suffixTypes := make([]int, n)
	computeBuckets(suffixes, n, bucketStart, bucketEnd)
	classifySuffixesLorS(suffixes, n, suffixTypes)
	lms := findLMSSuffixes(suffixes, n, suffixTypes)
	_ = lms // Avoid unused variable warning
	induceSort(sa, suffixes, suffixTypes, bucketStart, bucketEnd, n)

	return sa
}

func computeLCP(suffixArray []int, dnaSequence string) []int {
	n := len(suffixArray)
	lcp := make([]int, n-1)
	rank := make([]int, n)

	for i, suffix := range suffixArray {
		rank[suffix] = i
	}

	lcpLength := 0

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

func findLCS(dnaSequence string, suffixArray []int) []string {
	lcp := computeLCP(suffixArray, dnaSequence)

	maxLCP := 0
	var lcsList []string
	for i, l := range lcp {
		if l > maxLCP {
			maxLCP = l
			lcsList = []string{dnaSequence[suffixArray[i] : suffixArray[i]+maxLCP]}
		} else if l == maxLCP {
			lcsList = append(lcsList, dnaSequence[suffixArray[i]:suffixArray[i]+maxLCP])
		}
	}

	return lcsList
}

func main() {
	// Define the FASTA file flag
	fastaFile := flag.String("fasta", "", "Path to the FASTA file containing the DNA sequence")
	flag.Parse()

	// Ensure the FASTA file flag is provided
	if *fastaFile == "" {
		fmt.Println("Please provide the path to a FASTA file using the -fasta flag.")
		return
	}

	// Read DNA sequence from the FASTA file
	dna, err := readFastaFile(*fastaFile)
	if err != nil {
		fmt.Printf("Error reading FASTA file: %v\n", err)
		return
	}

	// Build Suffix Array and compute LCP
	sa := buildSuffixArray(dna)
	fmt.Println("Suffix Array:", sa)

	// Find LCS sequences
	lcs := findLCS(dna, sa)
	fmt.Println("Longest Common Substrings:")
	for _, seq := range lcs {
		fmt.Println(seq)
	}
}
