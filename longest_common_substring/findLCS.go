package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"sort"
	"strings"
)

// Boundary records the start and end indices of a sequence within the concatenated string.
type Boundary struct {
	id         string
	start, end int // [start, end) holds the sequence content (not including the separator)
}

// readFastaFile reads the FASTA file and returns a map of sequence IDs to sequences.
func readFastaFile(inputFilePath string) (map[string]string, error) {
	file, err := os.Open(inputFilePath)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	sequences := make(map[string]string)
	var currentID string
	var currentSequence strings.Builder
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			if currentID != "" {
				sequences[currentID] = currentSequence.String()
			}
			currentID = strings.TrimSpace(line[1:])
			currentSequence.Reset()
		} else {
			currentSequence.WriteString(strings.TrimSpace(line))
		}
	}
	if currentID != "" {
		sequences[currentID] = currentSequence.String()
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return sequences, nil
}

// concatenateInputSequences builds a concatenated string from the sequences,
// inserting a '$' separator after each sequence. It also returns the boundaries
// (start and end positions) of each sequence in the concatenated string.
func concatenateInputSequences(sequences map[string]string) (string, []Boundary) {
	// To have a consistent order (since map order is random), sort the keys.
	var ids []string
	for id := range sequences {
		ids = append(ids, id)
	}
	sort.Strings(ids)

	var concatenatedSeq strings.Builder
	var boundaries []Boundary

	for _, id := range ids {
		start := concatenatedSeq.Len()
		concatenatedSeq.WriteString(sequences[id])
		end := concatenatedSeq.Len() // end of the sequence (not including separator)
		// Append a separator after the sequence.
		concatenatedSeq.WriteString("$")
		boundaries = append(boundaries, Boundary{
			id:    id,
			start: start,
			end:   end,
		})
	}
	return concatenatedSeq.String(), boundaries
}

// buildSuffixArrayNaive builds a suffix array for the input string using a naive method.
func buildSuffixArrayNaive(s string) []int {
	n := len(s)
	// Build slice of all starting indices.
	suffixArray := make([]int, n)
	for i := 0; i < n; i++ {
		suffixArray[i] = i
	}
	// Use sort.Slice comparing the suffix strings.
	sort.Slice(suffixArray, func(i, j int) bool {
		return s[suffixArray[i]:] < s[suffixArray[j]:]
	})
	return suffixArray
}

// computeLCParray computes the LCP array for the string given the suffix array.
func computeLCParray(s string, suffixArray []int) []int {
	n := len(s)
	lcpArray := make([]int, n)
	rank := make([]int, n)
	for i := 0; i < n; i++ {
		rank[suffixArray[i]] = i
	}
	lcpLength := 0
	for i := 0; i < n; i++ {
		if rank[i] > 0 {
			j := suffixArray[rank[i]-1]
			for i+lcpLength < n && j+lcpLength < n && s[i+lcpLength] == s[j+lcpLength] {
				lcpLength++
			}
			lcpArray[rank[i]] = lcpLength
			if lcpLength > 0 {
				lcpLength--
			}
		}
	}
	return lcpArray
}

// getSequenceIndex returns the index (in boundaries slice) of the sequence that contains pos,
// and true if pos lies within a sequence (i.e. not in a separator). Otherwise, it returns false.
func getSequenceIndex(pos int, boundaries []Boundary) (int, bool) {
	for i, b := range boundaries {
		if pos >= b.start && pos < b.end {
			return i, true
		}
	}
	return -1, false
}

// minIntSlice returns the minimum value in a slice of ints.
func minIntSlice(a []int) int {
	if len(a) == 0 {
		return 0
	}
	minVal := a[0]
	for _, v := range a {
		if v < minVal {
			minVal = v
		}
	}
	return minVal
}

// findLCS finds the longest common substring among all sequences using the suffix and LCP arrays.
// It uses a sliding window over the suffix array. Only suffixes that start within sequence boundaries
// (i.e. not at a separator) are considered. Additionally, any candidate substring that contains the '$'
// separator is rejected.
func findLCS(s string, suffixArray, lcpArray []int, boundaries []Boundary) string {
	n := len(s)
	numSequences := len(boundaries)
	bestLen := 0
	bestStart := 0

	// frequency map: key is sequence index, value is count of suffixes in the current window
	freq := make(map[int]int)
	distinct := 0
	left := 0

	// sliding window over the suffix array indices
	for right := 0; right < n; right++ {
		// For the suffix starting at suffixArray[right], get its sequence index.
		if seqIdx, ok := getSequenceIndex(suffixArray[right], boundaries); ok {
			freq[seqIdx]++
			if freq[seqIdx] == 1 {
				distinct++
			}
		}

		// Try to shrink the window from the left while it still covers all sequences.
		for distinct == numSequences && left <= right {
			// We need at least two suffixes to have an LCP.
			if right > left {
				// Compute the minimum LCP in the window [left+1, right]
				candidateLCPs := lcpArray[left+1 : right+1]
				candidateLen := minIntSlice(candidateLCPs)
				// Candidate substring, taken from suffixArray[right] (could choose any in the window).
				candidate := s[suffixArray[right] : suffixArray[right]+candidateLen]
				// Ensure candidate does not contain the separator.
				if candidateLen > bestLen && !strings.Contains(candidate, "$") {
					bestLen = candidateLen
					bestStart = suffixArray[right]
				}
			}
			// Remove the leftmost suffix from the window.
			if seqIdx, ok := getSequenceIndex(suffixArray[left], boundaries); ok {
				freq[seqIdx]--
				if freq[seqIdx] == 0 {
					distinct--
				}
			}
			left++
		}
	}
	if bestLen == 0 {
		return ""
	}
	return s[bestStart : bestStart+bestLen]
}

func main() {
	fastaFile := flag.String("fasta", "", "Path to input FASTA file")
	flag.Parse()

	if *fastaFile == "" {
		fmt.Println("Error: Please provide a FASTA file path using the -fasta flag.")
		return
	}

	// Read FASTA file.
	sequences, err := readFastaFile(*fastaFile)
	if err != nil {
		fmt.Println("Error reading FASTA file:", err)
		return
	}

	// Concatenate sequences and record boundaries.
	concatenatedSeqs, boundaries := concatenateInputSequences(sequences)
	fmt.Println("Concatenated Sequences:")
	fmt.Println(concatenatedSeqs)
	fmt.Println("Number of Sequences:", len(boundaries))

	// Build suffix array.
	suffixArray := buildSuffixArrayNaive(concatenatedSeqs)
	fmt.Println("Suffix Array:", suffixArray)

	// Compute LCP array.
	lcpArray := computeLCParray(concatenatedSeqs, suffixArray)
	fmt.Println("LCP Array:", lcpArray)

	// Find the longest common substring.
	lcs := findLCS(concatenatedSeqs, suffixArray, lcpArray, boundaries)
	fmt.Println("Longest Common Substring:", lcs)
}
