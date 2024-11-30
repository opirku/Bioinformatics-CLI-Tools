#include <iostream>
#include <machine/limits.h>
#include <string>
#include <map>   
#include <fstream>

class ParseFasta {
public:

        std::map<std::string, std::string> dnaSequences;

        ParseFasta(const char* filePath) {
        // Open the file using the file path
        std::ifstream file(filePath);
        if (!file) {
            std::cerr << "Error opening file: " << filePath << std::endl;
            exit(1); // Or handle the error appropriately
        }

        // Read the file and store sequences
        std::string line;
        std::string currentHeader;
        while (std::getline(file, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {  // Header line
                currentHeader = line.substr(1);  // Remove the '>'
            } else {  // Sequence line
                dnaSequences[currentHeader] += line;
            }
        }

        file.close();
    }
    
};

class ReverseComplement {
public:
    static std::string reverse_Complement(const std::string& sequence) {
        // Reserve space for the reverse complement string to avoid reallocation
        std::string revComp;
        revComp.reserve(sequence.size());

        // Iterate through the sequence in reverse order
        for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
            // Switch case for reverse complement of each base
            switch (*it) {
                case 'A': revComp += 'T'; break;
                case 'T': revComp += 'A'; break;
                case 'C': revComp += 'G'; break;
                case 'G': revComp += 'C'; break;
                default: revComp += 'N'; break; // For any non-standard base
            }
        }
        return revComp;
    }
};

class PalindromeFinder {
private:
        ReverseComplement reverseComplement;
public:
    void findReversePalindromes(const std::string& header, const std::string& sequence) {
        int sequenceLength = sequence.length();
        for (int length = 4; length <= 12; ++length) {
            for (int start = 0; start <= sequenceLength - length; ++start) {
                std::string substring = sequence.substr(start, length);
                std::string rev_complement = reverseComplement.reverse_Complement(substring);

                if (substring == rev_complement) {
                    std::cout << start + 1 << "\t" << length << std::endl;
                }
            }
        }
    }
};


int main(int argc, char* argv[]){
        
        if (argc < 2) {
                std::cerr << "Usage: " << argv[0] << "<fasta file>" << std::endl;
                return 1;
        }
        
        ParseFasta parser(argv[1]);

    // Instantiate the palindrome finder
        PalindromeFinder palindromeFinder;

    // Loop through each sequence and find reverse palindromes
    for (const auto& [header, sequence] : parser.dnaSequences) {
        palindromeFinder.findReversePalindromes(header, sequence);
    }

    return 0;}

