//this script is the solution to the Rosalind Reverse Complement problem
#include <iostream>
#include <machine/limits.h>
#include <string>
#include <algorithm>

char Complement(char nucleotide){
switch (nucleotide) {
	case 'A' : return 'T';
	case 'G' : return 'C';
	case 'T' : return 'A';
	case 'C' : return 'G';
	default : return nucleotide;
	} 
}

std::string reverseComplement(const std::string& dna) {
	std::string revComp;
	revComp.reserve(dna.size());

	for (auto it=dna.rbegin(); it != dna.rend(); ++it) {
		revComp += Complement(*it);
	}
	return revComp;
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		std::cerr << "Usage:" << argv[0] << "<Dna Sequence>" << std::endl;
		return 1;
	}
	std::string dna = argv[1];
	std::string revCompDna = reverseComplement(dna);

	std::cout << "reverse complement: " << revCompDna << std::endl;
}
