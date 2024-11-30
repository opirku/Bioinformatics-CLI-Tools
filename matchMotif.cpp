#include <iostream>
#include <vector>
#include <string>

// Utilization of the Knuth-Morris-Pratt Algorithm to identify the positions of a DNA motif in a DNA sequence


void computeLPS(const std::string& dnaMotif, std::vector<int>& lps) {    
       
        int m = dnaMotif.length();
        lps[0]=0;
        int j = 0;

        for (int i = 1; i < m; ++i) { 
                while (j > 0 && dnaMotif[i] != dnaMotif[j]) {
                        j=lps[j-1];
                }
                if (dnaMotif[i] == dnaMotif[j]) {
                        ++j;
                }
                lps[i]=j;
        }
}

//KMP DNA motif matching algorithm 
void KMP(const std::string& dnaSeq, const std::string& dnaMotif) {

        int n = dnaSeq.length();
        int m = dnaMotif.length();

        //Compute LPS array
        std::vector<int> lps(m,0);
        computeLPS(dnaMotif, lps);
        
        int i = 0;
        int j = 0;
        std::vector<int> matches; 

        while (i < n) {
        
                if (dnaSeq[i]==dnaMotif[j]) {
                        i++;
                        j++;
                }
                if (j==m) {
                       // Match found
                        matches.push_back((i-m)+1);
                        j = lps[j - 1]; // Continue searching
                }

                else if ( i < n && dnaSeq[i] != dnaMotif[j]) {
                        
                        if (j > 0){

                                j=lps[j-1];
                        }
                        else {
                                i++;
                        }
                }
                
        }
        
        for (size_t k = 0; k < matches.size(); ++k) {
                if (k > 0) std::cout << '\t';
                std::cout << matches[k];
                }

        std::cout << std::endl;

}

int main(int argc, char* argv[]) {

        if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <text> <pattern>" << std::endl;
        return 1;
    }
        std::string dnaSeq = argv[1];
        std::string dnaMotif = argv[2];

        KMP(dnaSeq, dnaMotif);

        return 0;

}
