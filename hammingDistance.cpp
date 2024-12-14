#include <iostream>
#include <string>
#include <stdexcept> // For std::invalid_argument

class HammingDistance {
private:
    std::string sequenceOne;
    std::string sequenceTwo;

public:
    HammingDistance(const std::string& first, const std::string& second)
        : sequenceOne(first), sequenceTwo(second) {
        if (sequenceOne.length() != sequenceTwo.length()) {
            throw std::invalid_argument("Strings must be of equal length!");
        }
    }

    int computeHammingDistance() const {
        
        int hammingDistance = 0;
       
        for (size_t i = 0; i < sequenceOne.length(); ++i) {

                if (sequenceOne[i] != sequenceTwo[i]) {
                ++hammingDistance;
            }
        }     
        return hammingDistance;
    }
};

int main(int argc, char* argv[]) {
    try {
        if (argc != 3) {
            throw std::invalid_argument("Usage: hammingDistance <Sequence One> <Sequence Two>");
        }

        std::string sequenceOne = argv[1];
        std::string sequenceTwo = argv[2];

        HammingDistance calculator(sequenceOne, sequenceTwo);

        std::cout << "Hamming Distance: " << calculator.computeHammingDistance() << std::endl;

    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}  
