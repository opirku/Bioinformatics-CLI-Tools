//Solution to the rabbit pair problem in Rosalind
#include <iostream>
#include <vector>

class RabbitPairs {
  private:
    int n, k;
    std::vector<long long> rabbit_pairs;

  public:
    RabbitPairs(int months, int pairs_per_generation)
        : n(months), k(pairs_per_generation), rabbit_pairs(months + 1, 0) {
        rabbit_pairs[1] = 1;
    }

    void calculateRabbitPairs() {
        for (int i = 2; i <= n; ++i) {
            rabbit_pairs[i] = (i==2) ? 1 : rabbit_pairs[i - 1] + k * rabbit_pairs[i - 2];
        }
    }

    void printRabbitPairs() const {
        std::cout << rabbit_pairs[n] << std::endl;
    }
};


int main(int argc, char* argv[]) {
    // Check if the correct number of arguments are provided
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <months> <initial_pairs>" << std::endl;
        return 1;
    }

    // Convert arguments to integers
    int n = std::stoi(argv[1]);  // Number of months
    int k = std::stoi(argv[2]);  // Number of initial rabbit pairs

    // Check if the input values are within valid range
    if (n <= 0 || n > 40 || k <= 0 || k > 5) {
        std::cerr << "Invalid input: n must be between 1 and 40, k must be between 1 and 5." << std::endl;
        return 1;
    }

    // Create an object of RabbitPairs class
    RabbitPairs rabbits(n, k);
    rabbits.calculateRabbitPairs();
    rabbits.printRabbitPairs();

    return 0;
}
