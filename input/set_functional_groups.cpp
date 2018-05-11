#include <random>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

double uniform(std::mt19937_64& engine, int start, int stop) {
    std::uniform_int_distribution<int> dis(start, stop); 
    return dis(engine); 
} 

int main(int argc, char* argv[]) {
	unsigned NumberOfChains {};
    std::vector<unsigned> ChainLengths {};
    std::vector<unsigned> FunctionalGroups {};
    double FractionFunctional {};
    unsigned Seed {}; 
    std::string FunctionalFileName{};
    
    if (argc != 5) {
    	std::cout << "usage: ./set_functional_groups FUNCTIONALFILE CHAINFILE FRACTIONFUNCTIONAL SEED" << std::endl;
    	return EXIT_FAILURE;
    }
    
    FunctionalFileName = std::string(argv[1]);
    std::ifstream inputfile(argv[2], std::ios::in);
    FractionFunctional = std::atof(argv[3]);
    std::cout << FractionFunctional;
    Seed = atoi(argv[4]);
    
    if (!inputfile.is_open()) {
		std::cout << "could not open file '" << argv[1] << "' , exiting" << std::endl;
		return EXIT_FAILURE;
    }
    
    unsigned Number {}, Length {};
    while (inputfile >> Number >> Length) {
    	ChainLengths.push_back(Length);
    	FunctionalGroups.push_back((unsigned)(FractionFunctional*Length));
    }
    
    for (unsigned i = 0; i < ChainLengths.size(); i++) {
    	std::cout << ChainLengths[i] << " " << FunctionalGroups[i] << std::endl;
    }


    std::mt19937_64 generator(Seed);
    generator.discard(10000); 
    
    FILE* file_functional = fopen(FunctionalFileName.c_str(), "w");
    
    bool found {};
    for (unsigned i = 0; i < ChainLengths.size(); i++) {
    	std::vector<int> ReactiveMonos(FunctionalGroups[i],ChainLengths[i]+1);
    	for (unsigned j = 0; j < FunctionalGroups[i]; j++) {
    		found = false;
			while (!found) {
				int Mono = uniform(generator, 0, ChainLengths[i]-1);
				found = true;
				for (unsigned k = 0; k < j; k++) {
					if (Mono >= ReactiveMonos[k] -1 && Mono <= ReactiveMonos[k]+1) {
						found = false;
						break;
					}
				}
				ReactiveMonos[j] = Mono;
			}
    	    std::cout << ReactiveMonos[j] << " ";
    	}
    	std::cout << std::endl;
    	std::sort(ReactiveMonos.begin(), ReactiveMonos.end());
    	for (auto& mono : ReactiveMonos) {
    		fprintf(file_functional, "%d %d \n", i, mono);
    	}
    }
    return 0; 
}
