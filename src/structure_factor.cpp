#include <vector>
#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <math.h>
#include "Molecule.h"
#include "HelperFunctions.h"
#include "BoundaryConditions.h"

int main(int argc, char* argv[]) {
    std::string Directory{}, ConfigFile{}, ConfigFileStart{}, KVecFileListName {};
    int StartStep{}, Step{}, SamplingStep{};
    unsigned NumberOfMonomers {}, Replicas{};
    double BoxSize {}; 
    std::vector<Particle> Monomers{}; 
    std::vector<std::string> KVecFileList {}; 
    
    std::map<double, double> structure_factor {};
    std::map<double, int> q_count {}; 
    
    if (argc != 8) {
        std::cout << "usage: ./structure_factor RUNPARENTDIR REPLICAS NUMBEROFMONOMERS BOXSIZE STARTSTEP SAMPLINGSTEP KVECFILELIST" << std::endl;
        return EXIT_FAILURE; 
    } 
    
    Directory=argv[1];
	Replicas = std::stoi(argv[2]);
	NumberOfMonomers = std::stoi(argv[3]);
	BoxSize = std::stod(argv[4]); 
	StartStep = std::stoi(argv[5]);
	SamplingStep = std::stoi(argv[6]);
	KVecFileListName = argv[7]; 
	
    std::cout << "Directory: " << Directory << std::endl;
	std::cout << "StartStep: " << StartStep << "   SamplingStep: " << SamplingStep << std::endl;
	std::cout << "Replicas: " << Replicas << "  Number of monomers: " << NumberOfMonomers << "  BoxSize: " << BoxSize << std::endl; 
	std::cout << "KVec File List: " << KVecFileListName << std::endl; 
    
    for (unsigned i = 0; i < NumberOfMonomers; i++) {
        Monomers.push_back(Particle(i)); 
    }
    
    if (!fillConfigPool(KVecFileList, KVecFileListName)) {
        std::cout << "problem filling KVecFileList, exiting!" << std::endl; 
        return EXIT_FAILURE; 
    } 
    
    std::string name;
	name = Directory+"structure_factor";
	std::ifstream test(name);
	if (test.good()){
		std::cout << "file " << name << " already exists! Aborting..." << std::endl;
		return EXIT_FAILURE;
	}
	std::ofstream output(name, std::ios::out | std::ios::trunc);
	output.precision(8);

    double Kabs {}; 
    Vector3d KVec {}; 
    unsigned count_K {0}, count_configs{0};
    double cossq {}, sinsq {}, sum {}; ;  
    double constant {2.*M_PI/BoxSize}; 
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
    for (unsigned repl = 0; repl < Replicas; repl++) {
        Step = StartStep;
	    ConfigFileStart = Directory+"REPL-"+std::to_string(repl)+"/configs/config"; 
	    while (true) {
		    ConfigFile = ConfigFileStart+std::to_string(Step)+".pdb";
		    std::cout << ConfigFile << std::endl; 
		    if (!initializePositions(Monomers, ConfigFile)) {
		        std::cout << "problem with initializing monomers" << std::endl;
		        break; 
		    }
		    for (auto& mono : Monomers) {
		        wrapUniformNoShear(mono, BoxSize);
		    }
		    for (auto& KVecFileName : KVecFileList) {
		        std::ifstream KVecFile(KVecFileName, std::ios::in); 
		        KVecFile >> Kabs; 
		        Kabs *= constant;  
		        count_K = 0;
		        sum = 0.0;  
		        while (KVecFile >> KVec(0) >> KVec(1) >> KVec(2)) {
		            count_K++; 
		            cossq = 0.0; 
		            sinsq = 0.0; 
		           
		            for (auto& mono : Monomers) {
		                cossq += cos(KVec.dot(mono.Position)*constant); 
		                sinsq += sin(KVec.dot(mono.Position)*constant); 
		            }
		            cossq = pow(cossq, 2); 
		            sinsq = pow(sinsq, 2); 
		            sum += (cossq + sinsq)/NumberOfMonomers;   
		        }
		        structure_factor[Kabs] += sum/count_K; 		        
		    }
		    count_configs++;
		    Step += SamplingStep;
	    }    
        
    }
    
    gettimeofday(&end, NULL); 
    double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    std::cout << "total time: " << realTime << " , time per conformation: " << realTime/count_configs << std::endl;
    
    
    for (auto& w : structure_factor) {
        output << w.first << " " << w.second/count_configs << std::endl; 
    }
    output.close(); 
    
    
    
    return 0; 
}
