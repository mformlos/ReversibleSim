#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <math.h>
#include <sys/stat.h>
#include <sys/time.h>
#include "Molecule.h"
#include "HelperFunctions.h"
#include "BoundaryConditions.h"

int main(int argc, char* argv[]) {
    std::string Directory{}, ReplicaDir {}, ConfigFile{}, ConfigFileStart{}, StepSampleFile{}, OutputDir{}, OutputFileStart{}, OutputFileName {};
    double StartR{0.0}, DR{0.1};
    unsigned StartStep{}, SamplingStep{}, Replicas{}, Monomers{}, BoxSize{}, T0Step {}, T1Step{};    
    double DeltaT{}, bin{}; 
    
    std::vector<unsigned long long> SampleSteps {}; 
    
    
    if (argc != 10) {
            std::cout << "usage: ./Van_Hove DIRECTORY STEPFILE REPLICAS STARTSTEP SAMPLINGSTEP DELTAT Monomers BOXSIZE INTERVAL" << std::endl;  
            return EXIT_FAILURE; 
    }
    
    Directory=argv[1];
    StepSampleFile=argv[2]; 
	Replicas = std::stoi(argv[3]); 
	StartStep = std::stoi(argv[4]);
    SamplingStep = std::stoi(argv[5]); 
	DeltaT = std::stod(argv[6]); 
	Monomers = std::stoi(argv[7]); 
	BoxSize = std::stoi(argv[8]); 
	DR = std::stod(argv[9]); 
	
	std::cout << "Directory: " << Directory << std::endl;
	std::cout << "StartStep: " << StartStep << "   SamplingStep: " << SamplingStep << std::endl;
	std::cout << "Replicas: " << Replicas << " DR: " << DR << std::endl;    
	
	initializeStepVector(SampleSteps, StepSampleFile);  
	std::vector<Particle> Monomers_zero{};
	std::vector<Particle> Monomers_t{}; 
    for (unsigned i = 0; i < Monomers; i++) {
        Monomers_zero.push_back(Particle(i)); 
        Monomers_t.push_back(Particle(i)); 
    }
	
    Vector3d relPos {}; 
    double distance {}; 
    std::array<unsigned,3> Box {}; 
    Box[0] = Box[1] = Box[2] = BoxSize; 
    
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
    for (unsigned repl = 0; repl < Replicas; repl++){
        std::cout << "Replica "  << repl << std::endl; 
        ReplicaDir = Directory+"/REPL-"+std::to_string(repl);
        OutputDir = ReplicaDir+"/VanHove";
        ConfigFileStart = Directory+"/REPL-"+std::to_string(repl)+"/configs/config"; 
        std::cout << OutputDir << std::endl; 
        const int dir_err {mkdir(OutputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)};
        if (-1 == dir_err) {
            printf("Error creating directory!");
            return EXIT_FAILURE;
        }
        OutputFileStart = OutputDir + "/VanHove"; 
        std::map<unsigned,std::map<double,double>> VanHove {}; // for each t, there is a map (r, p(r))
        std::map<unsigned,unsigned> VanHove_count {}; // for each t, there is a map (r, p(r))
        //initialize the VanHoveFunctions 
	    for (unsigned i = 0; i < SampleSteps.size(); i++) {
	        VanHove[SampleSteps[i]] = std::map<double, double>(); 
	        VanHove_count[SampleSteps[i]] = 0; 
	    }
        T0Step = StartStep; 
        while(true) {
            ConfigFile = ConfigFileStart+std::to_string(T0Step)+".pdb";
		    std::cout << "ConfigFile T0: " << ConfigFile << std::endl; 
		    if (!initializePositions(Monomers_zero, ConfigFile)) {
		        std::cout << "problem with initializing monomers" << std::endl;
		        break; 
		    }
            for (auto& Step : SampleSteps) {
                ConfigFile = ConfigFileStart+std::to_string(T0Step+Step)+".pdb";
		        std::cout << "ConfigFile T1: " << ConfigFile << std::endl; 
		        if (!initializePositions(Monomers_t, ConfigFile)) {
		            std::cout << "problem with initializing monomers" << std::endl;
		            break; 
		        }
                for (unsigned i = 0; i < Monomers; i++) {
                    for (unsigned j = 0; j < Monomers; j++) {
                        relPos = relative(Monomers_zero[i], Monomers_t[j], Box, 0.0); 
                        distance = relPos.norm(); 
                        bin = unsigned(distance/DR)*DR; 
                        VanHove[Step][bin]++; 
                        VanHove_count[Step]++; 
                    }
                }
            }
            T0Step += SamplingStep; 
        }
        for (auto& Step : SampleSteps) {
            OutputFileName = OutputFileStart + std::to_string(Step*DeltaT); 
            std::ofstream OutputFile (OutputFileName, std::ios::out | std::ios::trunc); 
            for (auto& corr : VanHove[Step]) {
                OutputFile << corr.first << " " << corr.second/(Monomers*VanHove_count[Step]) << std::endl; 
            } 
            OutputFile.close();  
        }
        
    }
    gettimeofday(&end, NULL); 
    double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    std::cout << "total time: " << realTime << std::endl;
    
    return 0; 
}
