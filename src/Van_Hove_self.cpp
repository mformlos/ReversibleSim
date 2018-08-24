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
    std::string Directory{}, ReplicaDir {}, ConfigFile{}, ConfigFileStart{}, StepSampleFile{}, SelfOutputDir{}, SelfOutputFileStart{}, SelfOutputFileName {};
    double DR{0.1};
    unsigned StartStep{}, SamplingStep{}, Monomers{}, BoxSize{}, T0Step {};    
    double DeltaT{}, bin{}; 
    
    std::vector<unsigned long long> SampleSteps {}; 
    
    
    if (argc != 9) {
            std::cout << "usage: ./Van_Hove DIRECTORY STEPFILE STARTSTEP SAMPLINGSTEP DELTAT Monomers BOXSIZE INTERVAL" << std::endl;  
            return EXIT_FAILURE; 
    }
    
    Directory=argv[1];
    StepSampleFile=argv[2]; 
	StartStep = std::stoi(argv[3]);
    SamplingStep = std::stoi(argv[4]); 
	DeltaT = std::stod(argv[5]); 
	Monomers = std::stoi(argv[6]); 
	BoxSize = std::stoi(argv[7]); 
	DR = std::stod(argv[8]); 
	
	std::cout << "Directory: " << Directory << std::endl;
	std::cout << "StartStep: " << StartStep << "   SamplingStep: " << SamplingStep << std::endl;
	std::cout << " DR: " << DR << std::endl;    
	
	initializeStepVector(SampleSteps, StepSampleFile);  
	std::vector<Particle> Monomers_zero{};
	std::vector<Particle> Monomers_t{}; 
	std::vector<Particle> Monomers_tlast{}; 
	
    for (unsigned i = 0; i < Monomers; i++) {
        Monomers_zero.push_back(Particle(i)); 
        Monomers_t.push_back(Particle(i)); 
        Monomers_tlast.push_back(Particle(i));
    }
	
    Vector3d relPos {}; 
    double distance {}; 
    std::array<unsigned,3> Box {}; 
    Box[0] = Box[1] = Box[2] = BoxSize; 
    
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
  

    SelfOutputDir = Directory+"/VanHoveSelf";
    ConfigFileStart = Directory+"/configs/config"; 
    
    int dir_err {mkdir(SelfOutputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)};
    if (-1 == dir_err) {
        printf("Error creating directory!");
        return EXIT_FAILURE;
    }
    
    SelfOutputFileStart = SelfOutputDir + "/VanHoveSelf"; 
   
    
    std::map<unsigned,std::map<double,double>> VanHoveSelf {}; // for each t, there is a map (r, p(r))
    std::map<unsigned,unsigned> VanHove_count {}; 
    //initialize the VanHoveFunctions 
    for (unsigned i = 0; i < SampleSteps.size(); i++) {
        VanHoveSelf[SampleSteps[i]] = std::map<double, double>(); 
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
	    Monomers_tlast = Monomers_zero; 
        for (auto& Step : SampleSteps) {
            ConfigFile = ConfigFileStart+std::to_string(T0Step+Step)+".pdb";
	        //std::cout << "ConfigFile T1: " << ConfigFile << std::endl; 
	        if (!initializePositions(Monomers_t, ConfigFile)) {
	            std::cout << "problem with initializing monomers" << std::endl;
	            break; 
	        }
	        for (unsigned mono = 0; mono < Monomers; mono++) {
		        wrapBack(Monomers_tlast[mono], Monomers_t[mono], BoxSize); 
		    }
	        
            for (unsigned i = 0; i < Monomers; i++) { 
                relPos = Monomers_t[i].Position-Monomers_zero[i].Position; 
                distance = relPos.norm(); 
                bin = unsigned(distance/DR)*DR;
                VanHoveSelf[Step][bin]++;
            }
            Monomers_tlast = Monomers_t;
            VanHove_count[Step]++; 
        }
        T0Step += SamplingStep; 
    }
    for (auto& Step : SampleSteps) {
        SelfOutputFileName = SelfOutputFileStart + std::to_string(Step*DeltaT); 
        std::ofstream SelfOutputFile (SelfOutputFileName, std::ios::out | std::ios::trunc); 
        for (auto& corr : VanHoveSelf[Step]) {
            SelfOutputFile << corr.first << " " << corr.second/(Monomers*VanHove_count[Step]) << std::endl; 
        } 
        SelfOutputFile.close();
    }

    gettimeofday(&end, NULL); 
    double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    std::cout << "total time: " << realTime << std::endl;
    
    return 0; 
}
