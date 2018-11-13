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
    std::string Directory{}, ReplicaDir {}, ConfigFile{}, ConfigFileStart{}, StepSampleFile{}, SelfOutputDir{}, DistinctOutputDir {},  SelfOutputFileStart{}, DistinctOutputFileStart{}, SelfOutputFileName {}, DistinctOutputFileName {};
    double DR{0.1}, BoxSize{};
    unsigned StartStep{}, SamplingStep{}, Monomers{}, T0Step {}, T1Step{};    
    double DeltaT{}, bin{}; 
    
    std::vector<unsigned long long> SampleSteps {}; 
    
    
    if (argc != 9) {
            std::cout << "usage: ./Van_Hove_Distinct DIRECTORY STEPFILE STARTSTEP SAMPLINGSTEP DELTAT MONOMERS BOXSIZE INTERVAL" << std::endl;  
            return EXIT_FAILURE; 
    }
    
    Directory=argv[1];
    StepSampleFile=argv[2]; 
	StartStep = std::stoi(argv[3]);
    SamplingStep = std::stoi(argv[4]); 
	DeltaT = std::stod(argv[5]); 
	Monomers = std::stoi(argv[6]); 
	BoxSize = std::stod(argv[7]); 
	DR = std::stod(argv[8]); 
	
	std::cout << "Directory: " << Directory << std::endl;
	std::cout << "StartStep: " << StartStep << "   SamplingStep: " << SamplingStep << std::endl;
	std::cout << " DR: " << DR << std::endl;    
	
	initializeStepVector(SampleSteps, StepSampleFile);  
	std::vector<Particle> Monomers_zero{};
	std::vector<Particle> Monomers_t{}; 
	
    for (unsigned i = 0; i < Monomers; i++) {
        Monomers_zero.push_back(Particle(i)); 
        Monomers_t.push_back(Particle(i)); 
    }
	
    Vector3d relPos {}; 
    double distance {}; 
    std::array<double,3> Box {}; 
    Box[0] = Box[1] = Box[2] = BoxSize; 
    
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
  
    DistinctOutputDir = Directory+"/VanHoveDistinct";
    ConfigFileStart = Directory+"/configs/config"; 
    
   
    int dir_err {mkdir(DistinctOutputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)};
    if (-1 == dir_err) {
        printf("Error creating directory!");
        return EXIT_FAILURE;
    }
   
    DistinctOutputFileStart = DistinctOutputDir + "/VanHoveDistinct"; 
    
    double density{(Monomers*(Monomers-1)*0.5)/(pow((double)BoxSize,3)-4.*M_PI*pow(double(BoxSize)*0.5,3)/3.)};
    std::cout << "lost density: " << density << std::endl; 
   
    std::map<unsigned,std::map<double,double>> VanHoveDistinct {}; 
    std::map<unsigned,unsigned> VanHoveDistinct_count {}; 
    std::map<unsigned,unsigned> VanHove_count {}; 
    //initialize the VanHoveFunctions 
    for (unsigned i = 0; i < SampleSteps.size(); i++) {
        VanHoveDistinct[SampleSteps[i]] = std::map<double, double>();
        VanHoveDistinct_count[SampleSteps[i]] = 0; 
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
	        //std::cout << "ConfigFile T1: " << ConfigFile << std::endl; 
	        if (!initializePositions(Monomers_t, ConfigFile)) {
	            std::cout << "problem with initializing monomers" << std::endl;
	            break; 
	        }
            for (unsigned i = 0; i < Monomers; i++) {
                for (unsigned j = 0; j < Monomers; j++) {
                    if (i==j) continue; 
                    else {
                        relPos = relative(Monomers_zero[i], Monomers_t[j], Box, 0.0);
                        distance = relPos.norm(); 
                        if (distance <= (double)BoxSize*0.5) {
                            bin = round(distance/DR)*DR; 
                            VanHoveDistinct[Step][bin]++; 
                            VanHoveDistinct_count[Step]++;
                        }
                    }
                }
            }
            VanHove_count[Step]++; 
        }
        T0Step += SamplingStep; 
    }
    double normalization{}, rlast{0.0};  
    for (auto& Step : SampleSteps) {
        DistinctOutputFileName = DistinctOutputFileStart + std::to_string(Step*DeltaT); 
        std::ofstream DistinctOutputFile (DistinctOutputFileName, std::ios::out | std::ios::trunc);
        for (auto& corr : VanHoveDistinct[Step]) {
            normalization = 4.*M_PI*(pow(corr.first,3)-pow(rlast,3))/3.0;
            DistinctOutputFile << corr.first << " " << corr.second/(density*VanHove_count[Step]*normalization) << " " << corr.second/VanHoveDistinct_count[Step]<< std::endl; 
            rlast = corr.first;
        } 
        DistinctOutputFile.close();  
    }

    gettimeofday(&end, NULL); 
    double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    std::cout << "total time: " << realTime << std::endl;
    
    return 0; 
}
