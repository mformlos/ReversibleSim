#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <sys/time.h>
#include <map>
#include <math.h>
#include <sys/stat.h>
#include "Particle.h"
#include "Molecule.h"
#include "HelperFunctions.h"
#include "BoundaryConditions.h"



int main(int argc, char* argv[]) {
    std::string Directory{}, ConfigFile{}, ConfigFileStart{}, OutputFileNameStart{}, OutputFileName{}, StepFile{}; 
    double DeltaTSim{}, CurrentDt{}; 
    unsigned StartStep{}, CurrentDtStep{}, NumberOfMonomers{}, BoxLength {}, MaxStep{}, FileStep{};   
    
    std::vector<Particle> Monomers_zero{}; 
    std::vector<Particle> Monomers_t{}; 
    std::vector<Particle> Monomers_tlast{};
    std::vector<unsigned long long> StepVector{}; 
    std::vector<unsigned long long>::iterator StepVectorIterator{}; 
    
    if (argc != 7) {
        std::cout << "usage: ./msd_stepfile DIRECTORY STARTSTEP STEPFILE MONOMERS DELTAT BOXLENGTH" << std::endl;  
        return EXIT_FAILURE; 
    }      
    
    Directory = argv[1]; 
    StartStep = std::stoi(argv[2]); 
    StepFile = argv[3];
    NumberOfMonomers = std::stoi(argv[4]);
    DeltaTSim = std::stod(argv[5]); 
    BoxLength = std::stoi(argv[6]);
    
    std::cout << "Directory: " << Directory << " StepFile: " << StepFile << std::endl;
	std::cout << "StartStep: " << StartStep << std::endl;
    
    initializeStepVector(StepVector, StepFile); 
    
    for(auto& s : StepVector) s -= StartStep; 
    
   
    
    
    
    OutputFileNameStart = Directory+"/MSD/MSD";
    ConfigFileStart = Directory+"/configs/config";
    
    for (unsigned i = 0; i < NumberOfMonomers; i++) {
        Monomers_zero.push_back(Particle(i)); 
        Monomers_t.push_back(Particle(i));
        Monomers_tlast.push_back(Particle(i));
    }
    
    
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
    std::map<double,double> MSD; 
    std::map<double,double>::iterator MSD_iter; 
    
    std::cout << "Current Step 0: " << StartStep << std::endl;  
    StepVectorIterator = StepVector.begin(); 
    OutputFileName = OutputFileNameStart + std::to_string(StartStep); 
    std::ofstream OutputFile (OutputFileName, std::ios::out | std::ios::trunc);
    if (!OutputFile.is_open()) {
        std::cout << "Cannot open output file, directory 'MSD' might not exist. Exiting..." << std::endl; 
        return EXIT_FAILURE; 
    }
     
    ConfigFile = ConfigFileStart+std::to_string(StartStep)+".pdb";
    if (!initializePositions(Monomers_zero, ConfigFile)) {
	        std::cout << "problem with initializing monomers" << std::endl;
	        return EXIT_FAILURE;  
	}
	Vector3d COMPosZero {Vector3d::Zero()}; 
	for (auto& mono : Monomers_zero) {
	    COMPosZero += mono.Position; 
	}
	COMPosZero /= NumberOfMonomers; 
	
	Monomers_tlast = Monomers_zero; 
	
	while (StepVectorIterator != StepVector.end()) {
	    CurrentDtStep = *StepVectorIterator; 
	    CurrentDt = CurrentDtStep*DeltaTSim;
	    FileStep = StartStep+CurrentDtStep;
	    ConfigFile = ConfigFileStart+std::to_string(FileStep)+".pdb";
	    if (!initializePositions(Monomers_t, ConfigFile)) {
	        std::cout << "reached last step" << std::endl;
	        break; 
	    } 
	    
	    
	    for (unsigned mono = 0; mono < NumberOfMonomers; mono++) {
	        wrapBack(Monomers_tlast[mono], Monomers_t[mono], BoxLength); 
	    }
	    Vector3d COMPosT {Vector3d::Zero()}; 
	    for (auto& mono : Monomers_t) {
	        COMPosT += mono.Position; 
	    }
	    COMPosT /= NumberOfMonomers; 
	    
	    double msd_current {0.0}; 
	    for (unsigned j = 0; j < NumberOfMonomers; j++) {
	        msd_current += (Monomers_zero[j].Position - COMPosZero - Monomers_t[j].Position + COMPosT).squaredNorm(); 
	    }
	    MSD[CurrentDt] = msd_current/NumberOfMonomers;		    
	    StepVectorIterator++; 
        Monomers_tlast = Monomers_t; 
	}
	for (auto& m : MSD) {
        OutputFile << m.first << " " << m.second << std::endl; 
    }          
    OutputFile.close(); 

    
    gettimeofday(&end, NULL); 
    double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    std::cout << "total time: " << realTime << std::endl;
    
    return 0;  

}
