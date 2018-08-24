#include <vector>
#include <sys/time.h>
#include <sys/stat.h>
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
    std::string Directory{}, ConfigFile{}, ConfigFileStart{}, KVecFileListName {}, CoherentDir {}, CoherentFileNameStart{}, CoherentFileName{}, StepFile{};
    int StartStep{}, Step{}, T0Step{}, T1Step{};
    unsigned NumberOfMonomers {};
    double BoxSize {}, DeltaT {}, Time{}; 
    std::vector<Particle> Monomers_zero{}; 
    std::vector<Particle> Monomers_t{}; 
    std::vector<std::string> KVecFileList {}; 
    
    std::map<double, std::map<double, double>> coherent_scattering_function {}; // for each q, there is a map<t, F(q,t)>
    std::map<double, int> coherent_count {}; // for each dt, there is a count map<t, count(t, q)>  
    std::vector<unsigned long long> StepVector{}; 
    std::vector<unsigned long long>::iterator StepVectorIterator{}; 
    
    
    if (argc != 8) {
        std::cout << "usage: ./Coherent_Scattering_stepfile RUNPARENTDIR NUMBEROFMONOMERS BOXSIZE STARTSTEP STEPFILE KVECFILELIST DELTAT" << std::endl;
        return EXIT_FAILURE; 
    } 
    
    Directory=argv[1];
	NumberOfMonomers = std::stoi(argv[2]);
	BoxSize = std::stod(argv[3]); 
	StartStep = std::stoi(argv[4]);
	StepFile = argv[5];
	KVecFileListName = argv[6]; 
	DeltaT = std::stod(argv[7]);
	
	std::cout << "Directory: " << Directory << std::endl;
	std::cout << "StartStep: " << StartStep << "   StepFile: " << StepFile << std::endl;
	std::cout << "Number of monomers: " << NumberOfMonomers << "  BoxSize: " << BoxSize << std::endl; 
	std::cout << "KVec File List: " << KVecFileListName << std::endl; 
    
    for (unsigned i = 0; i < NumberOfMonomers; i++) {
        Monomers_zero.push_back(Particle(i)); 
        Monomers_t.push_back(Particle(i));
    }
    
    if (!fillConfigPool(KVecFileList, KVecFileListName)) {
        std::cout << "problem filling KVecFileList, exiting!" << std::endl; 
        return EXIT_FAILURE; 
    } 
    
    initializeStepVector(StepVector, StepFile); 
    for(auto& s : StepVector) s -= StartStep; 
    
	CoherentDir = Directory+"/Fqt";
	int dir_err {mkdir(CoherentDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)};
    if (-1 == dir_err) {
        printf("Error creating directory!");
        return EXIT_FAILURE;
    }
	
	CoherentFileNameStart=CoherentDir+"/Fqt_q";
	
	double Kabs {}, sum {0.0}, constant {2.*M_PI/BoxSize}; 
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
    for (auto& KVecFileName : KVecFileList) {
        std::vector<Vector3d> KVecs; 
        initializeKVectors(KVecs, Kabs, KVecFileName);
        Kabs *= constant; 
        std::cout << KVecFileName << " " << Kabs<< " " << KVecs.size() << std::endl; 
        Step = StartStep;
        ConfigFileStart = Directory+"/configs/config"; 
        
        T0Step = StartStep;     
        ConfigFile = ConfigFileStart+std::to_string(T0Step)+".pdb";    
        //std::cout << "T0Step: " << T0Step << std::endl; 
        if (!initializePositions(Monomers_zero, ConfigFile)) {
            std::cout << "problem with initializing monomers" << std::endl;
            break; 
        }
        for (auto& mono : Monomers_zero) {
            wrapUniformNoShear(mono, BoxSize);
        }
        std::vector<double> cos_zero(KVecs.size(),0.0); 
        std::vector<double> sin_zero(KVecs.size(), 0.0);
        for (unsigned i = 0; i < KVecs.size(); i++) {
            for (auto& mono: Monomers_zero) {
                cos_zero[i]+= cos(KVecs[i].dot(mono.Position))*constant; 
                sin_zero[i]+= sin(KVecs[i].dot(mono.Position))*constant; 
            }
        }
        StepVectorIterator = StepVector.begin(); 	        
        while (StepVectorIterator != StepVector.end()) {
            Step = *StepVectorIterator;
            Time = Step*DeltaT; 
            T1Step = T0Step + Step; 
            //std::cout << "T1Step: " << T1Step << std::endl; 
            ConfigFile = ConfigFileStart+std::to_string(T1Step)+".pdb";
            if (!initializePositions(Monomers_t, ConfigFile)) {
                std::cout << "problem with initializing monomers" << std::endl;
                break; 
            }
            for (auto& mono : Monomers_t) {
                wrapUniformNoShear(mono, BoxSize);
            }
            sum = 0.0; 
            for (unsigned i = 0; i < KVecs.size(); i++) {
                double cos_t{0.0}, sin_t{0.0}; 
                for (auto& mono: Monomers_t) {
                    cos_t+= cos(KVecs[i].dot(mono.Position))*constant; 
                    sin_t+= sin(KVecs[i].dot(mono.Position))*constant; 
                }
                sum += cos_zero[i]*cos_t+sin_zero[i]*sin_t; 
            }
            coherent_scattering_function[Kabs][Time] += sum/KVecs.size(); 
            coherent_count[Time]++; 
            StepVectorIterator++; 
        } 
        CoherentFileName = CoherentFileNameStart+std::to_string(Kabs); 
        std::ofstream OutputFile (CoherentFileName, std::ios::out | std::ios::trunc);
        for (auto& corr :  coherent_scattering_function[Kabs]) {
            OutputFile << corr.first << " " << corr.second / coherent_count[corr.first] << std::endl; 
    }
    OutputFile.close();
    }
    gettimeofday(&end, NULL); 
    double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    std::cout << "total time: " << realTime << std::endl;
    return 0; 
}
