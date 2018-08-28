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
    std::string Directory{}, ConfigFile{}, ConfigFileStart{}, KVecFileListName {}, CoherentDir {}, CoherentFileNameStart{}, CoherentFileName{};
    int StartStep{}, MaxStep{}, Step{}, DStep{}, SamplingStep{}, T0Step{}, T1Step{};
    unsigned NumberOfMonomers {}, count{};
    double BoxSize {}, DeltaT {}, Time{}, realTime{}; 
    std::vector<Particle> Monomers_zero{}; 
    std::vector<Particle> Monomers_t{}; 
    std::vector<std::string> KVecFileList {}; 
    
    std::map<double, std::map<double, double>> coherent_scattering_function {}; // for each q, there is a map<t, F(q,t)>
    std::map<double, int> coherent_count {}; // for each dt, there is a count map<t, count(t, q)>  
    
    if (argc != 10) {
        std::cout << "usage: ./coherent_scattering RUNPARENTDIR NUMBEROFMONOMERS BOXSIZE STARTSTEP MAXSTEP DELTASTEP SAMPLINGSTEP KVECFILELIST DELTAT" << std::endl;
        return EXIT_FAILURE; 
    } 
    
    Directory=argv[1];
	NumberOfMonomers = std::stoi(argv[2]);
	BoxSize = std::stod(argv[3]); 
	StartStep = std::stoi(argv[4]);
	MaxStep = std::stoi(argv[5]);
	DStep = std::stoi(argv[6]);
	SamplingStep = std::stoi(argv[7]);
	KVecFileListName = argv[8]; 
	DeltaT = std::stod(argv[9]);
	
	std::cout << "Directory: " << Directory << std::endl;
	std::cout << "StartStep: " << StartStep << "   SamplingStep: " << SamplingStep << std::endl;
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
    
    
	CoherentDir = Directory+"/Fqt";
	int dir_err {mkdir(CoherentDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)};
    if (-1 == dir_err) {
        printf("Error creating directory!");
        return EXIT_FAILURE;
    }
	double Kabs {}, sum {0.0}, constant {2.*M_PI/BoxSize}; 
	
	std::vector<std::vector<Vector3d>> QKVecs {}; // for each q, there is a list of KVectors 
	std::vector<double> QKabs {}; // for each q, there is a Kabs; 
	for (auto& KVecFileName : KVecFileList) {
	    QKVecs.push_back(std::vector<Vector3d>()); 
        initializeKVectors(QKVecs.back(),Kabs, KVecFileName); 
        Kabs *= constant;
        QKabs.push_back(Kabs); 
        //std::cout << KVecFileName << " " << Kabs<< " " << QKVecs.back().size() << std::endl; 
	
	}
	
	std::vector<std::vector<double>> cos_zero(QKVecs.size(), std::vector<double>()); 
	std::vector<std::vector<double>> sin_zero(QKVecs.size(), std::vector<double>()); 
	
	for (unsigned i = 0; i < QKVecs.size(); i++) {
	    cos_zero[i] = std::vector<double>(QKVecs[i].size(), 0.0);
	    sin_zero[i] = std::vector<double>(QKVecs[i].size(), 0.0); 
	}
	
	CoherentFileNameStart=CoherentDir+"/Fqt_q";
	
	
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
    Step = StartStep;
    ConfigFileStart = Directory+"/configs/config"; 
    
    T0Step = StartStep;   
     
    while (true) {
        ConfigFile = ConfigFileStart+std::to_string(T0Step)+".pdb";    
        std::cout << "T0Step: " << T0Step << std::endl; 
        if (!initializePositions(Monomers_zero, ConfigFile)) {
            std::cout << "problem with initializing monomers" << std::endl;
            break; 
        }
        for (auto& mono : Monomers_zero) {
            wrapUniformNoShear(mono, BoxSize);
        }
        for (unsigned i = 0; i < QKVecs.size(); i++) {
            std::fill(cos_zero[i].begin(), cos_zero[i].end(), 0.0); 
            std::fill(sin_zero[i].begin(), sin_zero[i].end(), 0.0); 
            for (unsigned j = 0; j < QKVecs[i].size(); j++) {
                for (auto& mono: Monomers_zero) {
                    cos_zero[i][j]+= cos(QKVecs[i][j].dot(mono.Position)*constant); 
                    sin_zero[i][j]+= sin(QKVecs[i][j].dot(mono.Position)*constant); 
                }
            }
        }
        	        
        T1Step = T0Step;
        Step = 0; 
        count = 0;
        //gettimeofday(&start, NULL); 
        while (Step <= MaxStep) {
            Time = Step*DeltaT; 
            ConfigFile = ConfigFileStart+std::to_string(T1Step)+".pdb";
            //std::cout << "T1Step: " << T1Step << std::cout << std::endl; 
            //gettimeofday(&start, NULL); 
            if (!initializePositions(Monomers_t, ConfigFile)) {
                std::cout << "problem with initializing monomers" << std::endl;
                break; 
            }
            for (auto& mono : Monomers_t) {
                wrapUniformNoShear(mono, BoxSize);
            }
            
            for (unsigned i = 0; i < QKVecs.size(); i++) {
                sum = 0.0; 
                std::vector<double> cos_t(QKVecs[i].size(), 0.0); 
                std::vector<double> sin_t(QKVecs[i].size(), 0.0);
                for (unsigned j = 0; j < QKVecs[i].size(); j++) {
                
                    for (auto& mono: Monomers_t) {
                        cos_t[j]+= cos(QKVecs[i][j].dot(mono.Position)*constant); 
                        sin_t[j]+= sin(QKVecs[i][j].dot(mono.Position)*constant); 
                    }
                    sum += cos_zero[i][j]*cos_t[j]+sin_zero[i][j]*sin_t[j]; 
                }
                
                coherent_scattering_function[QKabs[i]][Time] += sum/QKVecs[i].size(); 
            }
            coherent_count[Time]++; 
                
            T1Step += DStep; 
            Step += DStep;
            count++; 
            //gettimeofday(&end, NULL);
            //realTime =  ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 ;
            //std::cout << "Time between file openinings: " << realTime << std::endl;    
        } 
        T0Step += SamplingStep; 
        //gettimeofday(&end, NULL);
        //realTime =  ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 ;
        //std::cout << "Time between file openinings: " << realTime/count << std::endl;          
    }
    
    for (unsigned i = 0; i < QKVecs.size(); i++) {
        CoherentFileName = CoherentFileNameStart+std::to_string(QKabs[i]); 
        std::ofstream OutputFile (CoherentFileName, std::ios::out | std::ios::trunc);
        for (auto& corr :  coherent_scattering_function[QKabs[i]]) {
            OutputFile << corr.first << " " << corr.second / coherent_count[corr.first] << std::endl; 
        }
        OutputFile.close();
    }
    gettimeofday(&end, NULL); 
    realTime =  ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 ;
    std::cout << "total time: " << realTime << std::endl;
    return 0; 
}
