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
    std::string Directory{}, ConfigFile{}, ConfigFileStart{}, OutputFileNameStart{}, OutputFileName{}, StepFile{}, OutputQFileNameStart{}, OutputQFileName{},  WQFile {}, IncoherentDir{}, IncoherentFileName{}, IncoherentFileNameStart{}; 
    double DeltaTSim{}, CurrentDt{}; 
    unsigned StartStep{}, CurrentDtStep{}, NumberOfMonomers{}, BoxLength {}, FileStep{};   
    
    std::vector<Particle> Monomers_zero{}; 
    std::vector<Particle> Monomers_t{}; 
    std::vector<Particle> Monomers_tlast{};
    std::vector<unsigned long long> StepVector{}; 
    std::vector<unsigned long long>::iterator StepVectorIterator{}; 
    
    std::vector<double> q {}; 
    std::map<double, std::map<double, double>> incoherent_scattering_function {}; // for each q, there is a map<t, F(q,t)>
    std::map<double, int> incoherent_count {}; // for each dt, there is a count map<t, count(t, q)>  
   
    
    if (argc != 8) {
        std::cout << "usage: ./msd_stepfile DIRECTORY STARTSTEP STEPFILE MONOMERS DELTAT BOXLENGTH WQFILE" << std::endl;  
        return EXIT_FAILURE; 
    }      
    
    Directory = argv[1]; 
    StartStep = std::stoi(argv[2]); 
    StepFile = argv[3];
    NumberOfMonomers = std::stoi(argv[4]);
    DeltaTSim = std::stod(argv[5]); 
    BoxLength = std::stoi(argv[6]);
    WQFile = argv[7];
     
    std::cout << "Directory: " << Directory << " StepFile: " << StepFile << std::endl;
	std::cout << "StartStep: " << StartStep << std::endl;
    
    initializeStepVector(StepVector, StepFile); 
    
    for(auto& s : StepVector) s -= StartStep; 
    
    IncoherentDir = Directory+"/Sqt";
    int dir_err {mkdir(IncoherentDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)};
    if (-1 == dir_err) {
        printf("Error creating directory!");
        return EXIT_FAILURE;
    }
    IncoherentFileNameStart=IncoherentDir+"/Sqt"; 
    
    if (!initializeDoubleVector(q, WQFile)){
        std::cout << "Problem with WQFile, exiting!" << std::endl; 
        return EXIT_FAILURE; 
    } 
    
    
    OutputFileNameStart = Directory+"/MSD/MSD";
    OutputQFileNameStart = Directory+"/MQD/MQD";
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
    std::map<double,double> MQD; 
    
    std::cout << "Current Step 0: " << StartStep << std::endl;  
    StepVectorIterator = StepVector.begin(); 
    OutputFileName = OutputFileNameStart + std::to_string(StartStep); 
    OutputQFileName = OutputQFileNameStart + std::to_string(StartStep);
    std::ofstream OutputFile (OutputFileName, std::ios::out | std::ios::trunc);
    std::ofstream OutputQFile (OutputQFileName, std::ios::out | std::ios::trunc);
    if (!OutputFile.is_open()) {
        std::cout << "Cannot open output file, directory 'MSD' might not exist. Exiting..." << std::endl; 
        return EXIT_FAILURE; 
    }
    
    if (!OutputQFile.is_open()) {
        std::cout << "Cannot open output file, directory 'MQD' might not exist. Exiting..." << std::endl; 
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
	
	double distance{0.0}, qr {0.0}, bessel {0.0}; 
	
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
	    
	    double msd_current {0.0}, mqd_current {0.0}, square{0.0}; 
	    std::vector<double> wq_current (q.size(), 0.0);
	    Vector3d dist {Vector3d::Zero()};  
	    for (unsigned j = 0; j < NumberOfMonomers; j++) {
	        dist = Monomers_zero[j].Position - COMPosZero - Monomers_t[j].Position + COMPosT; 
	        square = dist.squaredNorm();
	        distance = sqrt(square);
            msd_current += square;
            square *= square; 
            mqd_current += square; 
            for (unsigned k = 0; k < q.size(); k++) {
                qr = q[k]*distance;
                bessel = sin(qr)/qr;  
                wq_current[k] += bessel;   
            }  
	        
	    }
	    for (unsigned k = 0; k < q.size(); k++){
            incoherent_scattering_function[k][CurrentDt] += wq_current[k]/NumberOfMonomers; 
        }
        incoherent_count[CurrentDt]++;
	    MSD[CurrentDt] = msd_current/NumberOfMonomers;		    
	    MQD[CurrentDt] = mqd_current/NumberOfMonomers;	
	    StepVectorIterator++; 
        Monomers_tlast = Monomers_t; 
	}
	for (auto& m : MSD) {
        OutputFile << m.first << " " << m.second << std::endl; 
    }          
    OutputFile.close(); 
    
    for (auto& m : MQD) {
        OutputQFile << m.first << " " << m.second << std::endl; 
        //OutputQFile << m.first << " " << (3./5.)*(m.second/pow(MSD[m.first],2)) -1. << std::endl; 
    }          
    OutputQFile.close(); 
    
    for (unsigned k = 0; k < q.size(); k++) {
	    IncoherentFileName = IncoherentFileNameStart+std::to_string(q[k]);
	    std::ofstream OutputFile (IncoherentFileName, std::ios::out | std::ios::trunc);
	    for (auto& corr :  incoherent_scattering_function[k]) {
	        if (corr.first == 0.0) continue; 
	        OutputFile << corr.first << " " << corr.second/incoherent_count[corr.first] << std::endl;
	    } 
	    OutputFile.close();
	}
    
    gettimeofday(&end, NULL); 
    double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    std::cout << "total time: " << realTime << std::endl;
    
    return 0;  

}
