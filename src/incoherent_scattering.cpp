/*
 * radial_distribution.cpp
 *
 *  Created on: Jul 4, 2018
 *      Author: maud
 */

#include <vector>
#include <sys/time.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <math.h>
#include "System.h"
#include "HelperFunctions.h"

int main(int argc, char* argv[]) {
    std::string Directory{}, ParameterFile {}, MoleculeFile{}, ConfigFile{}, ConfigFileStart{}, WQFile {}, IncoherentDir{}, IncoherentFileName{}, IncoherentFileNameStart{};
    int StartStep{}, MaxStep{}, Step{}, DStep{}, SamplingStep{}, T0Step {}, T1Step{};
    unsigned NumberOfMonomers {}, BoxSize {}; 
    std::vector<Particle> Monomers{}; 
    std::vector<std::string> KVecFileList {}; 
    bool ParameterInitialized {}; 
    double DeltaT{}, Time {}; 
    std::vector<double> q {}; 

    unsigned count {0}; 
    
    
    if (argc != 9) {
        std::cout << "usage: ./Incoherent_scattering RUNPARENTDIR PARAMETERFILE WQFILE STARTSTEP MAXSTEP DELTASTEP SAMPLINGSTEP DELTAT" << std::endl;  
        return EXIT_FAILURE; 
    }    
    
    Directory=argv[1];
	ParameterFile = argv[2];
	WQFile = argv[3];
	StartStep = std::stoi(argv[4]);
	MaxStep = std::stoi(argv[5]);
	DStep = std::stoi(argv[6]);
	SamplingStep = std::stoi(argv[7]);
	DeltaT = std::stod(argv[8]);
	
	std::ifstream inputfile(ParameterFile, std::ios::in);
    if (!inputfile.is_open()) {
        std::cout << "could not open file '" << argv[1] << "' , exiting" << std::endl; 
        return EXIT_FAILURE;  
    } 
    
    std::string search{"/"}; 
    std::size_t found{ParameterFile.find(search)};
    while (ParameterFile.find(search, found+1)!= std::string::npos) {
        found = ParameterFile.find(search, found+1); 
    }
    std::string ParameterDir = ParameterFile.substr(0, found+1); 
    
     
    
    
    BoxSize = extractParameter<unsigned>("BoxX", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    MoleculeFile = extractParameter<std::string>("MoleculeFile", inputfile, ParameterInitialized);      
    if (!ParameterInitialized) return EXIT_FAILURE;
    MoleculeFile = ParameterDir+MoleculeFile; 
    
    std::cout << "MoleculeFile : " << MoleculeFile << std::endl; 
    
    inputfile.close(); 
    
    if (!initializeDoubleVector(q, WQFile)){
        std::cout << "Problem with WQFile, exiting!" << std::endl; 
        return EXIT_FAILURE; 
    } 
    //std::vector<double> wq (q.size(), 0.0); 
    std::map<double, std::map<double, double>> incoherent_scattering_function {}; // for each q, there is a map<t, F(q,t)>
    std::map<double, int> incoherent_count {}; // for each dt, there is a count map<t, count(t, q)>  
    
    
    std::cout << "Directory: " << Directory << std::endl;
	std::cout << "StartStep: " << StartStep << "   SamplingStep: " << SamplingStep << " Number of qs sampled: " << q.size() << std::endl;
	std::cout << "  MoleculeFile: " << MoleculeFile << "  BoxSize: " << BoxSize << std::endl; 
	
	System SysZero(BoxSize, BoxSize, BoxSize, 0.0, 0.0, true);
	System SysT(BoxSize, BoxSize, BoxSize, 0.0, 0.0, true);
	
	if (!SysZero.addMolecules(MoleculeFile, 1.0)) {
        std::cout << "MoleculeFile does not exist!" << std::endl; 
        return EXIT_FAILURE;     
    }
    
    if (!SysT.addMolecules(MoleculeFile, 1.0)) {
        std::cout << "MoleculeFile does not exist!" << std::endl; 
        return EXIT_FAILURE;     
    }
    
    NumberOfMonomers = SysZero.NumberOfParticles(); 
    std::cout << "Total Monomers: " << NumberOfMonomers << std::endl; 
    
    IncoherentDir = Directory+"/Sqt";
    
    int dir_err {mkdir(IncoherentDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)};
    if (-1 == dir_err) {
        printf("Error creating directory!");
        return EXIT_FAILURE;
    }

	  
	
	timeval start {}, end {};
    gettimeofday(&start, NULL); 
	
	Vector3d relPos{}; 
	double distance{}, qr {}, bessel {}; 
    //Step = StartStep; 
    T0Step = StartStep;   
    ConfigFileStart = Directory+"/configs/config"; 
    while (true) {
        ConfigFile = ConfigFileStart+std::to_string(T0Step)+".pdb";
	    std::cout << "T0: " << ConfigFile << std::endl;
	    if (!SysZero.initializePositionsPDB(ConfigFile)) {
	        std::cout << "problem with initializing monomers" << std::endl;
	        break; 
	    } 
	    T1Step = T0Step;
        Step = 0; 
        count = 0;
        while (Step <= MaxStep) {
            Time = Step*DeltaT; 
            ConfigFile = ConfigFileStart+std::to_string(T1Step)+".pdb";
            std::cout << "T1: " << ConfigFile << std::endl;
            if (!SysT.initializePositionsPDB(ConfigFile)) {
	            std::cout << "problem with initializing monomers" << std::endl;
	            break; 
	        } 
	        
	        std::vector<double> wq_current (q.size(), 0.0);
	        for (unsigned mol = 0; mol < SysZero.Molecules.size(); mol++) {
	            for (unsigned i = 0; i < SysZero.Molecules[mol].NumberOfMonomers; i++) {
	                for (unsigned j = 0; j < SysT.Molecules[mol].NumberOfMonomers; j++) {
	                    relPos = SysT.Molecules[mol].Monomers[j].Position - SysZero.Molecules[mol].Monomers[i].Position; 
	                    distance = relPos.norm(); 
	                    for (unsigned k = 0; k < q.size(); k++) {
	                        qr = q[k]*distance;
	                        bessel = sin(qr)/qr; 
	                        wq_current[k] += bessel;      
	                    }    
	                }
	            }    
	        }
	        for (unsigned k = 0; k < q.size(); k++){
	            incoherent_scattering_function[k][Time] += wq_current[k]/NumberOfMonomers; 
	        }
	        incoherent_count[Time]++;
	        T1Step += DStep; 
            Step += DStep;
            count++; 
        }    
	    T0Step += SamplingStep; 
    }

	
	
	
	for (unsigned k = 0; k < q.size(); k++) {
	    IncoherentFileName = IncoherentFileNameStart+std::to_string(q[k]);
	    std::ofstream OutputFile (IncoherentFileName, std::ios::out | std::ios::trunc);
	    for (auto& corr :  incoherent_scattering_function[k]) {
	        OutputFile << corr.first << " " << 1.0+(2.0*corr.second/incoherent_count[corr.first]) << std::endl;
	        OutputFile.close();
	    } 
	}
   	gettimeofday(&end, NULL); 
    double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    std::cout << "total time: " << realTime << " , time per conformation: " << realTime/count << std::endl;
    return 0; 
    
}
