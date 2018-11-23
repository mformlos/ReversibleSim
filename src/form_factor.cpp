/*
 * radial_distribution.cpp
 *
 *  Created on: Jul 4, 2018
 *      Author: maud
 */

#include <vector>
#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <math.h>
#include "System.h"
#include "HelperFunctions.h"

int main(int argc, char* argv[]) {
    std::string Directory{}, ParameterFile {}, MoleculeFile{}, ConfigFile{}, ConfigFileStart{}, WQFile {};
    int StartStep{}, Step{}, SamplingStep{};
    unsigned NumberOfMonomers {}, Replicas{}, BoxSize {}; 
    std::vector<Particle> Monomers{}; 
    std::vector<std::string> KVecFileList {}; 
    bool ParameterInitialized {}; 
    
    std::vector<double> q {}; 

    unsigned count {0}; 
    
    
    if (argc != 8) {
        std::cout << "usage: ./Form_factor RUNPARENTDIR REPLICAS PARAMETERFILE MOLECULEFILE WQFILE STARTSTEP SAMPLINGSTEP" << std::endl;  
    }    
    
    Directory=argv[1];
	Replicas = std::stoi(argv[2]);
	ParameterFile = argv[3];
	MoleculeFile = argv[4];
	WQFile = argv[5];
	StartStep = std::stoi(argv[6]);
	SamplingStep = std::stoi(argv[7]);
	
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
    /*MoleculeFile = extractParameter<std::string>("MoleculeFile", inputfile, ParameterInitialized);      
    if (!ParameterInitialized) return EXIT_FAILURE;
    MoleculeFile = ParameterDir+MoleculeFile; 
    */
    std::cout << "MoleculeFile : " << MoleculeFile << std::endl; 
    
    inputfile.close(); 
    
    if (!initializeDoubleVector(q, WQFile)){
        std::cout << "Problem with WQFile, exiting!" << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::vector<double> wq (q.size(), 0.0); 
    
    std::cout << "Directory: " << Directory << std::endl;
	std::cout << "StartStep: " << StartStep << "   SamplingStep: " << SamplingStep << " Number of qs sampled: " << q.size() << std::endl;
	std::cout << "Replicas: " << Replicas << "  MoleculeFile: " << MoleculeFile << "  BoxSize: " << BoxSize << std::endl; 
	
	System Sys(BoxSize, BoxSize, BoxSize, 0.0, 0.0, true);
	
	if (!Sys.addMolecules(MoleculeFile, 1.0)) {
        std::cout << "MoleculeFile does not exist!" << std::endl; 
        return EXIT_FAILURE;     
    }
    
    NumberOfMonomers = Sys.NumberOfParticles(); 
    std::cout << "Total Monomers: " << NumberOfMonomers << std::endl; 
    
    std::string name;
	name = Directory+"form_factor";
	std::ifstream test(name);
	if (test.good()){
		std::cout << "file " << name << " already exists! Aborting..." << std::endl;
		return EXIT_FAILURE;
	}
	std::ofstream output(name, std::ios::out | std::ios::trunc);
	output.precision(8);   
	
	timeval start {}, end {};
    gettimeofday(&start, NULL); 
	
	Vector3d relPos{}; 
	double distance{}, qr {}, bessel {}; 
	for (unsigned repl = 0; repl < Replicas; repl++) {
	    Step = StartStep; 
	    ConfigFileStart = Directory+"REPL-"+std::to_string(repl)+"/configs/config"; 
	    while (true) {
	        ConfigFile = ConfigFileStart+std::to_string(Step)+".pdb";
		    std::cout << ConfigFile << std::endl;
		    if (!Sys.initializePositionsPDB(ConfigFile)) {
		        std::cout << "problem with initializing monomers" << std::endl;
		        break; 
		    } 
		    std::vector<double> wq_current (q.size(), 0.0);
		    for (auto& mol : Sys.Molecules) {
		        for (unsigned i = 0; i < mol.NumberOfMonomers-1; i++) {
		            for (unsigned j = i+1; j < mol.NumberOfMonomers; j++) {
		                relPos = mol.Monomers[i].Position - mol.Monomers[j].Position; 
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
		        wq[k] += wq_current[k]/NumberOfMonomers; 
		    }
		    count++;
		    Step += SamplingStep; 
	    }
	}
	
	gettimeofday(&end, NULL); 
    double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    std::cout << "total time: " << realTime << " , time per conformation: " << realTime/count << std::endl;
	
	for (unsigned k = 0; k < q.size(); k++) {
	    output << q[k] << " " << 1.0+(2.0*wq[k]/(double)count) << std::endl; 
	}
	
	output.close(); 
    	
    
}
