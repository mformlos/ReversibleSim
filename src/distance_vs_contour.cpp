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
    std::string Directory{}, ReplicaDir{}, OutputDir {}, MoleculeFileName{}, OutputFileName {}, ConfigFile{}, ConfigFileStart{}; 
    unsigned Replicas{}, StartStep{}, Step{}, SampleStep{}, Count{0}, first{}, second{}, NumberOfMolecules{0};
    
    std::map<unsigned, double> distance_contour; 
    
    if (argc != 6) {
        std::cout << "usage: ./distance_contour DIRECTORY MOLECULEFILE REPLICAS STARTSTEP SAMPLESTEP" << std::endl;  
        return EXIT_FAILURE; 
    }
    
    Directory = argv[1];
    MoleculeFileName = argv[2];
	Replicas = std::stoi(argv[3]);  
    StartStep = std::stoi(argv[4]);
    SampleStep = std::stoi(argv[5]); 
    
    std::cout << "Directory: " << Directory << " Replicas: " << Replicas << "  MoleculeFile: " << MoleculeFileName <<  std::endl;
	std::cout << "StartStep: " << StartStep << " SamplingStep: " << SampleStep  << std::endl;
	
	//setup Molecules
	System Sys(200, 200, 200, 0.0, 0.0, false);
	
	if (!Sys.addMolecules(MoleculeFileName, 1.0)) {
        std::cout << "MoleculeFile does not exist!" << std::endl; 
        return EXIT_FAILURE;     
    }
	/*std::ifstream MoleculeFile(MoleculeFileName, ios::in);
	std::string line;
    unsigned current_mol{0}, mol, monos, monos_total{0};
	if (!file.is_open()) return false;
	while (file >> mol >> monos) {
	    Molecules.push_back(Molecule(monos, mass, 0)); 
	    NumberOfMolecules++;   
	    monos_total += monos; 
	}
	std::cout << "initialized " << Molecules.size() << " molecules with a total of " << mono_total << " monomers." << std::endl;*/
	/////
	Vector3d relPos{};
	double dissq{};  
	
	for (unsigned repl = 0; repl < Replicas; repl++) {
        std::cout << "Replica "  << repl << std::endl; 
        ReplicaDir = Directory+"/REPL-"+std::to_string(repl);
        Step = StartStep; 
        ConfigFileStart = Directory+"/REPL-"+std::to_string(repl)+"/configs/config";
        while (true) {
	        ConfigFile = ConfigFileStart+std::to_string(Step)+".pdb";
		    std::cout << ConfigFile << std::endl;
		    if (!Sys.initializePositionsPDB(ConfigFile)) {
		        std::cout << "problem with initializing monomers" << std::endl;
		        break; 
		    } 
		    for (auto& mol : Sys.Molecules) {
		        for (unsigned i = 0; i < mol.NumberOfMonomers -1; i++) {
		            for (unsigned j = i+1; j < mol.NumberOfMonomers; j++) {
		                relPos = mol.Monomers[i].Position - mol.Monomers[j].Position;   
		                dissq = relPos.squaredNorm(); 
		                distance_contour[j-i] += dissq;  
		            }    
		        }
		    }
		    Count++; 
		    Step += SampleStep; 
		}
    }
    
    OutputFileName = Directory+"/distance_contour"; 
    std::ofstream OutputFile (OutputFileName, std::ios::out | std::ios::trunc); 
    std::cout << OutputFileName << std::endl; 
    
    for (auto& contour : distance_contour) {
        OutputFile << contour.first << " " << sqrt(contour.second/(Count*Sys.Molecules.size()*(Sys.Molecules[0].NumberOfMonomers-contour.first))) << std::endl; 
        //std::cout << contour.first << " " << contour.second/Count << std::endl; 
    }
    OutputFile.close();
    return 0; 	
}
