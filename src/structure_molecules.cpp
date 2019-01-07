#include <vector>
#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <math.h>
#include <sys/stat.h>
#include "System.h"
#include "HelperFunctions.h"


int main(int argc, char* argv[]) {
    std::string Directory{}, OutputDir {}, StepFile {}, MoleculeFileName{}, ConfigFile{}, ConfigFileStart{}; 
    unsigned NumberOfMolecules{0};
   
    std::vector<unsigned long long> StepVector{}; 
    std::vector<unsigned long long>::iterator StepVectorIterator{}; 
    
    if (argc != 4) {
        std::cout << "usage: ./structure_averages DIRECTORY MOLECULEFILE STEPFILE" << std::endl;  
        return EXIT_FAILURE; 
    }
    
    Directory = argv[1];
    MoleculeFileName = argv[2];
    StepFile = argv[3];
    
    std::cout << "Directory: " << Directory << "  MoleculeFile: " << MoleculeFileName <<  std::endl;
	std::cout << " StepFile: " << StepFile << std::endl;
	
	initializeStepVector(StepVector, StepFile); 
	System Sys(200, 200, 200, 0.0, 0.0, false);
	
	if (!Sys.addMolecules(MoleculeFileName, 1.0)) {
        std::cout << "MoleculeFile does not exist!" << std::endl; 
        return EXIT_FAILURE;     
    }
    
    NumberOfMolecules = Sys.Molecules.size(); 
    
    std::cout << "Number of Molecules to analyse: " << NumberOfMolecules << std::endl; 
    
    OutputDir = Directory + "/data";
    const int dir_err {mkdir(OutputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)};
    if (-1 == dir_err) {
        printf("Error creating directory!");
        return EXIT_FAILURE;
    }
    
    std::vector<std::string> OutputFileNames;
    std::vector<std::ofstream*> OutputFiles;   
    
    for (unsigned i = 0; i < NumberOfMolecules; i++) {
        OutputFileNames.push_back(OutputDir+"/stats-Mol-"+std::to_string(i)); 
        std::ofstream* f = new std::ofstream(OutputFileNames[i], std::ios::out | std::ios::trunc);
        OutputFiles.push_back(f);
    }
    
    //// Measurement Variables
    Vector3d COMPos {Vector3d::Zero()}; 
    Matrix3d GyrTensor {Matrix3d::Zero()};
    std::tuple<double, Matrix3d> GyrTuple {}; 
    
    ConfigFileStart = Directory+"/configs/config";
    
    for (auto& Step : StepVector) {
        ConfigFile = ConfigFileStart+std::to_string(Step)+".pdb";
        if (!Sys.initializePositionsPDB(ConfigFile)) {
	        std::cout << "problem with initializing monomers" << std::endl;
	        break; 
	    } 
	    std::cout << Step << std::endl; 
	    for (unsigned i = 0; i < NumberOfMolecules; i++) {
	        GyrTuple = Sys.Molecules[i].GyrationTensor(); 
	        GyrTensor = std::get<1>(GyrTuple); 
	        OutputFiles[i] -> precision(10); 
	        OutputFiles[i] -> width(16); 
	        *OutputFiles[i] << Step << " "; 
            OutputFiles[i] -> precision(6); 
            OutputFiles[i] -> width(14); 
            *OutputFiles[i] << std::get<0>(GyrTuple); 
            for (unsigned m = 0; m < 3; m++) {
                for (unsigned n = m; n < 3; n++) {
                    OutputFiles[i] -> width(14); 
                    *OutputFiles[i] << GyrTensor(m,n) << " ";
                }
            }
            *OutputFiles[i] << std::endl; 
	    }
		        
    }
	
	for (unsigned i = 0; i < NumberOfMolecules; i++) {
	    OutputFiles[i] -> close();
	    delete OutputFiles[i];
	}
    
    return 0; 
}
