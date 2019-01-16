#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <map>
#include <math.h>
#include <sys/stat.h>
#include "System.h"
#include "HelperFunctions.h"

int main(int argc, char* argv[]) {
    std::string Directory{}, ReplicaDir{},  OutputDir {}, OutputFileName{}, StepFile{}, MoleculeFile{}, ConfigFile{}, ConfigFileStart{}, makeDir{}; 
    unsigned StartStep{}, DStep {}, SampleStep{}, Replicas{}, Molecules{}, NumberOfMonomers{}; 
    double DeltaT{}, Time{}; 
    int dir_err{};
    
    std::vector<unsigned long long> StepVector{}; 
    std::vector<unsigned long long>::iterator StepVectorIterator{}; 
    
    std::map<double,double> MSD; 
    std::map<double,unsigned> count;
    std::map<double,double> MQD; 
    
    std::map<double,double> MSD_t; 
    std::map<double,unsigned> count_t;
    std::map<double,double> MQD_t; 
    
    std::ofstream OutputFile{}; 

    if (argc != 8) {
        std::cout << "usage: ./msd_new DIRECTORY MOLECULEFILE REPLICAS STARTSTEP SAMPLESTEP DELTAT STEPFILE " << std::endl;  
        return EXIT_FAILURE; 
    }
    
    Directory=argv[1];
    MoleculeFile=argv[2];
	Replicas = std::stoi(argv[3]); 
	StartStep = std::stoi(argv[4]);
    SampleStep = std::stoi(argv[5]); 
	DeltaT = std::stod(argv[6]); 
	StepFile = argv[7];
	
    std::cout << "Directory: " << Directory << " Replicas: " << Replicas << " MoleculeFile: " << MoleculeFile << std::endl;
	std::cout << "StartStep: " << StartStep << " SamplingStep: " << SampleStep << " StepFile: " << StepFile << std::endl;
	std::cout << "DeltaT: " << DeltaT << std::endl;

	
    initializeStepVector(StepVector, StepFile);
    System Sys0(200, 200, 200, 0.0, 0.0, false);
    System SysT(200, 200, 200, 0.0, 0.0, false);
    
    if (!Sys0.addMolecules(MoleculeFile, 1.0)) {
        std::cout << "MoleculeFile does not exist!" << std::endl; 
        return EXIT_FAILURE;     
    }
    if (!SysT.addMolecules(MoleculeFile, 1.0)) {
        std::cout << "MoleculeFile does not exist!" << std::endl; 
        return EXIT_FAILURE;     
    }
    Molecules = Sys0.Molecules.size(); 
    NumberOfMonomers = Sys0.NumberOfParticles(); 
    std::cout << "Total Monomers: " << NumberOfMonomers << std::endl; 
    
    
    int T0Step {}; 
    int T1Step {};
    Vector3d relPos {Vector3d::Zero()}; 
    Vector3d COMPos0 {Vector3d::Zero()};
    Vector3d COMPosT {Vector3d::Zero()};
    double square {}; 
    for (unsigned repl = 0; repl < Replicas; repl++) {
        std::cout << "Replica "  << repl << std::endl; 
        ReplicaDir = Directory+"/REPL-"+std::to_string(repl);
        ConfigFileStart = ReplicaDir+"/configs/config"; 
        T0Step = StartStep;
        /*makeDir=ReplicaDir+"/data";
        dir_err = mkdir(makeDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err) {
            printf("Error creating directory!n");
            exit(1);
        }*/
        makeDir=ReplicaDir+"/data/MSD";
        dir_err = mkdir(makeDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err) {
            printf("Error creating directory!n");
            exit(1);
        }
        makeDir=ReplicaDir+"/data/MQD";
        dir_err = mkdir(makeDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err) {
            printf("Error creating directory!n");
            exit(1);
        }
        /// loop over T0
        while (true) {
            std::cout <<"T0: " <<  T0Step << std::endl; 
            ConfigFile = ConfigFileStart+std::to_string(T0Step)+".pdb";
            //std::cout << ConfigFile << std::endl;
            if (!Sys0.initializePositionsPDB(ConfigFile)) {
		        std::cout << "problem with initializing monomers" << std::endl;
		        break; 
		    }
            MSD_t.clear();
            MQD_t.clear(); 
            count_t.clear();
            /// loop over t > 0
            StepVectorIterator = StepVector.begin(); 
            while (StepVectorIterator != StepVector.end()) {
                DStep = *StepVectorIterator; 
                T1Step = T0Step+DStep;     
                Time = DStep*DeltaT;
                ConfigFile = ConfigFileStart+std::to_string(T1Step)+".pdb";
                //std::cout << ConfigFile << std::endl;
                if (!SysT.initializePositionsPDB(ConfigFile)) {
		            std::cout << "problem with initializing monomers" << std::endl;
		            break; 
		        }
		        for (unsigned mol = 0; mol < Molecules; mol++) {
		            //COMPos0 = Sys0.Molecules[mol].centerOfMassPosition(); 
		            //COMPosT = SysT.Molecules[mol].centerOfMassPosition(); 
		            //std::cout << "mol: " << mol << std::endl; 
		            for (unsigned mono = 0; mono < Sys0.Molecules[mol].NumberOfMonomers; mono++) {
		                //std::cout << "mono: " << mono << std::endl; 
		                //relPos = Sys0.Molecules[mol].Monomers[mono].Position - COMPos0 - SysT.Molecules[mol].Monomers[mono].Position + COMPosT; 
		                relPos = Sys0.Molecules[mol].Monomers[mono].Position - SysT.Molecules[mol].Monomers[mono].Position;
		                square = relPos.squaredNorm(); 
		                MSD[Time] += square; 
		                MSD_t[Time] += square;
		                square *= square; 
		                MQD[Time] += square; 
		                MQD_t[Time] += square; 
		                count[Time]++;  
		                count_t[Time]++;   
		            }
		        }
		        
		        StepVectorIterator++;
            }
            //// Output for current T0
            OutputFileName = ReplicaDir+"/data/MSD/MSD-t-"+std::to_string(T0Step); 
            OutputFile.open(OutputFileName, std::ios::out | std::ios::trunc);
            for (auto& el : MSD_t) {
                OutputFile << el.first << " " << el.second/(count_t[el.first]) << std::endl; 
            }
            OutputFile.close(); 
            OutputFileName = ReplicaDir+"/data/MQD/MQD-t-"+std::to_string(T0Step); 
            OutputFile.open(OutputFileName, std::ios::out | std::ios::trunc);
            for (auto& el : MQD_t) {
                OutputFile << el.first << " " << el.second/(count_t[el.first]) << std::endl; 
            }
            OutputFile.close(); 
            T0Step += SampleStep; 
        }
    }
    
    OutputFileName = Directory+"/MSD"; 
    OutputFile.open(OutputFileName, std::ios::out | std::ios::trunc);
    for (auto& el : MSD) {
        OutputFile << el.first << " " << el.second/(count[el.first]) << std::endl; 
    }
    OutputFile.close(); 
    
    OutputFileName = Directory+"/MQD"; 
    OutputFile.open(OutputFileName, std::ios::out | std::ios::trunc);
    for (auto& el : MQD) {
        OutputFile << el.first << " " << el.second/(count[el.first]) << std::endl; 
    }
    OutputFile.close(); 
    
    return 0; 
}
