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
#include "HelperFunctions.h"

int main(int argc, char* argv[]) {
    std::string Directory{}, ReplicaDir{},  OutputDir {}, BondDir {}, BondFileStart{}, BondFileName{}, BondOutputFileStart{}, SpinOutputFileStart{}, BondOutputFileName{}, SpinOutputFileName{}, StepFile{}; 
    unsigned StartStep{}, MaxStep{}, Step{}, DStep {}, SampleStep{}, Replicas{}, FunctionalGroups{}, Molecules{}; 
    double DeltaT{}, Time{}; 
    std::map<double,double> Bond_corr_average;
    std::map<double, double>::iterator Bond_corr_average_iter; 
    std::map<double,unsigned> Bond_corr_average_count;
    std::map<double, unsigned>::iterator Bond_corr_average_count_iter; 
    
    std::map<double,double> Spin_corr_average;
    std::map<double, double>::iterator Spin_corr_average_iter; 
    std::map<double,unsigned> Spin_corr_average_count;
    std::map<double, unsigned>::iterator Spin_corr_average_count_iter; 
    
    std::vector<unsigned long long> StepVector{}; 
    std::vector<unsigned long long>::iterator StepVectorIterator{}; 
    
    if (argc != 8) {
        std::cout << "usage: ./cluster_corr DIRECTORY REPLICAS FUNCTIONALGROUPS STARTSTEP SAMPLESTEP DELTAT STEPFILE" << std::endl;  
        return EXIT_FAILURE; 
    }
    
    Directory=argv[1];
	Replicas = std::stoi(argv[2]); 
	FunctionalGroups = std::stoi(argv[3]);
	StartStep = std::stoi(argv[4]);
    SampleStep = std::stoi(argv[5]); 
	DeltaT = std::stod(argv[6]); 
	StepFile = argv[7];
	
	std::cout << "Directory: " << Directory << " Replicas: " << Replicas << " Functional groups: " << FunctionalGroups << std::endl;
	std::cout << "StartStep: " << StartStep << " SamplingStep: " << SampleStep << " StepFile: " << StepFile << std::endl;
	std::cout << "DeltaT: " << DeltaT << std::endl;
	
	initializeStepVector(StepVector, StepFile); 

    int T0Step {}; 
    int T1Step {};
    double Normalization {0.0}; 	
    unsigned count {0};
    for (unsigned repl = 0; repl < Replicas; repl++) {
        std::cout << "Replica "  << repl << std::endl; 
        ReplicaDir = Directory+"/REPL-"+std::to_string(repl);
        BondDir = ReplicaDir+"/Bonds"; 
        /*OutputDir = ReplicaDir+"/BondCorr"; 
        const int dir_err {mkdir(OutputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)};
        if (-1 == dir_err) {
            printf("Error creating directory!");
            return EXIT_FAILURE;
        }*/
        BondFileStart = BondDir+"/Bonds"; 
        //BondOutputFileStart = OutputDir+"/BondCorr"; 
        //SpinOutputFileStart = OutputDir+"/SpinCorr"; 
        /*BondOutputFileName = ReplicaDir+"/BondCorr";
        SpinOutputFileName = ReplicaDir+"/SpinCorr";
        std::ofstream BondOutputFile (OutputFileName, std::ios::out | std::ios::trunc);
        if (!BondOutputFile.is_open()) {
            std::cout << "Cannot open output file. Exiting..." << std::endl; 
            return EXIT_FAILURE; 
        }*/
        
        std::map<double,double> Spin_corr;
        std::map<double, double>::iterator Spin_corr_iter; 
        
        T0Step = StartStep;  
        
        
        while (true) {
            BondFileName = BondFileStart + std::to_string(T0Step);   
            std::cout <<"T0: " <<  T0Step << std::endl; 
            std::ifstream BondFile (BondFileName, std::ios::in);   
            if (!BondFile.is_open()) {
                std::cout << "reached last frame " << T0Step-SampleStep << ", continuing with next replica" << std::endl; 
                break;  
            }
            std::map<double,double> Bond_corr;
            std::map<double, double>::iterator Bond_corr_iter; 
            std::map<unsigned, unsigned> Bonds_zero{}; 
            std::map<unsigned, unsigned> Spin_t{};
            /// get input from Bondfile at t=0
            unsigned first {}, second{}; 
            while (BondFile >> first >> second) {
                if (first < second) {
                    Bonds_zero[first] = second; 
                    Normalization++; 
                }    
            }
            Spin_t = Bonds_zero; 
            
            /////////////////////////////////
            
            /// loop over t > 0
            StepVectorIterator = StepVector.begin(); 
            while (StepVectorIterator != StepVector.end()) {
                DStep = *StepVectorIterator; 
                T1Step = T0Step+DStep; 
                //std::cout << T1Step << std::endl; 
                Time = DStep*DeltaT;      
                
                // get input from BondFile at t
                BondFileName = BondFileStart + std::to_string(T1Step);   
                std::ifstream BondFile (BondFileName, std::ios::in);   
                if (!BondFile.is_open()) {
                    std::cout << "reached last frame " << T0Step-SampleStep << ", continuing with next replica" << std::endl; 
                    break;  
                }
                std::map<unsigned, unsigned> Bonds_t{}; 
                while (BondFile >> first >> second) {
                    if (first < second) {
                        Bonds_t[first] = second;  
                    }    
                }
                BondFile.close(); 
               
                
                //// calculating bond correlation
                
                for (auto& bond : Bonds_zero){
                    if (Bonds_t.find(bond.first) -> second == bond.second) {
                        Bond_corr[Time]++; 
                    }
                }
                
                /// calculating spin correlation
                for (auto it = Spin_t.cbegin(); it != Spin_t.cend();) {
                    if (!(Bonds_t.find(it -> first) -> second == it -> second)) { 
                        Spin_t.erase(it++);    // or "it = m.erase(it)" since C++11
                    }
                    else {
                        ++it;
                    }
                }
                
                Spin_corr[Time] += Spin_t.size(); 
                Bond_corr_average[Time] += Bond_corr[Time]; 
                Spin_corr_average[Time] += Spin_t.size(); 
                Bond_corr_average_count[Time]++; 
                
                
                StepVectorIterator++; 
            }
            count++; 
            T0Step += SampleStep; 
        }
    }   
    
    Normalization /= count;
    BondOutputFileName = Directory+"/BondCorr";
    SpinOutputFileName = Directory+"/SpinCorr";
    std::ofstream BondOutputFile (BondOutputFileName, std::ios::out | std::ios::trunc);
    std::ofstream SpinOutputFile (SpinOutputFileName, std::ios::out | std::ios::trunc);
    
    for (auto& corr : Bond_corr_average) {
        BondOutputFile << corr.first << " " << corr.second/(Normalization*Bond_corr_average_count[corr.first]) << std::endl; 
    }
    BondOutputFile.close();
    
    for (auto& corr : Spin_corr_average) {
        SpinOutputFile << corr.first << " " << corr.second/(Normalization*Bond_corr_average_count[corr.first]) << std::endl; 
    }
    SpinOutputFile.close();
    
    return 0;  
}
