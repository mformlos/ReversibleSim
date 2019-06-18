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
    std::string Directory{}, ReplicaDir{},  OutputDir {}, BondDir {}, BondFileStart{}, BondFileName{}, BondOutputFileStart{}, SpinOutputFileStart{}, BondOutputFileName{}, SpinOutputFileName{}, StepFile{}, IntervalFile{}; 
    unsigned StartStep{}, DStep {}, SampleStep{}, Replicas{}, FunctionalGroups{}, NIntervals {}; 
    double DeltaT{}, Time{}; 
    
    std::vector<unsigned> IntervalEnds; 
    std::vector<unsigned>::iterator IntervalEnds_iter;
    
    std::vector<unsigned long long> StepVector{}; 
    std::vector<unsigned long long>::iterator StepVectorIterator{}; 
    
    if (argc != 9) {
        std::cout << "usage: ./cluster_corr DIRECTORY REPLICAS FUNCTIONALGROUPS STARTSTEP SAMPLESTEP DELTAT STEPFILE INTERVALFILE" << std::endl;  
        return EXIT_FAILURE; 
    }
    
    Directory=argv[1];
	Replicas = std::stoi(argv[2]); 
	FunctionalGroups = std::stoi(argv[3]);
	StartStep = std::stoi(argv[4]);
    SampleStep = std::stoi(argv[5]); 
	DeltaT = std::stod(argv[6]); 
	StepFile = argv[7];
	IntervalFile = argv[8]; 
	
	std::cout << "Directory: " << Directory << " Replicas: " << Replicas << " Functional groups: " << FunctionalGroups << std::endl;
	std::cout << "StartStep: " << StartStep << " SamplingStep: " << SampleStep << " StepFile: " << StepFile << std::endl;
	std::cout << "DeltaT: " << DeltaT  << "IntervalFile: " << IntervalFile << std::endl;
	
	initializeStepVector(StepVector, StepFile); 
	initializeVector(IntervalEnds, IntervalFile); 
	
	NIntervals = IntervalEnds.size(); 
	std::cout << "Intervals: " << std::endl; 
	for (auto& interval : IntervalEnds) {
	    std::cout << interval << " "; 
	}
	
	std::vector<std::map<double,double>> Bond_corr_average (NIntervals, std::map<double,double>()); 
	std::map<double,unsigned> Bond_corr_average_count; 
    //std::vector<std::map<double,unsigned>> Bond_corr_average_count(NIntervals, std::map<double,unsigned>());
    
    std::vector<std::map<double,double>> Spin_corr_average(NIntervals, std::map<double,double>());
    std::map<double,unsigned> Spin_corr_average_count;
	

    int T0Step {}; 
    int T1Step {};
    double Normalization {0.0}; 	
    std::vector<double> Normalization_intervals(NIntervals, 0.0); 
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
        
        
        unsigned contour{}; 
        unsigned intervalindex {}; 
        
        T0Step = StartStep;  
        
        
        while (true) {
            BondFileName = BondFileStart + std::to_string(T0Step);   
            std::cout <<"T0: " <<  T0Step << std::endl; 
            std::ifstream BondFile (BondFileName, std::ios::in);   
            if (!BondFile.is_open()) {
                std::cout << "reached last frame " << T0Step-SampleStep << ", continuing with next replica" << std::endl; 
                break;  
            }
            
            std::vector<std::map<double,double>> Bond_corr (NIntervals, std::map<double,double>());
            std::vector<std::map<double,double>> Spin_corr (NIntervals, std::map<double,double>());
            std::map<unsigned, unsigned> Bonds_zero{}; 
            std::map<unsigned, unsigned> Spin_t{};
            /// get input from Bondfile at t=0
            unsigned first {}, second{}; 
            while (BondFile >> first >> second) {
                if (first < second) {
                    Bonds_zero[first] = second; 
                    contour = second-first; 
                    intervalindex = 0; 
                    while ( IntervalEnds[intervalindex] < contour) {
                        intervalindex++;     
                    }
                    Normalization_intervals[intervalindex]++; 
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
                        contour = bond.second-bond.first;
                        intervalindex = 0; 
                        while ( IntervalEnds[intervalindex] < contour) {
                            intervalindex++;     
                        }
                        Bond_corr[intervalindex][Time]++; 
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
                for (auto& bond : Spin_t) {
                    contour =  bond.second-bond.first;
                    intervalindex = 0; 
                    while ( IntervalEnds[intervalindex] < contour) {
                        intervalindex++;     
                    }
                    Spin_corr[intervalindex][Time]++; 
                }
                //Spin_corr[Time] += Spin_t.size(); 
                for (unsigned i = 0; i < NIntervals; i++) {
                    Bond_corr_average[i][Time] += Bond_corr[i][Time]; 
                    Spin_corr_average[i][Time] += Spin_corr[i][Time]; 
                }
                Bond_corr_average_count[Time]++; 
                
                
                StepVectorIterator++; 
            }
            count++; 
            T0Step += SampleStep; 
        }
    }   
    
    Normalization /= count;
    
    
    
    std::ofstream BondOutputFile{}; 
    unsigned IntervalStart {2}, IntervalEnd{}; 
    for (unsigned i = 0; i < NIntervals; i++) {
        IntervalEnd = IntervalEnds[i]; 
        BondOutputFileName = Directory+"/BondCorr-"+std::to_string(IntervalStart)+"-"+std::to_string(IntervalEnd);
        BondOutputFile.open(BondOutputFileName, std::ios::out | std::ios::trunc);
        for (auto& corr : Bond_corr_average[i]) {
            BondOutputFile << corr.first << " " << corr.second*count/(Normalization_intervals[i]*Bond_corr_average_count[corr.first]) << std::endl; 
        }
        BondOutputFile.close();
        IntervalStart = IntervalEnd+1;     
    }
    
    std::ofstream SpinOutputFile{}; 
    IntervalStart = 2; 
    for (unsigned i = 0; i < NIntervals; i++) {
        IntervalEnd = IntervalEnds[i]; 
        SpinOutputFileName = Directory+"/SpinCorr-"+std::to_string(IntervalStart)+"-"+std::to_string(IntervalEnd);
        SpinOutputFile.open(SpinOutputFileName, std::ios::out | std::ios::trunc);
    
        for (auto& corr : Spin_corr_average[i]) {
            SpinOutputFile << corr.first << " " << corr.second*count/(Normalization_intervals[i]*Bond_corr_average_count[corr.first]) << std::endl; 
        }
        SpinOutputFile.close();
        IntervalStart = IntervalEnd+1;
    }
    
    return 0;  
}
