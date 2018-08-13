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
    std::string Directory{}, ReplicaDir{},  OutputDir {}, ClusterDir {}, ClusterFileStart{}, ClusterFileName{}, OutputFileStart{}, OutputFileName{}, ClusterLine {}; 
    unsigned StartStep{}, MaxStep{}, Step{}, DStep {}, SampleStep{}, Replicas{}, MinClusterSize{}, Molecules{}; 
    double DeltaT{}, Time{}; 
    std::map<double,double> Cluster_corr_average;
    std::map<double, double>::iterator Cluster_corr_average_iter; 
    std::map<double,unsigned> Cluster_corr_average_count;
    std::map<double, unsigned>::iterator Cluster_corr_average_count_iter; 
    
    if (argc != 10) {
            std::cout << "usage: ./cluster_corr DIRECTORY REPLICAS STARTSTEP MAXSTEP STEP SAMPLESTEP DELTAT MINCLUSTERSIZE MOLECULES" << std::endl;  
            return EXIT_FAILURE; 
    }
    
    Directory=argv[1];
	Replicas = std::stoi(argv[2]); 
	StartStep = std::stoi(argv[3]);
	MaxStep = std::stoi(argv[4]); 
    Step = std::stoi(argv[5]); 
    SampleStep = std::stoi(argv[6]); 
	DeltaT = std::stod(argv[7]); 
	MinClusterSize = std::stoi(argv[8]);
	Molecules = std::stoi(argv[9]); 

    std::cout << "Directory: " << Directory << " Replicas: " << Replicas << " Molecules: " << Molecules << std::endl;
	std::cout << "StartStep: " << StartStep << " MaxStep: " << MaxStep << " Step: " << Step << "   SamplingStep: " << SampleStep << std::endl;
	std::cout << "DeltaT: " << DeltaT << " MinClusterSize: " << MinClusterSize << std::endl;    
    
    int T0Step {StartStep}; 
    int T1Step {};
    double Normalization {0.0}; 
    
    
    
    for (unsigned repl = 0; repl < Replicas; repl++) {
        std::cout << "Replica "  << repl << std::endl; 
        ReplicaDir = Directory+"/REPL-"+std::to_string(repl);
        ClusterDir = ReplicaDir+"/Clusters";
        OutputDir = ReplicaDir+"/ClusterCorr";
        std::cout << OutputDir; 
        const int dir_err {mkdir(OutputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)};
        if (-1 == dir_err) {
            printf("Error creating directory!");
            return EXIT_FAILURE;
        }
        ClusterFileStart = ClusterDir+"/Cluster"; 
        OutputFileStart = OutputDir+"/clustercorr"; 
        T0Step = StartStep; 
        while(true) {
            ClusterFileName = ClusterFileStart + std::to_string(T0Step); 
            std::cout <<"T0: " <<  T0Step << std::endl; 
            std::ifstream ClusterFile (ClusterFileName, std::ios::in); 
            if (!ClusterFile.is_open()) {
                std::cout << "reached last frame " << T0Step-SampleStep << ", continuing with next replica" << std::endl; 
                break;  
            }
            OutputFileName = OutputFileStart + std::to_string(T0Step); 
            std::ofstream OutputFile (OutputFileName, std::ios::out | std::ios::trunc);
            if (!OutputFile.is_open()) {
                std::cout << "Cannot open output file, directory 'ClusterCorr' might not exist. Exiting..." << std::endl; 
                return EXIT_FAILURE; 
            }
            std::map<double,double> Cluster_corr;
            std::map<double, double>::iterator Cluster_corr_iter; 
            std::vector<unsigned> S_zero(Molecules, 0);
            
            /// get input from Clusterfile at t=0
            getline(ClusterFile, ClusterLine); 
            std::istringstream buf1(ClusterLine); 
            std::istream_iterator<std::string> beg1(buf1), end1;
            std::vector<std::string> tokens1(beg1, end1);
            ClusterFile.close(); 
            /////////////////////////////////
            
            //// fill S_zero vector 
            unsigned mol {};
            for (auto mol_iter = tokens1.begin()+2; mol_iter != tokens1.end(); mol_iter++) {
                mol = std::stoi(*mol_iter); 
                S_zero[mol] = 1;     
            }
            ///////////////////////
            
            /// Normalization ////
            Normalization = 0.0; 
            for (auto& S : S_zero) {
                Normalization += (double)(S*S); 
            }
            Normalization /= Molecules; 
            std::cout << "Normalization: " << Normalization << std::endl; 
            //Cluster_corr[0.0] = Normalization; 
            
            std::vector<unsigned> S_t(S_zero.begin(), S_zero.end()); 
            T1Step = T0Step;
            DStep = 0; 
            
            /// loop over t > 0
            while (DStep <= MaxStep) {
                Time = DStep*DeltaT; 
                Cluster_corr[Time] = 0.0; 
                //// open Clusterfile for time t 
                ClusterFileName = ClusterFileStart + std::to_string(T1Step); 
                std::ifstream ClusterFile (ClusterFileName, std::ios::in); 
                if (!ClusterFile.is_open()) {
                    std::cout << "reached last frame " << T1Step-Step << ", continuing with next starting point" << std::endl; 
                    break;  
                }
                /// getline and split
               
                getline(ClusterFile, ClusterLine); 
                std::istringstream buf2(ClusterLine); 
                std::istream_iterator<std::string> beg2(buf2), end2;
                std::vector<std::string> tokens2(beg2, end2);
                ClusterFile.close(); 
               
                
                
                /// fill vector S_t
                std::vector<std::string>::iterator mol_cluster_iter = tokens2.begin()+2; 
                
                unsigned mol_cluster {}; 
                mol = 0; 
                 
                while (mol_cluster_iter != tokens2.end()) {
                    mol_cluster = std::stoi(*mol_cluster_iter); 
                   
                    while (mol_cluster != mol) {
                        S_t[mol] = 0; 
                        mol++; 
                    }
                    mol++; 
                    mol_cluster_iter++; 
                }
              
                //// calculating correlation
                for (unsigned i = 0; i < Molecules; i++) {
                    Cluster_corr[Time] += S_zero[i]*S_t[i]; 
                    
                }
                Cluster_corr[Time] /= (Molecules*Normalization);
                Cluster_corr_average[Time] += Cluster_corr[Time];  
                Cluster_corr_average_count[Time]++;      
                DStep += Step; 
                T1Step += Step; 
            }  
            //// write correlator
            for (auto& corr : Cluster_corr) {
                OutputFile << corr.first << " " << corr.second << std::endl; 
            }          
            OutputFile.close(); 
            T0Step += SampleStep;
        }
    }
    
    OutputFileName = Directory+"/cluster_corr"; 
    std::ofstream OutputFile (OutputFileName, std::ios::out | std::ios::trunc);
    
    for (auto& corr : Cluster_corr_average) {
        OutputFile << corr.first << " " << corr.second/Cluster_corr_average_count[corr.first] << std::endl; 
    }
    OutputFile.close(); 

    return 0;
}
