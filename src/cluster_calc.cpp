/*
 * cluster_calc.cpp
 *
 *  Created on: Jul 25, 2018
 *      Author: maud
 */

#include <vector>
#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <map>
#include <math.h>
#include "System.h"
#include "HelperFunctions.h"

bool size_compare(std::vector<unsigned>& a, std::vector<unsigned>& b) { return (a.size() > b.size()); }

int main(int argc, char* argv[]) {
    std::string Directory{}, BondFileStart{}, BondFile{}, OutputFileNameStart{}, OutputFileName{}; 
    unsigned StartStep{}, Step{}, SamplingStep{}, StepCount{}, Monomers{}, Chains{}; 

    
    if (argc != 6) {
        std::cout << "usage: ./cluster_calc DIRECTORY STARTSTEP SAMPLINGSTEP MONOMERSPERCHAIN NUMBEROFCHAINS" << std::endl;  
    }  
    
    Directory=argv[1];
    StartStep = std::stoi(argv[2]);
	SamplingStep = std::stoi(argv[3]);  
	Monomers = std::stoi(argv[4]);   
	Chains = std::stoi(argv[5]);
    Step = StartStep; 
    
    std::cout << "Directory: " << Directory << std::endl;
	std::cout << "StartStep: " << StartStep << "   SamplingStep: " << SamplingStep << std::endl;
    
    OutputFileNameStart = Directory+"/Clusters/Cluster";
    // Maybe add some precautions so as to not override data
    BondFileStart = Directory+"/Bonds/Bonds";
    
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
    unsigned MonoFirst{}, MonoSecond{}, ChainFirst{}, ChainSecond{}; 
    
    while (true) {
        std::vector<std::vector<unsigned>> ChainContacts(Chains,std::vector<unsigned>()); 
        std::vector<std::vector<unsigned>> Clusters (Chains, std::vector<unsigned>()); 
        BondFile = BondFileStart+std::to_string(Step);
        std::ifstream file (BondFile, std::ios::in);
        OutputFileName = OutputFileNameStart+std::to_string(Step); 
        
        /// read BondFile and fill IntermolBonds Matrix; 
        if (!file.is_open()) {
            std::cout << "File " << BondFile << "does not exist!"; 
            break;  
        } 
        while(file >> MonoFirst >> MonoSecond) {
            ChainFirst = MonoFirst/Monomers; 
            ChainSecond = MonoSecond/Monomers; 
            if (ChainFirst != ChainSecond) {
                ChainContacts[ChainFirst].push_back(ChainSecond);     
            } 
        } 
        for (unsigned i = 0; i < Chains; i++) Clusters[i].push_back(i); // each chain is one cluster
        
        //loop over all clusters
        for (unsigned i = 0; i < Chains -1; i++) {
            std::vector<unsigned>::iterator iter{}; 
            if (Clusters[i].empty()) continue; 
            iter = Clusters[i].begin(); 
            unsigned n{0}; 
            do {
                if (Clusters[i].empty()) continue;
                ChainFirst = *iter; // first chain in cluster
                for (unsigned j = i+1; j < Chains; j++) { // loop through all clusters j > i 
                    for (unsigned ChainSecond : Clusters[j]) {
                        if(std::find(ChainContacts[ChainFirst].begin(), ChainContacts[ChainFirst].end(), ChainSecond) != ChainContacts[ChainFirst].end()){
                            Clusters[i].insert(Clusters[i].end(), Clusters[j].begin(), Clusters[j].end()); // move all to Cluster[i]; 
                            Clusters[j].clear(); //delete all in j
                            if (n > 0) iter = Clusters[i].begin()+n; 
                            else iter = Clusters[i].begin(); 
                        }
                    }
                }
                n++;
                if (iter != Clusters[i].end()) iter++; 
            } while (iter != Clusters[i].end());
        }
        std::sort(Clusters.begin(), Clusters.end(), size_compare);  /// sort according to descending cluster size      
        std::ofstream output(OutputFileName, std::ios::out | std::ios::trunc);
        if (!output.is_open()) {
            std::cout << "Cannot open output file, directory 'Clusters' might not exist. Exiting..." << std::endl; 
            return EXIT_FAILURE; 
        }
        unsigned cluster_count {0};
        for (auto& Cluster : Clusters) {
            if (Cluster.size() > 1) {
                std::sort(Cluster.begin(), Cluster.end()); 
                output << cluster_count << " " << Cluster.size() ; 
                for (auto& Chain : Cluster) output << " " <<  Chain;
                output << std::endl;  
                cluster_count++;  
            }
        }  
        Step += SamplingStep; 
        StepCount++;         
    }
    gettimeofday(&end, NULL); 
    double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    std::cout << "total time: " << realTime << " , time per step: " << realTime/StepCount << std::endl;
    
    return 0; 
}

