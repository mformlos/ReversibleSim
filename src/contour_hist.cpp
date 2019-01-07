#include <vector>
#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <math.h>
#include "HelperFunctions.h"


int main(int argc, char* argv[]) {
    std::string Directory{}, ReplicaDir{}, OutputDir {}, OutputFileName {}, BondDir{}, BondFileStart{}, BondFileName{};
    unsigned Replicas{}, StartStep{}, Step{}, SampleStep{}, Count{0}, first{}, second{};
    
    std::map<unsigned, double> contour_hist; 
    std::map<unsigned, double>::iterator contour_hist_iter;
    
    
    if (argc != 5) {
        std::cout << "usage: ./cluster_corr DIRECTORY REPLICAS STARTSTEP SAMPLESTEP" << std::endl;  
        return EXIT_FAILURE; 
    }
     
    Directory=argv[1];
	Replicas = std::stoi(argv[2]);  
    StartStep = std::stoi(argv[3]);
    SampleStep = std::stoi(argv[4]); 
    
    std::cout << "Directory: " << Directory << " Replicas: " << Replicas << std::endl;
	std::cout << "StartStep: " << StartStep << " SamplingStep: " << SampleStep  << std::endl;
	
	for (unsigned repl = 0; repl < Replicas; repl++) {
        std::cout << "Replica "  << repl << std::endl; 
        ReplicaDir = Directory+"/REPL-"+std::to_string(repl);
        BondDir = ReplicaDir+"/Bonds"; 
        BondFileStart = BondDir+"/Bonds"; 
        Step = StartStep; 
        while (true) {
            BondFileName = BondFileStart + std::to_string(Step);   
            std::ifstream BondFile (BondFileName, std::ios::in);    
            if (!BondFile.is_open()) {
                std::cout << "reached last frame " << Step << ", continuing with next replica" << std::endl; 
                break;  
            }
            while (BondFile >> first >> second) {
                if (first < second) {
                    contour_hist[second-first] += 1; 
                    Count++; 
                }
            }
            Step += SampleStep; 
        }   
    }
    
    OutputFileName = Directory+"/contour_hist"; 
    std::ofstream OutputFile (OutputFileName, std::ios::out | std::ios::trunc); 
    std::cout << OutputFileName << std::endl; 
    
    for (auto& contour : contour_hist) {
        OutputFile << contour.first << " " << contour.second/Count << std::endl; 
        //std::cout << contour.first << " " << contour.second/Count << std::endl; 
    }
    OutputFile.close();
    return 0;    
}
