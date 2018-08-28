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
    unsigned long StartStep{}, CurrentStep{}, Step{}, Replicas{}, MinClusterSize{}, Molecules{}; 
    
    double DeltaT{}, Time_outside{}; 
    /*std::map<double,double> Time_outside_histo;
    std::map<double, double>::iterator Time_outside_histo_iter;
    unsigned count{0}; */
    
    
    if (argc != 8) {
            std::cout << "usage: ./cluster_corr DIRECTORY REPLICAS STARTSTEP STEP DELTAT MINCLUSTERSIZE MOLECULES" << std::endl;  
            return EXIT_FAILURE; 
    }
    
    Directory=argv[1];
	Replicas = std::stoi(argv[2]); 
	StartStep = std::stoi(argv[3]);
    Step = std::stoi(argv[4]); 
	DeltaT = std::stod(argv[5]); 
	MinClusterSize = std::stoi(argv[6]);
	Molecules = std::stoi(argv[7]); 
	
	std::cout << "Directory: " << Directory << " Replicas: " << Replicas << " Molecules: " << Molecules << std::endl;
	std::cout << "StartStep: " << StartStep << " Step: " << Step << std::endl;
	std::cout << "DeltaT: " << DeltaT << " MinClusterSize: " << MinClusterSize << std::endl;    
	
	std::vector<int> S_current(Molecules, 0); // vector saying whether a chain belongs to a cluster 
	std::vector<int> S_last(Molecules, 0); 
	
	std::vector<int long> Step_left(Molecules, 0); // vector saying at which step a chain left the cluster
	std::vector<int long> Step_entered(Molecules, 0);
	
	OutputFileName = Directory+"/TimeOutside"; 
	std::ofstream OutputFile (OutputFileName, std::ios::out | std::ios::trunc);
	
	for (unsigned repl = 0; repl < Replicas; repl++) {
	    std::cout << "Replica "  << repl << std::endl; 
        ReplicaDir = Directory+"/REPL-"+std::to_string(repl);
        ClusterDir = ReplicaDir+"/Clusters";
        ClusterFileStart = ClusterDir+"/Cluster"; 
        /// reset everything: 
        std::fill(S_current.begin(), S_current.end(), 0); 
        std::fill(S_last.begin(), S_last.end(), 0);
        std::fill(Step_left.begin(), Step_left.end(), 0);
        std::fill(Step_entered.begin(), Step_entered.end(), 0);
        
        CurrentStep = StartStep; 
        ////
        ClusterFileName = ClusterFileStart + std::to_string(CurrentStep); 
        std::ifstream ClusterFile (ClusterFileName, std::ios::in); 
        if (!ClusterFile.is_open()) {
            std::cout << "reached last frame " <<  CurrentStep << ", continuing with next replica" << std::endl; 
            break;  
        }
        // get input from Clusterfile
        getline(ClusterFile, ClusterLine); 
        std::istringstream buf1(ClusterLine); 
        std::istream_iterator<std::string> beg1(buf1), end1;
        std::vector<std::string> tokens1(beg1, end1);
        ClusterFile.close(); 
        /////////////////////////////////
        //// fill S_last vector 
        unsigned mol {};
        for (auto mol_iter = tokens1.begin()+2; mol_iter != tokens1.end(); mol_iter++) {
            mol = std::stoi(*mol_iter); 
            S_last[mol] = 1;     
        }
        ////////

    
        while (true) {
            ClusterFileName = ClusterFileStart + std::to_string(CurrentStep); 
            std::ifstream ClusterFile (ClusterFileName, std::ios::in); 
            if (!ClusterFile.is_open()) {
                std::cout << "reached last frame " << CurrentStep << ", continuing with next replica" << std::endl; 
                break;  
            }
            // get input from Clusterfile
            getline(ClusterFile, ClusterLine); 
            std::istringstream buf1(ClusterLine); 
            std::istream_iterator<std::string> beg1(buf1), end1;
            std::vector<std::string> tokens1(beg1, end1);
            ClusterFile.close(); 
            /////////////////////////////////
           
            ///// reset S_current vector
            std::fill(S_current.begin(), S_current.end(), 0); 
            
            //// fill S_current vector 
            unsigned mol {};
            for (auto mol_iter = tokens1.begin()+2; mol_iter != tokens1.end(); mol_iter++) {
                mol = std::stoi(*mol_iter); 
                S_current[mol] = 1;     
            }
            ////////
            
            //// compare S_current & S_last 
            for (unsigned i = 0; i < Molecules; i++) {
                if (S_current[i] - S_last[i] == -1){ // chain left cluster 
                    Step_left[i] = CurrentStep;     
                }
                else if (S_current[i] - S_last[i] == 1) { // chain enters cluster 
                    Step_entered[i] = CurrentStep;
                    if (Step_left[i] > 0) {
                        Time_outside = (Step_entered[i]-Step_left[i])*DeltaT; 
                        if (Time_outside < 0.0 ) {
                            std::cout << "left at Step: " << Step_left[i] << " , entered at Step: " << Step_entered[i] << std::endl; 
                        }
                        Step_left[i] = 0; 
                        Step_entered[i] = 0; 
                        OutputFile << Time_outside << std::endl; 
                    }
                    //Time_outside_histo[Time_outside]++; 
                    //count++;     
                }
            }
            //// copy S_current to S_last 
            S_last = S_current;
            CurrentStep += Step; 
        }
	}
	/*OutputFileName = Directory+"/TimeOutsideHisto"; 
    std::ofstream OutputFile (OutputFileName, std::ios::out | std::ios::trunc);
	for (auto& time : Time_outside_histo) {
	    OutputFile << time.first << " " << time.second/count << std::endl; 
	}*/
	OutputFile.close(); 
	return 0; 
}
