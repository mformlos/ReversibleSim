/*
 * radial_distribution.cpp
 *
 *  Created on: Jul 4, 2018
 *      Author: maud
 */

#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <math.h>
#include "Particle.h"
#include "HelperFunctions.h"


int main(int argc, char* argv[]) {
	std::string Directory{}, ConfigFile{}, ConfigFileStart{};
	double FractionReactive{}, StartR{0.0}, EndR{5.0}, DR{0.1};
	int StartStep{}, EndStep{}, Step{}, SamplingStep{}, N{};
	unsigned NumberOfMonomers {}, NumberReactive{};
	std::vector<Particle> ReactiveMonomers {};

	std::map<double, double> radial_dist {};
	std::map<double, int> radial_dist_count{};


	unsigned total_count {0};

	if (argc != 9) {
	            std::cout << "usage: ./radial_dist DIRECTORY MONOMERS REACTIVEFRACTION STARTSTEP ENDSTEP SAMPLINGSTEP INTERVAL ENDR" << std::endl;
	            return EXIT_FAILURE;
	}

	Directory=argv[1];
	NumberOfMonomers=std::stoi(argv[2]);
	FractionReactive = std::stod(argv[3]);
	StartStep = std::stoi(argv[4]);
	EndStep = std::stoi(argv[5]);
	SamplingStep = std::stoi(argv[6]);
	DR = std::stod(argv[7]);
	EndR  = std::stod(argv[8]); 

	std::cout << "Directory: " << Directory << std::endl;
	std::cout << "StartStep: " << StartStep << "   EndStep: " << EndStep << "   SamplingStep: " << SamplingStep << std::endl;
	std::cout << "StartR: " << StartR << " EndR: " << EndR << " DR: " << DR << std::endl;    

	N = ((EndStep-StartStep)/SamplingStep);

	NumberReactive = NumberOfMonomers*FractionReactive;

	for (unsigned i = 0; i < NumberReactive; i++)
	{
		ReactiveMonomers.push_back(Particle(i));
	}
	

	ConfigFileStart = Directory+"/configs/config";
	std::string name;
	name = Directory+"/radial_dist";
	std::ifstream test(name);
	if (test.good()){
		std::cout << "file " << name << " already exists! Aborting..." << std::endl;
		return EXIT_FAILURE;
	}
	std::ofstream output(name, std::ios::out | std::ios::trunc);
	output.precision(8);

	Step = StartStep;
	double distance {};
	double bin {};
	while (true) {
		ConfigFile = ConfigFileStart+std::to_string(Step)+".pdb";
		std::cout << ConfigFile << std::endl; 
		if (!initializeReactivePositions(ReactiveMonomers, ConfigFile)) {
		    std::cout << "problem with initializing monomers" << std::endl;
		    break; 
		}
		for (unsigned i = 0; i < NumberReactive; i++) {
			for (unsigned j = i + 1; j < NumberReactive; j++) {
				distance = (ReactiveMonomers[i].Position-ReactiveMonomers[j].Position).norm();
				if (distance <= EndR) {
				    bin = unsigned(distance/DR)*DR;
				    //std::cout << "d: " << distance << " bin: " << bin << std::endl;  
				    radial_dist_count[bin]++;
				    
				}
			}
		}
		total_count++;
		Step += SamplingStep;
	}
	double normalization {1./(2.*M_PI*DR*(NumberReactive-1)*total_count)};
	for (auto& r : radial_dist_count) {
		output << r.first+0.5*DR << " " << (double)r.second*normalization/pow((r.first+0.5*DR),2) << std::endl;
	}
	output.close();
}


