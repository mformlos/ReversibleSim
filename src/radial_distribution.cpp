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
	            std::cout << "usage: ./msd DIRECTORY MONOMERS REACTIVEFRACTION STARTSTEP ENDSTEP SAMPLINGSTEP MaxRadius Interval " << std::endl;
	            return EXIT_FAILURE;
	}

	Directory=argv[1];
	NumberOfMonomers=std::stoi(argv[2]);
	FractionReactive = std::stod(argv[3]);
	StartStep = std::stoi(argv[4]);
	EndStep = std::stoi(argv[5]);
	SamplingStep = std::stoi(argv[6]);
	EndR = std::stod(argv[7]);
	DR = std::stod(argv[8]);

	std::cout << "Directory: " << Directory << std::endl;
	std::cout << "StartStep: " << StartStep << "   EndStep: " << EndStep << "   SamplingStep: " << SamplingStep << std::endl;

	N = ((EndStep-StartStep)/SamplingStep);

	NumberReactive = NumberOfMonomers*FractionReactive;

	for (unsigned i = 0; i < NumberReactive; i++)
	{
		ReactiveMonomers.push_back(Particle(i));
	}
	for (double r = StartR; r <= EndR; r += DR) {
		radial_dist_count[r] = 0;
	}


	ConfigFileStart = Directory+"/configs/config";
	std::string name;
	name = Directory+"/data/radial_dist";
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
	while (Step <= EndStep) {
		ConfigFile = ConfigFileStart+std::to_string(Step)+".pdb";
		initializeReactivePositions(ReactiveMonomers, ConfigFile);
		for (unsigned i = 0; i < NumberReactive; i++) {
			for (unsigned j = i + 1; j < NumberReactive; j++) {
				distance = (ReactiveMonomers[i].Position-ReactiveMonomers[j].Position).norm();
				bin = unsigned(distance/DR)*DR;
				radial_dist_count.at(bin)++;
				total_count++;
			}
		}
		Step += SamplingStep;
	}
	for (auto& r : radial_dist_count) {
		output << r.first << " " << r.second/total_count << std::endl;
	}
	output.close();
}


