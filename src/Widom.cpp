/*
 * Widom.cpp
 *
 *  Created on: May 22, 2018
 *      Author: maud
 */

#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
#include <csignal>
#include <string>
#include <vector>
#include "Molecule.h"
#include "HelperFunctions.h"
#include "Rand.h"
#include "System.h"

using namespace Eigen;

int main(int argc, char* argv[]) {
	unsigned Seed {};
	bool ParameterInitialized {false};
	size_t NInsertions {}, NConfigs {}, NTopols {}, Count {0};
	double Rmax {50}, DeltaR {0.1}, Energy {}, BoltzmannFactor {};
	std::string ParameterFile {}, MoleculeFile {}, HistogramFile {};
	std::vector<std::string> ConfigPoolFiles{}, ConfigFiles1 {}, ConfigFiles2 {};
	std::map<double, double> RadialDistHist {};
	std::map<double, double>::iterator RadialDistHistIter {};

	if (argc != 2) {
		std::cout << "usage: ./Widom.cpp PARAMETER-INPUT-FILE " << std::endl;
		return EXIT_FAILURE;
	}

	//////////////////////////////

	////// PARAMETER READS ///////

	std::ifstream inputfile(argv[1], std::ios::in);
	if (!inputfile.is_open()) {
		std::cout << "could not open file '" << argv[1] << "' , exiting" << std::endl;
		return EXIT_FAILURE;
	}

	NInsertions = extractParameter<size_t>("Insertions", inputfile, ParameterInitialized);
	if (!ParameterInitialized) return EXIT_FAILURE;
	Rmax = extractParameter<double>("Rmax", inputfile, ParameterInitialized);
	if (!ParameterInitialized) return EXIT_FAILURE;
	DeltaR = extractParameter<double>("DeltaR", inputfile, ParameterInitialized);
	if (!ParameterInitialized) return EXIT_FAILURE;
	Seed = extractParameter<double>("Seed", inputfile, ParameterInitialized);
	if (!ParameterInitialized) return EXIT_FAILURE;
	MoleculeFile = extractParameter<std::string>("MoleculeFile", inputfile, ParameterInitialized);
	if (!ParameterInitialized) return EXIT_FAILURE;
	HistogramFile = extractParameter<std::string>("HistogramFile", inputfile, ParameterInitialized);
	if (!ParameterInitialized) return EXIT_FAILURE;

	fillConfigPoolVector(ConfigPoolFiles, inputfile);
	NTopols = ConfigPoolFiles.size();

	inputfile.close();

	/////////////////////////////////////////

	/////// PARAMETER CONTROL OUTPUT ////////
	std::cout << "Rmax is " << Rmax << std::endl;
	std::cout << "DeltaR is " << DeltaR << std::endl;
	std::cout << "Number of Insertions per Config is " << NInsertions << std::endl;
	std::cout << "RNG seed is " << Seed << std::endl;
	std::cout << "MoleculeFile is " << MoleculeFile << std::endl;
	std::cout << "HistogramFile is " << HistogramFile << std::endl;
	std::cout << "Number of distinct topologies is " << ConfigPoolFiles.size() << std::endl;

	std::cout << "TopologyConfigPoolFiles: " << std::endl;
	for (auto& filename : ConfigPoolFiles) std::cout << filename << std::endl;

	////// RANDOM ENGINE SEEDING & WARMUP //////
	Rand::seed(Seed);
	Rand::warmup(10000);

	/////////////
	//Initialization

	System SystemWidom((unsigned)(4*Rmax), (unsigned)(4*Rmax), (unsigned)(4*Rmax), 0.0, 0.0, false);
	SystemWidom.addMolecules(MoleculeFile);

	////////////

	//// initialize histogram of radial distribution function
	for (double Distance = 0.0; Distance < Rmax; Distance += DeltaR) {
		RadialDistHist[Distance] = 0.0;
	}

	timeval start {}, end {};
	gettimeofday(&start, NULL);


	//// loop over different topologies of molecule 1
	for (size_t Topol1 = 0; Topol1 < NTopols; Topol1++) {
		fillConfigPool(ConfigFiles1, ConfigPoolFiles[Topol1]);

		//// loop over different topologies of molecule 2
		for (size_t Topol2 = Topol1; Topol2 < NTopols; Topol2++) {
			fillConfigPool(ConfigFiles2, ConfigPoolFiles[Topol2]);

			//// loop over different configurations of molecule 1
			for (size_t config1 = 0; config1 < ConfigFiles1.size(); config1++)  {
				SystemWidom.Molecules[0].initializePositions(ConfigFiles1[config1]);
				SystemWidom.centerMolecule(0);
				std::cout << Count << std::endl;
				size_t config2start { Topol1 == Topol2 ? config1 : 0 };

				//// loop over different configurations of molecule 2
				for (size_t config2 = config2start; config2 < ConfigFiles2.size(); config2++) {
					SystemWidom.Molecules[1].initializePositions(ConfigFiles2[config2]);
					SystemWidom.centerMolecule(1);

					//// loop over different distances r_12
					for (double Distance = 0.0; Distance < Rmax; Distance += DeltaR) {
						BoltzmannFactor = 0.0;

						//// loop over different insertions at different angles of molecule 2
						for (size_t Insertion = 0; Insertion < NInsertions; Insertion++) {
							SystemWidom.centerMolecule(1);
							SystemWidom.Molecules[1].randomRotation();
							double phi {}, theta {};
							Vector3d Direction {};
							phi = 2.*M_PI*(Rand::real_uniform());
							theta = 2.*(Rand::real_uniform()-0.5);
							Direction(0) = sqrt(1-theta*theta)*cos(phi);
							Direction(1) = sqrt(1-theta*theta)*sin(phi);
							Direction(2) = theta;
							Direction *= Distance;
							SystemWidom.Molecules[1].translate(Direction);
							Energy = SystemWidom.calculateIntermolecularEnergy(0,1);
							BoltzmannFactor += exp(-Energy);
						}
						RadialDistHist.at(Distance) += BoltzmannFactor/NInsertions;
					}
					Count++;
				}
			}
		}
	}
	for (auto& Histvalue : RadialDistHist) {
		Histvalue.second /= Count;
	}

	std::ofstream Output(HistogramFile, std::ios::out | std::ios::trunc);
	size_t samples {0};
	for (double Distance = 0.0; Distance < Rmax; Distance += DeltaR) {
		Output << Distance << " " << RadialDistHist.at(Distance) << std::endl;
		samples++;
	}

	std::cout << "averaged over a total of " << Count << " configurations." << std::endl;

	gettimeofday(&end, NULL);
	double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    std::cout << "total time: " << realTime << " , time per insertion: " << realTime/(Count*NInsertions*samples) << std::endl;


	return 1;
}



