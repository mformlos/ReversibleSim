/*
 * Widom_parallel.cpp
 *
 *  Created on: May 24, 2018
 *      Author: maud
 */

/*
 * Widom.cpp
 *
 *  Created on: May 22, 2018
 *      Author: maud
 */

#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>
#include <iostream>
#include <csignal>
#include <string>
#include <vector>
#include <map>
#include "Molecule.h"
#include "HelperFunctions.h"
#include "Rand.h"

using namespace Eigen;

int main(int argc, char* argv[]) {
	unsigned Seed {}, NumberOfMonomers{};
	bool ParameterInitialized {false};
	size_t NInsertions {}, NConfigs {}, NTopols {}, Count {0}, NIntervals {};
	double Rmax {50}, DeltaR {0.1};
	std::string ParameterFile {}, HistogramFile {};
	std::vector<std::string> ConfigPoolFiles{}, ConfigFiles1 {}, ConfigFiles2 {};
	std::vector<std::pair<size_t, size_t>> TopologyPairs {};
	std::map<double, double> RadialDistHist {};
	std::map<double, double>::iterator RadialDistHistIter {};

	if (argc != 2) {
		std::cout << "usage: ./Widom.cpp PARAMETER-INPUT-FILE " << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "using " << omp_get_num_threads() << " threads " << std::endl;

	//////////////////////////////

	////// PARAMETER READS ///////

	std::ifstream inputfile(argv[1], std::ios::in);
	if (!inputfile.is_open()) {
		std::cout << "could not open file '" << argv[1] << "' , exiting" << std::endl;
		return EXIT_FAILURE;
	}

	NumberOfMonomers = extractParameter<unsigned>("NumberOfMonomers", inputfile, ParameterInitialized);
	if (!ParameterInitialized) return EXIT_FAILURE;
	NInsertions = extractParameter<size_t>("Insertions", inputfile, ParameterInitialized);
	if (!ParameterInitialized) return EXIT_FAILURE;
	Rmax = extractParameter<double>("Rmax", inputfile, ParameterInitialized);
	if (!ParameterInitialized) return EXIT_FAILURE;
	DeltaR = extractParameter<double>("DeltaR", inputfile, ParameterInitialized);
	if (!ParameterInitialized) return EXIT_FAILURE;
	Seed = extractParameter<double>("Seed", inputfile, ParameterInitialized);
	if (!ParameterInitialized) return EXIT_FAILURE;
	HistogramFile = extractParameter<std::string>("HistogramFile", inputfile, ParameterInitialized);
	if (!ParameterInitialized) return EXIT_FAILURE;

	fillConfigPoolVector(ConfigPoolFiles, inputfile);
	NTopols = ConfigPoolFiles.size();
	NIntervals = (size_t)(Rmax/DeltaR);
	inputfile.close();

	/////////////////////////////////////////

	/////// PARAMETER CONTROL OUTPUT ////////
	std::cout << "NumberOfMonomers is " << NumberOfMonomers << std::endl;
	std::cout << "Rmax is " << Rmax << std::endl;
	std::cout << "DeltaR is " << DeltaR << std::endl;
	std::cout << "Number of Insertions per Config is " << NInsertions << std::endl;
	std::cout << "RNG seed is " << Seed << std::endl;
	std::cout << "HistogramFile is " << HistogramFile << std::endl;
	std::cout << "Number of distinct topologies is " << ConfigPoolFiles.size() << std::endl;

	std::cout << "TopologyConfigPoolFiles: " << std::endl;
	for (auto& filename : ConfigPoolFiles) std::cout << filename << std::endl;

	////// RANDOM ENGINE SEEDING & WARMUP //////
	Rand::seed(Seed);
	Rand::warmup(10000);

	/////////////
	//Initialization

	////////////

	//// initialize histogram of radial distribution function
	for (double Distance = 0.0; Distance < Rmax; Distance += DeltaR) {
		RadialDistHist[Distance] = 0.0;
	}

	//// initialize TopologyPairs to loop over
	for (size_t Topol1 = 0; Topol1 < NTopols; Topol1++) {
		for (size_t Topol2 = Topol1; Topol2 < NTopols; Topol2++) {
			TopologyPairs.push_back(std::pair<size_t, size_t>(Topol1, Topol2));
		}
	}

	//// test TopologyPairs
	for (auto& pair : TopologyPairs) {
		std::cout << pair.first << " " << pair.second << std::endl;
	}


	timeval start {}, end {};
	gettimeofday(&start, NULL);

	#pragma omp parallel private(ConfigFiles1, ConfigFiles2)
	{
		if (omp_get_thread_num()==0) {
			std::cout << "using " << omp_get_num_threads() << " threads " << std::endl;
			std::cout << "on " << omp_get_num_procs() << " processors " << std::endl;
		}

		//// loop over different topologies of molecule 1 and 2
		#pragma omp for
		for (size_t TopolPair = 0; TopolPair < TopologyPairs.size(); TopolPair++) {
			size_t Topol1 { TopologyPairs[TopolPair].first };
			size_t Topol2 { TopologyPairs[TopolPair].second };

			fillConfigPool(ConfigFiles1, ConfigPoolFiles[Topol1]);
			fillConfigPool(ConfigFiles2, ConfigPoolFiles[Topol2]);
			std::map<double, double> RadialDistHist_local {};
			size_t Count_local {0};
			for (double Distance = 0.0; Distance < Rmax; Distance += DeltaR) {
				RadialDistHist_local[Distance] = 0.0;
			}
			//// loop over different configurations of molecule 1
			for (size_t config1 = 0; config1 < ConfigFiles1.size(); config1++)  {
				Molecule Molecule1(NumberOfMonomers, 0);
				Molecule1.initializePositions(ConfigFiles1[config1]);
				Vector3d COMPos1 {-Molecule1.centerOfMassPosition()};
				Molecule1.translate(COMPos1);
				size_t config2start { Topol1 == Topol2 ? config1 : 0 };

				//// loop over different configurations of molecule 2
				for (size_t config2 = config2start; config2 < ConfigFiles2.size(); config2++) {
					Molecule Molecule2(NumberOfMonomers, (int)NumberOfMonomers);
					Molecule2.initializePositions(ConfigFiles1[config1]);

					//// loop over different distances r_12
					for (size_t Interval = 0; Interval < NIntervals; Interval++) {
					//for (double Distance = 0.0; Distance < Rmax; Distance += DeltaR) {
						double Distance {Interval*DeltaR};
						double BoltzmannFactor {0.0};

						//// loop over different insertions at different angles of molecule 2
						for (size_t Insertion = 0; Insertion < NInsertions; Insertion++) {
							Vector3d COMPos2 {-Molecule2.centerOfMassPosition()};
							Molecule2.translate(COMPos2);
							Molecule2.randomRotation();
							double phi {}, theta {};
							Vector3d Direction {};
							phi = 2.*M_PI*(Rand::real_uniform());
							theta = 2.*(Rand::real_uniform()-0.5);
							Direction(0) = sqrt(1-theta*theta)*cos(phi);
							Direction(1) = sqrt(1-theta*theta)*sin(phi);
							Direction(2) = theta;
							Direction *= Distance;
							Molecule2.translate(Direction);
							double Energy {calculateIntermolecularEnergy(Molecule1, Molecule2)};
							BoltzmannFactor += exp(-Energy);
						}
						RadialDistHist_local.at(Distance) += BoltzmannFactor/NInsertions;
					}
					Count_local++;
				}
			}
			#pragma omp atomic
			Count += Count_local;
			for (auto& Histvalue : RadialDistHist_local) {
				#pragma omp atomic
				RadialDistHist.at(Histvalue.first) += Histvalue.second;
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







