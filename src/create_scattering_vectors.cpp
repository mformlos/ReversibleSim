/*
 * create_scattering_vectors.cpp
 *
 *  Created on: Jul 11, 2018
 *      Author: maud
 */

#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>


int main(int argc, char* argv[]) {
	double deltaK {}, Kmax{}, tolerance{}, Kabs{};

	if (argc != 4) {
		std::cout << "usage: ./create_scattering_vectors DELTAK KMAX TOLERANCE" << std::endl;
		return EXIT_FAILURE;
	}

	deltaK = std::stod(argv[1]);
	Kmax = std::stod(argv[2]);
	tolerance = std::stod(argv[3]);

	Kabs = 1.0;
	int nx{}, ny{}, nz{}, count{};
	double norm {};
	while (Kabs <= Kmax) {
		count = 0;
		for (nx = 0; nx <= Kabs; nx++) {
			for (ny = 0; ny <= Kabs; ny++) {
				for (nz = 0; nz <= Kabs; nz++) {
					norm = nx*nx+ny*ny+nz*nz;
					norm = sqrt(norm);
					if (abs(norm - Kabs) < tolerance) {
						std::cout << nx << " " << ny << " " << nz << std::endl;
						count++;
					}
				}
			}
		}
		std::cout << Kabs << " " << count << std::endl;
		Kabs += deltaK
	}

}
