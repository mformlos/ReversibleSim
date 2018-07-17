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
#include "HelperFunctions.h"

struct three_dim {
    int x; 
    int y; 
    int z; 
    three_dim(int i, int j, int k) :
        x{i}, 
        y{j},
        z{k} {}
};

int main(int argc, char* argv[]) {
	double deltaK {}, Kmax{}, tolerance{}, tolerance_change {};
	std::string KabsFile{}, OutputDir{}; 
	std::vector<double> KabsVec {}; 
	unsigned MaxVecs {}; 

	if (argc != 6) {
		std::cout << "usage: ./create_scattering_vectors KABSFILE TOLERANCE TOLERANCECHANGE MAXVECS OUTPUTDIR" << std::endl;
		return EXIT_FAILURE;
	}

    KabsFile = argv[1];
	tolerance = std::stod(argv[2]);
	tolerance_change = std::stod(argv[3]); 
	MaxVecs = std::stoi(argv[4]);
	OutputDir = argv[5]; 
	

    std::cout << "KabsFile: " << KabsFile << ", OutputDir: " << OutputDir << ", tolerance: " << tolerance << " , Maxvecs: " << MaxVecs << std::endl; 

    if (!initializeDoubleVector(KabsVec, KabsFile)){
        std::cout << "KabsFile doesn't exist, exiting!" << std::endl; 
        return EXIT_FAILURE; 
    }
    std::vector<three_dim> VecK{};
    
    

	int nx{}, ny{}, nz{}, count{}, filecount {0};
	double norm {};
    int mod_x {}, mod_y {}, mod_z {}; 
    std::string OutputName {}, OutputNameStart {OutputDir+"VecK"};
    std::ofstream Kcount (OutputDir+"Kcount", std::ios::out | std::ios::trunc); 
    std::ofstream KVecFileList (OutputDir+"KVecFileList", std::ios::out | std::ios::trunc); 
    	
	for (auto& Kabs : KabsVec) {
		count = 0;
		VecK.clear(); 
		mod_x = 0; 
		//std::cout << Kabs << std::endl; 
		
		for (nx = 0; nx <= Kabs; nx++) {
		    mod_y = 0;
            if (nx == 1) mod_x = -1; 
			for (ny = 0; ny <= Kabs; ny++) {
			    mod_z = 0; 
                if (ny == 1) mod_y = -1; 
				for (nz = 0; nz <= Kabs; nz++) {
                    if (nz == 1) mod_z = -1; 
					norm = nx*nx+ny*ny+nz*nz;
					norm = sqrt(norm);
					if (fabs(norm - Kabs) < tolerance) {   	
					    for (int i = 1; i >= mod_x; i -= 2) {
					        for (int j = 1; j >= mod_y; j -= 2) {
					            for (int k = 1; k >= mod_z; k -= 2) {
					                VecK.push_back(three_dim(i*nx, j*ny, k*nz));
					                //std::cout << i*nx << " " << j*ny << " " << k*nz << " , " << norm << std::endl;
					                count++;
					            }
					        }
					    } 
					}
				}
			}
		}
		std::cout << Kabs << " " << count << std::endl;
		
		if (count > 0) { 
		    Kcount << Kabs << " " << count << std::endl; 
		    OutputName = OutputNameStart+std::to_string(filecount); 
		    std::ofstream KFile(OutputName, std::ios::out | std::ios::trunc); 
		    KVecFileList << OutputName << std::endl; 
		    KFile << Kabs << std::endl;
		    //KFile.precision(6); 
		    for (auto& K : VecK) {
		        KFile.width(10); 
		        KFile << K.x << " "; 
		        KFile.width(10); 
		        KFile << K.y << " ";
		        KFile.width(10); 
		        KFile << K.z << std::endl;
		    }
		}
		if (count >= MaxVecs) {
		    tolerance *= tolerance_change; 
		    std::cout << "new tolerance: " << tolerance << std::endl; 
		}
		filecount++; 
	}

}
