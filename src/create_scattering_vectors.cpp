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
	double deltaK {}, Kmax{}, tolerance{}, Kabs{};

	if (argc != 4) {
		std::cout << "usage: ./create_scattering_vectors DELTAK KMAX TOLERANCE" << std::endl;
		return EXIT_FAILURE;
	}

	deltaK = std::stod(argv[1]);
	Kmax = std::stod(argv[2]);
	tolerance = std::stod(argv[3]);

    std::cout << "DeltaK: " << deltaK << ", Kmax: " << Kmax << ", tolerance: " << tolerance << std::endl; 

    std::vector<three_dim> VecK{};
    
    

	Kabs = 1.0;
	int nx{}, ny{}, nz{}, count{}, filecount {0};
	double norm {};
    int mod_x {}, mod_y {}, mod_z {}; 
    std::string OutputName {}, OutputNameStart {"/home/formanek/REVERSIBLE/input/structure_factor/VecK"};
    std::ofstream Kcount ("/home/formanek/REVERSIBLE/input/structure_factor/Kcount", std::ios::out | std::ios::trunc); 
    	
	while (Kabs <= Kmax) {
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
		if (count > 0) {
		    std::cout << Kabs << " " << count << std::endl;
		    Kcount << Kabs << " " << count << std::endl; 
		    OutputName = OutputNameStart+std::to_string(filecount); 
		    std::ofstream KFile(OutputName, std::ios::out | std::ios::trunc); 
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
		filecount++; 
		
		Kabs += deltaK;
	}

}
