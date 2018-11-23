
#include <sys/time.h>
#include "System.h"
#include "HelperFunctions.h"


int main() {
    System sys_test(200,200,200, true);
    std::vector<unsigned long long> OutputSteps {};
    std::vector<unsigned long long>::iterator OutputStepsIt{};
    unsigned long long TotalSteps {100000000};
    
    bool test {};
    test = sys_test.addMolecules("input/SCNPs/SCNP-PRECURSOR-2chains", 1.0);
    if (!test) {
    	std::cout << "problem with adding molecules" << std::endl;
    }
    /*for (auto& mol : sys_test.Molecules) {
    	std::cout << "Molecule: " << count<< std::endl;
    	for (auto& part : mol.Monomers) {
    		std::cout << part.Identifier << std::endl;
    	}
    	count++;
    }*/
    unsigned Index {}, monoIndex{13}, molIndex{2};
    sys_test.makeIndex(molIndex, monoIndex, Index);

    std::cout << Index << " " << molIndex << " " << monoIndex;

    //sys_test.addFunctional("input/SCNPs/SCNP-PRECURSOR-functional");

    for (auto& bonds : sys_test.ReversibleBonds) {
    	std::cout << bonds.first << " " << bonds.second << std::endl;
    }

    //sys_test.addLinks("/home/formanek/HYBRIDSIM/input/SCNPs/SCNP-0-bonds"); 
    test = sys_test.initializePositions("input/SCNPs/SCNP-PRECURSOR-config");
    if (!test) {
        	std::cout << "problem with adding positions" << std::endl;
    }
    sys_test.initializeVelocitiesRandom(1.0); 
    initializeStepVector(OutputSteps, "input/filetimeoutput");
    OutputStepsIt = OutputSteps.begin();    
    
    std::cout << "radius of gyration: " << sys_test.Molecules[0].radiusOfGyration() << std::endl; 
    

    FILE* pdb {}; 
    pdb = fopen("pdb_equil_langevin_test.pdb", "w");
    
    Vector3d NewCOM (100., 100., 100.); 
    sys_test.setMoleculeCOM(0, NewCOM);
     
   
    
      
    sys_test.printPDB(pdb, 0); 
    
    sys_test.calculateForcesBrute(true);
    sys_test.updateVerletLists();
    sys_test.breakBonds();
    sys_test.makeBonds();
    //sys_test.calculateForces(true);
    std::cout << sys_test.PotentialEnergy() << std::endl; 
    ofstream gyr{"./results/reversibletest.dat"};
    timeval start {}, end {};
    gettimeofday(&start, NULL); 

    
    for (unsigned long long i = 0; i < TotalSteps; i++) {
        if (i == *OutputStepsIt) {
            sys_test.propagateLangevin(true);
            //sys_test.propagate(0.001, true); 
            sys_test.printStatistics(gyr, i*0.001);
            sys_test.printPDB(pdb,i);
            std::cout << 0.01*i << " " << sys_test.PotentialEnergy() << std::endl;
            //sys_test.printBonds();
            OutputStepsIt++;
        }
        else sys_test.propagateLangevin(); //sys_test.propagate(0.001); //

        
    }
    
    fclose(pdb); 
    
    gettimeofday(&end, NULL); 
    double realTime = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << "total time: " << realTime << " , time per step: " << realTime/TotalSteps << std::endl;
    
    return 0; 
}
