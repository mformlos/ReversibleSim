#include "System.h"
#include "HelperFunctions.h"

int main() {
    System sys_test(40,40,40, true);
    std::vector<unsigned long long> OutputSteps; 
    std::vector<unsigned long long>::iterator OutputStepsIt{};
    
    bool add{};
    add = sys_test.addMolecules("input/SCNPs/SCNP-0-chain", 5.0);
    if (add == false) {
    	std::cout << "adding molecules failed" << std::endl;

    }
    //sys_test.addLinks("/home/formanek/HYBRIDSIM/input/SCNPs/SCNP-0-bonds");
    sys_test.initializePositions("input/SCNPs/SCNP-0-config");
    sys_test.initializeVelocitiesRandom(1.0); 
    initializeStepVector(OutputSteps, "input/teststeps");
    OutputStepsIt = OutputSteps.begin();    
    
    std::vector<Vector3d> COMPos {};
    for (auto& mol : sys_test.Molecules) {
        COMPos.push_back(mol.centerOfMassPosition()); 
    }
    
    for (unsigned i = 0; i < COMPos.size(); i++) {
        std::cout << i << " " << COMPos[i].transpose() << std::endl; 
        for (unsigned j = i+1; j < COMPos.size(); j++) {
            double distance {(COMPos[i] - COMPos[j]).norm()}; 
            if (distance < 100) std::cout << i << " " << j << " " << distance << std::endl;     
        } 
    }
    
    double rgyr_mean {}; 
    unsigned mol_count {}; 
    for (auto& mol : sys_test.Molecules) {
        double rgyr {mol.radiusOfGyration()}; 
        rgyr_mean += rgyr; 
        ++mol_count; 
        std::cout << rgyr << std::endl; 
    }
    rgyr_mean /= mol_count; 
    std::cout << "mean radius of gyration: " << rgyr_mean << std::endl; 
    
    FILE* pdb {}; 
    pdb = fopen("pdb_equil_brute.pdb", "w");
    
    wrapCOM(sys_test.Molecules.front(), sys_test.BoxSize, 0.0, 0.0);
    sys_test.printPDB(pdb, 0); 
    
   
    
    std::cout << "new center of mass: " << sys_test.Molecules.front().centerOfMassPosition().transpose() << std::endl; 
    
    //sys_test.updateVerletLists();
    //sys_test.calculateForces(true);

    sys_test.calculateForcesBrute(true);
    std::cout << sys_test.PotentialEnergy() << std::endl; 
    ofstream gyr{"./results/NH-SCNP-0-stats-2.dat"};

    
    for (unsigned long long i = 0; i < 10000000; i++) {
        if (i == *OutputStepsIt) {
            sys_test.propagate(0.001, true); 
            sys_test.printStatistics(gyr, i*0.001, i);
            std::cout << 0.001*i << std::endl; 
            OutputStepsIt++;
        }
        else sys_test.propagate(0.001); 
    }
    
    fclose(pdb); 
    
    
    return 0; 
}
