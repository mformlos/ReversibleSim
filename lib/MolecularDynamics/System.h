#ifndef LIB_SYSTEM_H
#define LIB_SYSTEM_H 

#include <string>
#include <fstream>
#include <map>
#include "Molecule.h"
#include "Potentials.h"
#include "Rand.h"
#include "BoundaryConditions.h"
#include "../Exceptions/LibraryException.h"


class System {
public: 
    double Cutoff; 
    double VerletRadius; 
    double VerletRadiusSq;
    double CaptureDistance;
    double K;  //constant for the reversible bond potential
    double Radius0; //equilibrium radius for the reversible bond potential
    bool PBC;
    
    std::array<unsigned,3> BoxSize; 
    std::array<unsigned,3> Cells; 
    std::array<double,3> CellSideLength; 
    
    std::vector<std::vector<std::vector<std::forward_list<Particle*>>>> CellList;
    
    std::vector<Molecule> Molecules; 
    
    std::map<int, int> ReversibleBonds;
    
    
    
    System(unsigned, unsigned, unsigned, double aK = 29.6, double aRadius0 = 1.448, bool PBCon = true); //initialize only boxsize;
    System(double, double, double, unsigned, unsigned, unsigned, double aK = 29.6, double aRadius0 = 1.448, bool PBCon = true); //initialize with cutoffs and capture distance
    
    void updateVerletLists(); 
    void checkVerletLists(); 
    
    void breakBonds();
    void makeBonds();

    void getIndex(unsigned, unsigned&, unsigned&);
    void makeIndex(unsigned, unsigned, unsigned&);

    void calculateForces(bool calcEpot=false); 
    void calculateForcesBrute(bool calcEpot=false); 
    
    void wrapMoleculesCOM(); 
    
    // Initialize Molecules 
    bool addMolecules(std::string, double mass = 1.0); 
    bool addLinks(std::string);  
    bool addFunctional(std::string);
    bool initializePositions(std::string); 
    bool initializeVelocities(std::string);
    void initializeVelocitiesRandom(double); 
    void setMoleculeCOM(unsigned, Vector3d); 
    void centerMolecule(unsigned); 

    
    void propagate(double dt, bool calcEpot=false); 
    void propagateLangevin(double dt, double Temperature, double gamma=0.05, bool calcEpot=false); 
    
    //getter:
    
    unsigned NumberOfParticles(); 
    unsigned NumberOfMolecules(); 
    unsigned NumberOfBonds();
    double KineticEnergy(); 
    double PotentialEnergy(); 
    std::tuple<double, Matrix3d> GyrationTensor(); 
    Vector3d RotationFrequency(); 
    std::vector<double> calculateExtension(unsigned); 
    
    
    void printPDB(FILE* pdb, int step, bool velocs = false); 
    void printStatistics(std::ofstream& os, double time); 

    void printBonds();

};


#endif
