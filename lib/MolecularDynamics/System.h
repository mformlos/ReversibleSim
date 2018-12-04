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
    double twoK; // 2*K
    double Radius0; //equilibrium radius for the reversible bond potential
    double DeltaT; //TimeStep
    double Temperature; 
    double Gamma; //Constant for Langevin thermostat
    double Epot; 
    /// Langevin Constants
    double velcexp;
    double velcsqrt; 
    double poscexp; 
    double poscsqrt; 
    double tau1; 
    double tau2; 
    double zc11;
    double zc21; 
    double zc22; 
    /////
    
    bool PBC;
    
    
    std::array<double,3> BoxSize; 
    std::array<unsigned,3> Cells; 
    std::array<double,3> CellSideLength; 
    
    std::vector<std::vector<std::vector<std::forward_list<Particle*>>>> CellList;
    std::array<std::array<int,3>,13> NeighbourDirections;
    
    std::vector<Molecule> Molecules; 
    
    std::map<int, int> ReversibleBonds;
    
    
    
    System(double, double, double, double aK = 29.6, double aRadius0 = 1.448, bool PBCon = true, double DeltaT = 0.01, double Temperature = 1.0, double Gamma = 0.05); //initialize only boxsize;
    System(double, double, double, double, double, double, double aK = 29.6, double aRadius0 = 1.448, bool PBCon = true, double DeltaT = 0.01, double Temperature = 1.0, double Gamma = 0.05); //initialize with cutoffs and capture distance
    
    void updateVerletLists(); 
    void checkVerletLists(); 
    void updateCellLists(); 
    
    void breakBonds();
    void makeBonds();
    void makeBondsCellList(); 

    void getIndex(unsigned, unsigned&, unsigned&);
    void makeIndex(unsigned, unsigned, unsigned&);

    void calculateForces(bool calcEpot=false); 
    void calculateForcesBrute(bool calcEpot=false); 
    void calculateForcesCellList(bool calcEpot=false); 
    double calculateIntermolecularEnergy(unsigned, unsigned);
    bool calculateOverlap(const Molecule&, const Molecule&);
    
    void wrapMoleculesCOM(); 
    
    // Initialize Molecules 
    bool addMolecules(std::string, double mass = 1.0); 
    bool addLinks(std::string);  
    bool addFunctional(std::string);
    bool initializePositions(std::string); 
    bool initializePositionsPDB(std::string); 
    bool initializeVelocities(std::string);
    void initializeVelocitiesRandom(double); 
    void setMoleculeCOM(unsigned, Vector3d); 
    void centerMolecule(unsigned); 
    bool arrangeMolecules();
    bool arrangeMoleculesFCC();
    bool setNeighbourDirections(std::string); 

    
    void propagate(bool calcEpot=false); 
    void propagateLangevin(bool calcEpot=false); 
    void scaleSystem(double);  //scale all directions equally
    void scaleSystem(double, double, double); //scale directions differently
    void updateCellListDimensions();  
    
    //getter:
    
    unsigned NumberOfParticles(); 
    unsigned NumberOfMolecules(); 
    std::tuple<unsigned, unsigned> NumberOfBonds();
    double KineticEnergy(); 
    double PotentialEnergy(); 
    std::tuple<double, Matrix3d> GyrationTensor(); 
    Vector3d RotationFrequency(); 
    std::vector<double> calculateExtension(unsigned); 
    
    
    void printPDB(FILE* pdb, int step, bool velocs = false); 
    void printStatistics(std::ostream& os, double time);

    void printBonds(std::ofstream& os);

};


#endif
