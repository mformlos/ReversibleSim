#ifndef MOLECULE_H_
#define MOLECULE_H_

#include "Particle.h"
#include <fstream>

class Molecule {
public: 
    
    unsigned NumberOfMonomers; 
    double Epot; 
    std::vector<Particle> Monomers;
    
    Molecule(unsigned, int StartParticleIndex = 0); //Standard Initialization of N Particles
    Molecule(unsigned, double, int StartParticleIndex = 0); //Initialization of N Particles of mass M
    
    Particle& operator[](unsigned); //random access
    const Particle& operator[](unsigned) const;
    void push_back(Particle&);
    
    void setChainBonds(); 
    void setLink(unsigned, unsigned); 
    void translate(Vector3d); 
    void removeAngularMomentum(); 
    bool initializePositions(std::string); 
    bool addLinks(std::string, unsigned);
    bool addFunctional(std::string, unsigned);
    
    void calculateInternalForces(); 
    void calculateSpringForces(); 
    
    //Property getter functions
    Vector3d centerOfMassPosition(); 
    Vector3d centerOfMassVelocity(); 
    double KineticEnergy(); 
    double PotentialEnergy() {return Epot;}
    double radiusOfGyration(); 
    std::tuple<double, Matrix3d> GyrationTensor(); 
    Vector3d RotationFrequency(); 
    Matrix3d StressTensor(); 
    
    
    void printForces(FILE* f, int step); 
    
}; 



#endif
