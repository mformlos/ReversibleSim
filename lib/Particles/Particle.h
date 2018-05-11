#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <iostream>
#include <forward_list>
#include <../Eigen/Eigen/Dense> 
#include "Potentials.h"
using namespace Eigen;


class Particle {
public: 
    //Members:
	int Identifier;
    double Mass;
    bool Functional;
    Vector3d Position; 
    Vector3d Velocity; 
    Vector3d Force;
    Vector3d VerletPosition;

    std::forward_list<Particle*> VerletList;
    std::forward_list<Particle*> Bonds;
    
    //Constructor:
    Particle(int);  //Initialize everything to 0, mass to 1.0
    Particle(int, double);  //initialize only mass
    Particle(int, double, Vector3d, Vector3d);
    Particle(int, double, bool, Vector3d, Vector3d); //initialize
    ~Particle() = default;
    
    void setBond(Particle&);
}; 


#endif


