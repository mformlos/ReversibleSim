#include "Particle.h" 

Particle::Particle(int Index) :
	Identifier {Index},
	Mass { 1.0 },
	Functional { false },
    Position { Vector3d::Zero() }, 
    Velocity { Vector3d::Zero() }, 
	Force { Vector3d::Zero() },
	VerletPosition { Vector3d::Zero() },
	VerletList { },
    Bonds { } { }

Particle::Particle(int Index, double aMass) :
	Identifier {Index},
	Mass { aMass },
	Functional { false },
	Position { Vector3d::Zero() },
	Velocity { Vector3d::Zero() },
	Force { Vector3d::Zero() },
	VerletPosition { Vector3d::Zero() },
	VerletList { },
	Bonds { } { }

Particle::Particle(int Index, double aMass, Vector3d aPosition, Vector3d aVelocity) :
	Identifier {Index},
	Mass { aMass },
	Functional { false },
	Position { aPosition },
	Velocity { aVelocity },
	Force { Vector3d::Zero() },
	VerletPosition { Vector3d::Zero() },
	VerletList { },
	Bonds { } { }

Particle::Particle(int Index, double aMass, bool aFunctional, Vector3d aPosition, Vector3d aVelocity) :
	Identifier {Index},
	Mass { aMass },
	Functional { aFunctional },
	Position { aPosition },
	Velocity { aVelocity },
	Force { Vector3d::Zero() },
	VerletPosition { Vector3d::Zero() },
	VerletList { },
	Bonds { } { }


void Particle::setBond(Particle& bonded) {
    Bonds.push_front(&bonded); 
}


