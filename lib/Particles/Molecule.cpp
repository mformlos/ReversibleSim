#include "Molecule.h" 

Molecule::Molecule(unsigned N, int IdZero) :
    NumberOfMonomers {N}, 
    Epot { } 
    {
        Monomers.reserve(N); 
        for (unsigned i = 0; i < NumberOfMonomers; i++) {
            Monomers.push_back(Particle(i+IdZero));
        }
    }

Molecule::Molecule(unsigned N, double Mass, int IdZero) : 
    NumberOfMonomers {N},
    Epot { }
    {
        Monomers.reserve(N); 
        for (unsigned i = 0; i < NumberOfMonomers; i++) {
            Monomers.push_back(Particle(i+IdZero, Mass));
        }
    }

Particle& Molecule::operator[](unsigned i) {return Monomers[i];}

const Particle& Molecule::operator[](unsigned i) const {return Monomers[i];}

void Molecule::push_back(Particle& part) {
    part.Identifier = NumberOfMonomers; 
    Monomers.push_back(part); 
    ++NumberOfMonomers; 
}

void Molecule::setChainBonds() {
    for (unsigned i = 0; i < NumberOfMonomers-1; i++) {
        Monomers[i].setBond(Monomers[i+1]); 
    }
}

void Molecule::setLink(unsigned first, unsigned second) {
    Monomers[first].setBond(Monomers[second]); 
}

bool Molecule::initializePositions(std::string filename) {
    std::ifstream file {filename};
    if (!file.is_open()) return false; 
    std::string dump; 
    double x, y, z;  
    unsigned count {0}; 
    if (file.is_open()) {  
        file >> dump >> dump; 
        for (auto& mono : Monomers) {
            if (file >> dump >> dump >> dump >> dump >> dump >> x >> y >> z >> dump >> dump >> dump) {
                mono.Position(0) = x; 
                mono.Position(1) = y; 
                mono.Position(2) = z; 
                count++; 
             }
             else {
                std::cout << "only " << count << " monomers were initialized" << std::endl; 
                return false; 
             }
        }
    }
    //std::cout << "all " << count << " monomers were initialized" << std::endl; 
    if (file >> dump >> dump >> dump >> dump >> dump >> x >> y >> z >> dump >> dump >> dump) std::cout << "...but there is more data..." << std::endl; 
    return true; 
}

bool Molecule::addLinks(std::string filename, unsigned molspecifier) {
    std::ifstream file {filename};
    if (!file.is_open()) return false; 
    unsigned mol, bond1, bond2, numberOfLines, count;  
    if (file.is_open()) {
        file >> numberOfLines;
        count = 0;
        while (file >> mol >> bond1 >> bond2) {
            if (mol == molspecifier){
                setLink(bond1-1, bond2-1);
                count++;  
            }
        }    
    } 
    std::cout << "set " << count << " out of " << numberOfLines << " links" << std::endl; 
    return true; 
}

bool Molecule::addFunctional(std::string filename, unsigned molspecifier) {
	std::ifstream file {filename};
	if (!file.is_open()) return false;
	unsigned mol, number, numberOfLines, count {0};
	if (file.is_open()) {
		file >> numberOfLines;
		while (file >> mol >> number) {
			if (mol == molspecifier) {
				Monomers[number].Functional = true;
				count++;
			}
		}
	}
	std::cout << "Molecule " <<  molspecifier << " has " << count << "functional groups." << std::endl;
	return true;
}


Vector3d Molecule::centerOfMassPosition() {
    Vector3d COMPos {Vector3d::Zero()}; 
    for (auto& mono : Monomers) {
        COMPos += mono.Position; 
    }
    COMPos /= (double)NumberOfMonomers; 
    return COMPos; 
}

Vector3d Molecule::centerOfMassVelocity() {
    Vector3d COMVel {Vector3d::Zero()}; 
    for (auto& mono : Monomers) {
        COMVel += mono.Velocity; 
    }
    COMVel /= NumberOfMonomers; 
    return COMVel; 
}

void Molecule::translate(Vector3d vec) {
    for (auto& mono : Monomers) {
        mono.Position += vec; 
    }
}

void Molecule::removeAngularMomentum() {
    Vector3d COMPos {centerOfMassPosition()}; 
    Vector3d omega {RotationFrequency()}; 
    for (auto& mono : Monomers) {
        mono.Velocity += (mono.Position-COMPos).cross(omega); 
    }
}

double Molecule::KineticEnergy() {
    double Ekin {0.0}; 
    for (auto& mono : Monomers) {
        Ekin += mono.Mass*mono.Velocity.squaredNorm(); 
    }
    Ekin /= 2; 
    return Ekin; 
}

double Molecule::radiusOfGyration() {
    double rgyr {}; 
    for (unsigned i = 0; i < NumberOfMonomers; i++) {
        for (unsigned j = i+1; j < NumberOfMonomers; j++) {
            Vector3d distance {Monomers[i].Position - Monomers[j].Position}; 
            rgyr += distance.squaredNorm(); 
        }
    }
    rgyr /= pow(NumberOfMonomers,2); 
    rgyr = sqrt(rgyr); 
    return rgyr; 
}

std::tuple<double, Matrix3d> Molecule::GyrationTensor() {
    Matrix3d gyrTensor {Matrix3d::Zero()}; 
    Vector3d COMPos {centerOfMassPosition()}; 
    double rgyr { }; 
    for (auto& mono : Monomers) {
        Vector3d relPos {mono.Position-COMPos}; 
        for (unsigned alpha = 0; alpha < 3; alpha++) {
            for (unsigned beta = alpha; beta < 3; beta++) {
                gyrTensor(alpha, beta) += relPos(alpha)*relPos(beta); 
            }
        }
    }
    gyrTensor /= NumberOfMonomers; 
    gyrTensor(1,0) = gyrTensor(0,1);
    gyrTensor(2,0) = gyrTensor(0,2);
    gyrTensor(2,1) = gyrTensor(1,2); 
    rgyr = gyrTensor(0,0) + gyrTensor(1,1) + gyrTensor(2,2);
    rgyr = sqrt(rgyr); 
    return std::make_tuple(rgyr, gyrTensor);  
}

Vector3d Molecule::RotationFrequency() {
    Vector3d omega {Vector3d::Zero()}; 
    Matrix3d inertiaTensor {Matrix3d::Zero()}; 
    Vector3d angularMomentum {Vector3d::Zero()};
    Vector3d COMPos {centerOfMassPosition()}; 
    Vector3d COMVel {centerOfMassVelocity()};  
    for (auto& mono : Monomers) {   
        Vector3d relPos {mono.Position - COMPos}; 
        Vector3d relVel {mono.Velocity - COMVel}; 
        double rsqr {relPos.squaredNorm()};
        for (unsigned alpha = 0; alpha < 3; alpha++) {
            inertiaTensor(alpha, alpha) += rsqr; 
            for (unsigned beta = alpha; beta < 3; beta++) {
                inertiaTensor(alpha, beta) -= relPos(alpha)*relPos(beta);       
            }
        }
        angularMomentum += relPos.cross(relVel); 
    }
    inertiaTensor(1,0) = inertiaTensor(0,1); 
    inertiaTensor(2,0) = inertiaTensor(0,2); 
    inertiaTensor(2,1) = inertiaTensor(1,2); 
    omega = inertiaTensor.ldlt().solve(angularMomentum); 
    return omega; 
}


void Molecule::printForces(FILE* f, int step) {
    fprintf(f, "MODEL     %d \n", step);
    for (auto& mono : Monomers) {
        fprintf(f, "%12.4f %12.4f %12.4f \n", mono.Force(0), mono.Force(1), mono.Force(2));
    }
    fprintf(f, "\n"); 
}
 

void Molecule::calculateInternalForces() {
    double force_abs {}; 
    Vector3d force {}; 
    for (auto& mono : Monomers) mono.Force = Vector3d::Zero(); 
    
    for (unsigned i = 0; i < NumberOfMonomers; i++) {
        for (unsigned j = i+1; j < NumberOfMonomers; j++) {
            Vector3d relPos;             
            relPos = Monomers[j].Position - Monomers[i].Position; 
            double radius2 {relPos.squaredNorm()}; 
            force_abs = RLJ_Force(radius2); 
            force = relPos*force_abs; 
            Monomers[i].Force -= force; 
            Monomers[j].Force += force;     
        }
        for (auto& bonded : Monomers[i].Bonds) {
            Vector3d relPos;
            relPos = bonded->Position - Monomers[i].Position;
            double radius2 {relPos.squaredNorm()}; 
            force_abs = FENE_Force(radius2); 
            force = relPos*force_abs; 
            Monomers[i].Force -= force; 
            bonded -> Force += force; 
        }
    }
}

void Molecule::calculateSpringForces() {
    double force_abs {}; 
    Vector3d force {}; 
    for (auto& mono : Monomers) mono.Force = Vector3d::Zero(); 
    
    for (unsigned i = 0; i < NumberOfMonomers; i++) {
        for (auto& bonded : Monomers[i].Bonds) {
            Vector3d relPos;
            relPos = bonded->Position - Monomers[i].Position;
            double radius2 {relPos.squaredNorm()}; 
            force_abs = FENE_Force(radius2); 
            force = relPos*force_abs; 
            Monomers[i].Force -= force; 
            bonded -> Force += force; 
        }
    }
}


Matrix3d Molecule::StressTensor() {
    Vector3d COMPos; 
    Vector3d RelPos;
    Matrix3d Stress{Matrix3d::Zero()}; 
    
    COMPos = centerOfMassPosition(); 
    for (auto& mono : Monomers) {
        RelPos = mono.Position - COMPos; 
        for (unsigned i = 0; i < 3; i++) {
            for (unsigned j = 0; j < 3; j++) {
                Stress(i, j) += RelPos(i)*mono.Force(j); 
            }
        }
    }
    
    //Stress(1,0) = Stress(0,1); 
    //Stress(2,0) = Stress(0,2); 
    //Stress(2,1) = Stress(1,2); 
    return Stress; 

}

