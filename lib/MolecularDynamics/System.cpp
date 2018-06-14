#include "System.h"

System::System(unsigned Lx, unsigned Ly, unsigned Lz, double aK, double aRadius0,  bool PBCon) :
    Cutoff {1.5}, 
    VerletRadius {2.0}, 
    VerletRadiusSq {4.0},
	CaptureDistance {1.3},
	K {aK},
	Radius0 {aRadius0},
    PBC{PBCon},
	ReversibleBonds {} {
        BoxSize[0] = Lx; 
        BoxSize[1] = Ly; 
        BoxSize[2] = Lz;
        Cells[0] = unsigned(BoxSize[0]/(sqrt(2)*VerletRadius)); 
        Cells[1] = unsigned(BoxSize[1]/(sqrt(2)*VerletRadius)); 
        Cells[2] = unsigned(BoxSize[2]/(sqrt(2)*VerletRadius));     
        CellSideLength[0] = BoxSize[0]/(double)Cells[0]; 
        CellSideLength[1] = BoxSize[1]/(double)Cells[1];
        CellSideLength[2] = BoxSize[2]/(double)Cells[2];
        CellList = std::vector<std::vector<std::vector<std::forward_list<Particle*>>>>(Cells[0], std::vector<std::vector<std::forward_list<Particle*>>>(Cells[1], std::vector<std::forward_list<Particle*>>(Cells[2], std::forward_list<Particle*>())));
    }
    

System::System(double aCutoff, double aVerletRadius, double aCapture, unsigned Lx, unsigned Ly, unsigned Lz, double aK, double aRadius0, bool PBCon) :
    Cutoff {aCutoff}, 
    VerletRadius {aVerletRadius}, 
    VerletRadiusSq {aVerletRadius*aVerletRadius},
	CaptureDistance {aCapture},
	K {aK},
	Radius0 {aRadius0},
    PBC {PBCon},
	ReversibleBonds {}{
        BoxSize[0] = Lx; 
        BoxSize[1] = Ly; 
        BoxSize[2] = Lz;
        Cells[0] = unsigned(BoxSize[0]/(sqrt(2)*VerletRadius)); 
        Cells[1] = unsigned(BoxSize[1]/(sqrt(2)*VerletRadius)); 
        Cells[2] = unsigned(BoxSize[2]/(sqrt(2)*VerletRadius));      
        CellSideLength[0] = BoxSize[0]/(double)Cells[0]; 
        CellSideLength[1] = BoxSize[1]/(double)Cells[1];
        CellSideLength[2] = BoxSize[2]/(double)Cells[2];
        CellList = std::vector<std::vector<std::vector<std::forward_list<Particle*>>>>(Cells[0], std::vector<std::vector<std::forward_list<Particle*>>>(Cells[1], std::vector<std::forward_list<Particle*>>(Cells[2], std::forward_list<Particle*>())));
    }
        
bool System::addMolecules(std::string filename, double mass) {
        std::ifstream file(filename, ios::in);
        if (!file.is_open()) return false;
        std::string line;
        unsigned current_mol{0}, mol, monos, mono_start{0};

        if (file.is_open()) {
            std::cout << "file " << filename << " successfully opened" << std::endl;
            while (file >> mol >> monos) {
            	if (mol != current_mol) {
            		std::cout << "order of molecules is wrong" << std::endl;
            		return false;
            	}
            	Molecules.push_back(Molecule(monos, mass, mono_start));
            	Molecules[mol].setChainBonds();
                mono_start += monos;
                current_mol++;
            }
        }
        std::cout << "initialized " << Molecules.size() << " molecules with a total of " << mono_start << " monomers." << std::endl;


        /*for (auto& mol : Molecules) {
            for (auto& mono : mol.Monomers) {
                std::cout << mono.Identifier << " ";
            }
            std::cout << std::endl;
        }*/

        return true;
    }


bool System::addLinks(std::string filename) {
        std::ifstream file {filename};
        if (!file.is_open()) return false;
        unsigned mol, bond1, bond2, numberOfLines, count;
        if (file.is_open()) {
            file >> numberOfLines;
            count = 0;
            while (file >> mol >> bond1 >> bond2) {
                Molecules[mol-1].setLink(bond1-1, bond2-1);
                count++;
            }
        }
        std::cout << "set " << count << " out of " << numberOfLines << " links" << std::endl;
        return true;
    }

bool System::addFunctional(std::string filename) {
    	std::ifstream file {filename};
    	if (!file.is_open()) return false;
    	else std::cout << "file " << filename << " successfully opened" << std::endl;
    	unsigned mol, mono, count {0};
    	while(file >> mol >> mono) {
    		if (mol >= Molecules.size()) {
    			std::cout << "Molecule #" << mol << " out of range!" << std::endl;
    			return false;
    		}
    		if (mono >= Molecules[mol].NumberOfMonomers) {
    			std::cout << "Index " << mono << "of Molecule #" << mol << " out of range!" << std::endl;
    			return false;
    		}
    		else {
    			unsigned index;
    			makeIndex(mol, mono, index);
    			Molecules[mol].Monomers[mono].Functional=true;
    			ReversibleBonds[index] = -1;
    			count++;
    		}
    	}
    	std::cout << "set " << count << " functional groups" << std::endl;
    	return true;
    }

bool System::initializePositions(std::string filename) {
        std::ifstream file {filename};
        if (!file.is_open()) return false;
        double x, y, z;
        unsigned count {0};
        if (file.is_open()) {
            for (auto& mol : Molecules) {
                for (auto& mono : mol.Monomers) {
                    if (file >> x >> y >> z) {
                        mono.Position(0) = x;
                        mono.Position(1) = y;
                        mono.Position(2) = z;
                        count++;
                    }
                    else {
                        std::cout << "only " << count << " monomer positions were initialized" << std::endl;
                        return false;
                    }
                }
            }
        }
        std::cout << "all " << count << " monomer positions were initialized" << std::endl;
        if (file >> x >> y >> z) std::cout << "...but there is more data..." << std::endl;
        return true;
    }

bool System::initializeVelocities(std::string filename) {
        std::ifstream file {filename};
        if (!file.is_open()) return false;
        double x, y, z;
        unsigned count {0};
        if (file.is_open()) {
            for (auto& mol : Molecules) {
                for (auto& mono : mol.Monomers) {
                    if (file >> x >> y >> z) {
                        mono.Velocity(0) = x;
                        mono.Velocity(1) = y;
                        mono.Velocity(2) = z;
                        count++;
                    }
                    else {
                        std::cout << "only " << count << " monomer velocities were initialized" << std::endl;
                        return false;
                    }
                }
            }
        }
        std::cout << "all " << count << " monomer velocities were initialized" << std::endl;
        if (file >> x >> y >> z) std::cout << "...but there is more data..." << std::endl;
        return true;
    }


void System::initializeVelocitiesRandom(double Temperature) {
        Vector3d COMVel {Vector3d::Zero()};
        unsigned totalMonomers {};
        double EKin {};
        for (auto& mol : Molecules) {
            totalMonomers += mol.NumberOfMonomers;
            for (auto& mono : mol.Monomers) {
                for (unsigned i = 0; i < 3; i++) {
                    mono.Velocity(i) = Rand::real_uniform() - 0.5;
                }
            }
            mol.removeAngularMomentum();
        }
        for (auto& mol : Molecules) {
            COMVel = mol.centerOfMassVelocity();
            for (auto& mono : mol.Monomers) {
                mono.Velocity -= COMVel;
                EKin += mono.Mass*mono.Velocity.squaredNorm();
            }
        }
        double scaling = sqrt(3.*totalMonomers*Temperature/EKin);
        for (auto& mol : Molecules) {
            for (auto& mono : mol.Monomers) {
                mono.Velocity *= scaling;
            }
        }
    }

bool System::arrangeMolecules() {
	bool overlap {true};
	unsigned count {};
	for (unsigned i = 0; i < Molecules.size(); i++) {
		count = 0;
		Vector3d COMPos {};
		overlap = true;
		while (overlap) {
			COMPos(0) = Rand::real_uniform(0.0, BoxSize[0]);
			COMPos(1) = Rand::real_uniform(0.0, BoxSize[1]);
			COMPos(2) = Rand::real_uniform(0.0, BoxSize[2]);
			setMoleculeCOM(i, COMPos);
			overlap = false;
            for (unsigned j = 0; j < i; j++) {
            	overlap = calculateOverlap(Molecules[i], Molecules[j]);
            	if (overlap) break;
            }
			count++;
			if (count%100 == 0) {
				std::cout << "already " << count << "trial insertions for molecule " << i << std::endl;
				if (count > 100000) {
					return false;
				}
			}
		}
	}
	std::cout << "successfully arranged all molecules without overlap" << std::endl;
	return !overlap;
}

bool System::arrangeMoleculesFCC() {
	bool overlap {true};
	unsigned N {};
	unsigned magic {};
	for (unsigned i = 0; i < 10; i++) {
		N = 4*pow(i,3);
		if (N >= Molecules.size()){
			magic = i;
			break;
		}
	}
	//magic = (unsigned)(pow(Molecules.size()/4, 1./3.));
	unsigned L = *std::min_element(BoxSize.begin(), BoxSize.end());
	double dL = double(L)/(2.*magic);
	unsigned count {0};
	Vector3d COMPos {};
	std::cout << "trying to initialize fcc lattice... \n magic: " << magic << ", dL: " << dL << std::endl;
	for (unsigned i = 0; i < 2*magic; i++) {
		for (unsigned j = 0; j < 2*magic; j++) {
			for (unsigned k = 0; k < 2*magic; k++) {
				if ((i+j+k)%2 == 0) {
					COMPos(0) = i*dL;
					COMPos(1) = j*dL;
					COMPos(2) = k*dL;
					setMoleculeCOM(count, COMPos);
					overlap = true;
					while (overlap) {
						Molecules[count].randomRotation();
						overlap = false;
						for (unsigned l = 0; l < count; l++) {
							overlap = calculateOverlap(Molecules[count], Molecules[l]);
							if (overlap) break;
						}
					}
					count++;

				}
				if (count >= Molecules.size()) {
					for (unsigned l = 0; l < Molecules.size(); l++) {
						for (unsigned m = l + 1; m < Molecules.size(); m++) {
							overlap = calculateOverlap(Molecules[l], Molecules[m]);
							if (overlap) {
								std::cout << "arranging molecules in an fcc lattice created an overlap between the following molecules: " << l << " , " << m << std::endl;
								return false;
							}
						}
					}
					std::cout << "placed all molecules." << std::endl;
					return true;
				}
			}
		}
	}
	if (count < Molecules.size() ) {
		std::cout << "only " << count << " molecules were placed!" << std::endl;
		return false;
	}
	return true;
}

void System::updateVerletLists() {
    //std::cout << "rebuilding Verlet list... " << std::endl;  
    for (auto& sheet : CellList) {
        for (auto& row : sheet) {
            for (auto& list : row) {
                list.clear(); 
            }    
        }
    }
    
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.VerletList.clear(); 
        }
    }
    
    std::array<int, 3> CellNumber{}; 
    
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            Vector3d imPos; 
            if (PBC) imPos = image(mono, BoxSize, 0.0);
            else imPos = mono.Position; 
            for (unsigned i = 0; i < 3; i++) {
                CellNumber[i] = (int)(imPos(i)/CellSideLength[i]); 
            }
            mono.VerletPosition = mono.Position; 
            try {
                CellList[CellNumber[0]][CellNumber[1]][CellNumber[2]].push_front(&mono); 
            }
            catch (std::exception& e) {
                throw CellAllocationException(mono, CellNumber);  
            }
        }
    }
    int l, m, n; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            Vector3d imPos; 
            if (PBC) imPos = image(mono, BoxSize, 0.0);
            else imPos = mono.Position; 
            for (unsigned i = 0; i < 3; i++) { 
                CellNumber[i] = (int)(imPos(i)/CellSideLength[i]); 
            }
            for (int i = CellNumber[0]-2; i < CellNumber[0]+3; i++) {
                
                for (int j = CellNumber[1]-2; j < CellNumber[1]+3; j++) {
                    
                    for (int k = CellNumber[2]-2; k < CellNumber[2]+3; k++) {
                        l = i - floor((double)i/Cells[0])*Cells[0];
                        m = j - floor((double)j/Cells[1])*Cells[1];
                        n = k - floor((double)k/Cells[2])*Cells[2]; 
                        for (auto& other : CellList[l][m][n]) {
                            if (other == &mono) continue; 
                            Vector3d relPos; 
                            if (PBC) relPos = relative(mono, *other, BoxSize, 0.0);
                            else relPos = other -> Position - mono.Position; 
                            double distance {relPos.squaredNorm()}; 
                            if (distance <= VerletRadiusSq) {
                                mono.VerletList.push_front(other); 
                            }
                        }
                    }  
                }
            }    
        }
    }
}    

void System::checkVerletLists() {
    Vector3d displacement {Vector3d::Zero()}; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            displacement = mono.Position - mono.VerletPosition; 
            if (displacement.norm() > (VerletRadius - Cutoff)*0.5) {
                updateVerletLists(); 
                return; 
            }
        }
    }
}

void System::breakBonds() {
	for (auto& bondedpair : ReversibleBonds) {
		if (bondedpair.second >= 0) {
			unsigned molIndex1, molIndex2, monoIndex1, monoIndex2;
			getIndex(bondedpair.first, molIndex1, monoIndex1);
			getIndex(bondedpair.second, molIndex2, monoIndex2);
			Vector3d relPos;
			if (PBC) relPos = relative(Molecules[molIndex1].Monomers[monoIndex1], Molecules[molIndex2].Monomers[monoIndex2], BoxSize, 0.0);
			else relPos = Molecules[molIndex2].Monomers[monoIndex2].Position - Molecules[molIndex1].Monomers[monoIndex1].Position;
			double radius {relPos.norm()};
			if (radius >= CaptureDistance) {
				ReversibleBonds.at(bondedpair.second) = -1;
				bondedpair.second = -1;
			}
		}
	}
}

void System::makeBonds() {
	for (auto& bondedpair : ReversibleBonds) {
		if (bondedpair.second < 0) {
			unsigned molIndex1 {}, monoIndex1 {};
			std::vector<unsigned> PossibleBonds {};
			getIndex(bondedpair.first, molIndex1, monoIndex1);
			Particle &mono {Molecules[molIndex1].Monomers[monoIndex1]};
			for (auto& other : mono.VerletList) {
				if (other -> Functional && ReversibleBonds.at(other -> Identifier) < 0) {
					Vector3d relPos;
					if (PBC) relPos = relative(mono, *other, BoxSize, 0.0);
					else relPos = other -> Position - mono.Position;
					double radius {relPos.norm()};
					if (radius < CaptureDistance) {
						PossibleBonds.push_back(other -> Identifier); //put possible candidates in list
						//bondedpair.second = other -> Identifier;
						//ReversibleBonds.at(bondedpair.second) = mono.Identifier;
					}
				}
			}
			// now choose one candidate
			if (PossibleBonds.size() > 1) {
				int ChosenIndex {Rand::int_uniform(0, PossibleBonds.size()-1)};
				bondedpair.second = PossibleBonds[ChosenIndex];
				ReversibleBonds.at(bondedpair.second) = mono.Identifier;
			}
			else if (PossibleBonds.size() == 1) {
				bondedpair.second = PossibleBonds[0];
				ReversibleBonds.at(bondedpair.second) = mono.Identifier;
			}
		}
	}
}

void System::getIndex(unsigned Index, unsigned& molIndex, unsigned& monoIndex) {
	unsigned currentmol {0}, currentmono {Molecules[0].NumberOfMonomers};
	while (Index >= currentmono) {
		currentmol++;
		currentmono += Molecules[currentmol].NumberOfMonomers;
	}
	monoIndex = Index - (currentmono - Molecules[currentmol].NumberOfMonomers);
	molIndex = currentmol;
}

void System::makeIndex(unsigned molIndex, unsigned monoIndex, unsigned& Index) {
	unsigned currentmono{0};
	for (unsigned i = 0; i < molIndex; i++) {
		currentmono += Molecules[i].NumberOfMonomers;
	}
	Index = currentmono+monoIndex;
}

void System::calculateForces(bool calcEpot) {
    //unsigned count_bonds {0}; 
    double force_abs {}; 
    Vector3d force {}; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) mono.Force = Vector3d::Zero(); 
        if (calcEpot) mol.Epot = 0.0;    
    }
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
        	int ReversiblyBonded {-1};
        	if (ReversibleBonds.count(mono.Identifier) > 0) ReversiblyBonded = ReversibleBonds.at(mono.Identifier);

            for (auto& other : mono.VerletList) {
                Vector3d relPos;
                if (PBC) relPos = relative(mono, *other, BoxSize, 0.0);
                else relPos = other -> Position - mono.Position;  
                double radius2 {relPos.squaredNorm()};
                if (calcEpot) {
                    mol.Epot += 0.5*RLJ_Potential(radius2);
                }
                force_abs = RLJ_Force(radius2); 
                if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                    throw(RLJException(mono.Identifier, mono.Position, mono.Velocity, other -> Identifier, other -> Position, other -> Velocity, force_abs)); 
                }
                force = relPos*force_abs; 
                mono.Force -= force;  
                if (ReversiblyBonded == other -> Identifier) {
                	double radius = sqrt(radius2);
                	if (calcEpot) {
                		mol.Epot += 0.5*Reversible_Bond_Potential(radius, Radius0, K);
                	}
                	force_abs = Reversible_Bond_Force(radius, Radius0, K);
                	if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
						throw(RLJException(mono.Identifier, mono.Position, mono.Velocity, other -> Identifier, other -> Position, other -> Velocity, force_abs));
					}
                	force = relPos*force_abs/radius;
                	mono.Force -= force;
                }
            }
            
            for (auto& bonded : mono.Bonds) {
                Vector3d relPos;
                if (PBC) relPos = relative(mono, *bonded, BoxSize, 0.0);
                else relPos = bonded -> Position - mono.Position;  
                double radius2 {relPos.squaredNorm()}; 
                if (calcEpot) {
                    mol.Epot += FENE_Potential(radius2);
                    
                }    
                force_abs = FENE_Force(radius2); 
                if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                    throw(FENEException(mono.Identifier, bonded->Identifier, force_abs)); 
                }
                force = relPos*force_abs; 
                mono.Force -= force; 
                bonded -> Force += force; 
            }
        }
    }
    /*std::map<int, bool> calculated;
    for (auto& pair : ReversibleBonds) {
    	calculated[pair.first] = false;
    }
    for (auto& pair : ReversibleBonds) {
    	if (pair.second >= 0) {
    		unsigned molIndex1, monoIndex1, molIndex2, monoIndex2;
    		getIndex(pair.first, molIndex1, monoIndex1);
    		getIndex(pair.second, molIndex2, monoIndex2);
    		Vector3d relPos
    	}

    }*/

}

void System::calculateForcesBrute(bool calcEpot) {
    double force_abs {}; 
    Vector3d force {}; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) mono.Force = Vector3d::Zero(); 
        if (calcEpot) mol.Epot = 0.0;    
    }
    for (auto& mol : Molecules) {
        for (unsigned i = 0; i < mol.NumberOfMonomers; i++) {
            for (unsigned j = i+1; j < mol.NumberOfMonomers; j++) {
                Vector3d relPos;             
                if (PBC) relPos = relative(mol.Monomers[i], mol.Monomers[j], BoxSize, 0.0);
                else relPos = mol.Monomers[j].Position - mol.Monomers[i].Position; 
                double radius2 {relPos.squaredNorm()}; 
                if (calcEpot) {
                    mol.Epot += RLJ_Potential(radius2); 
                }
                force_abs = RLJ_Force(radius2); 
                if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                    Particle* mono {&mol.Monomers[i]};
                    Particle* other {&mol.Monomers[j]};
                    throw(RLJException(mono -> Identifier, mono -> Position, mono -> Velocity, other -> Identifier, other -> Position, other -> Velocity, force_abs)); 
                }
                //std::cout << i << " " << j << " " << force_abs << std::endl;
                force = relPos*force_abs; 
                mol.Monomers[i].Force -= force; 
                mol.Monomers[j].Force += force;     
            }
            for (auto& bonded : mol.Monomers[i].Bonds) {
                Vector3d relPos;
                if (PBC) relPos = relative(mol.Monomers[i], *bonded, BoxSize, 0.0);
                else relPos = bonded->Position - mol.Monomers[i].Position;
                double radius2 {relPos.squaredNorm()}; 
                if (calcEpot) {
                    mol.Epot += FENE_Potential(radius2); 
                }
                force_abs = FENE_Force(radius2); 
                if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                    throw(FENEException(mol.Monomers[i].Identifier, bonded->Identifier, force_abs)); 
                }
                force = relPos*force_abs; 
                mol.Monomers[i].Force -= force; 
                bonded -> Force += force; 
            }

        }
    }

}

double System::calculateIntermolecularEnergy(unsigned firstIndex, unsigned secondIndex) {
	double Energy {0.0};
	Vector3d relPos {};
	double radius2 {};
	for (auto& mono_first : Molecules[firstIndex].Monomers) {
		for (auto& mono_second : Molecules[secondIndex].Monomers) {
			if (PBC) relPos = relative(mono_first, mono_second, BoxSize, 0.0);
			else relPos = mono_second.Position - mono_first.Position;
			radius2 = relPos.squaredNorm();
			Energy += RLJ_Potential(radius2);
		}
	}
	return Energy;
}

bool System::calculateOverlap(const Molecule& first, const Molecule& second) {
	Vector3d relPos {};
	double radius2 {};
	double Force_abs {0.0};
	for (auto& mono_first : first.Monomers) {
		for (auto& mono_second : second.Monomers) {
			relPos = relative(mono_first, mono_second, BoxSize, 0.0);
			radius2 = relPos.squaredNorm();
			Force_abs = RLJ_Force(radius2);
			if (fabs(Force_abs) > 1e4 || std::isinf(Force_abs) || std::isnan(Force_abs)) {
				return true;
			}
		}
	}
	return false;
}

void System::setMoleculeCOM(unsigned molIndex, Vector3d newCOM) {
    if (molIndex >= Molecules.size()) {
        std::cout << "Molecule number " << molIndex << " does not exist." << std::endl; 
        return; 
    }
    Vector3d currentCOM {Molecules[molIndex].centerOfMassPosition()}; 
    newCOM -= currentCOM; 
    for (auto& mono : Molecules[molIndex].Monomers) {
        mono.Position += newCOM; 
    }
}

void System::centerMolecule(unsigned molIndex) {
    Vector3d BoxCenter(BoxSize[0]*0.5, BoxSize[1]*0.5, BoxSize[2]*0.5); 
    setMoleculeCOM(molIndex, BoxCenter); 
}

void System::wrapMoleculesCOM() {
    if (!PBC) {
        std::cout << "periodic boundary conditions not activated!" << std::endl; 
        return;
    }
    for (auto& mol : Molecules) {
        wrapCOM(mol, BoxSize, 0.0, 0.0);
    }
}


void System::propagate(double dt, bool calcEpot) {
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity += (mono.Force/mono.Mass)*dt*0.5; 
            mono.Position += mono.Velocity*dt; 
            //TODO: boundary 
        }
    }
    checkVerletLists();
    breakBonds();
    makeBonds();
    calculateForces(calcEpot);
    //calculateForcesBrute(calcEpot);
    
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity += (mono.Force/mono.Mass)*dt*0.5;  
        }
    }
}

void System::propagateLangevin(double dt, double Temperature, double gamma, bool calcEpot) {
    double velcexp {exp(-gamma*dt)}; 
    double velcsqrt {sqrt(2.*gamma*Temperature)}; 
    double poscexp {(1.-velcexp)/(gamma)}; 
    double poscsqrt {velcsqrt/gamma}; 
    double tau1 {(1.-exp(-gamma*dt))/gamma}; 
    double tau2 {(1.-exp(-2.*gamma*dt))/(2.*gamma)}; 
    double zc11 {sqrt(tau2)};
    double zc21 {(tau1-tau2)/zc11}; 
    double zc22 {sqrt(dt - tau1*tau1/tau2)}; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            Vector3d Z1; 
            Vector3d Z2; 
            Vector3d O1 {Rand::real_normal(), Rand::real_normal(), Rand::real_normal()}; 
            Vector3d O2 {Rand::real_normal(), Rand::real_normal(), Rand::real_normal()}; 
            Z1 = zc11*O1;
            Z2 = zc21*O1 + zc22*O2; 
            Vector3d veltilde {mono.Velocity + mono.Force*dt/(2.*mono.Mass)}; 
            mono.Velocity = velcexp*veltilde + velcsqrt*Z1/sqrt(mono.Mass); 
            mono.Position += poscexp*veltilde + poscsqrt*Z2/sqrt(mono.Mass); 
        }
    }
    checkVerletLists();
    breakBonds();
    makeBonds();
    calculateForces(calcEpot);
    //calculateForcesBrute(calcEpot);
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity += (mono.Force/mono.Mass)*dt*0.5;  
        }
    }
}

/*void System::propagateLangevin(double dt, double Temperature, double gamma, bool calcEpot) { //Angel's version with gamma -> gamma*m 

    double velcsqrt {sqrt(2.*gamma*Temperature)}; 
    
    double poscsqrt {velcsqrt/gamma}; 
    double tau1 {(1-exp(-gamma*dt))/gamma}; 
    double tau2 {(1-exp(-2.*gamma*dt))/(2.*gamma)}; 
    double zc11 {sqrt(tau2)};
    double zc21 {(tau1-tau2)/zc11}; 
    double zc22 {sqrt(dt - tau1*tau1/tau2)}; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            double velcexp {exp(-gamma*dt/mono.Mass)};
            double poscexp {(1-velcexp)*mono.Mass/(gamma)};  
            Vector3d Z1; 
            Vector3d Z2; 
            Vector3d O1 {Rand::real_normal(), Rand::real_normal(), Rand::real_normal()}; 
            Vector3d O2 {Rand::real_normal(), Rand::real_normal(), Rand::real_normal()}; 
            Z1 = zc11*O1;
            Z2 = zc21*O1 + zc22*O2; 
            Vector3d veltilde {mono.Velocity + mono.Force*dt/(2.*mono.Mass)}; 
            mono.Velocity = velcexp*veltilde + velcsqrt*Z1/mono.Mass; 
            mono.Position += poscexp*veltilde + poscsqrt*Z2; 
        }
    }
    calculateForcesBrute(calcEpot); 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity += (mono.Force/mono.Mass)*dt*0.5;  
        }
    }
}*/


unsigned System::NumberOfParticles() {
    unsigned count {}; 
    for (auto& mol : Molecules) {
        count += mol.NumberOfMonomers; 
    }
    return count; 
}

unsigned System::NumberOfMolecules() {return Molecules.size();}

std::tuple<unsigned, unsigned> System::NumberOfBonds() {
	unsigned Bonds {0}, IntramolecularBonds {0};
	unsigned mol1, mol2, mono1, mono2;
	for (auto& bond : ReversibleBonds) {
		if (bond.second >= 0) {
			Bonds++;
			getIndex(bond.first, mol1, mono1);
			getIndex(bond.second, mol2, mono2);
			if (mol1==mol2) IntramolecularBonds++;
		}
	}
	if (Bonds%2 != 0) {
		std::cout << "uneven number of bonds!" << std::endl;
	}
	else Bonds /= 2;
	IntramolecularBonds /= 2;
	return std::make_tuple(Bonds, IntramolecularBonds);
}

double System::KineticEnergy() {
    double Ekin = 0.0; 
    for (auto& mol : Molecules) {
        Ekin += mol.KineticEnergy(); 
    }
    return Ekin; 
}


double System::PotentialEnergy() {
    double Epot {0.0}; 
    for (auto& mol : Molecules) {
        Epot += mol.PotentialEnergy(); 
    }
    return Epot; 
}

std::tuple<double, Matrix3d> System::GyrationTensor() {
    double rgyr {}; 
    Matrix3d gyrTensor {Matrix3d::Zero()}; 
    std::tuple<double, Matrix3d> gyrMol {}; 
    for (auto& mol : Molecules) {
        gyrMol = mol.GyrationTensor(); 
        gyrTensor += std::get<1>(gyrMol); 
        rgyr += std::get<0>(gyrMol); 
    }
    gyrTensor /= Molecules.size(); 
    rgyr /= Molecules.size(); 
    return std::make_tuple(rgyr, gyrTensor);    
}

Vector3d System::RotationFrequency() {
    Vector3d omega {Vector3d::Zero()}; 
    for (auto& mol : Molecules) {
        omega += mol.RotationFrequency(); 
    }
    omega /= Molecules.size(); 
    return omega; 
}

std::vector<double> System::calculateExtension(unsigned dim) {
    std::vector<double> Extension; 
    for (auto& mol : Molecules) {
        Extension.push_back(abs(mol.Monomers.front().Position(dim) - mol.Monomers.back().Position(dim))); 
        Extension.push_back((mol.Monomers.front().Position-mol.Monomers.back().Position).norm());
    }
    return Extension; 
}

void System::printPDB(FILE* pdb, int step, bool velocs) {
    int mol_count{0}; 
    fprintf(pdb, "MODEL     %d \n", step);
    for (auto& mol : Molecules) {
    	wrapCOM(mol, BoxSize, 0.0, 0.0);
		mol_count++;
		int mono_count{0}; 
		for (auto& mono : mol.Monomers) {
			char type {mono.Functional ? 'O' : 'C'};
		    if (velocs) {
		        fprintf(pdb, "ATOM %6d  %c   GLY   %3d     %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f \n", mono_count+1, type, mol_count, mono.Position(0), mono.Position(1), mono.Position(2), mono.Velocity(0), mono.Velocity(1), mono.Velocity(2));
		    }
			else {
			    fprintf(pdb, "ATOM %6d  %c   GLY   %3d     %7.3f %7.3f %7.3f \n", mono_count+1, type, mol_count, mono.Position(0), mono.Position(1), mono.Position(2));
		    }
		    mono_count++; 
		}	
        fprintf(pdb, "TER \n");
    }
    fprintf(pdb, "ENDMDL \n");   
    fflush(pdb);       
}

void System::printStatistics(std::ostream& os, double time) {
    double Ekin {KineticEnergy()}, Epot {PotentialEnergy()}; 
    std::tuple<unsigned, unsigned> Bonds {NumberOfBonds()};
    std::tuple<double, Matrix3d> GyrTuple {GyrationTensor()};
    Matrix3d GyrTensor {std::get<1>(GyrTuple)}; 
    os.precision(10); 
    os.width(16); 
    os << time << " "; 
    os.precision(6); 
    os.width(14); 
    os << Ekin << " ";
    os.width(14); 
    os << Epot << " "; 
    os.width(14);
    os << std::get<0>(Bonds) << " ";
    os.width(14);
    os << std::get<1>(Bonds) << " ";
    os.width(14);
    os << std::get<0>(GyrTuple); 
    for (unsigned i = 0; i < 3; i++) {
        for (unsigned j = 0; j < 3; j++) {
            os.width(14); 
            os << GyrTensor(i,j) << " ";
        } 
    }
    os << std::endl;  
}

void System::printBonds(std::ofstream& os) {
	for (auto& bond : ReversibleBonds) {
		if (bond.second >= 0) os << bond.first << " " << bond.second << std::endl;
	}
}

    
