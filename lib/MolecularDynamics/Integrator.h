#ifndef LIB_INTEGRATOR_H
#define LIB_INTEGRATOR_H 


class Integrator {
public: 
    double TimeStep; 
    double CutOff; 
    double VerletRadius; 
    
    std::array<unsigned,3> BoxSize; 
    std::array<unsigned,3> Cells; 
    std::array<double,3> CellSideLength; 
    
    std::vector<std::forward_list<MDParticle*>> CellList; 
    
    Integrator(double, double, unsigned, unsigned, unsigned);
    
    void propagate(Molecule& ) 


};

#endif
