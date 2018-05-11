#ifndef LIB_BOUNDARY_H_
#define LIB_BOUNDARY_H_

#include <array>

inline void wrap(Particle& part, const std::array<unsigned, 3>& BoxSize, const double& Shear, const double& delrx) {
    double cy {floor(part.Position(1)/BoxSize[1])}; 
    part.Position(0) -= cy*delrx; 
    part.Position(0) -= BoxSize[0]*floor(part.Position(0)/BoxSize[0]); 
    part.Position(1) -= BoxSize[1]*cy; 
    part.Position(2) -= BoxSize[2]*floor(part.Position(2)/BoxSize[2]); 
    part.Velocity(0) -= cy*Shear*BoxSize[1]; 
}

inline void wrapVelocityBack(Particle& part, Particle& image, const std::array<unsigned, 3>& BoxSize, const double& Shear, const double& delrx) {
    double cy {floor(part.Position(1)/BoxSize[1])};
    //if (cy != 0.0) std::cout << "vel before: " << part.Velocity.transpose() << " ";
    part.Velocity = image.Velocity;  
    part.Velocity(0) += cy*Shear*BoxSize[1]; 
    //if (cy != 0.0) std::cout << "vel after: " << part.Velocity.transpose() << std::endl;
}


inline Vector3d image(const Particle& part, const std::array<unsigned, 3>& BoxSize, const double& delrx) {
    double cy {floor(part.Position(1)/BoxSize[1])}; 
    Vector3d pos {Vector3d::Zero()}; 
    pos(0) = part.Position(0) - cy*delrx; 
    pos(0) -= BoxSize[0]*floor(pos(0)/BoxSize[0]); 
    pos(1) = part.Position(1) - BoxSize[1]*cy; 
    pos(2) = part.Position(2) - BoxSize[2]*floor(part.Position(2)/BoxSize[2]); 
    //part.Velocity(0) -= cy*Shear*BoxSize[1]; 
    return pos; 
}




inline Vector3d relative(const Particle& one, const Particle& two, const std::array<unsigned, 3>& BoxSize, const double& delrx) {
    Vector3d dist {two.Position - one.Position}; 
    double cy {round(dist(1)/BoxSize[1])}; 
    dist(0) -= cy*delrx; 
    dist(0) -= BoxSize[0]*round(dist(0)/BoxSize[0]); 
    dist(1) -= BoxSize[1]*cy; 
    dist(2) -= BoxSize[2]*round(dist(2)/BoxSize[2]); 
    return dist;
} 

inline Vector3d relative(const Vector3d& one, const Vector3d& two, const std::array<unsigned, 3>& BoxSize, const double& delrx) {
    Vector3d dist {two - one}; 
    double cy {round(dist(1)/BoxSize[1])}; 
    dist(0) -= cy*delrx; 
    dist(0) -= BoxSize[0]*round(dist(0)/BoxSize[0]); 
    dist(1) -= BoxSize[1]*cy; 
    dist(2) -= BoxSize[2]*round(dist(2)/BoxSize[2]); 
    return dist;
} 

inline void wrapCOM(Molecule& mol, const std::array<unsigned, 3>& BoxSize, const double& Shear, const double& delrx) {
    Vector3d COMPos {mol.centerOfMassPosition()}; 
    Vector3d Translate {Vector3d::Zero()}; 
    double cy {floor(COMPos(1)/BoxSize[1])};
    double velchange {};
    Translate(0) = COMPos(0)-cy*delrx; 
    Translate(0) = cy*delrx + BoxSize[0]*floor(Translate(0)/BoxSize[0]); 
    Translate(1) = BoxSize[1]*cy; 
    Translate(2) = BoxSize[2]*floor(COMPos(2)/BoxSize[2]); 
    velchange = cy*Shear*BoxSize[1]; 
    for (auto& mono : mol.Monomers) {
        mono.Position -= Translate; 
        mono.Velocity(0) -= velchange; 
    }
}

#endif
