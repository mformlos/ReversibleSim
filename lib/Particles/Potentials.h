#ifndef LIB_POTENTIALS_H_
#define LIB_POTENTIALS_H_

#include <cmath>

inline double FENE_Potential(double r2) {
    double potential {0.0}; 
    //potential = -15.0*2.25*log(1.0 - r2/2.25);
    potential = -33.75*log(1.0 - r2/2.25);
    return potential; 
}

inline double FENE_Force(double r2) {
    double force {0.0}; 
    force = -30.0/(1.0 - (r2/2.25));
    return force; 
}

inline double RLJ_Potential(double r2) {
    double potential {0.0}; 
    if (r2 < 1.25992105) {
		double rm6 { 1.0 / r2 };
		rm6 *= rm6*rm6;
		potential = 4.0*rm6*(rm6 - 1.0) + 1.0;
	}
	return potential;
} 

inline double RLJ_Force(double r2) {
    double force{0.0};
	if (r2 <1.25992105) {
		double rm2 { 1.0 / r2 };
		double rm6 { rm2*rm2*rm2 };
		force = 24.*rm2*rm6*(2.*rm6 - 1);
	}
	return force;
}

inline double LJSurf_Potential(double r2, double eps = 1.0) {
    double potential {0.0}; 
    if (r2 < 6.25) {
        double rm6 {1.0/r2}; 
        rm6 *= rm6*rm6; 
        potential = 4.0*eps*rm6*(rm6-1.0); 
    }
    return potential; 
}

inline double LJSurf_Force(double r2, double eps = 1.0) {
    double force{0.0};
	if (r2 < 6.25) {
		double rm2 { 1.0 / r2 };
		double rm6 { rm2*rm2*rm2 };
		force = 24.*eps*rm2*rm6*(2.*rm6 - 1);
	}
	return force;
}

inline double Harmonic_Potential(double r2, double k) {
    double potential {k*r2/2.0}; 
    return potential; 
}

inline double Reversible_Bond_Potential(double r, double r0, double K = 29.6) {
	double potential {0.0};
	potential = K*(2*exp(r0-r)-exp(2*(r0-r)));
	return potential;
}

inline double Reversible_Bond_Force(double r, double r0, double twoK) {
	double force {0.0};
	//force = 2.0*K*(exp(r0-r)-exp(2.0*(r0-r)));
	force = twoK*(exp(r0-r)-exp(2.0*(r0-r)));
	return force;
}



#endif
