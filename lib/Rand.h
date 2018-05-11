#ifndef NEWRAND_H_
#define NEWRAND_H_

#include <omp.h>
#include <thread>
#include <random>
#include <chrono>
#include <algorithm>
#include <iostream>



using namespace std;

class Rand{
public:
	static thread_local mt19937_64 generator;
	static thread_local minstd_rand0 lc_generator; 
    static thread_local normal_distribution<double> dis_normal;
	static thread_local uniform_real_distribution<double> dis_uniform;
	static thread_local uniform_int_distribution<int> dis_intuniform;
	static thread_local chi_squared_distribution<double> dis_chisquared;
	static thread_local gamma_distribution<double> dis_gamma;


	static double real_normal(double, double);
	static double real_normal();
	static double real_uniform();
	static double real_uniform(double);
	static double real_uniform(double, double);
	static double real_chisquared(unsigned n);
	//static double real_gamma(double shape);
	static double real_gamma(double, double);

	static void seed(int);
	static void seed(std::seed_seq); 
	static void warmup(unsigned); 

};

#endif /* RAND_H_ */
