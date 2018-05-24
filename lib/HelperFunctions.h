#include <fstream>
#include <istream>
#include <iostream>
#include <string>
#include <sstream>
#include <../Eigen/Eigen/Dense>

using namespace Eigen;

template <typename Type>
Type extractParameter(std::string key, std::ifstream& inputfile, bool& found) {
    inputfile.clear(); 
    inputfile.seekg(0, std::ios::beg);
    std::string line{}, fkey{}, equal{};
    Type param {};   
    found = false; 
    while (getline(inputfile, line)){
        if (line.find(key) != std::string::npos) {
            std::istringstream iss(line); 
            if (iss >> fkey >> equal >> param) {
                found = true; 
                return param; 
            }
        }
    }
    std::cout << "keyword " << key << " not found" << std::endl;
    return param; 
}

bool initializeStepVector(std::vector<unsigned long long>& vec, std::string filename) {
	vec.clear();
    std::ifstream file (filename, std::ios::in);
    if (!file.is_open()) {
        return false; 
    }
    unsigned long long Step{}; 
    while(file >> Step) {
        vec.push_back(Step); 
    } 
    return true; 
} 

struct ForceUpdate {
    unsigned long long Step; 
    Vector3d Force; 
    ForceUpdate(unsigned long long s, Vector3d f) : 
        Step {s}, 
        Force {f} {}
}; 

struct ConstraintUpdate {
    unsigned long long Step; 
    double Constraint; 
    ConstraintUpdate(unsigned long long s, double c) : 
        Step {s}, 
        Constraint {c} {}
}; 

bool initializeForceUpdateVector(std::vector<ForceUpdate>& vec, std::string filename) {
    std::ifstream file (filename, std::ios::in);
    vec.clear();
    if (!file.is_open()) {
        return false; 
    }
    unsigned long long Step{};
    double fx{}, fy{}, fz{}; 
    while(file >> Step >> fx >> fy >> fz) {
        Vector3d Force(fx, fy, fz); 
        vec.push_back(ForceUpdate(Step, Force)); 
    } 
    return true; 
} 

bool initializeConstraintUpdateVector(std::vector<ConstraintUpdate>& vec, std::string filename) {
    std::ifstream file (filename, std::ios::in);
    vec.clear();
    if (!file.is_open()) {
        return false; 
    }
    unsigned long long Step{};
    double c {};  
    while(file >> Step >> c) {
        vec.push_back(ConstraintUpdate(Step, c)); 
    } 
    return true; 
} 

struct TopologyAndConfigPool {
	std::string TopologyFile;
	std::string ConfigPoolFile;
	TopologyAndConfigPool(std::string Topology, std::string ConfigPool) :
		TopologyFile {Topology},
		ConfigPoolFile {ConfigPool} {}
};


bool fillConfigPool(std::vector<std::string>& vec, std::string filename) {
	vec.clear();
	std::ifstream file (filename, std::ios::in);
	if (!file.is_open()) return false;
	std::string configfile {};
	while(file >> configfile) vec.push_back(configfile);
	return true;
}

bool fillConfigPoolVector(std::vector<std::string>& vec, std::ifstream& inputfile) {
	vec.clear();
    inputfile.clear();
    inputfile.seekg(0, std::ios::beg);
    std::string line{}, key{"TopologyConfigPool:"};
    while (getline(inputfile, line)){
        if (line.find(key) != std::string::npos) {
            break;
        }
    }
	std::string configfile {};
	while(inputfile >> configfile) vec.push_back(configfile);
	return true;
}
