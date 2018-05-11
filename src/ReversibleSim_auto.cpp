#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
#include <omp.h>
#include <csignal>
#include <string>
#include <vector>
#include "HelperFunctions.h"
#include "System.h"


using namespace Eigen;

int SignalCaught {}; 
void signalHandler(int signum) 
{
    SignalCaught = signum; 
    std::cout << "signal " << signum << " caught!" << std::endl; 
}

int main(int argc, char* argv[]) {
    SignalCaught = 0;
    signal(SIGINT, signalHandler); 
    unsigned Lx{}, Ly{}, Lz{}, Seed {};  
    unsigned long long TotalSteps{}, n {}, m {}; 
    double MDStep{}, SimTime{}, EquilTime{}, Temperature{}, Time{}, Gamma {}, ConstantK {}, ConstantR0 {};
    bool ParameterInitialized{}, Equilibrated{false};
    std::string OutputStepFile{}, MoleculeFile{}, FunctionalFile{}, ConfigFile{}, StatisticsFile{}, ConfigOutFile{};
    std::vector<unsigned long long> OutputSteps {};
    std::vector<unsigned long long>::iterator OutputStepsIt{};
    
    
    if (argc != 2) {
        std::cout << "usage: ./ReversibleSim.cpp PARAMETER-INPUT-FILE " << std::endl;
        return EXIT_FAILURE; 
    }
    
    //////////////////////////////
    
    ////// PARAMETER READS ///////
    
    std::ifstream inputfile(argv[1], std::ios::in);
    if (!inputfile.is_open()) {
        std::cout << "could not open file '" << argv[1] << "' , exiting" << std::endl; 
        return EXIT_FAILURE;  
    } 
    
    Lx = extractParameter<unsigned>("BoxX", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    Ly = extractParameter<unsigned>("BoxY", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    Lz = extractParameter<unsigned>("BoxZ", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    Temperature = extractParameter<double>("Temperature", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE; 
    Gamma = extractParameter<double>("Gamma", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE; 
    ConstantK = extractParameter<double>("ConstantK", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    ConstantR0 = extractParameter<double>("ConstantR0", inputfile, ParameterInitialized);
	if (!ParameterInitialized) return EXIT_FAILURE;
    MDStep = extractParameter<double>("MDStep", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    Time = extractParameter<double>("StartTime", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    SimTime = extractParameter<double>("SimTime", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    EquilTime = extractParameter<double>("EquilTime", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    Seed = extractParameter<double>("Seed", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    OutputStepFile = extractParameter<std::string>("OutputStepFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    MoleculeFile = extractParameter<std::string>("MoleculeFile", inputfile, ParameterInitialized);      
    if (!ParameterInitialized) return EXIT_FAILURE;
    FunctionalFile = extractParameter<std::string>("FunctionalFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE; 
    ConfigFile = extractParameter<std::string>("ConfigFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    StatisticsFile = extractParameter<std::string>("StatisticsFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    ConfigOutFile = extractParameter<std::string>("ConfigOutFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    
    inputfile.close(); 
    
    TotalSteps = (unsigned long long) (SimTime/MDStep); 
    std::cout << "Total number of MD steps: " << TotalSteps << std::endl;  
    
    ///////////////////////////////////////// 
     
    /////// PARAMETER CONTROL OUTPUT //////// 
    std::cout << "Lx is " << Lx << std::endl;
    std::cout << "Ly is " << Ly << std::endl;
    std::cout << "Lz is " << Lz << std::endl;
    std::cout << "Temperature is " << Temperature << std::endl;
    std::cout << "Friction Coefficient is " << Gamma << std::endl;
    std::cout << "Constant K is " << ConstantK << std::endl;
    std::cout << "Constant R0 is " << ConstantR0 << std::endl;
    std::cout << "MDStep is " << MDStep << std::endl;
    std::cout << "RNG seed is " << Seed << std::endl; 
    std::cout << "Equiltime is " << EquilTime << std::endl;
    std::cout << "OutputStepFile is " << OutputStepFile << std::endl;
    std::cout << "MoleculeFile is " << MoleculeFile << std::endl;
    std::cout << "FunctionalFile is " << FunctionalFile << std::endl;
    std::cout << "ConfigFile is " << ConfigFile << std::endl;

    /////////////////////////////////////
    
    /////// SYSTEM INITIALIZATION ///////
    System Sys(Lx, Ly, Lz, ConstantK, ConstantR0, true);
    
    if (!Sys.addMolecules(MoleculeFile, 1.0)) {
        std::cout << "MoleculeFile does not exist!" << std::endl; 
        return EXIT_FAILURE;     
    }
    
    if (!Sys.addFunctional(FunctionalFile)) {
        std::cout << "LinkFile does not exist!" << std::endl; 
        return EXIT_FAILURE; 
    } 
    
    if (!Sys.initializePositions(ConfigFile)) {
        std::cout << "ConfigFile does not exist or contains too little lines!" << std::endl; 
        return EXIT_FAILURE; 
    } 
 
    Sys.initializeVelocitiesRandom(Temperature); 
    
    if (Sys.NumberOfMolecules() == 1)  {
        Sys.centerMolecule(0);
    }
    else if (Sys.NumberOfMolecules() > 1){
        std::cout << "Code is supposed to be used with only 1 Molecule" << std::endl; 
        return EXIT_FAILURE;  
    }
    
    try {
    	Sys.breakBonds();
    	Sys.makeBonds();
        Sys.updateVerletLists();
        Sys.calculateForces();
        //Sys.calculateForcesBrute();
    }
    catch (const LibraryException &ex) {
        std::cout << ex.what() << std::endl; 
        std::cout << "bad initial configuration! Terminating program." << std::endl;   
        return EXIT_FAILURE; 
    }
    
    
    /////////////////////////////////////
    
    /////// OUTPUT INITIALIZATION ///////
    
    if (!initializeStepVector(OutputSteps, OutputStepFile)) {
        std::cout << "OutputStepFile does not exist!" << std::endl; 
        return EXIT_FAILURE; 
    }
    std::cout << "Output will be done " << OutputSteps.size() << " times. " << std::endl; 
    OutputStepsIt = OutputSteps.begin(); 
    
    
    std::ofstream StatisticsStream(StatisticsFile, std::ios::out | std::ios::trunc);
    FILE* PDBout{}; 
    
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
    ////// RANDOM ENGINE SEEDING & WARMUP //////
    Rand::seed(Seed); 
    Rand::warmup(10000); 
    
    ///////////////////////////////////
    ////// MAIN SIMULATION LOOP ///////
    m = -TotalSteps; 
    for (n = 0; n <= TotalSteps; n++, m++) {
        Time += MDStep; 
        try {
            if (m == *OutputStepsIt) Sys.propagateLangevin(MDStep, Temperature, Gamma, true);
            else Sys.propagateLangevin(MDStep, Temperature, Gamma, false);
        }
        catch (const LibraryException &ex) {
            PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "w");
            Sys.printPDB(PDBout, n, 1);
            fclose(PDBout);    
            std::cout << "Terminating..." << std::endl;   
            std::terminate(); 
        }
        if (!Equilibrated && Time >= EquilTime) {
            std::cout << "System equilibrated, starting production run" << std::endl;
            Equilibrated = true; 
            m = 0;
            Time = 0.0;
        }
        
        
        
        ////////// OUTPUT ////////////
        if (SignalCaught) {
            std::cout << "writing job data..." << std::endl; 
            PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "w");
            Sys.printPDB(PDBout, n, 1);
            fclose(PDBout);
            fclose(PDBout);  
            std::cout << "Terminating..." << std::endl;   
            std::terminate(); 
        }
        
        if (m == *OutputStepsIt) {
            Sys.printStatistics(StatisticsStream, Time); 
            PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "w");
            Sys.printPDB(PDBout, n, 1); 
            fclose(PDBout);
            std::cout << Time <<  std::endl; 
            OutputStepsIt++;        
        }
    }    
    
    //////////////////////////////////////
    ////// MAIN SIMULATION LOOP END///////
    
    gettimeofday(&end, NULL); 
    
    double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    std::cout << "total time: " << realTime << " , time per particle and step: " << realTime/Sys.Molecules[0].NumberOfMonomers/TotalSteps << std::endl;
    
    return EXIT_SUCCESS;

}
    

