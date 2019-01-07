#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
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
    unsigned Seed {}; 
    double Lx{}, Ly{}, Lz{};  
    unsigned long long TotalSteps{}, n {}, m {}; 
    double MDStep{}, SimTime{}, EquilTime{}, Temperature{}, Time{}, Density {}, Gamma {}, ConstantK {}, ConstantR0 {};
    bool ParameterInitialized{}, Equilibrated{false}, Arranged {};
    std::string OutputStepFile{}, MoleculeFile{}, FunctionalFile{}, ConfigFile{}, VelocFile{}, StatisticsFile{}, ConfigOutFile{}, NeighbourCellFile{}, ArrangeMolecules{};
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
    
    Lx = extractParameter<double>("BoxX", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    Ly = extractParameter<double>("BoxY", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    Lz = extractParameter<double>("BoxZ", inputfile, ParameterInitialized); 
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
    VelocFile = extractParameter<std::string>("VelocFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    StatisticsFile = extractParameter<std::string>("StatisticsFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    ConfigOutFile = extractParameter<std::string>("ConfigOutFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    ArrangeMolecules = extractParameter<std::string>("ArrangeMolecules", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    NeighbourCellFile = extractParameter<std::string>("NeighbourCellFile", inputfile, ParameterInitialized);
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
    std::cout << "Starttime is " << Time << std::endl;
    std::cout << "Totaltime is " << SimTime << std::endl;
    std::cout << "Equiltime is " << EquilTime << std::endl;
    std::cout << "OutputStepFile is " << OutputStepFile << std::endl;
    std::cout << "MoleculeFile is " << MoleculeFile << std::endl;
    std::cout << "VelocityFile is " << VelocFile << std::endl;
    std::cout << "FunctionalFile is " << FunctionalFile << std::endl;
    std::cout << "ConfigFile is " << ConfigFile << std::endl;
    std::cout << "NeighbourCellFile is " << NeighbourCellFile << std::endl;

    ////// RANDOM ENGINE SEEDING & WARMUP //////
    Rand::seed(Seed);
    Rand::warmup(10000);

    /////////////////////////////////////
    
    /////// SYSTEM INITIALIZATION ///////

    System Sys(1.3, 1.5, 1.3, Lx, Ly, Lz, ConstantK, ConstantR0, true, MDStep, Temperature, Gamma);
    
    std::cout << "Neighbourcells in each direction: " << Sys.Cells[0] << " " <<  Sys.Cells[1] << " " << Sys.Cells[2] << std::endl; 
    std::cout << "Cell side lengths: " << Sys.CellSideLength[0] << " " <<  Sys.CellSideLength[1] << " " << Sys.CellSideLength[2] << std::endl;
    
    if (!Sys.setNeighbourDirections(NeighbourCellFile)){
        return EXIT_FAILURE;
    }
    
    if (!Sys.addMolecules(MoleculeFile, 1.0)) {
        std::cout << "MoleculeFile does not exist!" << std::endl; 
        return EXIT_FAILURE;     
    }

    Density = double(Sys.NumberOfParticles())/(Lx*Ly*Lz);
    std::cout << "Number of chains: " << Sys.NumberOfMolecules() << ", Number of monomers: " << Sys.NumberOfParticles() << std::endl;
    std::cout << "Number density of the system: " << Density << std::endl;
    if (!Sys.addFunctional(FunctionalFile)) {
        std::cout << "LinkFile does not exist!" << std::endl; 
        return EXIT_FAILURE; 
    } 
    
    if (!Sys.initializePositions(ConfigFile)) {
        std::cout << "ConfigFile does not exist or contains too little lines!" << std::endl; 
        return EXIT_FAILURE; 
    } 
    
    if (VelocFile == "RANDOM") Sys.initializeVelocitiesRandom(Temperature);
    else {
    	if (!Sys.initializeVelocities(VelocFile)) {
    		std::cout << "VelocFile does not exist or contains too little lines!" << std::endl;
    		return EXIT_FAILURE;
    	}
    }
    
    unsigned count {0};
    if (ArrangeMolecules == "yes") {
        Arranged = false; 
        while (!Arranged) {
            Arranged = Sys.arrangeMoleculesFCC();
            try {
                Sys.updateCellLists(); 
                Sys.calculateForcesCellList(); 
		        /*Sys.updateVerletLists();
		        Sys.calculateForces(true);*/
		    }
		    catch (const LibraryException &ex) {
		        Arranged = false; 
		        std::cout << "bad placement, trying again..." << std::endl; 
		    }
		    if (count >= 50) {
		        std::cout << "not able to place molecules! Terminating program." << std::endl;
			    return EXIT_FAILURE;
		    }
	    }
    }

    std::ofstream StatisticsStream(StatisticsFile, std::ios::out | std::ios::app);
	FILE* PDBout{};

    /////////////////////////////////////
	/// print PDB to check intitial placement
	if (Time < EquilTime) {
		PDBout = fopen("initial.pdb", "w");
		Sys.printPDB(PDBout, 0, 1);
		fclose(PDBout);
		Sys.printStatistics(StatisticsStream, -EquilTime);
	}
	////////////////////////////////////


	try {
		Sys.breakBonds();
		Sys.updateCellLists(); 
		Sys.makeBondsCellList(); 
		Sys.calculateForcesCellList(); 
		/*Sys.makeBonds();
		Sys.updateVerletLists();
		Sys.calculateForces(true);*/
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


    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
    ///////////////////////////////////
    ////// MAIN SIMULATION LOOP ///////
    if (Time < EquilTime) {
    	m = TotalSteps;
    	n = 0;
    }
    else {
    	n = (unsigned long long) (Time/MDStep);
    	m = n;
    	Equilibrated = true;
    }
    std::cout << " Time: " << Time << ", m: " << m << ", n: " << n << ", next output at m = " << *OutputStepsIt << std::endl;
    for (; n <= TotalSteps; n++, m++) {
        Time += MDStep; 
        try {
            if (m == *OutputStepsIt) Sys.propagateLangevin(true);
            else Sys.propagateLangevin(false);
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
            std::cout << "Terminating..." << std::endl;   
            std::terminate(); 
        }
        
        if (m == *OutputStepsIt) {
            PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "w");
            Sys.printPDB(PDBout, n, 1); 
            fclose(PDBout);
            Sys.printStatistics(StatisticsStream, Time);
        	std::ofstream Bonds ("Bonds/Bonds"+std::to_string(m),std::ios::out | std::ios::trunc);
        	Sys.printBonds(Bonds);
        	Bonds.close();
            std::cout << Time <<  std::endl; 
            OutputStepsIt++;        
        }
    }    
    
    //////////////////////////////////////
    ////// MAIN SIMULATION LOOP END///////
    
    gettimeofday(&end, NULL); 
    
    double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    unsigned TotalMonomers {Sys.NumberOfParticles()};
    std::cout << "total time: " << realTime << " , time per particle and step: " << realTime/TotalMonomers/TotalSteps << std::endl;
    
    return EXIT_SUCCESS;

}
    

