/* Simulation to equilibrate system towards a given box size by 
scaling it by "Scaling" every "ScalingStep" steps. No functional groups are included. 
*/

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
    double Lx{}, Ly{}, Lz{}, StartLx {}, StartLy {}, StartLz {};  
    unsigned long long n {}, ScalingStep {}, SamplingStep {}; 
    double MDStep{}, Temperature{}, Time{0.0}, Density {}, Gamma {}, ConstantK {}, ConstantR0 {}, ScalingX {0.99}, ScalingY{sqrt(ScalingX)}, ScalingZ{sqrt(ScalingX)};
    bool ParameterInitialized{}, BoxSizeReached{false}, Arranged {};
    std::string MoleculeFile{}, ConfigFile{}, VelocFile{}, StatisticsFile{}, ConfigOutFile{}, LinkFile{},  NeighbourCellFile{}, ArrangeMolecules{"no"};
    
    
    if (argc != 2) {
        std::cout << "usage: ./Equilibrate_Box PARAMETER-INPUT-FILE " << std::endl;
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
    StartLx = extractParameter<double>("StartBoxX", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    StartLy = extractParameter<double>("StartBoxY", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    StartLz = extractParameter<double>("StartBoxZ", inputfile, ParameterInitialized); 
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
    ScalingX = extractParameter<double>("ScalingX", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    ScalingY = extractParameter<double>("ScalingY", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    ScalingZ = extractParameter<double>("ScalingZ", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    ScalingStep = extractParameter<double>("ScalingStep", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    SamplingStep = extractParameter<double>("SamplingStep", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    Seed = extractParameter<double>("Seed", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    MoleculeFile = extractParameter<std::string>("MoleculeFile", inputfile, ParameterInitialized);      
    if (!ParameterInitialized) return EXIT_FAILURE;
    LinkFile = extractParameter<std::string>("LinkFile", inputfile, ParameterInitialized);
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

    ///////////////////////////////////////// 
     
    /////// PARAMETER CONTROL OUTPUT //////// 
    std::cout << "Lx is " << Lx << std::endl;
    std::cout << "Ly is " << Ly << std::endl;
    std::cout << "Lz is " << Lz << std::endl;
    std::cout << "StartLx is " << StartLx << std::endl;
    std::cout << "StartLy is " << StartLy << std::endl;
    std::cout << "StartLz is " << StartLz << std::endl;
    std::cout << "Temperature is " << Temperature << std::endl;
    std::cout << "Friction Coefficient is " << Gamma << std::endl;
    std::cout << "Constant K is " << ConstantK << std::endl;
    std::cout << "Constant R0 is " << ConstantR0 << std::endl;
    std::cout << "MDStep is " << MDStep << std::endl;
    std::cout << "RNG seed is " << Seed << std::endl; 
    std::cout << "MoleculeFile is " << MoleculeFile << std::endl;
    std::cout << "VelocityFile is " << VelocFile << std::endl;
    std::cout << "ConfigFile is " << ConfigFile << std::endl;
    std::cout << "Arrange molecules in FCC lattice? " << ArrangeMolecules << std::endl; 
    std::cout << "NeighbourCellFile is " << NeighbourCellFile << std::endl;
    std::cout << "ScalingX is " << ScalingX << std::endl;
    std::cout << "ScalingY is " << ScalingY << std::endl;
    std::cout << "ScalingZ is " << ScalingZ << std::endl;

    
    ////// RANDOM ENGINE SEEDING & WARMUP //////
    Rand::seed(Seed);
    Rand::warmup(10000);
    
    /////////////////////////////////////
    
    /////// SYSTEM INITIALIZATION ///////
    /*if (StartLx <= Lx || StartLy <= Ly || StartLz <= Lz) {
        std::cout << "Starting box size must be bigger than finishing box size! " << std::cout; 
        return EXIT_FAILURE; 
    }*/
    
    System Sys(1.3, 1.5, 1.3, StartLx, StartLy, StartLz, ConstantK, ConstantR0, true, MDStep, Temperature, Gamma);
    //System Sys(StartLx, StartLy, StartLz, ConstantK, ConstantR0, true, MDStep, Temperature, Gamma);
    
    std::cout << "Neighbourcells in each direction: " << Sys.Cells[0] << " " <<  Sys.Cells[1] << " " << Sys.Cells[2] << std::endl; 
    std::cout << "Cell side lengths: " << Sys.CellSideLength[0] << " " <<  Sys.CellSideLength[1] << " " << Sys.CellSideLength[2] << std::endl;
    
    if (!Sys.setNeighbourDirections(NeighbourCellFile)){
        return EXIT_FAILURE;
    }
    
    if (!Sys.addMolecules(MoleculeFile, 1.0)) {
        std::cout << "MoleculeFile does not exist!" << std::endl; 
        return EXIT_FAILURE;     
    }
    
    if (!Sys.addLinks(LinkFile)) {
        std::cout << "LinkFile does not exist!" << std::endl; 
        return EXIT_FAILURE; 
    } 
    
    if (!Sys.initializePositionsPDB(ConfigFile)) {
        std::cout << "ConfigFile does not exist or contains too little lines!" << std::endl; 
        return EXIT_FAILURE; 
    } 
    
    if (VelocFile == "RANDOM") Sys.initializeVelocitiesRandom(Temperature);
    
    /*else {
    	if (!Sys.initializeVelocities(VelocFile)) {
    		std::cout << "VelocFile does not exist or contains too little lines!" << std::endl;
    		return EXIT_FAILURE;
    	}
    }*/
    
    if (ArrangeMolecules == "yes") {
        Arranged = false; 
        while (!Arranged) {
            Arranged = Sys.arrangeMoleculesFCC();
            try {
		        Sys.updateVerletLists();
		        Sys.calculateForces(true);
		    }
		    catch (const LibraryException &ex) {
		        Arranged = false; 
		        std::cout << "bad placement, trying again..." << std::endl; 
		    }
	    }
    }
    
    Density = double(Sys.NumberOfParticles())/(Lx*Ly*Lz);
    std::cout << "Number of chains: " << Sys.NumberOfMolecules() << ", Number of monomers: " << Sys.NumberOfParticles() << std::endl;
    std::cout << "Number density of the system: " << Density << std::endl;
    
    std::ofstream StatisticsStream(StatisticsFile, std::ios::out | std::ios::app);
	FILE* PDBout{};
    
    PDBout = fopen("initial.pdb", "w");
    Sys.printPDB(PDBout, 0, 1);
	fclose(PDBout);
	Sys.printStatistics(StatisticsStream, Time);
    
    ////////////////////////////////////
    ////// TIME MEASUREMENT START //////
    
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
    ///////////////////////////////////
    ////// MAIN SIMULATION LOOP ///////
    
    n = 1; 
    while (!BoxSizeReached) {
        try {   
            if (n % SamplingStep == 0) Sys.propagateLangevin(true);
            else Sys.propagateLangevin(false);
        }
        catch (const LibraryException &ex) {
            PDBout = fopen((ConfigOutFile+std::to_string(n)+".pdb").c_str(), "w");
            Sys.printPDB(PDBout, n, 1);
            fclose(PDBout);    
            std::cout << "Terminating..." << std::endl;   
            std::terminate(); 
        }
        if (n % SamplingStep == 0) {
            Sys.printStatistics(StatisticsStream, Time);  
        }
        if (n % ScalingStep == 0) {
            PDBout = fopen((ConfigOutFile+"-L"+std::to_string(Sys.BoxSize[0])+".pdb").c_str(), "w");
            Sys.printPDB(PDBout, n, 1);
            fclose(PDBout);
            Sys.scaleSystem(ScalingX, ScalingY, ScalingZ); 
            
            std::cout << "BoxSize currently: " << Sys.BoxSize[0] << " , " << Sys.BoxSize[1] << " , " << Sys.BoxSize[2] << std::endl; 
            std::cout << "Cell side lengths: " << Sys.CellSideLength[0] << " " <<  Sys.CellSideLength[1] << " " << Sys.CellSideLength[2] << std::endl;
            if (Sys.BoxSize[0] >= Lx) {
            //if (abs(Sys.BoxSize[0] - Lx) < 0.1 && abs(Sys.BoxSize[1] - Ly) < 0.1 && abs(Sys.BoxSize[2] - Lz) < 0.1) {
                BoxSizeReached = true;
                n = 0;  
                Sys.wrapMoleculesCOM(); 
                std::cout << "BoxSize reached" << std::endl; 
            }    
        }
        
        if (SignalCaught) {
            std::cout << "writing job data..." << std::endl; 
            PDBout = fopen((ConfigOutFile+std::to_string(n)+".pdb").c_str(), "w");
            Sys.printPDB(PDBout, n, 1);
            fclose(PDBout);
            std::cout << "Terminating..." << std::endl;   
            std::terminate(); 
        }
        
        Time += MDStep; 
        n++; 
    }
    
    std::cout << "End of simulation. Writing job data..." << std::endl; 
    Sys.printStatistics(StatisticsStream, Time);
    PDBout = fopen((ConfigOutFile+std::to_string(n)+".pdb").c_str(), "w");
    Sys.printPDB(PDBout, n, 1);
    fclose(PDBout);
    StatisticsStream.close(); 
    
    //////////////////////////////////////
    ////// MAIN SIMULATION LOOP END///////
    
    gettimeofday(&end, NULL); 
    
    double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    unsigned TotalMonomers {Sys.NumberOfParticles()};
    std::cout << "total time: " << realTime << " , time per particle and step: " << realTime/TotalMonomers/n << std::endl;
    
    return EXIT_SUCCESS;
   
}
