import sys
import glob 
import os 
import random
import shutil
import subprocess 
from collections import namedtuple 


random.seed()
RunsPerParameter = 16
Parameter = namedtuple("Parameters", "Box, Temperature, Gamma, BondPotential, FunctionalFraction, MDStep, SimTime, EquilTime")

BondPotential = namedtuple("BondPotential", ["ConstantK", "ConstantR0"])
Box = namedtuple("Box", ["Lx", "Ly", "Lz"])

ParameterSets = []

ParameterSets.append(Parameter(Box(200, 200, 200), 1.0, 0.05, [BondPotential(24.6, 1.38), BondPotential(29.6, 1.448), BondPotential(33.7, 1.477)], 0.1, 0.005, 5000000.0, 10000.))

ParameterSets.append(Parameter(Box(200, 200, 200), 1.0, 0.05, [BondPotential(24.6, 1.38)], 0.3, 0.005, 5000000.0, 10000.))

homepath="/home/formanek/REVERSIBLE/"
clusterpath="/scratch-new/formanek/REVERSIBLE/"


submit_files =open("runs_to_submit.dat", "w")
if not os.path.exists(clusterpath): 
    os.makedirs(clusterpath)

    
for paramSet in ParameterSets: 
    Lx = paramSet.Box.Lx
    Ly = paramSet.Box.Ly
    Lz = paramSet.Box.Lz
    Temperature = paramSet.Temperature
    Gamma = paramSet.Gamma
    MDStep = paramSet.MDStep
    SimTime = paramSet.SimTime
    EquilTime = paramSet.EquilTime
    FunctionalFraction = paramSet.FunctionalFraction
    runparentdir = clusterpath+"runs-f-"+str(FunctionalFraction)+"-rho-0.0/"
    if not os.path.exists(runparentdir):
        os.makedirs(runparentdir)
    else: 
        printf("directory "+runparentdir+"already exists, exiting")
        sys.exit()
    for Potential in paramSet.BondPotential: 
        ConstantK = Potential.ConstantK
        ConstantR0 = Potential.ConstantR0
        rundir=runparentdir+"K-"+str(ConstantK)+"-R0-"+str(ConstantR0)+"/"
        if not os.path.exists(rundir):
            os.makedirs(rundir)  
        for i in range(RunsPerParameter): 
            Seed = random.randint(0,999999)
            run_name = "REPL-"+str(i)
            directory = rundir+run_name
            if os.path.exists(directory):
                print("path "+directory+" already exists!") 
                break 
            os.makedirs(directory) 
            os.makedirs(directory+"/configs")
            os.makedirs(directory+"/Bonds")

            shutil.copyfile(homepath+"input/SCNPs/Reversible-config", directory+"/config")

            shutil.copyfile(homepath+"input/SCNPs/Reversible-chain", directory+"/chain")
            command = "./set_functional_groups "+directory+"/functional "+directory+"/chain "+str(FunctionalFraction)+" "+str(random.randrange(1,9999999))
            print(command) 
            subprocess.call(command, shell=True)
            #### change later
            shutil.copyfile(homepath+"input/time", directory+"/steps")
            ####
            submit_files.write(directory+"\n") 
            with open("parameter_template_reversible.dat", "r") as inF:
                lines = inF.readlines()
            with open(directory+"/parameters.dat", "w") as outF:
                outF.write("BoxX = "+str(Lx)+"\n")
                outF.write("BoxY = "+str(Ly)+"\n")
                outF.write("BoxZ = "+str(Lz)+"\n\n")   
                outF.write("Temperature = "+str(Temperature)+"\n")
                outF.write("Gamma = "+str(Gamma)+"\n\n") 
                outF.write("ConstantK = "+str(ConstantK)+"\n")
                outF.write("ConstantR0 = "+str(ConstantR0)+"\n\n")	
                outF.write("MDStep = "+str(MDStep)+"\n\n") 
                outF.write("StartTime = "+str(0.0)+"\n")
                outF.write("SimTime = "+str(SimTime)+"\n")
                outF.write("EquilTime = "+str(EquilTime)+"\n\n")
                outF.write("Seed = "+str(Seed)+"\n\n")
                
                lineInd = 0
                for currLine in lines: 
                    if lineInd >= 18 : 
                        outF.write(currLine)
                    lineInd += 1

