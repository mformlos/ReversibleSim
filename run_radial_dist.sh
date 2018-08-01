#!/bin/bash 

for i in $( seq 0 15 ); do 
    ./bin/Radial_dist /scratch-new/formanek/REVERSIBLE/runs-f-0.1-rho-0.0/K-29.6-R0-1.448/REPL-$i 200 0.1 120000 800160000 40000 0.01
done
