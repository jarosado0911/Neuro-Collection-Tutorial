#!/bin/sh 
#PBS -l walltime=02:00:00 
#PBS -N sampleRyR-SERCAjob0
#PBS -q normal 
#PBS -l nodes=1:ppn=28 
#PBS -e error.txt 
#PBS -o output.txt 
cd $PBS_O_WORKDIR

mpirun -np 28 ugshell -ex /home/tuh36651/work/Example-Spine-Dend/reconstructed_spine_wER.lua -grid /home/tuh36651/work/Example-Spine-Dend/Spine5_wER.ugx -numRefs 0 -caInflux 0.0086 -tstep 2.5e-6 -endTime 0.015 -setting ryrserca -ryrDensity 0.009 -minDef 1e-12 -outName output -solver GS -vtk -pstep 1e-4