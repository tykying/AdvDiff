#!/bin/bash          
echo Running AdvDiff

# Compile Particle Advection Code
#( cd ./lib/lagrangian_particles; make clean; make)
#ln -sf /home/s1046972/opt/FORTRAN_LIB/fftw/lib/libfftw3.a ./lib/libfftw3.a

export LD_LIBRARY_PATH=$PWD/lib

echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH

#sleep 5

make clean
make

date
#nohup ./advdiff > mon_advdiff64.txt &
./advdiff

date 

#make clean

#python ./python/animate_field.py

