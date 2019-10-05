#!/bin/bash          
echo Running AdvDiff

# Compile Particle Advection Code
#( cd ./lib/lagrangian_particles; make clean; make)
#ln -sf /home/s1046972/opt/FORTRAN_LIB/fftw/lib/libfftw3.a ./lib/libfftw3.a
#ln -sf /home/s1046972/opt/FORTRAN_LIB/fftw/include/fftw3.f03 ./lib/fftw3.f03

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/intel/2013/composer_xe_2013.0.079/compiler/lib/intel64

#echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH

#sleep 5

make clean
make

date
SECONDS=0
#nohup ./advdiff > mon_advdiff_bilinear_h32d_prior2_L2.txt &
#nohup ./advdiff > mon_advdiff_bilinear_h32d_canon_lmt.txt &
#nohup ./advdiff > mon_advdiff_TTG.txt &
#nohup ./advdiff > mon_advdiff_h64d_bilinear.txt &
#nohup ./advdiff > mon_advdiff_h32d_unstructured2.txt &
#nohup ./advdiff > mon_advdiff_h32d_smallK.txt &
#nohup ./advdiff > mon_advdiff_h32d_highres80.txt &
#nohup ./advdiff > mon_advdiff_h64d_L2_sigma_64K.txt &
#nohup ./advdiff > mon_advdiff_h32d_L1_sigma.txt &
#nohup ./advdiff > LW_convtest2.txt &
#nohup ./advdiff > MC_convtest2.txt &
#./advdiff
#mpirun -np 32 ./advdiff
mpirun -np 4 ./advdiff
elapsedseconds=$SECONDS
echo Time taken = $elapsedseconds
date 

#make clean

#python ./python/animate_field.py

# Memory debugger
#valgrind --track-origins=yes --leak-check=full  --show-reachable=yes ./advdiff


