#!/bin/bash          
echo Running AdvDiff

# Compile Particle Advection Code
#( cd ./lib/lagrangian_particles; make clean; make)
#ln -sf /home/s1046972/opt/FORTRAN_LIB/fftpack/fftpack5.1/lib/libfftpack.a ./lib/libfftpack.a

export LD_LIBRARY_PATH=$PWD/lib

echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH

#sleep 5

make clean
make

#nohup ./qg > ./qg_scr_$COUNTER.txt &
#nohup mpirun -n 1 ./advdiff > mon_advdiff.txt &
#nohup ./advdiff > mon_advdiff32.txt &
#./advdiff
#mpirun -n 48  ./advdiff
#mpirun -n 32  ./advdiff


#nohup mpirun -n 1  ./advdiff > temporal_error_LaxWendoff_dx.txt &
#nohup mpirun -n 1  ./advdiff > ./unittest/convergence/log/temporal_error_LW2D2nd_dx.txt &
#nohup mpirun -n 1  ./advdiff > ./unittest/convergence/log/temporal_error_MC2D2nd_dx.txt &
#nohup mpirun -n 1  ./advdiff > ./unittest/convergence/log/temporal_error_LW2D_dx.txt &
#nohup mpirun -n 1  ./advdiff > ./unittest/convergence/log/temporal_error_MC2D_dx.txt &
#nohup mpirun -n 1  ./advdiff > ./unittest/convergence/log/temporal_error_MM2D_dx.txt &
#nohup mpirun -n 1  ./advdiff > ./unittest/convergence/log/temporal_error_SB2D_dx.txt &
#nohup mpirun -n 1  ./advdiff > temporal_error_LaxWendoff_dxsq.txt &
#nohup mpirun -n 1  ./advdiff > temporal_error_MC_dx.txt &
#nohup mpirun -n 1  ./advdiff > temporal_error_2DMC_dx.txt &
#nohup mpirun -n 1  ./advdiff > temporal_error_MM_dx.txt &
#nohup mpirun -n 1  ./advdiff > temporal_error_2DMM_dx.txt &
#nohup mpirun -n 1  ./advdiff > temporal_error_Sweby_dx.txt &
#nohup mpirun -n 1  ./advdiff > temporal_error_2DSweby_dx.txt &
#nohup mpirun -n 1  ./advdiff > temporal_error_MC_dxsq.txt &
#nohup mpirun -n 1  ./advdiff > spatial_error_MC_dx.txt &
#nohup mpirun -n 1  ./advdiff > spatial_error_LaxWendoff_dx.txt &
mpirun -n 1  ./advdiff 


 

#make clean

#python ./python/animate_field.py

