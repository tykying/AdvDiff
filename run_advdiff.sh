#!/bin/bash          
echo Running AdvDiff

# Compile Particle Advection Code
#( cd ./lib/lagrangian_particles; make clean; make)
#ln -sf /home/s1046972/opt/FORTRAN_LIB/fftw/lib/libfftw3.a ./lib/libfftw3.a
#ln -sf /home/s1046972/opt/FORTRAN_LIB/fftw/include/fftw3.f03 ./lib/fftw3.f03

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib
export LD_RUN_PATH=$LD_RUN_PATH:$PWD/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/intel/2013/composer_xe_2013.0.079/compiler/lib/intel64
#echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH

#sleep 5

make clean
make

date
SECONDS=0

Td=32
layer=2
NPart=676
for Phase in {2..3}
do
  Seed_ID=1
  
  # QG
  output_fld="./output/N4096_D256_I256/QGM2_L"$layer"_NPART"$NPart""
  
  # TTG
  output_fld="./output/N4096_D64_I64/TTG_sinusoidal"
  Td=30
  
  output_dir=$output_fld"/h"$Td"d/Seed"$Seed_ID
  mkdir -p $output_dir"/SpinUp"
  mkdir -p $output_dir"/Tuned"
  mkdir -p $output_dir"/Data"
  mkdir -p $output_dir"/screen"
  nohup mpirun -np 16 ./advdiff $Td $layer $NPart $Phase $Seed_ID > $output_dir"/screen/Phase"$Phase".txt"
  #mpirun -np 16 ./advdiff $Td $layer $NPart $Phase $Seed_ID
  #./advdiff $Td $layer $NPart $Phase $Seed_ID
done

elapsedseconds=$SECONDS
echo Time taken = $elapsedseconds
date 

#make clean

#python ./python/animate_field.py

# Memory debugger
#valgrind --track-origins=yes --leak-check=full  --show-reachable=yes ./advdiff


