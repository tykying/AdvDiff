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
for Phase in {1..1}
do
  Seed_ID=1
  
  case $HOSTNAME in
    *"maths.ed.ac.uk")
      output_prefix=".";;
    *"ecdf.ed.ac.uk")
      output_prefix="/exports/eddie/scratch/s1046972"
  esac


  # QG
  output_fld=$output_prefix"/output/N4096_D256_I256/QGM2_L"$layer"_NPART"$NPart
  echo $output_fld
  
  ## TTG
  #output_fld=$output_prefix"/output/N4096_D64_I64/TTG_sinusoidal"
  #Td=30
  
  output_dir=$output_fld"/h"$Td"d/Seed"$Seed_ID
  mkdir -p $output_dir"/SpinUp"
  mkdir -p $output_dir"/Tuned"
  mkdir -p $output_dir"/Data"
  mkdir -p $output_dir"/screen"
  #./advdiff $Td $layer $NPart $Phase $Seed_ID $output_dir
  #nohup mpirun -np 16 ./advdiff $Td $layer $NPart $Phase $Seed_ID $output_dir  > $output_dir"/screen/Phase"$Phase".txt"
  mpirun -np 32 ./advdiff $Td $layer $NPart $Phase $Seed_ID $output_dir
  #mpirun -np 10 ./advdiff $Td $layer $NPart $Phase $Seed_ID $output_dir
  echo Finish!
  #./advdiff $Td $layer $NPart $Phase $Seed_ID
done

elapsedseconds=$SECONDS
echo Time taken = $elapsedseconds
date 

#make clean

#python ./python/animate_field.py

# Memory debugger
#valgrind --track-origins=yes --leak-check=full  --show-reachable=yes ./advdiff $Td $layer $NPart $Phase $Seed_ID $output_dir 
