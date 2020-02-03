#FC = ifort
FC = mpiifort

#FCFLAGS += -Ofast -xHost -ipo -finline-functions -real-size 64 -extend-source 132
#FCFLAGS += -openmp

#FC = gfortran
FC = mpif90
FCFLAGS += -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -std=f2003 -Wall
FCFLAGS += -O3 -ffast-math -march=native -DNDEBUG
#FCFLAGS += -fopenmp

#FCFLAGS += -g

LDFLAGS += -llapack -lblas

