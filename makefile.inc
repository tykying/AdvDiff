FC = ifort

FCFLAGS += -I$(PWD)/lib
FCFLAGS += -Ofast -xHost -ipo -finline-functions
FCFLAGS += -openmp

#FC = gfortran
#FCFLAGS += -I$(PWD)/lib
#FCFLAGS += -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -std=f2003 -Wall
#FCFLAGS += -O3 -ffast-math -march=native -DNDEBUG
#FCFLAGS += -g



# Self-teaching notes
#-L<*> = Search in directory *
#-l<*> = Use library lib* (fftpack is linked manually)
#-I<*> = Include files in directory *

#LDFLAGS += -llapack -L. -lfftpack -L$(PWD)/lagrangian_particles -lparticles
LDFLAGS += -llapack -lblas -lfftw3 -L$(PWD)/lib

# REMARK:
# I put only symbolic links in the folder './lib' only
