FC = ifort
FC = gfortran

FCFLAGS += -I$(PWD)/lib
#FCFLAGS += -Ofast -xHost -ipo -finline-functions
#FCFLAGS += -openmp


# Self-teaching notes
#-L<*> = Search in directory *
#-l<*> = Use library lib* (fftpack is linked manually)
#-I<*> = Include files in directory *

#LDFLAGS += -llapack -L. -lfftpack -L$(PWD)/lagrangian_particles -lparticles
LDFLAGS += -llapack -lblas -lfftw3 -L$(PWD)/lib

# REMARK:
# I put only symbolic links in the folder './lib' only
