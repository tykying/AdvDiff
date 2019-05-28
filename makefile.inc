FC = mpif90
#FC = gfortran
#FC = mpifort
#FC = ifort
#FCFLAGS += -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -std=f2003 -Wall -I$(PWD)/lagrangian_particles
#FCFLAGS += -O0 -fbounds-check -g -finit-real=nan -finit-integer=-1000 -DDDEBUG
FCFLAGS += -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -std=f2003 -Wall -I$(PWD)/lib
FCFLAGS += -O3 -ffast-math -march=native -DNDEBUG
#FCFLAGS += -fopenmp
#FCFLAGS += -O3 -ffast-math -march=native

# Self-teaching notes
#-L<*> = Search in directory *
#-l<*> = Use library lib* (fftpack is linked manually)
#LDFLAGS += -llapack -L. -lfftpack -L$(PWD)/lagrangian_particles -lparticles
LDFLAGS += -llapack -lblas -L$(PWD)/lib

# REMARK:
# I put only symbolic links in the folder './lib' only
