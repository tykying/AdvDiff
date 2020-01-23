include ../makefile.inc

OBJS = $(addsuffix .o, $(basename $(wildcard *.F90)))

default: $(OBJS)

.PHONY: clean

%.o: %.F90
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c -o $@ $<

advdiff.o: advdiff_field.o advdiff_precision.o advdiff_timestep.o \
  advdiff_debug.o advdiff_timing.o advdiff_io.o advdiff_trajdata.o \
  advdiff_inference.o 

advdiff_field.o: advdiff_precision.o \
  advdiff_debug.o

advdiff_timestep.o: advdiff_field.o advdiff_precision.o advdiff_timing.o \
  advdiff_debug.o

advdiff_timing.o: advdiff_precision.o

advdiff_io.o: advdiff_precision.o advdiff_debug.o advdiff_field.o \
  advdiff_timing.o advdiff_trajdata.o advdiff_inference.o 

advdiff_trajdata.o: advdiff_precision.o advdiff_debug.o advdiff_timing.o \
  advdiff_field.o advdiff_complib.o

advdiff_complib.o: advdiff_precision.o

advdiff_debug.o: advdiff_precision.o

advdiff_inference.o: advdiff_field.o advdiff_precision.o advdiff_timestep.o \
  advdiff_debug.o advdiff_timing.o advdiff_trajdata.o 

advdiff_unittest.o: advdiff_field.o advdiff_precision.o advdiff_timestep.o \
  advdiff_debug.o advdiff_timing.o advdiff_trajdata.o advdiff_inference.o \
  advdiff_io.o

clean:
	rm -fv $(OBJS)
	rm -fv $(addsuffix .mod, $(basename $(OBJS)))
