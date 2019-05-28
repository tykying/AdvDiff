include makefile.inc

default: advdiff

.PHONY: objs test-bin test clean clean-test
.IGNORE: test-bin clean-test

OBJS = $(filter-out src/advdiff.o, $(addsuffix .o, $(basename $(wildcard src/*.F90))))

advdiff: objs
	$(FC) $(OBJS) src/advdiff.o -o $@ $(LDFLAGS)

objs:
	make -C src

test-bin: advdiff
	make $@ -C tests 

test: clean-test advdiff
	make $@ -C tests

clean: clean-test
	make $@ -C src
	rm -fv advdiff

clean-test:
	make clean -C tests
