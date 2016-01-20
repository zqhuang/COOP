SHELL = /bin/sh

TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background
FIRSTORDER=../firstorder
NONLINEAR=../nonlinear
FORECAST=../forecast

include ../configure.in
include ../compile_rules.in

OBJ = lattice.o 
$(LIB): libcooplattice.a

default: HLattice

$(LIB):	$(OBJ) 
	ar -r $@ $?

Test: $(OBJ) test.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

HLattice: $(OBJ) runlattice.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) runlattice.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~  HLattice

runlattice.o: lattice.o

