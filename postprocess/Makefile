SHELL = /bin/sh
TYPEDEF=../typedef
UTILS=../utils
include ../configure.in
include ../compile_rules.in

LIB = libcooppostprocess.a
OBJ = latex.o statchains.o postprocess.o 

default: $(LIB)

$(LIB): $(OBJ) 
	ar -r $@ $?


Test: $(OBJ) test.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

GetDist: $(OBJ) getdist.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) getdist.o $(LINK)


clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test GetDist
