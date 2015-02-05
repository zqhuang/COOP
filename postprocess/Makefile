SHELL=/bin/sh
COOPLIB=../lib
include ../configure.in
include ../compile_rules.in

LIB = libcooppostprocess.a
OBJ = latex.o statchains.o 

default: $(LIB) GD

$(LIB): $(OBJ) 
	ar -r $@ $?


Test: 	$(OBJ) test.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

GD: 	$(OBJ) getdist.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) getdist.o $(LINK)


clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test GD
