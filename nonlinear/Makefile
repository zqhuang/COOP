SHELL = /bin/sh

TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background
FIRSTORDER=../firstorder

include ../configure.in
include ../compile_rules.in

LIB = libcoopnonlinear.a
OBJ = nonlinear.o 

default: $(LIB)

$(LIB):	$(OBJ) 
	ar -r $@ $?

FNL: $(OBJ) testfnl.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) testfnl.o $(LINK)

Test: $(OBJ) test.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ 


test.o: $(OBJ)

fnl.o: $(OBJ)
