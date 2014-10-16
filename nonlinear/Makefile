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

Test: $(LIB) test.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(LIB) test.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ 

