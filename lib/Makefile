SHELL = /bin/sh

TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background
FIRSTORDER=../firstorder
NONLINEAR=../nonlinear
FORECAST=../forecast

include ../configure.in
include ../compile_rules.in
LIB = libcoopwrapper.a
OBJ = cambpp.o wrapper.o 


default: $(LIB)

$(LIB):	$(OBJ) 
	ar -r $@ $?

Test: $(OBJ) test.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ 

