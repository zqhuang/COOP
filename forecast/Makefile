SHELL = /bin/sh

TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background
FIRSTORDER=../firstorder
NONLINEAR=../nonlinear

include ../configure.in
include ../compile_rules.in

LIB = libcoopforecast.a
OBJ = forecast.o 


default: $(LIB)

$(LIB):	$(OBJ) 
	ar -r $@ $?

Test:	test.o 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)



clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ 

