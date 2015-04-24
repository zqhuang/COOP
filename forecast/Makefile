SHELL = /bin/sh

TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background
FIRSTORDER=../firstorder
NONLINEAR=../nonlinear

include ../configure.in
include ../compile_rules.in

LIB = libcoopforecast.a
OBJ = SNlike_JLA.o HSTlike.o cliklike.o forecast.o 

default: $(LIB) DOCLIK

$(LIB):	$(OBJ) 
	ar -r $@ $?

Test:	test.o $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK) 


DOSN:	testSN.o $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) testSN.o $(LINK) 


DOCLIK:	testclik.o $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) testclik.o $(LINK) 

STD6:	teststd6.o $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) teststd6.o $(LINK) 


clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~  DOCLIK STD6 DOSN Test

test.o: $(OBJ)
testSN.o: $(OBJ)
testclik.o: $(OBJ)
teststd6.o: $(OBJ)
