SHELL = /bin/sh

TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background
FIRSTORDER=../firstorder
NONLINEAR=../nonlinear

include ../configure.in
include ../compile_rules.in

LIB = libcoopforecast.a
OBJ = SNlike_JLA.o bao.o HSTlike.o wl.o cliklike.o forecast.o 

default: $(LIB) MCMC

$(LIB):	$(OBJ) 
	ar -r $@ $?

Test:	test.o $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK) 


PRINTCLIK:	print_clik.o $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) print_clik.o $(LINK) 


DOSN:	testSN.o $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) testSN.o $(LINK) 


MCMC:	runmcmc.o $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) runmcmc.o $(LINK)

STD6:	teststd6.o $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) teststd6.o $(LINK) 


clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~  MCMC Test STD6 DOSN PRINTCLIK

runmcmc.o: $(OBJ)
forecast.o: cliklike.o bao.o wl.o SNlike_JLA.o HSTlike.o
