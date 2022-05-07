SHELL = /bin/sh

TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background
FIRSTORDER=../firstorder
NONLINEAR=../nonlinear

include ../configure.in
include ../compile_rules.in

LIB = libcoopmcmc.a
OBJ = gobs.o SNlike_JLA.o bao.o sltd.o HSTlike.o chronometer.o wl.o cliklike.o mcmc.o fisher.o

default: $(LIB) MCMC FISHER PLOTF

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

FISHER:	runfisher.o $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) runfisher.o $(LINK)

PLOTF:	plotfisher.o $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) plotfisher.o $(LINK)


STD6:	teststd6.o $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) teststd6.o $(LINK) 


clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~  FISHER MCMC Test PRINTCLIK PLOTF

runmcmc.o: $(OBJ)
mcmc.o: cliklike.o bao.o wl.o SNlike_JLA.o HSTlike.o
runfisher: $(OBJ)
