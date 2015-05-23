SHELL = /bin/sh

TYPEDEF = ../typedef/

include ../configure.in
include ../compile_rules.in

LIB = libcooputils.a
OBJ = MPIutils.o specfunc.o list.o random.o cholesky.o matrix.o interpolation.o integrate.o ode.o file_io.o ndprob.o asy.o  fft.o jlbrute.o  inifile.o gaussian_peak.o  wrapper.o 
#jllpm.o  threejs.o

default: $(LIB)

$(LIB): $(OBJ) 
	ar -r $@ $?

Test: $(OBJ) test.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

MPLOT: $(OBJ) makeplot.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) makeplot.o $(LINK)


Enc: $(OBJ) encrypt.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) encrypt.o $(LINK)

Dec: $(OBJ) decrypt.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) decrypt.o $(LINK)


clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test MPLOT Dec Enc


test.o: $(LIB)

ode.o:	../include/dverk.h

integrate.o:	../include/qromb.h

$(LIB): wrapper.o

wrapper.o:  MPIutils.o specfunc.o list.o random.o cholesky.o matrix.o interpolation.o integrate.o ode.o file_io.o ndprob.o asy.o  fft.o jlbrute.o  inifile.o gaussian_peak.o 
