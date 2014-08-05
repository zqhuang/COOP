SHELL = /bin/sh

TYPEDEF = ../typedef/

include ../configure.in
include ../compile_rules.in

LIB = libcooputils.a
OBJ = MPIutils.o specfunc.o list.o random.o matrix.o interpolation.o integrate.o ode.o file_io.o asy.o  fft.o jlbrute.o  solvepower.o zqwrap.o inifile.o  wrapper.o 
#jllpm.o  threejs.o

default: $(LIB)

$(LIB): $(OBJ) 
	ar -r $@ $?

Test: $(OBJ) test.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

Enc: $(OBJ) encrypt.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) encrypt.o $(LINK)

Dec: $(OBJ) decrypt.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) decrypt.o $(LINK)


clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test


test.o: $(LIB)
