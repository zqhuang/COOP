SHELL = /bin/sh

TYPEDEF = ../typedef/

include ../configure.in
include ../compile_rules.in

LIB = libcooputils.a
OBJ = cwrapper.o  MPIutils.o sort.o list.o random.o specfunc.o matrix.o interpolation.o integrate.o ode.o file_io.o asy.o  wrapper.o
#jlbrute.o sphericalbessels.o threejs.o jllpm.o  fft.o wraputils.o solvepower.o zqwrap.o

default: $(LIB)

$(LIB): $(OBJ) 
	ar -r $@ $?

Test: $(OBJ) test.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test
