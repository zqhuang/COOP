SHELL = /bin/sh

TYPEDEF = ../typedef/

include ../configure.in
include ../compile_rules.in

LIB = libcooputils.a
OBJ = MPIutils.o evalstr.o  specfunc.o list.o random.o cholesky.o matrix.o interpolation.o integrate.o ode.o file_io.o ndprob.o asy.o  fft.o jlbrute.o  hankel.o nn.o inifile.o gaussian_peak.o sphere.o coop_fitsio.o wrapper.o parabolaU.o

#threejs.o 

default: $(LIB) MPLOT

$(LIB): $(OBJ) 
	ar -r $@ $?

Test: $(OBJ) test.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

Test0: $(OBJ) test0.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test0.o $(LINK)

Test1: $(OBJ) test1.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test1.o $(LINK)

Test2: $(OBJ) test2.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test2.o $(LINK)

Test3: $(OBJ) test3.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test3.o $(LINK)

Test4: $(OBJ) test4.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test4.o $(LINK)

Test5: $(OBJ) test5.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test5.o $(LINK)


Szeq: 	$(OBJ) search_zeq.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) search_zeq.o $(LINK)


LSQ: $(OBJ) leastsq.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) leastsq.o $(LINK)

IDP: $(OBJ) IDPossible.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) IDPossible.o $(LINK)

MPLOT: $(OBJ) makeplot.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) makeplot.o $(LINK)

SMOOTH:	$(OBJ) smooth_data.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) smooth_data.o $(LINK)

EBin: $(OBJ) encbin.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) encbin.o $(LINK)

Enc: $(OBJ) encrypt.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) encrypt.o $(LINK)

Dec: $(OBJ) decrypt.o $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) decrypt.o $(LINK)


clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test MPLOT Dec Enc SMOOTH EBin 


test.o: $(LIB)

rungas.o: gasmodel.o

ode.o:	../include/dverk.h

integrate.o:	../include/qromb.h

$(LIB): wrapper.o

wrapper.o:  MPIutils.o specfunc.o list.o random.o cholesky.o matrix.o interpolation.o integrate.o ode.o file_io.o ndprob.o asy.o  fft.o jlbrute.o  inifile.o gaussian_peak.o evalstr.o


