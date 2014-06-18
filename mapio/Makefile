SHELL = /bin/sh
TYPEDEF=../typedef
UTILS=../utils
include ../configure.in
include ../compile_rules.in

LIB = libfits.a
OBJ = sphere.o readfits.o fitswrap.o healgeo.o

default: $(LIB)

$(LIB): $(OBJ) 
	ar -r $@ $?

Test: $(OBJ) test.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

GetSpots: $(OBJ) exportspots.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) exportspots.o $(LINK)


Stack: $(OBJ) dostack.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) dostack.o $(LINK)


clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test


test.o: $(LIB)
