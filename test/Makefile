SHELL = /bin/sh
COOPLIB=../lib

include ../configure.in
include ../compile_rules.in

Ex1: example1.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ example1.o $(LINK)

Ex2: example2.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ example2.o $(LINK)

Ex3: example3.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ example3.o $(LINK)

Test:	test.o 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ test.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test Ex1 Ex2 Ex3
