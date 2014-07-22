SHELL = /bin/sh
COOPLIB=../lib

include ../configure.in
include ../compile_rules.in

Example1: example1.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ example1.o $(LINK)

Example2: example2.o
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ example2.o $(LINK)


Test:	test.o 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ test.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test Example1 Example2 *.txt *.pdf
