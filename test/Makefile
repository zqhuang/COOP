SHELL = /bin/sh
COOPLIB=../lib

include ../configure.in
include ../compile_rules.in


Test: test.o 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ test.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test
