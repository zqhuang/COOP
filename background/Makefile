SHELL = /bin/sh

TYPEDEF = ../typedef/

include ../configure.in
include ../compile_rules.in

BGLIB = libcoopbackground.a
BGOBJ = background.o wrapper.o

default: $(BGLIB)

$(BGLIB): $(BGOBJ) 
	ar -r $@ $?

Test: $(BGOBJ) test.o ../include/background.h $(TYPEDEF)/libcooptypedef.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(BGOBJ) test.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test