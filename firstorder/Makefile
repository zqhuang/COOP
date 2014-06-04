SHELL = /bin/sh

TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background

include ../configure.in
include ../compile_rules.in

BGLIB = libcoopfirstorder.a
BGOBJ = firstorder.o wrapper.o

default: $(BGLIB)

$(BGLIB): $(BGOBJ)
	ar -r $@ $?

Test: $(BGOBJ) test.o $(BACKGROUND)/libcoopbackground.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(BGOBJ) test.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test
