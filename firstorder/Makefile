SHELL = /bin/sh

BACKGROUND=../background

include ../configure.in
include ../compile_rules.in

BGLIB = libcoopfirstorder.a
BGOBJ = firstorder.o wrapper.o

default: $(BGLIB)

$(BGLIB): $(BGOBJ) 
	ar -r $@ $?

Test: $(BGOBJ) test.o ../background/libcoopbackground.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(BGOBJ) test.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test
