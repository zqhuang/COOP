SHELL = /bin/sh
include ../configure.in
include ../compile_rules.in
TYPEDEFLIB = libcooptypedef.a
TYPEDEFOBJ = constants.o basictype.o basicutils.o typedef.o utils.o

default: $(TYPEDEFLIB)

$(TYPEDEFLIB): $(TYPEDEFOBJ) ../headfile/constants.h 
	ar -r $@ $?

Test: $(TYPEDEFOBJ) test.o ../headfile/constants.h
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(TYPEDEFOBJ) test.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test
