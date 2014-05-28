SHELL = /bin/sh
include ../configure.in
include ../compile_rules.in
TYPEDEFLIB = libcooptypedef.a
TYPEDEFOBJ = constants.o svd.o basicutils.o string.o arguments.o function.o  species.o cosmology.o wrapper.o

default: $(TYPEDEFLIB)

$(TYPEDEFLIB): $(TYPEDEFOBJ) ../headfile/constants.h 
	ar -r $@ $?

Test: $(TYPEDEFOBJ) test.o ../headfile/constants.h
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(TYPEDEFOBJ) test.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test


wrapper.o 	:typedef.o
typedef.o	:basictype.o
basictype.o	:basicutils.o
basicutils.o	:svd.o
svd.o		:constants.o
