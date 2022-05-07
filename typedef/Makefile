SHELL = /bin/sh
include ../configure.in
include ../compile_rules.in
TYPEDEFLIB = libcooptypedef.a
TYPEDEFOBJ = page.o cwrapper.o constants.o svd.o basicutils.o sort.o sortrev.o string.o arguments.o function.o  particle.o species.o cosmology.o wrapper.o

default: $(TYPEDEFLIB)

$(TYPEDEFLIB): $(TYPEDEFOBJ) 
	ar -r $@ $?

Test: $(TYPEDEFOBJ) test.o 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(TYPEDEFOBJ) test.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test


test.o		:wrapper.o
wrapper.o 	:cosmology.o 
cosmology.o	:species.o page.o
species.o	:particle.o
particle.o	:function.o 
function.o	:arguments.o
arguments.o	:string.o
string.o	:basicutils.o
basicutils.o	:svd.o
svd.o		:constants.o
constants.o	: ../include/constants.h
