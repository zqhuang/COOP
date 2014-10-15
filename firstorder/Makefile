SHELL = /bin/sh

TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background

include ../configure.in
include ../compile_rules.in

OFLIB = libcoopfirstorder.a
OFOBJ = pertobj.o firstorder.o lensing.o Cls_post_process.o wrapper.o

default: $(OFLIB)

$(OFLIB): $(OFOBJ)
	ar -r $@ $?

DOCLS: $(OFOBJ) testcls.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OFOBJ) testcls.o $(LINK)

Test: $(OFOBJ) test.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OFOBJ) test.o $(LINK)

Z3D: $(OFOBJ) zeta.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OFOBJ) zeta.o $(LINK)

testZ: $(OFOBJ) testzeta.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OFOBJ) testzeta.o $(LINK)


clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test


firstorder.o:	pertobj.o  ../include/firstorder_equations.h ../include/firstorder_source.h ../include/recfast_source.h ../include/recfast_head.h
wrapper.o: 	firstorder.o lensing.o Cls_post_process.o
test.o:		wrapper.o
Cls_post_process.o: firstorder.o lensing.o
zeta.o:	wrapper.o
testzeta.o:	wrapper.o
