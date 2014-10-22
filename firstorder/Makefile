SHELL = /bin/sh

TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background

include ../configure.in
include ../compile_rules.in

LIB = libcoopfirstorder.a
OBJ = pertobj.o firstorder.o lensing.o Cls_post_process.o wrapper.o

default: $(LIB)

$(LIB): $(OBJ)
	ar -r $@ $?

FNL: 	$(OBJ) testfnl.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) testfnl.o $(LINK)


DOCLS: $(OBJ) testcls.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) testcls.o $(LINK)

Test: $(OBJ) test.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

Z3D: $(OBJ) zeta.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) zeta.o $(LINK)

testZ: $(OBJ) testzeta.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) testzeta.o $(LINK)


clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test


firstorder.o:	pertobj.o  ../include/firstorder_equations.h ../include/firstorder_source.h ../include/recfast_source.h ../include/recfast_head.h
wrapper.o: 	firstorder.o lensing.o Cls_post_process.o
test.o:		wrapper.o
Cls_post_process.o: firstorder.o lensing.o
zeta.o:	wrapper.o
testzeta.o:	wrapper.o
