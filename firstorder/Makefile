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


SolvePert: $(OBJ) solve_pert.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) solve_pert.o $(LINK)

CalcCls: $(OBJ) calc_cls.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) calc_cls.o $(LINK)


CalcMP:  $(OBJ) calc_mp.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) calc_mp.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test CalcCls SolvePert CalcMP


firstorder.o:	pertobj.o  ../include/firstorder_equations.h ../include/firstorder_source.h ../include/recfast_source.h ../include/recfast_head.h ../include/firstorder_basic_utils.h ../include/firstorder_ic.h ../include/firstorder_ic_EFT.h ../include/firstorder_equations_EFT.h
wrapper.o: 	firstorder.o lensing.o Cls_post_process.o
test.o:		wrapper.o
test_one_mode.o: wrapper.o
Cls_post_process.o: firstorder.o lensing.o
zeta.o:	wrapper.o
testzeta.o:	wrapper.o
testcls.o: 	$(OBJ)
