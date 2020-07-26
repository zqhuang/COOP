SHELL = /bin/sh

TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background

include ../configure.in
include ../compile_rules.in

LIB = libcoopfirstorder.a
OBJ = pertobj.o cls_indices.o lensing.o firstorder.o Cls_post_process.o wrapper.o

default: $(LIB) CalcCls CalcMP SolvePert

$(LIB): $(OBJ)
	ar -r $@ $?

Test: $(OBJ) test.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

TestZ: $(OBJ) testzeta.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) testzeta.o $(LINK)

NGMAP: $(OBJ) nonGmaps.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) nonGmaps.o $(LINK)


SolvePert: $(OBJ) solve_pert.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) solve_pert.o $(LINK)

CalcCls: $(OBJ) calc_cls.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) calc_cls.o $(LINK)


CalcMP:  $(OBJ) calc_mp.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) calc_mp.o $(LINK)

ZetaCorr: $(OBJ) do_zeta_corr.o $(BACKGROUND)/libcoopbackground.a 
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) do_zeta_corr.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test CalcCls SolvePert CalcMP


firstorder.o:	pertobj.o  ../include/firstorder_equations.h ../include/firstorder_source.h ../include/recfast_source.h ../include/recfast_head.h ../include/firstorder_basic_utils.h ../include/firstorder_ic.h ../include/firstorder_ic_EFT.h ../include/firstorder_equations_EFT.h

wrapper.o: 	firstorder.o lensing.o Cls_post_process.o
test.o:		$(OBJ)
calc_cls.o: $(OBJ)
Cls_post_process.o: firstorder.o lensing.o
zeta.o: $(OBJ)
testzeta.o:	wrapper.o
solve_pert.o: 	$(OBJ)
calc_mp.o:	$(OBJ)
