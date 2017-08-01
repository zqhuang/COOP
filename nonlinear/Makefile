SHELL = /bin/sh

TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background
FIRSTORDER=../firstorder

include ../configure.in
include ../compile_rules.in

LIB = libcoopnonlinear.a
OBJ = nonlinear.o ellipse_collapse.o coupledDE_collapse.o fR1Dobj.o 

default: $(LIB) GetZ1 Solve TabZ1 CDSolve CDExport

$(LIB):	$(OBJ) 
	ar -r $@ $?

FR1D: $(OBJ) fR_collapse_1Dtoy.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) fR_collapse_1Dtoy.o $(LINK)

SCC: solve_chameleon_collapse.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) solve_chameleon_collapse.o $(LINK)


Test: $(OBJ) test.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)

TestSimp: $(OBJ) testsimple.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) testsimple.o $(LINK)

CDExport: $(OBJ) coupledDE_export_cosmology.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) coupledDE_export_cosmology.o $(LINK)

CDSolve: $(OBJ) coupledDE_solve_collapse.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) coupledDE_solve_collapse.o $(LINK)


Solve: $(OBJ) solve_collapse.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) solve_collapse.o $(LINK)


TabZ1: $(OBJ) make_zvir1_table.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) make_zvir1_table.o $(LINK)


GetZ1: $(OBJ) get_zvir1.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) get_zvir1.o $(LINK)


clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ GetZ1 Solve TestSimp Test TabZ1 CDSolve CDExport FR1D SCC


test.o: $(OBJ)

make_zvir1_table.o:	$(OBJ)

get_zvir1.o: 	$(OBJ)

solve_coupledDE_collapse.o: $(OBJ)

solve_collapse.o: $(OBJ)

coupledDE_export_cosmology.o: $(OBJ)

fR_collapse_1Dtoy.o: fR1Dobj.o
fR1Dobj.o: ../include/updatephi_QS.h
