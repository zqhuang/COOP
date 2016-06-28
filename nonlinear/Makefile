SHELL = /bin/sh

TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background
FIRSTORDER=../firstorder

include ../configure.in
include ../compile_rules.in

LIB = libcoopnonlinear.a
OBJ = nonlinear.o ellipse_collapse.o coupledDE_collapse.o

default: $(LIB) GetZ1 Solve TabZ1 CDSolve

$(LIB):	$(OBJ) 
	ar -r $@ $?

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
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ GetZ1 Solve TestSimp Test TabZ1 CDSolve CDExport


test.o: $(OBJ)

fnl.o: $(OBJ)

solve_coupledDE_collapse.o: $(OBJ)

solve_collapse.o: $(OBJ)