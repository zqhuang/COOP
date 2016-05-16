SHELL = /bin/sh

TYPEDEF=../typedef
UTILS=../utils
BACKGROUND=../background
FIRSTORDER=../firstorder

include ../configure.in
include ../compile_rules.in

LIB = libcoopnonlinear.a
OBJ = nonlinear.o ellipse_collapse.o

default: $(LIB) GetZ1 Solve TabZ1

$(LIB):	$(OBJ) 
	ar -r $@ $?

Test: $(OBJ) test.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) test.o $(LINK)


TestSimp: $(OBJ) testsimple.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) testsimple.o $(LINK)

Solve: $(OBJ) solve_collapse.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) solve_collapse.o $(LINK)


TabZ1: $(OBJ) make_zvir1_table.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) make_zvir1_table.o $(LINK)


GetZ1: $(OBJ) get_zvir1.o $(FIRSTORDER)/libcoopfirstorder.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(OBJ) get_zvir1.o $(LINK)


clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ GetZ1 Solve TestSimp Test TabZ1


test.o: $(OBJ)

fnl.o: $(OBJ)
