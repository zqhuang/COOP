SHELL = /bin/sh

TYPEDEF = ../typedef
UTILS = ../utils/

include ../configure.in
include ../compile_rules.in

BGLIB = libcoopbackground.a
BGOBJ = background.o wrapper.o

default: $(BGLIB)

$(BGLIB): $(BGOBJ) 
	ar -r $@ $?

Test: $(BGOBJ) test.o  ../include/background.h $(UTILS)/libcooputils.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(BGOBJ) test.o $(LINK)

TestEFT: $(BGOBJ) testeft.o  ../include/background.h $(UTILS)/libcooputils.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(BGOBJ) testeft.o $(LINK)


TestCPLDE: $(BGOBJ) test_cpl_de.o  ../include/background.h $(UTILS)/libcooputils.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(BGOBJ) test_cpl_de.o $(LINK)

FEPS: $(BGOBJ) fiteps.o  ../include/background.h $(UTILS)/libcooputils.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(BGOBJ) fiteps.o $(LINK)

TEPS: $(BGOBJ) epstbl.o  ../include/background.h $(UTILS)/libcooputils.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(BGOBJ) epstbl.o $(LINK)


PlotData: $(BGOBJ) plotdata.o ../include/background.h $(UTILS)/libcooputils.a
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(BGOBJ) plotdata.o $(LINK)

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ Test TestCPLDE TEPS FEPS

test.o:   $(BGOBJ)
testde.o:   $(BGOBJ)
plotdata.o: $(BGOBJ)
