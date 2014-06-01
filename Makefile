SHELL=/bin/bash

include configure.in
include compile_rules.in

INCLUDE=-I./include

TestRun:	test/test.o
	$(FC) $(FFLAGS) -o test/Test test/test.o -L./lib/ -lcoop

clean:
	rm -f *.*~ Makefile~ include/*.*~ data/*.*~ test/*.o test/*.mod
