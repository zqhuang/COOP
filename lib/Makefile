SHELL = /bin/sh

FIRSTORDER=../firstorder

include ../configure.in
include ../compile_rules.in
LIB = libcoopwrapper.a
OBJ = wrapper.o 


default: $(LIB)

$(LIB):	$(OBJ) 
	ar -r $@ $?

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ 

