SHELL = /bin/sh
TYPEDEF = ../typedef
BACKGROUND = ../background
include ../configure.in
include ../compile_rules.in
LIB = libcoop.a
OBJ = wrapper.o 


default: $(LIB)

$(LIB):	$(OBJ)
	ar -r $@ $?

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ 

