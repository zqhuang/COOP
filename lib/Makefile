SHELL = /bin/sh

FIRSTORDER=../firstorder

include ../configure.in
include ../compile_rules.in
LIB = libcoop.a
OBJ = wrapper.o 


default: $(LIB)

$(LIB):	$(OBJ) $(FIRSTORDER)/libcoopfirstorder.a
	ar -r $@ $?

clean:
	rm -f *.o *.*~ *.mod \#* *.a Makefile~ 

