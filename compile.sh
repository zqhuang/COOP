echo "#define MAINPATH \""`pwd`"\"" > headfiles/paths.h
cd typedef
make -w
cd ../
