echo "#define MAINPATH \""`pwd`"\"" > headfiles/paths.h
cd typedef
make -w
cd ../
cd background
make -w
cd ../lib
make -w
cp coop_wrapper.mod ../include/
cd ..
