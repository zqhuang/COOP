echo "#define MAINPATH \""`pwd`"\"" > include/paths.h
cd typedef
make -w
cd ../utils
make -w
cd ../background
make -w
cd ../firstorder
make -w
cd ../lib
make -w
cd ..

