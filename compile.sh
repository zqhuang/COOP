echo "#define MAINPATH \""`pwd`"\"" > include/paths.h
cd typedef
make -w
cp *.mod ../include/

cd ../utils
make -w
cp *.mod ../include/

cd ../background
make -w
cp *.mod ../include/

cd ../firstorder
make -w
cp *.mod ../include/

cd ../nonlinear
make -w
cp *.mod ../include/

cd ../forecast
make -w
cp *.mod ../include/

cd ../mapio
make -w
cp *.mod ../include/

cd ../lib
make -w
cp *.mod ../include/

cd ..
make -w

cd postprocess
make -w
make -w GD
cp *.mod ../include/
cd ..
