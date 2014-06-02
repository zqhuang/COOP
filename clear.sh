#! /bin/bash
make clean
cd typedef 
make clean
cd ../background
make clean
cd ../firstorder
make clean
cd ../lib
make clean
cd ../include
rm -f *.mod \#*\# *.*~
cd ..
make clean
