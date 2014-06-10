#! /bin/bash
mv configure.in configure_local.in
mv test/test.f90 test/test_local.f90
mv test/example1.f90 test/example1_local.f90
mv test/example2.f90 test/example2_local.f90
git pull origin master
mv configure_local.in configure.in
./clear.sh
./compile.sh
