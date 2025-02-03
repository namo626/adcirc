#!/usr/bin/env sh

cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS_DEBUG="-Wall -Wextra -ffpe-trap=invalid,zero,overflow -g -fbounds-check"

cmake .. -DBUILD_ADCIRC=ON -DBUILD_PADCIRC=OFF -DBUILD_ADCPREP=OFF -DBUILD_PADCSWAN=OFF
