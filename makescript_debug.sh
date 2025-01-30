#!/usr/bin/env sh

cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS_DEBUG="   -g " -DENABLE_OUTPUT_NETCDF=ON

cmake .. -DBUILD_ADCIRC=ON -DBUILD_PADCIRC=ON -DBUILD_ADCPREP=OFF -DBUILD_PADCSWAN=OFF -DENABLE_WARN_ELEV_DEBUG=ON
