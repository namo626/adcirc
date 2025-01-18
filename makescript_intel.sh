export NETCDFHOME=$TACC_NETCDF_DIR

cmake .. -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCMAKE_Fortran_COMPILER=ifx -DENABLE_WARN_ELEV_DEBUG=ON -DCMAKE_Fortran_FLAGS='-xHost'

cmake .. -DBUILD_ADCIRC=ON -DBUILD_PADCIRC=ON -DBUILD_ADCPREP=ON -DBUILD_PADCSWAN=OFF -DENABLE_OUTPUT_NETCDF=ON
