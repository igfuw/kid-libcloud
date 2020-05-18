#!/usr/bin/env sh
# one argument - python command
set -ex
gcc -fPIC -shared ptrutil.c -o ptrutil.so
tar xvzf kid_a_setup-20180125.tar.gz
cp mphys_libcloud_lgr.f90 kid_a_setup/src/mphys_libcloud_lgr.f90 
cp kida_icmwSC_2D_libcloud_lgr.nml kid_a_setup/namelists/
cp kida_icmw1D_libcloud_lgr.nml kid_a_setup/namelists/
cp ICMW_SC_input.nml kid_a_setup/namelists/
cd kid_a_setup
patch -p1 < ../kid_a_setup_20180125.diff
make SHELL=/bin/bash CASE=1D COMPILER=gfortran NCPATH=/usr all

# run 1D case with w=3m/s
sed -i 's/wctrl(1)=2/wctrl(1)=3/g' namelists/kida_icmw1D_libcloud_lgr.nml
FILEOUT=output LD_LIBRARY_PATH=..:bin:$LD_LIBRARY_PATH LD_PRELOAD=ptrutil.so $1 ../kid_1D.py --backend=OpenMP
# test the results
cd output
ncdiff -O -v liquid_water_path 1D_out.nc ../../refdata/w3_N50_NORAIN.nc diff.nc
# check if the difference is small enough
ncap2 -A -s "diff_flag=liquid_water_path<2e-3 && liquid_water_path>-2e-3" diff.nc
ncra -O diff.nc diff_mean.nc
ncks -V -v diff_flag -C -H diff_mean.nc
ncks -V -v diff_flag -C -H diff_mean.nc | grep '1'
cd ..

# run 1D case with w=2m/s
sed -i 's/wctrl(1)=3/wctrl(1)=2/g' namelists/kida_icmw1D_libcloud_lgr.nml
FILEOUT=output LD_LIBRARY_PATH=..:bin:$LD_LIBRARY_PATH LD_PRELOAD=ptrutil.so $1 ../kid_1D.py --backend=OpenMP
# test the results
cd output
ncdiff -O -v liquid_water_path 1D_out.nc ../../refdata/w2_N50_NORAIN.nc diff.nc
# check if the difference is small enough
ncap2 -A -s "diff_flag=liquid_water_path<2e-3 && liquid_water_path>-2e-3" diff.nc
ncra -O diff.nc diff_mean.nc
ncks -V -v diff_flag -C -H diff_mean.nc
ncks -V -v diff_flag -C -H diff_mean.nc | grep '1'
cd ..

set +ex # see https://github.com/travis-ci/travis-ci/issues/6522

