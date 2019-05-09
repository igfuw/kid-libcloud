#/bin/bash
#diff -ruN \
#  --exclude=mphys_libcloud_lgr.f90 \
#  --exclude=kida_SC_2D_libcloud_lgr.nml \
#  --exclude='*~' \
#  --exclude='*.mod' \
#  --exclude='*.o' \
#  --exclude=bin \
#  --exclude=4A_code \
#  --exclude=UM_source \
#  --exclude=.includes \
#  --exclude=case.used \
#  --exclude=output \
#  --exclude=gmon.out \
#  kid_a_setup.orig kid_a_setup \
#  > kid_a_setup_20180125.diff

diff -ruN \
  --exclude=mphys_libcloud_lgr.f90 \
  --exclude=kida_icmwSC_2D_libcloud_lgr.nml \
  --exclude=SC_2D_libcloud_lgr.nml \
  --exclude='*~' \
  --exclude='*.mod' \
  --exclude='*.o' \
  --exclude=bin \
  --exclude=4A_code \
  --exclude=UM_source \
  --exclude=.includes \
  --exclude=case.used \
  --exclude='output*' \
  --exclude=gmon.out \
  --exclude='*.orig*' \
  --exclude='*rej*' \
  --exclude='*.out' \
  --exclude=kida_icmw1D_libcloud_lgr.nml \
  --exclude=ICMW_SC_input.nml \
  --exclude='.*' \
  --exclude='out' \
  kid_a_original_20180125_setup/ kid_a_setup > kid_a_setup_20180125.diff
