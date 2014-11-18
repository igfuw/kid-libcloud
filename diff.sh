#/bin/bash
diff -ruN \
  --exclude=mphys_libcloud_lgr.f90 \
  --exclude=*~ \
  --exclude=*.mod \
  --exclude=*.o \
  --exclude=bin \
  --exclude=4A_code \
  --exclude=UM_source \
  --exclude=.includes \
  --exclude=case.used \
  --exclude=output \
  --exclude=gmon.out \
  kid_a_setup.orig kid_a_setup \
  > kid_a_setup.diff
