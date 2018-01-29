kid-libcloud
============

This is work in progress related to KiD-A comparizon project 
(http://appconv.metoffice.com/kid_a_intercomparison/kid_a/home.html).

The lagrangian scheme from libcloudph++ library 
(https://github.com/igfuw/libcloudphxx) is joined to the Kinematic Driver 
model (KiD) using python bindings to the library (http://arxiv.org/abs/1504.01161).
The KiD model is started by python code kid.py.   

The KiD-A code stored in this repository was downloaded from:
http://appconv.metoffice.com/kid_a_intercomparison/kid_a/kid_a_setup.tar.gz

All changes needed to the KiD code are in kid_a_setup.diff.

In order to run the model with the lagrangian scheme (having already installed 
the libcloudph++ library) you need to follow these steps: 

  - cd kid-libcloud
  - gcc -fPIC -shared ptrutil.c -o ptrutil.so
  - tar xvzf kid_a_setup.tar.gz
  - cp mphys_libcloud_lgr.f90 kid_a_setup/src/mphys_libcloud_lgr.f90 
  - cp SC_2D_libcloud_lgr.nml kid_a_setup/namelists/SC_2D_libcloud_lgr.nml 
  - cp kida_SC_2D_libcloud_lgr.nml kid_a_setup/namelists/kida_SC_2D_libcloud_lgr.nml 
  - cd kid_a_setup

when using the new repo from github:
  - patch -p1 < ../kid_a_setup_20180125.diff
  - make SHELL=/bin/bash CASE=SC_2D COMPILER=gfortran NCPATH=/usr all
when using the original repo:
  - patch -p1 < ../kid_a_setup.diff
  - make SHELL=/bin/bash CASE=WMO_CASE1 COMPILER=gfortran NCPATH=/usr all

  - LD_LIBRARY_PATH=..:bin LD_PRELOAD=ptrutil.so python ../kid.py 
