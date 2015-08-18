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
  - cd kid_a_setup
  - patch -p1 < ../kid_a_setup.diff
  - make SHELL=/bin/bash CASE=WMO_CASE1 COMPILER=gfortran NCPATH=/usr all
  - LD_LIBRARY_PATH=..:bin LD_PRELOAD=ptrutil.so python ../kid.py 
