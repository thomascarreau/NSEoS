NSEoS
=====

Installation
------------

    git clone https://github.com/thomascarreau/NSEoS
    cd NSEoS/source/libs/libnseos/testing
    make

Usage
-----

In the testing directory:

    ./nseos compo.out eos.out

The first file gives you the crust composition (mass of the cluster, global asymmetry in the cluster, number of charges, cluster density, gas density, radius of the cell, and energy density in the cell) as a function of the baryonic density. The second file give you the EoS that is the pressure in the cell as a function of the baryonic density. 
You can change the parameter set to be use and the nuclear modeling in `choices.h`.

To do
-----

* Replace ELFc for ELFd modeling

* Fit the surface energy parameters from AME2012

* Write outer-crust EoS calculation code

* Write a single code for the crust and core calculation

* Work out the crust-core transition at constant pressure
