NSEoS
=====

NSEoS (Neutron Star Equation of State) is a library that aims to provide the useful tools to calculate 
the composition of the crust and the equation of state of neutrons stars according to different nuclear models.

Requirements
------------

* [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)

Installation
------------

    git clone https://github.com/thomascarreau/NSEoS
    cd NSEoS/source/testing
    make

Usage
-----

In the testing directory:

    ./nseos compo.out eos.out

The first file gives you the crust composition (mass of the cluster, global asymmetry in the cluster, number of charges, 
cluster density, gas density, radius of the cell, and energy density in the cell) as a function of the baryonic density. 
The second file gives you the EoS that is the pressure in the cell as a function of the baryonic density. 
You can change the parameter set to be use and the nuclear modeling options in `source/nseos/modeling.h`.

### Available sets of empirical parameters

#### Skyrme

* BSk14

* BSk16

* BSk17

* NRAPR

* RATP

* SkO

* SLy230a

* SLy230b

* SLy4

#### Relativistic

* NL3

* TM1

* DD-ME1

* DD-ME2

To do
-----

* Study the influence of the Taylor expansion order on the crust-core transition: the higher the expansion order, 
    the higher the transition density probably because of the low density neutron gas description.

* Why does the crust-core transition density depend on the initial baryonic density?: the density step must be 
    small enough.

* Relative uncertainty of b is large: why?: lack of data.

* Replace ELFc for ELFd modeling: not pertinent to take a high density point as it is done is the paper, ELFc 
    modeling give better results.

* Work out the crust-core transition at constant pressure

* Write outer-crust EoS calculation code

* Write a single code for the crust and core calculation (current situation: inner crust + core)
