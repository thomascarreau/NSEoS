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

    ./nseos set.in crust.out core.out eos.out

The first output file gives you the crust composition (mass of the cluster, global asymmetry in the cluster, number of charges, 
cluster density, gas density, radius of the cell, and energy density in the cell) as a function of the baryonic density. 
The second output file gives you the core composition (fraction of protons, of electrons, and of muons). 
Finally, the third output file gives you the EoS that is the pressure in the cell as a function of the baryonic density.

### Available sets of empirical parameters

#### Skyrme

* BSk14: `bsk14.in`

* BSk16: `bsk16.in`

* BSk17: `bsk17.in`

* NRAPR: `nrapr.in`

* RATP: `ratp.in`

* SkO: `sko.in`

* SLy230a: `sly230a.in`

* SLy230b: `sly230b.in`

* SLy4: `sly4.in`

#### Relativistic

* NL3: `nl3.in`

* TM1: `tm1.in`

* DD-ME1: `ddme1.in`

* DD-ME2: `ddme2.in`

To do
-----

* Study the influence of the Taylor expansion order on the crust-core transition: the higher the expansion order, 
    the higher the transition density probably because of the low density neutron gas description.

* Why does the crust-core transition density depend on the initial baryonic density?: the density step must be 
    small enough.

* Relative uncertainty of b is large: why?: lack of data.

* Work out the crust-core transition at constant pressure
