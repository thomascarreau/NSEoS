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
    cd NSEoS/source/apps/nseos
    make

Usage
-----

In `NSEoS/source/apps/nseos`:

    ./nseos set.in crust.out core.out eos.out tov.out

The first output file gives you the crust composition (number density, mass of the cluster, global asymmetry in the cluster, 
number of charges, cluster density, gas density, and radius of the cell). 
The second output file gives you the core composition (number density, fraction of protons, electrons, and muons). 
The third output file gives you the EoS (mass density, and pressure). Finally, the last output file gives you the TOV solution (central density, central pressure, radius, mass, core radius, core mass, normalized moment of inertia, and fraction of moment of inertia residing in the crust).

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

You can add new sets in `NSEoS/source/input/satdata`.
