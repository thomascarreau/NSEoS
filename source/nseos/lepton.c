#include <math.h>
#include <complex.h>

#include "mathconst.h"
#include "phyconst.h"
#include "lepton.h"
 
double calc_egas_energy_density(double ne_)
{
    double t;
    double egas_energy_density;

    t = HBARC*cpow(3.0*PI2*ne_,1./3.)/MEL;
    egas_energy_density = cpow(MEL,4.)/8./PI2/cpow(HBARC,3.)*((2.*t*t+1.)*t*cpow(t*t+1.,1./2.) - log(t + cpow(t*t+1.,1./2.)));

    return egas_energy_density;
}

double calc_egas_chemical_potential(double ne_)
{
    double t;
    double egas_chemical_potential;

    t = HBARC*cpow(3.0*PI2*ne_,1./3.)/MEL;
    egas_chemical_potential = cpow(MEL,3.)/8./cpow(3.*ne_*PI2,2./3.)/cpow(HBARC,2.)
        *(cpow(t*t+1.,1./2.)*(1.+6.*t*t) + t*t*(2.*t*t+1.)/cpow(t*t+1.,1./2.)
                - 1./cpow(t*t+1.,1./2.));

    return egas_chemical_potential;
}

double calc_egas_pressure(double ne_)
{
    double egas_energy_density;
    double egas_chemical_potential;
    double egas_pressure;

    egas_energy_density = calc_egas_energy_density(ne_);
    egas_chemical_potential = calc_egas_chemical_potential(ne_);
    egas_pressure = ne_*egas_chemical_potential - egas_energy_density;

    return egas_pressure;
}

double calc_ugas_energy_density(double nu_)
{
    double t;
    double ugas_energy_density;

    t = HBARC*cpow(3.0*PI2*nu_,1./3.)/MMU;
    ugas_energy_density = cpow(MMU,4.)/8./PI2/cpow(HBARC,3.)*((2.*t*t+1.)*t*cpow(t*t+1.,1./2.) - log(t + cpow(t*t+1.,1./2.)));

    return ugas_energy_density;
}

double calc_ugas_chemical_potential(double nu_)
{
    double t;
    double ugas_chemical_potential;

    t = HBARC*cpow(3.0*PI2*nu_,1./3.)/MMU;
    ugas_chemical_potential = cpow(MMU,3.)/8./cpow(3.*nu_*PI2,2./3.)/cpow(HBARC,2.)
        *(cpow(t*t+1.,1./2.)*(1.+6.*t*t) + t*t*(2.*t*t+1.)/cpow(t*t+1.,1./2.)
                - 1./cpow(t*t+1.,1./2.));

    return ugas_chemical_potential;
}

double calc_ugas_pressure(double nu_)
{
    double ugas_energy_density;
    double ugas_chemical_potential;
    double ugas_pressure;

    ugas_energy_density = calc_ugas_energy_density(nu_);
    ugas_chemical_potential = calc_ugas_chemical_potential(nu_);
    ugas_pressure = nu_*ugas_chemical_potential - ugas_energy_density;

    if (nu_ == 0)
        return 0.;
    else
        return ugas_pressure;
}
