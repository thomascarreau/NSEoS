#include <math.h>
#include <complex.h>

#include "mathconst.h"
#include "phyconst.h"
#include "lepton.h"

double calc_egas_free_energy_density(double ne_, double tt_)
{
    // see: eq. (2.65) of "Neutron Stars 1: Equation of State and Structure"
    double pr = pow(MEL,4.)/pow(HBARC,3.);
    double tr = tt_/MEL;
    double xr = HBARC*pow(3.*PI2*ne_,1./3.)/MEL;
    double gammar = pow(1.+xr*xr,0.5);
    double egas_energy_density;

    egas_energy_density = cpow(MEL,4.)/8./PI2/cpow(HBARC,3.)
        *((2.*xr*xr+1.)*xr*gammar - log(xr + gammar));

    return egas_energy_density - pr/6.*tr*tr*xr*gammar;
}

double calc_egas_chemical_potential(double ne_, double tt_)
{
    double tr = tt_/MEL;
    double xr = HBARC*pow(3.*PI2*ne_,1./3.)/MEL;
    double gammar = pow(1.+xr*xr,0.5);
    double egas_fermi_energy;

    egas_fermi_energy = 
        cpow(MEL,3.)/8./cpow(3.*ne_*PI2,2./3.)/cpow(HBARC,2.)
        *(gammar*(1.+6.*xr*xr) + xr*xr*(2.*xr*xr+1.)/gammar - 1./gammar);

    return egas_fermi_energy
        - pow(MEL,3.)/18./pow(HBARC,2.)*pow(3.*PI2,1./3.)*pow(ne_,-2./3.)
        *(xr*xr/gammar + gammar)*tr*tr;
}

double calc_egas_pressure(double ne_, double tt_)
{
    return ne_*calc_egas_chemical_potential(ne_, tt_) 
        - calc_egas_free_energy_density(ne_, tt_);
}

double calc_ugas_energy_density(double nu_)
{
    double t;
    double ugas_energy_density;

    t = HBARC*cpow(3.0*PI2*nu_,1./3.)/MMU;
    ugas_energy_density = cpow(MMU,4.)/8./PI2/cpow(HBARC,3.)
        *((2.*t*t+1.)*t*cpow(t*t+1.,1./2.) - log(t + cpow(t*t+1.,1./2.)));

    return ugas_energy_density;
}

double calc_ugas_chemical_potential(double nu_)
{
    double t;
    double ugas_chemical_potential;

    t = HBARC*cpow(3.0*PI2*nu_,1./3.)/MMU;
    ugas_chemical_potential = 
        cpow(MMU,3.)/8./cpow(3.*nu_*PI2,2./3.)/cpow(HBARC,2.)
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
