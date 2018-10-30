#include <math.h>

#include "mathconst.h"
#include "phyconst.h"
#include "coulomb.h"

double calc_coulomb_en(struct parameters satdata, double aa_, double ii_)
{
    double rsat;
    double ac;
    double zz;

    rsat = pow(3./4./PI/satdata.rhosat0,1./3.);
    ac = 3./5.*ALPHAFS*HBARC/rsat;
    zz = aa_*(1.-ii_)/2.;

    return ac*zz*zz*pow(aa_,-1./3.);
}

double calc_lattice_en(struct parameters satdata, double aa_, double ii_, double n0_, double np_)
{
    double rsat;
    double ac;
    double rpt;
    double fws; // screening function
    double zz;

    rsat = pow(3./4./PI/satdata.rhosat0,1./3.);
    ac = 3./5.*ALPHAFS*HBARC/rsat;
    rpt = 2.*np_/(1.-ii_)/n0_;
    fws = -1.5*pow(rpt,1./3.) + 0.5*rpt;
    zz = aa_*(1.-ii_)/2.;

    return ac*zz*zz*pow(aa_,-1./3.)*fws;
}

double calc_lattice_derivative(struct parameters satdata, double aa_, double ii_, double n0_, double np_)
{
    double rsat;
    double ac;
    double rpt;

    rsat = pow(3./4./PI/satdata.rhosat0,1./3.);
    ac = 3./5.*ALPHAFS*HBARC/rsat;
    rpt = 2.*np_/(1.-ii_)/n0_;

    return 0.5*ac*pow(aa_,2./3.)*np_/n0_*(pow(rpt,-2./3.) - 1.);
}

double calc_lattice_pressure(struct parameters satdata, double aa_, double ii_, double n0_, double np_)
{
    double rsat;
    double ac;
    double rpt;
    double lattice_pressure;

    rsat = pow(3./4./PI/satdata.rhosat0,1./3.);
    ac = 3./5.*ALPHAFS*HBARC/rsat;
    rpt = 2.*np_/(1.-ii_)/n0_;
    lattice_pressure = 0.5*ac*pow(aa_,2./3.)*np_*np_/n0_
        *(1.-pow(rpt,-2/3.));

    return lattice_pressure;
}
