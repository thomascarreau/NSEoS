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

double calc_lattice_en(struct parameters satdata, 
        double aa_, double ii_, double n0_, double np_)
{
    double rsat;
    double cm;
    double rpt;
    double zz;

    rsat = pow(3./4./PI/satdata.rhosat0,1./3.);
    cm = 0.895929255682; // bcc lattice; Table 2.4 of Haensel book
    rpt = 2.*np_/(1.-ii_)/n0_;
    zz = aa_*(1.-ii_)/2.;

    return ALPHAFS*HBARC/rsat*zz*zz*pow(aa_,-1./3.)
        *(-cm*pow(rpt,1./3.) + 3./10.*rpt);
}

double calc_lattice_pressure(struct parameters satdata, 
        double aa_, double ii_, double n0_, double np_)
{
    double rsat;
    double cm;
    double rpt;
    double zz;

    rsat = pow(3./4./PI/satdata.rhosat0,1./3.);
    cm = 0.895929255682; // bcc lattice; Table 2.4 of Haensel book
    rpt = 2.*np_/(1.-ii_)/n0_;
    zz = aa_*(1.-ii_)/2.;

    return ALPHAFS*HBARC/rsat*np_*zz*pow(aa_,-1./3.)
        *(-cm/3.*pow(rpt,1./3.) + 3./10.*rpt);
}
