#include <math.h>

#include "mathconst.h"
#include "phyconst.h"
#include "nuclear_matter.h"
#include "coulomb_en.h"
#include "observables.h"

double calc_crust_ws_cell_energy_density(double aa_, double n0_, double np_, double ng_, double enuc, double epsg,
        double nb_)
{
    double vws;
    double epseltot;
    double epsws; 

    vws = aa_*(1. - ng_/n0_)/(nb_ - ng_);
    epseltot = calc_egas_energy_density(np_);
    epsws = enuc/vws + epseltot + epsg*(1.-aa_/n0_/vws)
        + np_*(RMP-RMN) + nb_*(RMN-AMU);

    return epsws;
}

double calc_core_ws_cell_energy_density(double ii_, struct hnm nucmat, double nb_)
{
    double np;
    double epseltot;
    double epsws;

    np = nb_*(1.-ii_)/2.;
    epseltot = calc_egas_energy_density(np);
    epsws = nb_*nucmat.enpernuc + epseltot + np*(RMP-RMN) + nb_*(RMN-AMU);

    return epsws;
}

double calc_egas_pressure(double np_)
{
    double egas_energy_density;
    double egas_chemical_potential;
    double egas_pressure;

    egas_energy_density = calc_egas_energy_density(np_);
    egas_chemical_potential = calc_egas_chemical_potential(np_);
    egas_pressure = np_*egas_chemical_potential - egas_energy_density;

    return egas_pressure;
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

double calc_ngas_pressure(double ng_, double epsg, double mug)
{
    return ng_*mug - epsg;
}

double calc_crust_ws_cell_pressure(struct parameters satdata, double aa_, double ii_, double n0_, double np_,
        double ng_, double epsg, double mug)
{
    double egas_pressure;
    double lattice_pressure;
    double ngas_pressure;
    double ws_cell_pressure;

    egas_pressure = calc_egas_pressure(np_);
    lattice_pressure = calc_lattice_pressure(satdata, aa_, ii_, n0_, np_);
    ngas_pressure = calc_ngas_pressure(ng_, epsg, mug);
    ws_cell_pressure = egas_pressure + lattice_pressure + ngas_pressure;

    return ws_cell_pressure;
}

double calc_core_ws_cell_pressure(double ii_, struct hnm nucmat, double nb_)
{
    double np;
    double egas_pressure;
    double ws_cell_pressure;

    np = nb_*(1.-ii_)/2.;
    egas_pressure = calc_egas_pressure(np);
    ws_cell_pressure = nucmat.p + egas_pressure; 

    return ws_cell_pressure;
}
