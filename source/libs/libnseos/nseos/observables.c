#include "phyconst.h"
#include "coulomb_en.h"
#include "observables.h"

double calc_ws_cell_energy_density(double aa_, double n0_, double np_, double ng_, double enuc, double epsg,
        double nb_)
{

    double vws;
    double epseltot;
    double epsws; 

    vws = aa_*(1. - ng_/n0_)/(nb_ - ng_);
    epseltot = calc_egas_energy_density(np_);

    epsws = enuc/vws + epseltot + epsg*(1.-aa_/n0_/vws)
        + np_*(rmp-rmn) + nb_*(rmn-amu);

    return epsws;
}
