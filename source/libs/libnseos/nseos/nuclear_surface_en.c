#include <math.h>

#include "nuclear_matter.h"
#include "nuclear_surface_en.h"

double calc_ldm_surf_en(double aa_)
{
    float as;

    as = 18.24; // SLy4 value; in MeV

    return as*pow(aa_,2./3.);
}

double calc_ldm_wbbpcorr_surf_en(struct skyrme_parameters coeff,
        double aa_, double ii_, double n0_, double ng_)
{
    double surf_en_vac;
    struct hnm bulk;
    struct hnm ngas;
    double surf_en_ns;

    surf_en_vac = calc_ldm_surf_en(aa_);
    bulk = calc_skyrme_nuclear_matter(coeff, n0_, ii_);
    ngas = calc_skyrme_nuclear_matter(coeff, ng_, 1.);
    surf_en_ns = surf_en_vac*(1.-ngas.enpernuc/bulk.enpernuc)*pow(1.-ng_/n0_,2./3.);

    return surf_en_ns;
}
