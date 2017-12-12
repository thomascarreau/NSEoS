#include <math.h>

#include "nuclear_surface_en.h"

double calc_ldm_surf_en(double aa_)
{
    float as;

    as = 18.24; // SLy4 value; in MeV

    return as*pow(aa_,2./3.);
}
