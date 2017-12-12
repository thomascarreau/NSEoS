#include <stdio.h>

#include "../nseos/nuclear_en.h"

int main(void)
{
    struct parameters satdata0;
    struct skyrme_parameters coeff0;
    double nn;
    double zz;
    double aa;
    double del;
    double rho0;
    double evac2;
    double evac3;
    double evac4;
    double evacsly;

    zz = 28.;
    satdata0 = assign_param(satdata0);
    coeff0 = assign_skyrme_param(coeff0);

    for(nn = 25.; nn < 46.; nn += 1.)
    {
        aa = nn + zz;
        del = (nn - zz)/aa;
        rho0 = satdata0.rhosat0*(1.-3.*satdata0.lsym0*del*del/(satdata0.ksat0+satdata0.ksym0*del*del));
        evac2 = calc_cldm_meta_model_nuclear_en(satdata0, 2, aa, del, rho0, 0.);
        evac3 = calc_cldm_meta_model_nuclear_en(satdata0, 3, aa, del, rho0, 0.);
        evac4 = calc_cldm_meta_model_nuclear_en(satdata0, 4, aa, del, rho0, 0.);
        evacsly = calc_cldm_skyrme_based_nuclear_en(satdata0, coeff0, aa, del, rho0, 0.);
        printf("%g %g %g %g %g\n", nn, evac2/aa, evac3/aa, evac4/aa, evacsly/aa);
    }

    return 0;
}
