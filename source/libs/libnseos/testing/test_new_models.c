#include <stdio.h>
#include <math.h>

#include "../nseos/nuclear_en.h"

int main(void)
{
    struct parameters satdata0;
    double zz;
    int inn;
    double nn;
    double aa;
    double ii;
    double n0;
    double els;
    double edl;

    satdata0 = assign_param(satdata0);
    zz = 28.;

    for(inn = 25; inn < 46; inn += 1)
    {
        nn = inn;
        aa = nn+zz;
        ii = 1.-2.*zz/aa;
        n0 = satdata0.rhosat0*(1.-3.*satdata0.lsym0*ii*ii/(satdata0.ksat0+satdata0.ksym0*ii*ii));
        els = calc_ls_etf_meta_model_nuclear_en(satdata0, 2, aa, ii, n0, 0.);
        edl = calc_dl_etf_meta_model_nuclear_en(satdata0, 2, aa, ii, n0, 0.);
        printf("%g %g %g\n", nn, els/aa, edl/aa);
    }

    return 0;
}
