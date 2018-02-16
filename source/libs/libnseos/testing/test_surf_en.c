#include <stdio.h>
#include <math.h>

#include "../nseos/empirical.h"
#include "../nseos/nuclear_en.h"

int main(void)
{
    int aa;
    int iii;
    double ii;
    struct parameters satdata_ref;
    double n0;
    double evac;
    struct hnm meta0;
    double eb0, eb2;
    double ebulk;
    double ecoul;
    double es;
    double es_normalized;

    aa = 200;
    satdata_ref = assign_param_ref(satdata_ref);

    for(iii = -80; iii < 81; iii++)
    {
        ii = iii/100.;
        n0 = satdata_ref.rhosat0*(1.-3.*satdata_ref.lsym0*ii*ii/(satdata_ref.ksat0 + satdata_ref.ksym0*ii*ii));
        /* n0 = satdata_ref.rhosat0; */

        evac = calc_dl_etf_meta_model_nuclear_en(satdata_ref, 2, aa, ii, n0, 0.);

        meta0 = calc_meta_model_nuclear_matter(satdata_ref, 2, n0, 0);
        eb0 = meta0.enpernuc*aa;
        eb2 = meta0.jsym*aa;
        ebulk = eb0 + eb2*ii*ii; // quadratic expansion

        ecoul = calc_coulomb_en(satdata_ref, aa, ii, n0, 0.);

        es = evac - ebulk - ecoul;
        es_normalized = es/pow(aa,2./3.);

        printf("%g %g %g\n", ii, es_normalized, eb2*(1./(1.+eb2/11.5514/pow(aa,4./3.))-1.));
    }

    return 0;
}
