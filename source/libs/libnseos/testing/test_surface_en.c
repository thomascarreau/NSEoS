#include <stdio.h>
#include <math.h>

#include "../nseos/empirical.h"
#include "../nseos/nuclear_matter.h"
#include "../nseos/nuclear_surface_en.h"

int main(void)
{
    //========================================================== IN VACUUM
    /* struct parameters satdata0; */
    /* double aa; */
    /* int iii; */
    /* double ii; */
    /* double n0; */
    /* double esurf_cldm; */
    /* double esurf_etf; */
    /* double esurf_ls; */

    /* satdata0 = assign_param(satdata0); */
    /* aa = 150.; */

    /* for(iii = -90; iii < 91; iii += 1) */
    /* { */
    /*     ii = iii/100.; */
    /*     n0 = satdata0.rhosat0*(1.-3.*satdata0.lsym0*ii*ii/(satdata0.ksat0+satdata0.ksym0*ii*ii)); */
    /*     esurf_cldm = calc_ldm_surface_en(aa); */
    /*     esurf_etf = calc_etf_surface_en(satdata0, aa, ii, n0); */
    /*     esurf_ls = calc_ls_surface_en(satdata0, aa, ii); */
    /*     printf("%g %g %g %g\n", ii, esurf_cldm/pow(aa,2./3.), esurf_etf/pow(aa,2./3.), esurf_ls/pow(aa,2./3.)); */
    /* } */

    //========================================================== IN MEDIUM
    struct parameters satdata0;
    double aa;
    int iii;
    double ii;
    double n0;
    double esurf_cldm, esurf_cldm_g;
    double esurf_etf, esurf_etf_g;
    double esurf_ls, esurf_ls_g;
    struct hnm test;

    satdata0 = assign_param(satdata0);
    aa = 300.;

    for(iii = 20; iii < 91; iii += 1)
    {
        ii = iii/100.;
        n0 = satdata0.rhosat0*(1.-3.*satdata0.lsym0*ii*ii/(satdata0.ksat0+satdata0.ksym0*ii*ii));
        esurf_cldm = calc_ldm_wbbpcorr_surface_en(satdata0, 4, aa, ii, n0, 0.);
        esurf_etf = calc_etf_wbbpcorr_surface_en(satdata0, 4, aa, ii, n0, 0.);
        esurf_ls = calc_ls_wbbpcorr_surface_en(satdata0, 2, aa, ii, n0, 0.);
        esurf_cldm_g = calc_ldm_wbbpcorr_surface_en(satdata0, 4, aa, ii, n0, 0.4*n0);
        esurf_etf_g = calc_etf_wbbpcorr_surface_en(satdata0, 4, aa, ii, n0, 0.4*n0);
        esurf_ls_g = calc_ls_wbbpcorr_surface_en(satdata0, 2, aa, ii, n0, 0.4*n0);
        test = calc_meta_model_nuclear_matter(satdata0, 4, n0, ii);
        printf("%g %g %g %g %g %g %g %g\n", ii, esurf_cldm/pow(aa,2./3.), esurf_etf/pow(aa,2./3.), esurf_ls/pow(aa,2./3.),
                esurf_cldm_g/pow(aa,2./3.), esurf_etf_g/pow(aa,2./3.), esurf_ls_g/pow(aa,2./3.), test.enpernuc);
    }

    return 0;
}
