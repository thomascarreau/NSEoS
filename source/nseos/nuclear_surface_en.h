#ifndef H_NUCLEAR_SURFACE_EN
#define H_NUCLEAR_SURFACE_EN

#include <gsl/gsl_vector.h>

#include "empirical.h"

// LDM
//====
double calc_ldm_surface_en(struct parameters satdata, double aa_);

// analytical ETF, see: Aymard et al. (2016)
//==========================================
struct f_params
{
    int k;
    double gamma;
}; // integration of Fermi functions
double my_integrand(double x, void *params_ptr);
double eta_function(int a, double b);
double calc_etf_ana_surface_en(struct parameters satdata, double aa_, double ii_, double n0_);

// LS corrected surface energy, see: Lattimer and Swesty (1991); arXiv:1110.4043
//==============================================================================
struct data
{
    int *zz;
    int *aa;
    double *be;
    struct parameters satdata;
};
int be_f (const gsl_vector * x, void *data, gsl_vector * f);
struct sf_params fit_sf_params(struct parameters satdata);
struct sf_params
{
    double sigma0;
    double b;
};
double calc_ls_surface_en(struct sf_params params, double aa_, double ii_, double n0_);

#endif // H_NUCLEAR_SURFACE_EN
