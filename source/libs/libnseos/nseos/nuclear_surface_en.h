#ifndef H_NUCLEAR_SURFACE_EN
#define H_NUCLEAR_SURFACE_EN

// LDM
double calc_ldm_surf_en(double aa_);
double calc_ldm_wbbpcorr_surf_en(struct skyrme_parameters coeff,
        double aa_, double ii_, double n0_, double ng_); // see: Baym et al., 1971

// ETF (see: Aymard et al., 2016)
struct f_params
{
    int k;
    double gamma;
}; // integration of Fermi functions
double my_integrand(double x, void *params_ptr);
double eta_function(int a, double b);
double calc_etf_surface_energy(struct parameters satdata, double aa_, double ii_, double n0_);

#endif // H_NUCLEAR_SURFACE_EN
