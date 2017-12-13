#ifndef H_NUCLEAR_SURFACE_EN
#define H_NUCLEAR_SURFACE_EN

double calc_ldm_surf_en(double aa_);
double calc_ldm_wbbpcorr_surf_en(struct skyrme_parameters coeff,
        double aa_, double ii_, double n0_, double ng_); // see: Baym et al., 1971

#endif // H_NUCLEAR_SURFACE_EN
