#ifndef H_COULOMB_EN
#define H_COULOMB_EN

#include "empirical.h"
#include "nuclear_surface_en.h"

double calc_coulomb_en(struct parameters satdata, double aa_, double ii_);

// Coulomb crystal
double calc_lattice_en(struct parameters satdata, 
        double aa_, double ii_, double n0_, double np_); // E_L
double calc_zp_en(struct parameters satdata, struct sf_params sparams, 
        double aa_, double ii_, double n0_, double np_); // E_{zp}
double calc_harmonic_contrib(
        struct parameters satdata, struct sf_params sparams, 
        double aa_, double del_, double n0_, double np_, 
        double tt_); // F_{harm}

// Coulomb liquid of ions
double calc_translational_free_en(
        struct parameters satdata, struct sf_params sparams, 
        double aa_, double del_, double n0_, double np_, 
        double tt_); // F_i^{non-int}
double calc_total_coulomb_contrib(
        double zz_, double tt_, double np_); // F_{ii,liq}

#endif // H_COULOMB_EN
