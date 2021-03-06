#ifndef H_COULOMB_EN
#define H_COULOMB_EN

#include "empirical.h"
#include "nuclear_surface_en.h"

// Coulomb energy for a spherical nucleus
double calc_coulomb_en(double aa_, double ii_, double n0_);

// Pasta
double calc_fd(double u_, int d_);

// Coulomb crystal
double calc_lattice_en(double aa_, double ii_, double n0_, double np_); // E_L
double calc_lattice_en_for_tm(double zz_, double np_);                  // E_L
double calc_finite_size_contrib(double aa_, double ii_, double n0_, double np_);
double calc_zp_en(double zz_, double np_, double mi_); // E_{zp}
double calc_harmonic_contrib(double zz_, double np_, double mi_,
    double tt_); // F_{harm}
double calc_anharmonic_contrib(double zz_, double np_,
    double tt_); // F_{anharm}

// Coulomb liquid of ions
double calc_translational_free_en(double zz_, double np_, double mi_,
    double tt_); // F_i^{non-int}
double calc_total_coulomb_contrib(
    double zz_, double np_, double tt_); // F_{ii,liq}

#endif // H_COULOMB_EN
