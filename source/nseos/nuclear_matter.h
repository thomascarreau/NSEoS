#ifndef H_NUCLEAR_MATTER
#define H_NUCLEAR_MATTER

#include "empirical.h"
#include "mathconst.h"
#include "phyconst.h"

struct hnm {
  double enpernuc;  // energy per nucleon
  double spernuc;   // entropy per nucleon
  double mun;       // neutron chemical potential
  double mup;       // proton chemical potential
  double fpernuc;   // free energy per nucleon
  double p;         // pressure
  double jsym;      // symmetry energy
  double lsym;      // first derivative of symmetry energy
  double ksym;      // second derivative of symmetry energy
};

// Meta-model, see: Margueron et al. (2018)
//=========================================
double calc_meta_model_low_density_correction(
    struct parameters satdata, int max_order, int order, double xx_);
double calc_meta_model_low_density_correction_derivative(
    struct parameters satdata, int max_order, int order, double xx_);
double calc_meta_model_low_density_correction_second_derivative(
    struct parameters satdata, int max_order, int order, double xx_);
struct hnm calc_meta_model_nuclear_matter(struct parameters satdata,
    int max_order, double nn_, double ii_, double tt_);
double     calc_squared_nucleon_sound_velocity( // (vs/c)^2
        struct parameters satdata, int max_order, double nn_, double ii_,
        double tt_);

// Skyrme functional, see: Chabanat et al. (1997)
// !! only for T=0 !!
//===============================================
double     calc_asymmetry_factor(double m_, double ii_);
double     calc_asymmetry_factor_derivative(double m_, double ii_);
struct hnm calc_skyrme_nuclear_matter(
    struct skyrme_parameters coeff, double nn_, double ii_);

#endif // H_NUCLEAR_MATTER
