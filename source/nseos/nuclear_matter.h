#ifndef H_NUCLEAR_MATTER
#define H_NUCLEAR_MATTER

#include "mathconst.h"
#include "phyconst.h"
#include "empirical.h"

struct hnm
{
    double enpernuc;
    double p; // pressure
    double mun;
    double jsym;
    double lsym;
    double ksym;
};

// Meta-model, see: Margueron et al. (2018)
//=========================================
double calc_meta_model_low_density_correction(struct parameters satdata, int max_order, int order, double xx_);
double calc_meta_model_low_density_correction_derivative(struct parameters satdata, int max_order, int order, double xx_);
double calc_meta_model_low_density_correction_second_derivative(struct parameters satdata, int max_order, int order, double xx_);
struct hnm calc_meta_model_nuclear_matter(struct parameters satdata, int max_order, double nn_, double ii_); // ELFc meta-modeling

// Skyrme functional, see: Chabanat et al. (1997)
//===============================================
double calc_asymmetry_factor(double m_, double ii_);
double calc_asymmetry_factor_derivative(double m_, double ii_);
struct hnm calc_skyrme_nuclear_matter(struct skyrme_parameters coeff, double nn_, double ii_);

#endif // H_NUCLEAR_MATTER
