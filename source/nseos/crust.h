#ifndef H_ICRUST
#define H_ICRUST

#include <gsl/gsl_vector.h>

#include "empirical.h"
#include "nuclear_surface_en.h"

struct crust_fun_4d
{
    double f_stability;
    double f_beta;
    double f_muneq;
    double f_presseq;
};

struct compo
{
    double aa;
    double del;
    double n0;
    double ng;
}; 

// interpolation of the proton shell energy 
// from ETFSI calculation with a BSk functional;
// see: https://academic.oup.com/mnras/article/481/3/2994/5090419
double get_shell_energy_per_nucleon(float nb_, float zz_);
double calc_ion_free_en_sol(
        struct parameters satdata, struct sf_params sparams,
        double aa_, double del_, double n0_, double np_, double tt_);
double calc_ion_free_en_liq(
        struct parameters satdata, struct sf_params sparams,
        double aa_, double del_, double n0_, double np_, double tt_);
struct crust_fun_4d calc_crust_fun_4d(
        struct parameters satdata, struct sf_params sparams, 
        double aa_, double del_, double n0_, double np_, double ng_,
        double tt_, char phase[]);
struct crust_fun_4d calc_crust_fun_zz_fixed(struct parameters satdata, 
        struct sf_params sparams, 
        double aa_, double del_, double n0_, double np_, double ng_);

struct rparams_crust
{
    double np;
    double tt;
    int zz;

    char phase[3];

    struct parameters satdata;
    struct sf_params sparams;
};

// outer-crust composition w/o quantum effects
int assign_ocrust_fun_3d(const gsl_vector * x, void *params, gsl_vector * f);
struct compo calc_ocrust3d_composition(double nb_, double tt_, 
        char phase[], double *guess,
        struct parameters satdata, struct sf_params sparams);

// inner-crust composition w/o quantum effects
int assign_icrust_fun_4d(const gsl_vector * x, void *params, gsl_vector * f);
struct compo calc_icrust4d_composition(double nb_, double tt_, 
        char phase[], double *guess,
        struct parameters satdata, struct sf_params sparams);

// inner-crust composition w/ perturbative addition of proton sh. corr.
int assign_icrust_fun_zz_fixed(
        const gsl_vector * x, void *params, gsl_vector * f);
struct compo calc_icrust_composition_zz_fixed(
        double nb_, int zz_, double *guess,
        struct parameters satdata, struct sf_params sparams);
struct compo calc_icrust_composition_w_shl(double nb_,
        struct parameters satdata, struct sf_params sparams);

double calc_muncl(struct parameters satdata, struct sf_params sparams, 
        struct compo eq, double nb_);
double calc_mupcl(struct parameters satdata, struct sf_params sparams, 
        struct compo eq, double nb_);
double calc_crust_ws_cell_free_energy_density(
        struct parameters satdata, struct sf_params sparams, 
        struct compo eq, 
        double nb_, double tt_,
        char phase[]);
double calc_ngas_pressure(struct parameters satdata, double ng_);
double calc_ion_pressure(struct parameters satdata, struct sf_params sparams, 
        double aa_, double del_, double n0_, double np_,
        double tt_,
        char phase[]);
double calc_crust_ws_cell_pressure(struct parameters satdata, 
        struct sf_params sparams, 
        struct compo eq, 
        double nb_, double tt_, 
        char phase[]);
void print_state_crust(struct parameters satdata, struct sf_params sparams, 
        struct compo eq, 
        double nb_, double tt_, char phase[], 
        FILE *compo, FILE *eos);

#endif // H_ICRUST
