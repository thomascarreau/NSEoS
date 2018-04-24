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
    double rho0;
    double rhog;
};
struct crust_fun_4d calc_crust_fun_4d(struct parameters satdata, struct sf_params sparams, 
        double aa_, double del_, double rho0_, double rhop_, double rhog_);
struct rparams_crust
{
    double rhop;

    struct parameters satdata;
    struct sf_params sparams;
};


int assign_ocrust_fun_3d(const gsl_vector * x, void *params, gsl_vector * f);
struct compo calc_ocrust3d_composition(double rhob_, double *guess,
        struct parameters satdata, struct sf_params sparams);

int assign_icrust_fun_4d(const gsl_vector * x, void *params, gsl_vector * f);
struct compo calc_icrust4d_composition(double rhob_, double *guess,
        struct parameters satdata, struct sf_params sparams);

double calc_muncl(struct parameters satdata, struct sf_params sparams, struct compo eq, double rhob_);
double calc_crust_ws_cell_energy_density(struct parameters satdata, struct sf_params sparams, struct compo eq,
        double rhob_); // total energy density in the inner-crust WS cell -amu*nB in MeV/fm^3
double calc_ngas_pressure(struct parameters satdata, double rhog_); // in MeV/fm^3
double calc_crust_ws_cell_pressure(struct parameters satdata, struct compo eq,
        double rhob_); // pressure in the inner-crust WS cell in MeV/fm^3
void print_state_crust(struct parameters satdata, struct sf_params sparams, struct compo eq,
        double rhob_, FILE *compo, FILE *eos);

#endif // H_ICRUST
