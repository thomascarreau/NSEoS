#ifndef H_CORE
#define H_CORE

#include <gsl/gsl_vector.h>

#include "empirical.h"

struct core_fun
{
    double f_beta;
    double f_mueq;
};
struct core_compo
{
    double del;
    double rhou;
};
struct core_fun calc_core_fun(struct parameters satdata, double del_, double rhou_, double rhob_);
struct rparams_core
{
    double rhob;
    struct parameters satdata;
};

int assign_npecore_fun(const gsl_vector * x, void *params, gsl_vector * f);
struct core_compo calc_npecore_composition(double rhob_, double *guess, struct parameters satdata);

int assign_npeucore_fun(const gsl_vector * x, void *params, gsl_vector * f);
struct core_compo calc_npeucore_composition(double rhob_, double *guess, struct parameters satdata);

double calc_core_ws_cell_energy_density(struct parameters satdata, struct core_compo eq,
        double nb_); // total energy density in the core WS cell -amu*nB in MeV/fm^3
double calc_core_ws_cell_pressure(struct parameters satdata, struct core_compo eq,
        double nb_); // pressure in the inner-crust WS cell in MeV/fm^3
void print_state_core(struct parameters satdata, struct core_compo eq, 
        double rhob_, FILE *core, FILE *eos);

#endif // H_CORE
