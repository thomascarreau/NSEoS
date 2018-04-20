#ifndef H_CORE
#define H_CORE

#include <gsl/gsl_vector.h>

#include "empirical.h"

double calc_core_fun(struct parameters satdata, double del_, double rhob_);
struct rparams_core
{
    double rhob;
    struct parameters satdata;
};
int assign_core_fun(const gsl_vector * x, void *params, gsl_vector * f);
double calc_core_eq_asym(double rhob_, double *guess, struct parameters satdata);
double calc_core_ws_cell_energy_density(struct parameters satdata, double del_eq_,
        double nb_); // total energy density in the core WS cell -amu*nB in MeV/fm^3
double calc_core_ws_cell_pressure(struct parameters satdata, double del_eq_,
        double nb_); // pressure in the inner-crust WS cell in MeV/fm^3
void print_state_core(struct parameters satdata, double del_eq_, double rhob_, FILE *eos);

#endif // H_CORE
