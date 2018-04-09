#ifndef H_F_CORE
#define H_F_CORE

#include <gsl/gsl_vector.h>

double calc_core_fun(double del_, double rhob_);
struct rparams_core
{
    double rhob;
};
int assign_core_fun(const gsl_vector * x, void *params, gsl_vector * f);
double calc_core_eq_asym(double rhob_, double *guess);
void print_state_core(double del_eq_, double rhob_, FILE *eos);

#endif // H_F_CORE
