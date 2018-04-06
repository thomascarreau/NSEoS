#ifndef H_F_CORE
#define H_F_CORE

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

/* #define ASSIGN_PARAM (assign_param_sly4) */
/* static const int taylor_exp_order = 2; */

double calc_core_fun(double del_, double rhob_);
struct rparams_core
{
    double rhob;
};
int assign_core_fun(const gsl_vector * x, void *params, gsl_vector * f);
double calc_core_eq_asym(double rhob_, double *guess);
void print_state_core(double del_eq_, double rhob_);

#endif // H_F_CORE
