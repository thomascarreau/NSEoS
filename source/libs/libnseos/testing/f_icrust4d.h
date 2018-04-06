#ifndef H_F_ICRUST4D
#define H_F_ICRUST4D

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

/* #define CALC_NUCLEAR_EN (calc_ls_meta_model_nuclear_en) */
/* #define ASSIGN_PARAM (assign_param_sly4) */
/* static const int p_surf_tension = 3; */
/* static const int taylor_exp_order = 2; */

struct icrust_fun_4d
{
    double f_stability;
    double f_beta;
    double f_muneq;
    double f_presseq;
};
struct ic_compo
{
    double aa;
    double del;
    double rho0;
    double rhog;
};
struct icrust_fun_4d calc_icrust_fun_4d(double aa_, double del_, double rho0_, double rhop_, double rhog_);
struct rparams_crust
{
    double rhop;
};
int assign_icrust_fun_4d(const gsl_vector * x, void *params, gsl_vector * f);
struct ic_compo calc_icrust4d_composition(double rhob_, double *guess);
void print_state_icrust(struct ic_compo eq, double rhob_);

#endif // H_F_ICRUST4D
