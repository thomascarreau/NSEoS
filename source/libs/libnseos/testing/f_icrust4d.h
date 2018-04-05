#ifndef H_F_ICRUST4D
#define H_F_ICRUST4D

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

struct icrust_fun_4d
{
    double f_stability;
    double f_beta;
    double f_muneq;
    double f_presseq;
};
struct icrust_fun_4d calc_icrust_fun_4d(double aa_, double del_, double rho0_, double rhop_, double rhog_);
struct rparams
{
    double rhop;
};
int assign_icrust_fun_4d(const gsl_vector * x, void *params, gsl_vector * f);
void print_state_icrust(gsl_multiroot_fsolver * s, double rhob_);
void calc_icrust4d(double rhob_, double *guess);

#endif // H_F_ICRUST4D
