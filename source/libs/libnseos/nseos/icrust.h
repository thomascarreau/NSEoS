#ifndef H_ICRUST
#define H_ICRUST

#include <gsl/gsl_vector.h>

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
struct icrust_fun_4d calc_icrust_fun_4d(double aa_, double del_, double rho0_, double rhop_, double rhog_,
        struct sf_params sparams);
struct rparams_crust
{
    double rhop;

    double sigma0; // needed to avoid the fit at each iteration
    double b;
};
int assign_icrust_fun_4d(const gsl_vector * x, void *params, gsl_vector * f);
struct ic_compo calc_icrust4d_composition(double rhob_, double *guess,
        struct sf_params sparams);
double calc_crust_ws_cell_energy_density(struct parameters satdata, struct sf_params sparams, struct ic_compo eq,
        double rhob_); // total energy density in the inner-crust WS cell -amu*nB in MeV/fm^3
double calc_ngas_pressure(struct parameters satdata, double rhog_); // in MeV/fm^3
double calc_crust_ws_cell_pressure(struct parameters satdata, struct ic_compo eq,
        double rhob_); // pressure in the inner-crust WS cell in MeV/fm^3
void print_state_icrust(struct ic_compo eq, struct sf_params sparams, double rhob_, FILE *compo, FILE *eos);

#endif // H_ICRUST
