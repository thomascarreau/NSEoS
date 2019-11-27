#ifndef H_TOV
#define H_TOV

#include <stdio.h>

#include <gsl/gsl_multiroots.h>

double calc_dm(double rho_, double r_, double dr_);
double calc_dp(double rho_, double p_, double r_, double dr_, double m_);
double calc_dw(
    double rho_, double p_, double r_, double dr_, double m_, double w_);
double calc_dy(double rho_, double p_, double r_, double dr_, double m_,
    double y_, double vs2_);

double calc_moment_of_inertia(double r_, double w_);
double calc_normalized_moment_of_inertia_approx(double r_, double m_);
double calc_normalized_crustal_moment_of_inertia_approx(double r_, double m_,
    double i_over_mr2, double epst, double pt, double rcore);

double calc_tidal_love_number(double r_, double m_, double y_);
double calc_dimensionless_tidal_deformability(double r_, double m_, double k2_);

struct tov_solution {
  double mmax;
  double rhoc;
  double pc;
  double r;
  double i_over_mr2;
  double rcore;
  double mcore;
  double icrust_over_mr2;
  double k2;
  double lambda_dimless;
};

double get_observable_for_a_given_mass(
    double m, double mm, double mp, double om, double op);
double solve_tov_equation(int lines, double pt, FILE *eos,
    struct tov_solution *tovs_m, double fixed_m, FILE *tov);

int eval_observables_from_glitch_activity(int lines, double pt, double epst,
    FILE *eos, double g, double sigma_g, FILE *f_pulsar);

// Lambda1-Lambda2 relation ===================================================
struct rparams_tov {
  double m1;
  double mchirp;
};
struct Lambda {
  double s1;
  double s2;
};
// define the equation to solve to get m2
int equation_for_m2(const gsl_vector *x, void *params, gsl_vector *f);
// calculate m2 for given m1 and mchirp
double m2_for_m1_mchirp(const double mchirp, const double m1, double *guess);
// calculate dimensionless lambda1 and lambda2 for given m1 and mchirp
void dimensionless_lambda1_lambda2_for_m1_mchirp(const int lines,
    const char *eos, const double mchirp, const double m1, double *m2,
    struct Lambda *Lambdai);
// write in a file the dimensionless lamda1-lambda2 relation for a given of
// empirical parameters and mchirp
void dimensionless_lambda1_lambda2_relation(
    const int lines, const char *eos, const double mchirp, const char *outfile);

#endif // H_TOV
