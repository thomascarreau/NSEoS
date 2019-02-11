#ifndef H_TOV
#define H_TOV

#include <stdio.h>

double calc_dm(double rho_, double r_, double dr_);
double calc_dp(double rho_, double p_, double r_, double dr_, double m_);
double calc_dy(double rho_, double p_, double r_, double dr_, double m_, 
        double y_, double vs2_);
double calc_tidal_love_number(double r_, double m_, double y_);
double calc_dimensionless_tidal_deformability(double r_, double m_, 
        double k2_);
struct tov_solution
{
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
double calc_normalized_moment_of_inertia(double r_, double m_);
double calc_normalized_crustal_moment_of_inertia(double r_, double m_, 
        double i_over_mr2, double epst, double pt, double rcore);
double get_observable_for_a_given_mass(double m, double mm, double mp, 
        double om, double op);
double solve_tov_equation(int lines, double pt, double epst, FILE *eos, 
        struct tov_solution *tovs_m, double fixed_m, FILE *tov);

int eval_observables_from_glitch_activity(
        int lines, double pt, double epst, FILE *eos, double g, double sigma_g,
        FILE *f_pulsar);

#endif // H_TOV
