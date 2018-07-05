#ifndef H_TOV
#define H_TOV

#include <stdio.h>

double calc_dm(double rho_, double r_, double dr_);
double calc_dp(double rho_, double p_, double r_, double dr_, double m_);
double solve_tov_equation(int lines, double pt, double epst, FILE *eos, FILE *tov);

#endif // H_TOV
