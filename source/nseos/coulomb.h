#ifndef H_COULOMB_EN
#define H_COULOMB_EN

#include "empirical.h"

double calc_coulomb_en(struct parameters satdata, double aa_, double ii_, double n0_, double np_); // including Coulomb screening
double calc_screening_derivative(struct parameters satdata, double aa_, double ii_, double n0_, double np_);

double calc_egas_energy_density(double ne_);
double calc_egas_chemical_potential(double ne_);
double calc_egas_pressure(double ne_); // in MeV/fm^3

double calc_ugas_energy_density(double nu_);
double calc_ugas_chemical_potential(double nu_);
double calc_ugas_pressure(double nu_); // in MeV/fm^3

double calc_lattice_pressure(struct parameters satdata, double aa_, double ii_, double n0_, double np_);

#endif // H_COULOMB_EN
