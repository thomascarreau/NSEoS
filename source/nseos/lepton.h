#ifndef H_LEPTON
#define H_LEPTON

double calc_egas_free_energy_density(double ne_, double tt_);
double calc_egas_chemical_potential(double ne_, double tt_);
double calc_egas_pressure(double ne_, double tt_); // in MeV/fm^3

double calc_ugas_energy_density(double nu_);
double calc_ugas_chemical_potential(double nu_);
double calc_ugas_pressure(double nu_); // in MeV/fm^3

#endif // H_LEPTON
