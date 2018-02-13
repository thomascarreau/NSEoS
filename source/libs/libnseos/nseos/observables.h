#ifndef H_OBSERVABLES
#define H_OBSERVABLES

// ENERGY
double calc_ws_cell_energy_density(double aa_, double n0_, double np_, double ng_, double enuc, double epsg,
        double nb_); // total energy density in the WS cell -amu*nB in MeV/fm^3

// PRESSURE
double calc_egas_pressure(double np_); // in MeV/fm^3
double calc_lattice_pressure(struct parameters satdata, double aa_, double ii_, double n0_, double np_); // in MeV/fm^3
double calc_ngas_pressure(double ng_, double epsg, double mug); // in MeV/fm^3
double calc_ws_cell_pressure(struct parameters satdata, double aa_, double ii_, double n0_, double np_,
        double ng_, double epsg, double mug); // in MeV/fm^3

#endif // H_OBSERVABLES
