#include <stdio.h>

#include "../nseos/binding_en.h"

int main(int argc, char* argv[])
{
    FILE *fout_plot;
    struct parameters satdata0;
    double aa;
    double del;
    int rho0_;
    double rho0;
    double evac;
    double ebulk;
    double esurf;
    double ecoul;

    satdata0 = assign_param(satdata0);
    aa = 70.;
    del = 0.;

    fout_plot = fopen(argv[1], "a+");

    for(rho0_ = 90; rho0_ < 154; rho0_++)
    {
        rho0 = rho0_/1000.;
        printf("%g\n", rho0);
        evac = calc_evac(satdata0, aa, del, rho0);
        ebulk = calc_bulk_energy(satdata0, aa, del, rho0);
        esurf = calc_surface_energy(satdata0, aa, del, rho0);
        ecoul = calc_coulomb_energy(aa, del, rho0);
        fprintf(fout_plot, "%g %g %g %g %g\n", rho0, evac/aa, ebulk/aa, esurf/aa, ecoul/aa);
        /* printf("%g %g %g %g %g\n", rho0, evac/aa, ebulk/aa, esurf/aa, ecoul/aa); */
    }

    fclose(fout_plot);

    return 0;
}
