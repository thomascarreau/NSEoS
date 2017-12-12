#include <stdio.h>
#include <math.h>

#include "../nseos/ws_cell.h"

int main(void)
{
    struct parameters satdata0;
    double rhob;
    /* double iaa; */
    double aa;
    double idel;
    double del;
    double rho0;
    double rhop;
    double enuc;
    double evac;
    double ebulk;
    double esurf;
    double ecoul;

    satdata0 = assign_param(satdata0);
    rhob = 1.e-10;
    aa = 200.;
    /* del = 0.2; */

    /* for(iaa = 40; iaa < 200; iaa++) */
    for(idel = -90; idel < 91; idel++)
    {
        del = idel/100.;
        /* aa = iaa; */
        rho0 = satdata0.rhosat0;
        rhop = (1.-del)/2.*rhob;
        enuc = calc_enuc(satdata0, aa, del, rho0, rhop, 0.);
        evac = calc_evac(satdata0, aa, del, rho0, 0.);
        ebulk = calc_bulk_energy(satdata0, aa, del, rho0, 0.);
        esurf = calc_surface_energy(satdata0, aa, del, rho0, 0.);
        ecoul = calc_coulomb_energy(aa, del, rho0);
        printf("%g %g %g %g %g %g %g\n", aa, del, enuc/aa, 
                evac/aa, ebulk/aa, esurf/pow(aa,2./3.), ecoul/aa);
    }

    return 0;
}
