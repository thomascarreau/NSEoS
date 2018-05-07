#include <math.h>

#include "../../nseos/nuclear_en.h"
#include "../../nseos/crust.h"
#include "../../nseos/core.h"
#include "../../nseos/modeling.h"

int main(int argc, char* argv[])
{
    if (argc != 5)
    {
        fprintf(stderr, "ERROR: syntax is './nseos set.in crust.out core.out eos.out'\n");
        return 1;
    }

    FILE *mycrust;
    FILE *mycore;
    FILE *myeos;

    mycrust = fopen(argv[2],"w+");
    mycore = fopen(argv[3],"w+");
    myeos = fopen(argv[4],"w+");

    struct parameters satdata = assign_param(argv[1]);
    print_parameters(satdata);
    fprintf(stderr, "p = %d\n\n", P_SURF_TENSION);
    fprintf(stderr, "==============================================\n\n");
    struct sf_params sparams = fit_sf_params(satdata);

    struct compo comp;
    double muncl = -1.; // sign of muncl is negative is the outer crust
    double guess_oc[3] = {60., 0.15, 0.1595}; // initial guess for the outer crust

    double nb = 1e-10;

    while(1)
    {
        comp = calc_ocrust3d_composition(nb, guess_oc, satdata, sparams);
        if (guess_oc[0] != guess_oc[0]) // break if nan
        {
            fprintf(stderr, "%g\n", nb); // debug
            break;
        }

        muncl = calc_muncl(satdata, sparams, comp, nb);
        if (muncl > 0.)
            break;

        print_state_crust(satdata, sparams, comp, nb, mycrust, myeos);

        nb += nb/50.;
    }

    fprintf(stderr, "n_d = %g /fm^3\n\n", nb);
    fprintf(stderr, "==============================================\n\n");

    double epsws_ic;
    double epsws_core;
    int transition = 0;
    double guess_ic[4] = {guess_oc[0], guess_oc[1], guess_oc[2], 1.e-4}; // initial guess for the inner crust
    struct core_compo ccomp;
    double guess_npecore = 0.7; // initial guess for the core

    while(1)
    {
        comp = calc_icrust4d_composition(nb, guess_ic, satdata, sparams);
        if (guess_ic[0] != guess_ic[0]) // break if nan
            break;

        if (nb > 0.04)
        {
            // calculation of the energy density in the cell in the inner crust
            epsws_ic = calc_crust_ws_cell_energy_density(satdata, sparams, comp, nb);

            ccomp = calc_npecore_composition(nb, &guess_npecore, satdata);
            if (guess_npecore != guess_npecore) // break if nan
                break;

            // calculation of the energy density in the cell in the core
            epsws_core = calc_core_ws_cell_energy_density(satdata, ccomp, nb);

            if (epsws_core < epsws_ic) // crust-core transition
            {
                fprintf(stderr, "e_core < e_crust\n");
                transition = 1; // transition occurs
                break;
            }
        }

        print_state_crust(satdata, sparams, comp, nb, mycrust, myeos);

        nb += 0.0001;
    }

    if (transition == 0)
        fprintf(stderr, "e_core - e_crust = %g MeV/fm^3\n", epsws_core - epsws_ic);

    fprintf(stderr, "n_t = %g /fm^3\n", nb);
    fprintf(stderr, "P_t = %g MeV/fm^3\n\n", calc_core_ws_cell_pressure(satdata, ccomp, nb));
    fprintf(stderr, "==============================================\n\n");

    double mueltot = 0.; // initializing

    while(1)
    {
        ccomp = calc_npecore_composition(nb, &guess_npecore, satdata);
        if (guess_npecore != guess_npecore) // break if nan
            break;

        mueltot = calc_egas_chemical_potential(nb*(1.-ccomp.del)/2.);
        if (mueltot - MMU > 0.)
            break;

        print_state_core(satdata, ccomp, nb, mycore, myeos);

        nb += 0.001;
    }

    fprintf(stderr, "muons appear at %g /fm^3\n\n", nb);
    fprintf(stderr, "==============================================\n\n");

    double guess_npeucore[2] = {guess_npecore, 1.e-5};

    while (nb < 4.*satdata.rhosat0)
    {
        ccomp = calc_npeucore_composition(nb, guess_npeucore, satdata);
        if (guess_npeucore[0] != guess_npeucore[0]) // break if nan
            break;

        print_state_core(satdata, ccomp, nb, mycore, myeos);

        nb += 0.005;
    }

    fclose(mycrust);
    fclose(mycore);
    fclose(myeos);

    fprintf(stderr, "\\o/\n");

    return 0;
}
