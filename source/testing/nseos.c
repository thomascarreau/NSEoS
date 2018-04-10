#include "../nseos/nuclear_en.h"
#include "../nseos/icrust.h"
#include "../nseos/core.h"
#include "../nseos/modeling.h"

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "ERROR: syntax is './nseos compo.out eos.out'\n");
        exit(0);
    }

    FILE *mycompo;
    FILE *myeos;

    mycompo = fopen(argv[1],"w+");
    myeos = fopen(argv[2],"w+");

    int irhob;
    double rhob;

    struct ic_compo comp;
    double guess_ic[4] = {100.,0.25,0.15,1.e-4}; // initial guess for the inner crust
    double del_eq;
    double guess_core = 0.7; // initial guess for the core

    struct parameters satdata = ASSIGN_PARAM(satdata);
    print_parameters(satdata);
    fprintf(stderr, "p = %d\n\n", P_SURF_TENSION);
    fprintf(stderr, "==============================================\n\n");

    struct sf_params sparams = fit_sf_params();

    double epsws_ic;
    double epsws_core;

    for(irhob = 3; irhob < 1001; irhob ++)
    {
        rhob = irhob/10000.;

        comp = calc_icrust4d_composition(rhob, guess_ic, sparams);
        if (guess_ic[0] != guess_ic[0]) // break if nan
            break;

        if (irhob > 399)
        {
            // calculation of the energy density in the cell in the inner crust
            epsws_ic = calc_crust_ws_cell_energy_density(satdata, sparams, comp, rhob);

            del_eq = calc_core_eq_asym(rhob, &guess_core);
            if (guess_core != guess_core) // break if nan
                break;

            // calculation of the energy density in the cell in the core
            epsws_core = calc_core_ws_cell_energy_density(satdata, del_eq, rhob);

            if (epsws_core < epsws_ic) // crust-core transition
            {
                fprintf(stderr, "e_core < e_crust\n");
                break;
            }
        }

        print_state_icrust(comp, sparams, rhob, mycompo, myeos);
    }

    fprintf(stderr, "n_t = %g /fm^3\n\n", rhob);
    fprintf(stderr, "==============================================\n\n");

    do
    {
        del_eq = calc_core_eq_asym(rhob, &guess_core);
        if (guess_core != guess_core) // break if nan
            break;

        print_state_core(del_eq, rhob, myeos);

        rhob += 0.001;
    }
    while (rhob < 0.5);

    fclose(mycompo);
    fclose(myeos);

    fprintf(stderr, "\\o/\n");

    return 0;
}
