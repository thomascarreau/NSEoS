#include "../nseos/nuclear_en.h"
#include "../nseos/icrust.h"
#include "../nseos/core.h"
#include "../nseos/modeling.h"

int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        fprintf(stderr, "ERROR: syntax is './nseos set.in compo.out eos.out'\n");
        return 1;
    }

    FILE *mycompo;
    FILE *myeos;

    mycompo = fopen(argv[2],"w+");
    myeos = fopen(argv[3],"w+");

    double rhob;

    struct ic_compo comp;
    double guess_ic[4] = {100.,0.25,0.15,1.e-4}; // initial guess for the inner crust
    double del_eq;
    double guess_core = 0.7; // initial guess for the core

    struct parameters satdata = assign_param(argv[1]);
    print_parameters(argv[1], satdata);
    fprintf(stderr, "p = %d\n\n", P_SURF_TENSION);
    fprintf(stderr, "==============================================\n\n");

    struct sf_params sparams = fit_sf_params(satdata);

    double epsws_ic;
    double epsws_core;
    int transition = 0;

    rhob = 0.0003;

    while(1)
    {
        comp = calc_icrust4d_composition(rhob, guess_ic, satdata, sparams);
        if (guess_ic[0] != guess_ic[0]) // break if nan
            break;

        if (rhob > 0.04)
        {
            // calculation of the energy density in the cell in the inner crust
            epsws_ic = calc_crust_ws_cell_energy_density(satdata, sparams, comp, rhob);

            del_eq = calc_core_eq_asym(rhob, &guess_core, satdata);
            if (guess_core != guess_core) // break if nan
                break;

            // calculation of the energy density in the cell in the core
            epsws_core = calc_core_ws_cell_energy_density(satdata, del_eq, rhob);

            if (epsws_core < epsws_ic) // crust-core transition
            {
                fprintf(stderr, "e_core < e_crust\n");
                transition = 1; // transition occurs
                break;
            }
        }

        print_state_icrust(comp, satdata, sparams, rhob, mycompo, myeos);

        rhob += 0.0001;
    }

    if (transition == 0)
        fprintf(stderr, "e_core - e_crust = %g MeV/fm^3\n", epsws_core - epsws_ic);

    fprintf(stderr, "n_t = %g /fm^3\n", rhob);
    fprintf(stderr, "P_t = %g MeV/fm^3\n\n", calc_core_ws_cell_pressure(satdata, del_eq, rhob));
    fprintf(stderr, "==============================================\n\n");

    do
    {
        del_eq = calc_core_eq_asym(rhob, &guess_core, satdata);
        if (guess_core != guess_core) // break if nan
            break;

        print_state_core(satdata, del_eq, rhob, myeos);

        rhob += 0.005;
    }
    while (rhob < 6.*satdata.rhosat0);

    fclose(mycompo);
    fclose(myeos);

    fprintf(stderr, "\\o/\n");

    return 0;
}
