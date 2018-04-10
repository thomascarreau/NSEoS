#include "../nseos/nuclear_en.h"
#include "../nseos/observables.h"
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
    double guess_ic[4] = {100.,0.25,0.15,1.e-4}; // initial guess for rhob = 0.0003 /fm3

    double del_eq;
    double guess_core = 0.7;

    // parts of the energy density
    struct parameters satdata = ASSIGN_PARAM(satdata);
    struct sf_params sparams;

    fprintf(stderr, "%.2f %.2f %.2f\n", satdata.rhosat0, satdata.lasat0, satdata.ksat0);
    fprintf(stderr, "%.2f %.2f %.2f\n", satdata.jsym0, satdata.lsym0, satdata.ksym0);
    fprintf(stderr, "%.2f %.2f\n", satdata.barm, satdata.bardel);
    fprintf(stderr, "p = %d\n\n", P_SURF_TENSION);
    fprintf(stderr, "==============================================\n\n");

    sparams = fit_sf_params();
    double rhop;
    double enuc;
    struct hnm ngas;
    double epsg;
    double epsws_ic;
    struct hnm meta;
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
            rhop = (rhob-comp.rhog)*(1.-comp.del)/2./(1.-comp.rhog/comp.rho0);
            enuc = CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, comp.aa, comp.del, comp.rho0, rhop);
            ngas = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, comp.rhog, 1.);
            epsg = comp.rhog*ngas.enpernuc;
            epsws_ic = calc_crust_ws_cell_energy_density(comp.aa, comp.rho0, rhop, comp.rhog, enuc, epsg, rhob);

            del_eq = calc_core_eq_asym(rhob, &guess_core);
            if (guess_core != guess_core) // break if nan
                break;

            // calculation of the energy density in the cell in the core
            meta = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, rhob, del_eq);
            epsws_core = calc_core_ws_cell_energy_density(del_eq, meta, rhob);

            /* fprintf(stderr, "%g %g %g\n", rhob, epsws_ic, epsws_core); */
            if (epsws_core < epsws_ic) // crust-core transition
            {
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
