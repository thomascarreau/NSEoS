#include "../nseos/nuclear_en.h"
#include "../nseos/observables.h"

#include "choices.h"
#include "f_icrust4d.h"
#include "f_core.h"

int main(int argc, char* argv[])
{
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

        comp = calc_icrust4d_composition(rhob, guess_ic);
        if (guess_ic[0] != guess_ic[0]) // break if nan
            break;

        if (irhob > 399)
        {
            // calculation of the energy density in the cell in the inner crust
            rhop = (rhob-comp.rhog)*(1.-comp.del)/2./(1.-comp.rhog/comp.rho0);
            enuc = CALC_NUCLEAR_EN(satdata, p_surf_tension, taylor_exp_order, comp.aa, comp.del, comp.rho0, rhop);
            ngas = calc_meta_model_nuclear_matter(satdata, taylor_exp_order, comp.rhog, 1.);
            epsg = comp.rhog*ngas.enpernuc;
            epsws_ic = calc_crust_ws_cell_energy_density(comp.aa, comp.rho0, rhop, comp.rhog, enuc, epsg, rhob);

            del_eq = calc_core_eq_asym(rhob, &guess_core);
            if (guess_core != guess_core) // break if nan
                break;

            // calculation of the energy density in the cell in the core
            meta = calc_meta_model_nuclear_matter(satdata, taylor_exp_order, rhob, del_eq);
            epsws_core = calc_core_ws_cell_energy_density(del_eq, meta, rhob);

            if (epsws_core < epsws_ic) // crust-core transition
            {
                break;
            }
        }

        print_state_icrust(comp, rhob, mycompo, myeos);
    }

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

    printf("\\o/\n");

    return 0;
}
