#include <math.h>

#include "phyconst.h"
#include "nuclear_matter.h"
#include "coulomb.h"
#include "crust.h"
#include "core.h"
#include "modeling.h"
#include "eos.h"

int calc_equation_of_state(struct parameters satdata, double p, 
        struct transition_qtt *tqtt, double *epst,
        FILE *crust, FILE *core, FILE *eos)
{
    int lines = 0; 

    print_parameters(satdata);
    fprintf(stderr, "p = %g\n\n", p);
    fprintf(stderr, "==============================================\n\n");
    struct sf_params sparams = fit_sf_params(satdata, p);

    struct compo comp;
    double muncl = -1.; // sign of muncl is negative is the outer crust
    double guess_oc[3] = {60., 0.15, 0.1595}; // initial guess for the outer crust

    double nb = 1e-10;

    while(1)
    {
        comp = calc_ocrust3d_composition(nb, guess_oc, satdata, sparams);
        if (guess_oc[0] != guess_oc[0]) // exit if nan
            return lines;

        muncl = calc_muncl(satdata, sparams, comp, nb);
        if (muncl > 0.) // neutron drip -> transtion to inner crust
            break;

        print_state_crust(satdata, sparams, comp, nb, crust, eos);

        nb += nb/50.;

        lines += 1;
    }

    fprintf(stderr, "n_d = %g /fm^3\n\n", nb);
    fprintf(stderr, "==============================================\n\n");

    double epsws_ic;
    double epsws_core;
    double guess_ic[4] = {guess_oc[0], guess_oc[1], guess_oc[2], 1.e-4}; // initial guess for the inner crust
    struct core_compo ccomp;
    double guess_npecore = 0.7; // initial guess for the core

    while(1)
    {
        comp = calc_icrust4d_composition(nb, guess_ic, satdata, sparams);
        if (guess_ic[0] != guess_ic[0]) // exit if nan
        {
            if (epsws_core - epsws_ic < 1.e-3)
            {
                fprintf(stderr, "e_core - e_crust = %g\n", epsws_core - epsws_ic);
                break;
            }
            else
                return lines;
        }

        if (nb > 0.001)
        {
            // calculation of the energy density in the cell in the inner crust
            epsws_ic = calc_crust_ws_cell_energy_density(satdata, sparams, comp, nb);

            ccomp = calc_npecore_composition(nb, &guess_npecore, satdata);
            if (guess_npecore != guess_npecore) // exit if nan
            {
                if (epsws_core - epsws_ic < 1.e-3)
                {
                    fprintf(stderr, "e_core - e_crust = %g\n", epsws_core - epsws_ic);
                    break;
                }
                else
                    return lines;
            }

            // calculation of the energy density in the cell in the core
            epsws_core = calc_core_ws_cell_energy_density(satdata, ccomp, nb);

            if (epsws_core - epsws_ic < 0.) // crust-core transition
            {
                fprintf(stderr, "e_core < e_crust\n");
                break;
            }
        }

        print_state_crust(satdata, sparams, comp, nb, crust, eos);

        nb += 0.0001;

        lines += 1;
    }

    tqtt->nt = nb;
    tqtt->pt = calc_core_ws_cell_pressure(satdata, ccomp, nb);
    *epst = calc_core_ws_cell_energy_density(satdata, ccomp, nb); // needed to calculate Icrust

    fprintf(stderr, "n_t = %g /fm^3\n", tqtt->nt);
    fprintf(stderr, "P_t = %g MeV/fm^3\n\n", tqtt->pt);
    fprintf(stderr, "==============================================\n\n");

    double mueltot = 0.; // initializing

    while(1)
    {
        ccomp = calc_npecore_composition(nb, &guess_npecore, satdata);
        if (guess_npecore != guess_npecore) // exit if nan
            return lines;

        mueltot = calc_egas_chemical_potential(nb*(1.-ccomp.del)/2.);
        if (mueltot - MMU > 0.) // transition to npeu matter
            break;

        print_state_core(satdata, ccomp, nb, core, eos);

        nb += 0.001;

        lines += 1;
    }

    fprintf(stderr, "muons appear at %g /fm^3\n\n", nb);

    double guess_npeucore[2] = {guess_npecore, 1.e-5};
    struct hnm test_causality;
    double vs2 = 0.5;

    while (vs2 > 0. && vs2 < 1.)
    {
        ccomp = calc_npeucore_composition(nb, guess_npeucore, satdata);
        if (guess_npeucore[0] != guess_npeucore[0]) // break if nan; see q&d part in core.c
            break;

        print_state_core(satdata, ccomp, nb, core, eos);

        test_causality = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, nb, ccomp.del);
        vs2 = test_causality.vs2;

        nb += 0.01;

        lines += 1;
    }

    return lines;
}

void eval_transition_qtt(struct parameters satdata, double p,
        struct transition_qtt *tqtt, double *epst)
{
    struct sf_params sparams = fit_sf_params(satdata, p);
    struct compo comp;
    double guess_oc[3] = {60., 0.15, 0.9*satdata.rhosat0}; // initial guess for the outer crust
    double muncl = -1.; // sign of muncl is negative is the outer crust
    double nb = 1.e-6;

    while(1)
    {
        comp = calc_ocrust3d_composition(nb, guess_oc, satdata, sparams);
        if (guess_oc[0] != guess_oc[0]) // exit if nan
            return;

        muncl = calc_muncl(satdata, sparams, comp, nb);
        if (muncl > 0.) // neutron drip -> transtion to inner crust
            break;

        nb += nb/50.;
    }

    double guess_ic[4] = {guess_oc[0], guess_oc[1], guess_oc[2], 1.e-4}; // initial guess for the inner crust
    struct core_compo ccomp;
    double guess_npecore = 0.7; // initial guess for the core
    double epsws_ic;
    double epsws_core;

    while(1)
    {
        comp = calc_icrust4d_composition(nb, guess_ic, satdata, sparams);
        if (guess_ic[0] != guess_ic[0]) // break if nan
        {
            if (epsws_core - epsws_ic < 1.e-3)
            {
                fprintf(stderr, "e_core - e_crust = %g\n", epsws_core - epsws_ic);
                break;
            }
            else
                return;
        }

        // calculation of the energy density in the cell in the inner crust
        epsws_ic = calc_crust_ws_cell_energy_density(satdata, sparams, comp, nb);

        ccomp = calc_npecore_composition(nb, &guess_npecore, satdata);
        if (guess_npecore != guess_npecore) // break if nan
        {
            if (epsws_core - epsws_ic < 1.e-3)
            {
                fprintf(stderr, "e_core - e_crust = %g\n", epsws_core - epsws_ic);
                break;
            }
            else
                return;
        }

        // calculation of the energy density in the cell in the core
        epsws_core = calc_core_ws_cell_energy_density(satdata, ccomp, nb);

        if (epsws_core < epsws_ic) // crust-core transition
        {
            fprintf(stderr, "e_core < e_crust\n");
            break;
        }

        nb += 0.0001;
    }

    fprintf(stderr, "\n==============================================\n");
    fprintf(stderr, "==============================================\n\n");

    tqtt->nt = nb; 
    tqtt->pt = calc_core_ws_cell_pressure(satdata, ccomp, nb);
    *epst = calc_core_ws_cell_energy_density(satdata, ccomp, nb); // needed to calculate Icrust
}
