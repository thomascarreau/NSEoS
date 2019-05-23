#include <math.h>

#include "phyconst.h"
#include "nuclear_matter.h"
#include "lepton.h"
#include "crust.h"
#include "core.h"
#include "modeling.h"
#include "eos.h"

int calc_zero_temperature_equation_of_state(
        struct parameters satdata, double p, 
        struct transition_qtt *tqtt, double *epst, int *hd_checker,
        FILE *crust, FILE *core, FILE *eos)
{
    int lines = 0; 

    print_parameters(satdata);
    fprintf(stderr, "p = %g\n\n", p);
    fprintf(stderr, "==============================================\n\n");
    struct sf_params sparams = fit_sf_params(satdata, p);

    double nb = 1e-10;
    double tt = 0.;
    char phase[] = "sol";
    double pressure;
    double pressure_sav = 0.;
    double vs2;
    struct hnm test_hd;

    // ============================= OUTER CRUST =============================

    struct compo comp;
    double muncl = -1.; // sign of muncl is negative is the outer crust
    double guess_oc[3] = {60., 0.15, 0.1595}; // initial guess for the ocrust

    while(1)
    {
        comp = calc_ocrust3d_composition(nb, tt, 
                phase, guess_oc, satdata, sparams);
        if (guess_oc[0] != guess_oc[0]) // exit if nan
        {
            *hd_checker = 1;
            return lines;
        }

        muncl = calc_muncl(satdata, sparams, comp, nb, tt, phase);
        if (muncl > 0.) // neutron drip -> transtion to inner crust
            break;

        pressure = calc_crust_ws_cell_pressure(satdata, sparams, comp, 
                nb, tt, phase);

        if (*hd_checker == 0 && pressure < pressure_sav)
        {
            fprintf(stderr, 
                    "HD CHECKER: dP/dnB < 0 (outer crust) ; nB = %g /fm^3\n\n", 
                    nb);
            *hd_checker = 1;
        }

        print_state_crust(satdata, sparams, comp, nb, tt, phase, crust, eos);

        pressure_sav = pressure;

        nb += nb/50.;

        lines += 1;
    }

    fprintf(stderr, "n_d = %g /fm^3\n\n", nb);
    fprintf(stderr, "==============================================\n\n");

    // ============================= INNER CRUST =============================

    double epsws_ic;
    double epsws_core;
    // initial guess for the icrust
    double guess_ic[4] = {guess_oc[0], guess_oc[1], guess_oc[2], 1.e-4}; 
    struct core_compo ccomp;
    double guess_npecore = 0.7; // initial guess for the core

    while(1)
    {
        comp = calc_icrust4d_composition(nb, tt, phase, 
                guess_ic, satdata, sparams);
        if (guess_ic[0] != guess_ic[0]) // exit if nan
        {
            if (nb > 0.001 
                    && epsws_core - epsws_ic < 1.e-3) // crust-core transition
            {
                fprintf(stderr, "e_core - e_crust = %g\n", 
                        epsws_core - epsws_ic);
                break;
            }
            else
            {
                *hd_checker = 1;
                return lines;
            }
        }

        if (nb > 0.001)
        {
            // calculation of the energy density in the cell in the inner crust
            epsws_ic = calc_crust_ws_cell_free_energy_density(
                    satdata, sparams, comp, nb, tt, phase);

            ccomp = calc_npecore_composition(nb, tt, &guess_npecore, satdata);
            if (guess_npecore != guess_npecore) // exit if nan
            {
                if (epsws_core - epsws_ic < 1.e-3) // crust-core transition
                {
                    fprintf(stderr, "e_core - e_crust = %g\n", 
                            epsws_core - epsws_ic);
                    break;
                }
                else
                {
                    *hd_checker = 1;
                    return lines;
                }
            }

            // calculation of the energy density in the cell in the core
            epsws_core = calc_core_ws_cell_free_energy_density(satdata, ccomp, 
                    nb, tt);

            if (epsws_core - epsws_ic < 0.) // crust-core transition
            {
                fprintf(stderr, "e_core < e_crust\n");
                break;
            }
        }

        pressure = calc_crust_ws_cell_pressure(satdata, sparams, 
                comp, nb, tt, phase);

        if (*hd_checker == 0 && pressure < pressure_sav)
        {
            fprintf(stderr, 
                    "HD CHECKER: dP/dnB < 0 (inner crust) ; nB = %g /fm^3\n\n",
                    nb);
            *hd_checker = 1;
        }

        print_state_crust(satdata, sparams, comp, nb, tt, phase, crust, eos);

        pressure_sav = pressure;

        nb += 0.0001;

        lines += 1;
    }

    tqtt->nt = nb;
    tqtt->pt = calc_core_ws_cell_pressure(satdata, ccomp, nb, tt);
    *epst = calc_core_ws_cell_free_energy_density(satdata, ccomp, nb, tt);

    fprintf(stderr, "n_t = %g /fm^3\n", tqtt->nt);
    fprintf(stderr, "P_t = %g MeV/fm^3\n\n", tqtt->pt);
    fprintf(stderr, "==============================================\n\n");

    // ============================== npe CORE ==============================

    double mueltot = 0.; // initializing

    while(1)
    {
        ccomp = calc_npecore_composition(nb, tt, &guess_npecore, satdata);
        if (guess_npecore != guess_npecore) // exit if nan
        {
            if (nb < 3.*satdata.rhosat0)
                *hd_checker = 1;
            return lines;
        }

        mueltot = calc_egas_chemical_potential(nb*(1.-ccomp.del)/2., 0.);
        if (mueltot - MMU > 0.) // transition to npeu matter
            break;

        vs2 = calc_squared_nucleon_sound_velocity(satdata, TAYLOR_EXP_ORDER, 
                nb, ccomp.del, tt);
        test_hd = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER,
                nb, ccomp.del, tt);

        if (vs2 < 0. || vs2 > 1.)
        {
            fprintf(stderr, 
                    "HD CHECKER: vs/c<0 or >1 (npe core) ; nB = %g /fm^3\n\n",
                    nb);
            if (nb < 3.*satdata.rhosat0)
                *hd_checker = 1;
            return lines;
        }

        pressure = calc_core_ws_cell_pressure(satdata, ccomp, nb, tt);

        if (*hd_checker == 0 && nb != tqtt->nt && pressure < pressure_sav)
        {
            fprintf(stderr, 
                    "HD CHECKER: dP/dnB < 0 (npe core) ; nB = %g /fm^3\n\n",
                    nb);
            *hd_checker = 1;
        }

        if (*hd_checker == 0 && test_hd.jsym < 0)
        {
            fprintf(stderr, 
                    "HD CHECKER: Jsym < 0 (npe core) ; nB = %g /fm^3\n\n",
                    nb);
            *hd_checker = 1;
        }

        print_state_core(satdata, ccomp, nb, tt, core, eos);

        pressure_sav = pressure;

        nb += 0.001;

        lines += 1;
    }

    fprintf(stderr, "muons appear at %g /fm^3\n\n", nb);

    // ============================== npeu CORE ==============================

    double guess_npeucore[2] = {guess_npecore, 1.e-5};
    vs2 = 0.5;

    while (vs2 > 0. && vs2 < 1.)
    {
        ccomp = calc_npeucore_composition(nb, tt, guess_npeucore, satdata);
        if (guess_npeucore[0] != guess_npeucore[0]) // exit if nan
        {
            if (nb < 3.*satdata.rhosat0)
                *hd_checker = 1;
            return lines;
        }

        pressure = calc_core_ws_cell_pressure(satdata, ccomp, nb, tt);

        if (*hd_checker == 0 && pressure < pressure_sav)
        {
            fprintf(stderr,
                    "HD CHECKER: dP/dnB < 0 (core npeu) ; nB = %g /fm^3\n\n",
                    nb);
            *hd_checker = 1;
        }

        test_hd = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, 
                nb, ccomp.del, tt);

        if (*hd_checker == 0 && test_hd.jsym < 0)
        {
            fprintf(stderr,
                    "HD CHECKER: Jsym < 0 (core npeu) ; nB = %g /fm^3\n\n",
                    nb);
            *hd_checker = 1;
        }

        print_state_core(satdata, ccomp, nb, tt, core, eos);

        pressure_sav = pressure;

        vs2 = calc_squared_nucleon_sound_velocity(satdata, TAYLOR_EXP_ORDER,
                nb, ccomp.del, tt);

        nb += 0.005;

        lines += 1;
    }

    return lines;
}

void eval_transition_qtt(struct parameters satdata, double p,
        struct transition_qtt *tqtt, double *epst)
{
    struct sf_params sparams = fit_sf_params(satdata, p);
    struct compo comp;
    // initial guess for the ocrust
    double guess_oc[3] = {60., 0.15, 0.9*satdata.rhosat0};
    double muncl = -1.; // sign of muncl is negative is the outer crust
    double nb = 1.e-6;
    double tt = 0.;
    char phase[] = "sol";

    while(1)
    {
        comp = calc_ocrust3d_composition(nb, tt, phase, 
                guess_oc, satdata, sparams);
        if (guess_oc[0] != guess_oc[0]) // exit if nan
            return;

        muncl = calc_muncl(satdata, sparams, comp, nb, tt, phase);
        if (muncl > 0.) // neutron drip -> transtion to inner crust
            break;

        nb += nb/50.;
    }

    // initial guess for the icrust
    double guess_ic[4] = {guess_oc[0], guess_oc[1], guess_oc[2], 1.e-4};
    struct core_compo ccomp;
    double guess_npecore = 0.7; // initial guess for the core
    double epsws_ic;
    double epsws_core;

    while(1)
    {
        comp = calc_icrust4d_composition(nb, tt, phase, 
                guess_ic, satdata, sparams);
        if (guess_ic[0] != guess_ic[0]) // exit if nan
        {
            if (nb > 0.001 
                    && epsws_core - epsws_ic < 1.e-3) // crust-core transition
            {
                fprintf(stderr, "e_core - e_crust = %g\n", 
                        epsws_core - epsws_ic);
                break;
            }
            else
            {
                return;
            }
        }

        // calculation of the energy density in the cell in the inner crust
        epsws_ic = calc_crust_ws_cell_free_energy_density(satdata, sparams, 
                comp, nb, tt, phase);

        ccomp = calc_npecore_composition(nb, tt, &guess_npecore, satdata);
        if (guess_npecore != guess_npecore) // break if nan
        {
            if (nb > 0.001
                    && epsws_core - epsws_ic < 1.e-3)
            {
                fprintf(stderr, "e_core - e_crust = %g\n", 
                        epsws_core - epsws_ic);
                break;
            }
            else
            {
                return;
            }
        }

        // calculation of the energy density in the cell in the core
        epsws_core = calc_core_ws_cell_free_energy_density(satdata, ccomp, 
                nb, tt);

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
    tqtt->pt = calc_core_ws_cell_pressure(satdata, ccomp, nb, tt);
    *epst = calc_core_ws_cell_free_energy_density(satdata, ccomp, nb, tt);
}
