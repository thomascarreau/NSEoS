#include <stdio.h>

#include "../nseos/nuclear_en.h"
#include "../nseos/observables.h"
#include "f_icrust4d.h"

#define CALC_NUCLEAR_EN (calc_ls_meta_model_nuclear_en)
#define ASSIGN_PARAM (assign_param_sly4)

static const int p_surf_tension = 3;
static const int taylor_exp_order = 2;

struct icrust_fun_4d calc_icrust_fun_4d(double aa_, double del_, double rho0_, double rhop_, double rhog_)
{
    struct icrust_fun_4d result;
    struct parameters satdata;
    double enuc;
    double epsa;
    double epsb;
    double epsr;
    double enuc_ap, enuc_am;
    double enuc_bp, enuc_bm;
    double enuc_rp, enuc_rm;
    double denucdaa;
    double denucddel;
    double denucdrho0;
    double muel;
    double dmu;
    struct hnm ngas;

    satdata = ASSIGN_PARAM(satdata);

    enuc = CALC_NUCLEAR_EN(satdata, p_surf_tension, taylor_exp_order, aa_, del_, rho0_, rhop_);
    epsa = 0.001;
    epsb = 0.0001;
    epsr = 0.0001;
    enuc_ap = CALC_NUCLEAR_EN(satdata, p_surf_tension, taylor_exp_order, aa_+epsa, del_, rho0_, rhop_);
    enuc_am = CALC_NUCLEAR_EN(satdata, p_surf_tension, taylor_exp_order, aa_-epsa, del_, rho0_, rhop_);
    enuc_bp = CALC_NUCLEAR_EN(satdata, p_surf_tension, taylor_exp_order, aa_, del_+epsb, rho0_, rhop_);
    enuc_bm = CALC_NUCLEAR_EN(satdata, p_surf_tension, taylor_exp_order, aa_, del_-epsb, rho0_, rhop_);
    enuc_rp = CALC_NUCLEAR_EN(satdata, p_surf_tension, taylor_exp_order, aa_, del_, rho0_+epsr, rhop_);
    enuc_rm = CALC_NUCLEAR_EN(satdata, p_surf_tension, taylor_exp_order, aa_, del_, rho0_-epsr, rhop_);
    denucdaa = (enuc_ap - enuc_am)/2./epsa; // 2 points derivatives
    denucddel = (enuc_bp - enuc_bm)/2./epsb;
    denucdrho0 = (enuc_rp - enuc_rm)/2./epsr;

    muel = calc_egas_chemical_potential(rhop_);
    dmu = calc_screening_derivative(satdata, aa_, del_, rho0_, rhop_);

    ngas = calc_meta_model_nuclear_matter(satdata, taylor_exp_order, rhog_, 1.);

    result.f_stability = denucdaa/aa_ - enuc/aa_/aa_;
    result.f_beta = denucddel*2./aa_ - muel - dmu - RMP + RMN;
    result.f_muneq = enuc/aa_ - (ngas.mun)*(1.-rhog_/rho0_) 
        + (1.-del_)/2.*(muel + dmu + RMP - RMN) - rhog_*(ngas.enpernuc)/rho0_ ;
    result.f_presseq = rho0_*rho0_*denucdrho0/aa_ - rhog_*ngas.mun + rhog_*ngas.enpernuc;

    return result;
}

int assign_icrust_fun_4d(const gsl_vector * x, void *params, gsl_vector * f)
{
    double rhop = ((struct rparams *) params)->rhop;

    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);
    const double x2 = gsl_vector_get (x, 2);
    const double x3 = gsl_vector_get (x, 3);

    struct parameters satdata;
    satdata = ASSIGN_PARAM(satdata);

    rhop = (rhop-x3)*(1.-x1)/2./(1.-x3/x2);

    struct icrust_fun_4d functs;
    functs = calc_icrust_fun_4d(x0, x1, x2, rhop, x3);

    const double y0 = functs.f_stability;
    const double y1 = functs.f_beta;
    const double y2 = functs.f_muneq;
    const double y3 = functs.f_presseq;

    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);
    gsl_vector_set(f, 2, y2);
    gsl_vector_set(f, 3, y3);

    return GSL_SUCCESS;
}

void print_state_icrust(gsl_multiroot_fsolver * s, double rhob_)
{
    double aa_eq, del_eq, rho0_eq, rhog_eq;
    double rhop_eq;
    double vws, rws;
    struct parameters satdata;
    double enuc;
    struct hnm ngas;
    double epsg;
    double epsws;
    double pressws;

    aa_eq = gsl_vector_get(s->x, 0);
    del_eq = gsl_vector_get(s->x, 1);
    rho0_eq = gsl_vector_get(s->x, 2);
    rhog_eq = gsl_vector_get(s->x, 3);
    rhop_eq = (rhob_-rhog_eq)*(1.-del_eq)/2./(1.-rhog_eq/rho0_eq);
    vws = aa_eq/(rhob_-rhog_eq)*(1.-rhog_eq/rho0_eq);
    rws = pow(3.*vws/4./PI,1./3.);

    satdata = ASSIGN_PARAM(satdata);
    enuc = CALC_NUCLEAR_EN(satdata, p_surf_tension, taylor_exp_order, aa_eq, del_eq, rho0_eq, rhop_eq);
    ngas = calc_meta_model_nuclear_matter(satdata, taylor_exp_order, rhog_eq, 1.);
    epsg = rhog_eq*ngas.enpernuc;
    epsws = calc_ws_cell_energy_density(aa_eq, rho0_eq, rhop_eq, rhog_eq, enuc, epsg, rhob_);
    pressws = calc_ws_cell_pressure(satdata, aa_eq, del_eq, rho0_eq, rhop_eq, rhog_eq, epsg, ngas.mun);

    printf ("%g %g %g %g %g %g %g %g %g\n", rhob_, aa_eq, del_eq, aa_eq*(1.-del_eq)/2., rho0_eq, rhog_eq, rws,
            epsws, pressws);
    /* printf("%g %g\n", rhob_, pressws); */
}

void calc_icrust4d(double rhob_, double *guess)
{
    double aa_new_guess;
    double del_new_guess;
    double rho0_new_guess;
    double rhog_new_guess;

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t iter = 0;

    double aa_old, aa_new, astep;
    double basym_old, basym_new, bstep;
    double rho0_old, rho0_new, rstep;
    double rhog_old, rhog_new, gstep;

    struct rparams p;
    p.rhop = rhob_; // modification in assign_icrust_fun_4d (DIRTY!)

    const size_t n = 4;
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector_set(x, 0, guess[0]);
    gsl_vector_set(x, 1, guess[1]);
    gsl_vector_set(x, 2, guess[2]);
    gsl_vector_set(x, 3, guess[3]);

    gsl_multiroot_function f = {&assign_icrust_fun_4d, n, &p};

    T = gsl_multiroot_fsolver_broyden;
    s = gsl_multiroot_fsolver_alloc(T,4);
    gsl_multiroot_fsolver_set(s, &f, x);

    do
    {
        /* print_state_icrust(s, rhob); */

        aa_old = gsl_vector_get (s->x, 0);
        basym_old = gsl_vector_get (s->x, 1);
        rho0_old = gsl_vector_get (s->x, 2);
        rhog_old = gsl_vector_get (s->x, 3);

        iter++;

        status = gsl_multiroot_fsolver_iterate(s);

        aa_new = gsl_vector_get (s->x, 0);
        basym_new = gsl_vector_get (s->x, 1);
        rho0_new = gsl_vector_get (s->x, 2);
        rhog_new = gsl_vector_get (s->x, 3);
        astep = aa_new - aa_old;
        bstep = basym_new - basym_old;
        rstep = rho0_new - rho0_old;
        gstep = rhog_new - rhog_old;

        if (status)
            break;

        // dirty backstepping
        while (gsl_vector_get (s->x, 0) < 0.) {
            astep = astep/4.;
            aa_new = aa_old + astep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, rho0_new);
            gsl_vector_set (x, 3, rhog_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            aa_new = gsl_vector_get (s->x, 0.);
        }
        while (gsl_vector_get (s->x, 1) < 0.1 || gsl_vector_get (s->x, 1) > 1.0) {
            bstep = bstep/4.;
            basym_new = basym_old + bstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, rho0_new);
            gsl_vector_set (x, 3, rhog_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            basym_new = gsl_vector_get (s->x, 1);
        }
        while (gsl_vector_get (s->x, 2) < 0.) {
            rstep = rstep/4.;
            rho0_new = rho0_old + rstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, rho0_new);
            gsl_vector_set (x, 3, rhog_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            rho0_new = gsl_vector_get (s->x, 2);
        }
        while (gsl_vector_get (s->x, 3) < 0. || gsl_vector_get (s->x, 3) > rhob_) {
            gstep = gstep/4.;
            rhog_new = rhog_old + gstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, rho0_new);
            gsl_vector_set (x, 3, rhog_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            rhog_new = gsl_vector_get (s->x, 3);
        }

        status = gsl_multiroot_test_residual (s->f, 9e-9);
    }

    while (status == GSL_CONTINUE && iter < 5000);

    if (iter < 4995)
        print_state_icrust(s, rhob_);

    aa_new_guess = gsl_vector_get(s->x, 0);
    del_new_guess = gsl_vector_get(s->x, 1);
    rho0_new_guess = gsl_vector_get(s->x, 2);
    rhog_new_guess = gsl_vector_get(s->x, 3);

    guess[0] = aa_new_guess;
    guess[1] = del_new_guess;
    guess[2] = rho0_new_guess;
    guess[3] = rhog_new_guess;

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
}
