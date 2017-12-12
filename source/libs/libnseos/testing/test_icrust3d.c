#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "../nseos/ws_cell.h"

// ==================== FUNCTIONS ====================

struct icrust_fun_3d
{
    double f_stability;
    double f_beta;
    double f_muneq;
};

struct icrust_fun_3d calc_icrust_fun_3d(double aa_, double del_, double rhop_, double rhog_);
struct icrust_fun_3d calc_icrust_fun_3d(double aa_, double del_, double rhop_, double rhog_)
{
    struct icrust_fun_3d result;
    struct parameters satdata;
    double rho0;
    double drho0ddel;
    double enuc;
    double epsa;
    double epsb;
    double enuc_ap, enuc_am;
    double enuc_bp, enuc_bm;
    double denucdaa;
    double denucddel;
    struct gas egas;
    struct hnm ngas;

    // One needs to calculate rho0 in this case...
    satdata = assign_param(satdata);
    rho0 = calc_anm_saturation_density(satdata, del_);
    drho0ddel = -6.*satdata.rhosat0
        *satdata.ksat0*satdata.lsym0*del_/pow(satdata.ksat0 + satdata.ksym0*del_*del_,2.);

    enuc = calc_enuc(satdata, aa_, del_, rho0, rhop_);
    epsa = 0.001;
    epsb = 0.0001;
    enuc_ap = calc_enuc(satdata, aa_+epsa, del_, rho0, rhop_);
    enuc_am = calc_enuc(satdata, aa_-epsa, del_, rho0, rhop_);
    enuc_bp = calc_enuc(satdata, aa_, del_+epsb, rho0, rhop_);
    enuc_bm = calc_enuc(satdata, aa_, del_-epsb, rho0, rhop_);
    denucdaa = (enuc_ap - enuc_am)/2./epsa; // 2 points derivatives
    denucddel = (enuc_bp - enuc_bm)/2./epsb;

    egas = calc_egas(rhop_);
    ngas = calc_hnm(satdata, rhog_, 1.);

    result.f_stability = denucdaa/aa_ - enuc/aa_/aa_;
    result.f_beta = denucddel*2./aa_ - egas.mutot - rmp + rmn
        + 2./rho0/rho0*drho0ddel*(rhog_*ngas.enpernuc - rhog_*ngas.mug);
    result.f_muneq = denucdaa - ngas.mug*(1.-rhog_/rho0) 
        + (1.-del_)/2.*(egas.mutot + rmp - rmn) - rhog_*ngas.enpernuc/rho0;

    return result;
}

struct rparams
{
    double rhop;
};

int assign_icrust_fun_3d(const gsl_vector * x, void *params, gsl_vector * f);
int assign_icrust_fun_3d(const gsl_vector * x, void *params, gsl_vector * f)
{
    double rhop = ((struct rparams *) params)->rhop;

    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);
    const double x2 = gsl_vector_get (x, 2);

    struct parameters satdata;
    double rho0;
    satdata = assign_param(satdata);
    rho0 = calc_anm_saturation_density(satdata, x1);
    rhop = (rhop-x2)*(1.-x1)/2./(1.-x2/rho0);

    struct icrust_fun_3d functs;
    functs = calc_icrust_fun_3d(x0, x1, rhop, x2);

    const double y0 = functs.f_stability;
    const double y1 = functs.f_beta;
    const double y2 = functs.f_muneq;

    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);
    gsl_vector_set(f, 2, y2);

    return GSL_SUCCESS;
}

void print_state_icrust(gsl_multiroot_fsolver * s, double rhob_);
void print_state_icrust(gsl_multiroot_fsolver * s, double rhob_)
{
    double aa_eq, del_eq, rhop_eq, rhog_eq;
    double zz_eq;
    struct parameters satdata;
    double rho0_eq;
    double epsws;

    aa_eq = gsl_vector_get(s->x, 0);
    del_eq = gsl_vector_get(s->x, 1);
    rhop_eq = rhob_*(1.-del_eq)/2.;
    rhog_eq = gsl_vector_get(s->x, 2);
    zz_eq = aa_eq*(1.-del_eq)/2.;
    satdata = assign_param(satdata);
    rho0_eq = calc_anm_saturation_density(satdata, del_eq);
    epsws = calc_epsws(satdata, aa_eq, del_eq, rho0_eq, rhog_eq, rhop_eq, rhob_);
    printf ("%g %g %g %g %g %g %g\n", aa_eq, del_eq, zz_eq, rhop_eq/rhob_, 
            rho0_eq, rhog_eq/rhob_, epsws);
}

// ==================== MAIN ====================

int main(void)
{
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t iter = 0;

    double aa_old, aa_new, astep;
    double basym_old, basym_new, bstep;
    double rhog_old, rhog_new, gstep;

    double rhob;
    rhob = 2.e-3;

    struct rparams p;
    p.rhop = rhob; // modificiation in assign_icrust_fun_3d (DIRTY!)

    double x_init[3] = {100., 0.3, 1.e-6};
    const size_t n = 3;
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector_set(x, 0, x_init[0]);
    gsl_vector_set(x, 1, x_init[1]);
    gsl_vector_set(x, 2, x_init[2]);

    gsl_multiroot_function f = {&assign_icrust_fun_3d, n, &p};

    T = gsl_multiroot_fsolver_dnewton;
    s = gsl_multiroot_fsolver_alloc(T,3);
    gsl_multiroot_fsolver_set(s, &f, x);

    do
    {
        print_state_icrust(s, rhob);

        aa_old = gsl_vector_get (s->x, 0);
        basym_old = gsl_vector_get (s->x, 1);
        rhog_old = gsl_vector_get (s->x, 2);

        iter++;

        status = gsl_multiroot_fsolver_iterate(s);

        aa_new = gsl_vector_get (s->x, 0);
        basym_new = gsl_vector_get (s->x, 1);
        rhog_new = gsl_vector_get (s->x, 2);
        astep = aa_new - aa_old;
        bstep = basym_new - basym_old;
        gstep = rhog_new - rhog_old;

        if (status)
            break;

        // dirty backstepping
        while (gsl_vector_get (s->x, 0) < 40. || gsl_vector_get (s->x, 0) > 1000.) {
            astep = astep/4.;
            aa_new = aa_old + astep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, rhog_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            aa_new = gsl_vector_get (s->x, 0.);
        }
        while (gsl_vector_get (s->x, 1) < 0.0 || gsl_vector_get (s->x, 1) > 0.9) {
            bstep = bstep/4.;
            basym_new = basym_old + bstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, rhog_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            basym_new = gsl_vector_get (s->x, 1);
        }
        while (gsl_vector_get (s->x, 2) < 0. || gsl_vector_get (s->x, 2) > rhob) {
            gstep = gstep/4.;
            rhog_new = rhog_old + gstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, rhog_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            rhog_new = gsl_vector_get (s->x, 2);
        }

        status = gsl_multiroot_test_residual (s->f, 9e-9);
    }

    while (status == GSL_CONTINUE && iter < 5000);

    print_state_icrust(s, rhob);

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return 0;
}
