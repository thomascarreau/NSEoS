#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "../nseos/ws_cell.h"

// ==================== FUNCTIONS ====================

struct ocrust_fun_3d
{
    double f_stability;
    double f_beta;
    double f_press;
    double mu_n;

    double denucddel;
    double mueltot;
    double dm;
};

struct ocrust_fun_3d calc_ocrust_fun_3d(double aa_, double del_, double rho0_, double rhop_);
struct ocrust_fun_3d calc_ocrust_fun_3d(double aa_, double del_, double rho0_, double rhop_)
{
    struct ocrust_fun_3d result;
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
    struct gas egas;
    struct coulomb_energy_shift coul_shift;
    double dmu;
    double mueltot;

    satdata = assign_param(satdata);
    enuc = calc_enuc(satdata, aa_, del_, rho0_, rhop_, 0.);
    epsa = 0.001;
    epsb = 0.0001;
    epsr = 0.0001;
    enuc_ap = calc_enuc(satdata, aa_+epsa, del_, rho0_, rhop_, 0.);
    enuc_am = calc_enuc(satdata, aa_-epsa, del_, rho0_, rhop_, 0.);
    enuc_bp = calc_enuc(satdata, aa_, del_+epsb, rho0_, rhop_, 0.);
    enuc_bm = calc_enuc(satdata, aa_, del_-epsb, rho0_, rhop_, 0.);
    enuc_rp = calc_enuc(satdata, aa_, del_, rho0_+epsr, rhop_, 0.);
    enuc_rm = calc_enuc(satdata, aa_, del_, rho0_-epsr, rhop_, 0.);
    denucdaa = (enuc_ap - enuc_am)/2./epsa; // 2 points derivatives
    denucddel = (enuc_bp - enuc_bm)/2./epsb;
    denucdrho0 = (enuc_rp - enuc_rm)/2./epsr;
    egas = calc_egas(rhop_);
    coul_shift = calc_coulomb_energy_shift(satdata, aa_, del_, rho0_, rhop_, 0.);
    dmu = coul_shift.derivative;
    mueltot = egas.mutot;

    result.f_stability = denucdaa/aa_ - enuc/aa_/aa_;
    /* double es; */
    /* double eca, ecb; */
    /* struct coulomb_energy_shift screen; */
    /* screen = calc_coulomb_energy_shift(aa_, del_, rho0_, rhop_); */
    /* es = calc_surface_energy(satdata, aa_, del_, rho0_); */
    /* eca = calc_coulomb_energy(aa_, del_, rho0_); */
    /* ecb = screen.energy; */
    /* result.f_stability = es - 2.*(eca + ecb); */ 
    result.f_beta = denucddel*2./aa_ - mueltot - dmu - rmp + rmn;
    result.f_press = rho0_*rho0_*denucdrho0/aa_;
    result.mu_n = denucdaa + (1.-del_)/aa_*denucddel;

    result.denucddel = denucddel;
    result.mueltot = mueltot;
    result.dm = rmp - rmn;

    return result;
}

struct rparams
{
    double rhop;
};

int assign_ocrust_fun_3d(const gsl_vector * x, void *params, gsl_vector * f);
int assign_ocrust_fun_3d(const gsl_vector * x, void *params, gsl_vector * f)
{
    double rhop = ((struct rparams *) params)->rhop;
    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);
    const double x2 = gsl_vector_get (x, 2);

    double bdel;
    double r1, r2;
    double ac, qq;
    double r0;
    r0 = pow(3./4./pi/x2,1./3.);
    ac = 3./5.*alphafs*hbarc/r0;
    qq = 144.5*33.43/(77.92 + 55.5);
    r1 = 3.*ac/32./qq;
    r2 = 9.*33.43/4./qq;
    bdel = (x1 + r1*(1.-x1)*(1.-x1)*pow(x0,1./3.))/(1. + r2/pow(x0,1./3.));

    bdel = x1;
    rhop = rhop*(1.-bdel)/2.;

    struct ocrust_fun_3d functs;
    functs = calc_ocrust_fun_3d(x0, x1, x2, rhop);
    const double y0 = functs.f_stability;
    const double y1 = functs.f_beta;
    const double y2 = functs.f_press;
    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);
    gsl_vector_set(f, 2, y2);

    return GSL_SUCCESS;
}

void print_state_ocrust(gsl_multiroot_fsolver * s, double rhob_);
void print_state_ocrust(gsl_multiroot_fsolver * s, double rhob_)
{
    double aa_eq, del_eq, rhop_eq;
    double zz_eq;
    struct parameters satdata;
    double rho0_eq;
    double vws;
    double epsws;
    struct ocrust_fun_3d fun_eq;
    double f0, f1, f2;

    aa_eq = gsl_vector_get(s->x, 0);
    del_eq = gsl_vector_get(s->x, 1);
    rhop_eq = rhob_*(1.-del_eq)/2.;
    zz_eq = aa_eq*(1.-del_eq)/2.;
    rho0_eq = gsl_vector_get(s->x, 2);
    vws = aa_eq/rhob_;
    satdata = assign_param(satdata);
    epsws = calc_epsws(satdata, aa_eq, del_eq, rho0_eq, 0., rhop_eq, rhob_);
    fun_eq = calc_ocrust_fun_3d(aa_eq, del_eq, rho0_eq, rhop_eq);
    f0 = gsl_vector_get(s->f, 0);
    f1 = gsl_vector_get(s->f, 1);
    f2 = gsl_vector_get(s->f, 2);
    printf ("%g %g %g %g %g %g %g %g %g %g %g %g\n", rhob_, aa_eq, del_eq, zz_eq, rhop_eq/rhob_, rho0_eq, 
            vws, epsws, f0, f1, f2, fun_eq.mu_n);
}

// ==================== MAIN ====================

int main(void)
{
    double rhob;
    rhob = 2.4e-5;

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t iter = 0;

    double aa_old, aa_new, astep;
    double basym_old, basym_new, bstep;
    double rho0_old, rho0_new, rstep;

    struct rparams p;
    p.rhop = rhob; // modification in assign_ocrust_fun_3d (DIRTY!)

    double x_init[3] = {137.759, 0.501242, 0.139738};
    const size_t n = 3;
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector_set(x, 0, x_init[0]);
    gsl_vector_set(x, 1, x_init[1]);
    gsl_vector_set(x, 2, x_init[2]);

    gsl_multiroot_function f = {&assign_ocrust_fun_3d, n, &p};

    T = gsl_multiroot_fsolver_dnewton;
    s = gsl_multiroot_fsolver_alloc(T,3);
    gsl_multiroot_fsolver_set(s, &f, x);

    do
    {
        print_state_ocrust(s, rhob);

        aa_old = gsl_vector_get (s->x, 0);
        basym_old = gsl_vector_get (s->x, 1);
        rho0_old = gsl_vector_get (s->x, 2);

        iter++;

        status = gsl_multiroot_fsolver_iterate(s);

        aa_new = gsl_vector_get (s->x, 0);
        basym_new = gsl_vector_get (s->x, 1);
        rho0_new = gsl_vector_get (s->x, 2);
        astep = aa_new - aa_old;
        bstep = basym_new - basym_old;
        rstep = rho0_new - rho0_old;

        if (status)
            break;

        // dirty backstepping
        while (gsl_vector_get (s->x, 0) < 40. || gsl_vector_get (s->x, 0) > 500.) {
            astep = astep/4.;
            aa_new = aa_old + astep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, rho0_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            aa_new = gsl_vector_get (s->x, 0.);
        }
        while (gsl_vector_get (s->x, 1) < 0.1 || gsl_vector_get (s->x, 1) > 0.7) {
            bstep = bstep/4.;
            basym_new = basym_old + bstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, rho0_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            basym_new = gsl_vector_get (s->x, 1);
        }
        while (gsl_vector_get (s->x, 2) < 0.07 || gsl_vector_get (s->x, 2) > 0.1540) {
            rstep = rstep/4.;
            rho0_new = rho0_old + rstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, rho0_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            rho0_new = gsl_vector_get (s->x, 2);
        }

        status = gsl_multiroot_test_residual (s->f, 9e-9);
    }

    while (status == GSL_CONTINUE && iter < 5000);

    print_state_ocrust(s, rhob);

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return 0;
}
