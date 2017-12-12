#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "../nseos/ws_cell.h"

// ==================== FUNCTIONS ====================

struct threeconstraints_fun
{
    double f_stability;
    double f_press;
};

struct threeconstraints_fun eval_threeconstraints_fun(double aa_, double del_, double rho0_, double rhop_, double rhog_);
struct threeconstraints_fun eval_threeconstraints_fun(double aa_, double del_, double rho0_, double rhop_, double rhog_)
{
    struct threeconstraints_fun result;
    struct parameters satdata;
    double enuc;
    double epsa;
    double epsr;
    double enuc_ap, enuc_am;
    double enuc_rp, enuc_rm;
    double denucdaa;
    double denucdrho0;
    struct hnm ngas;

    satdata = assign_param(satdata);

    enuc = calc_enuc(satdata, aa_, del_, rho0_, rhop_);
    epsa = 0.001;
    epsr = 0.0001;
    enuc_ap = calc_enuc(satdata, aa_+epsa, del_, rho0_, rhop_);
    enuc_am = calc_enuc(satdata, aa_-epsa, del_, rho0_, rhop_);
    enuc_rp = calc_enuc(satdata, aa_, del_, rho0_+epsr, rhop_);
    enuc_rm = calc_enuc(satdata, aa_, del_, rho0_-epsr, rhop_);
    denucdaa = (enuc_ap - enuc_am)/2./epsa; // 2 points derivatives
    denucdrho0 = (enuc_rp - enuc_rm)/2./epsr;

    ngas = calc_hnm(satdata, rhog_, 1.);

    result.f_stability = denucdaa/aa_ - enuc/aa_/aa_;
    result.f_press = rho0_*rho0_*denucdrho0/aa_ - rhog_*ngas.mug + rhog_*ngas.enpernuc;

    return result;
}

struct rparams
{
    double del;
    double rhop;
    double rhog;
};

int assign_threeconstraints_fun(const gsl_vector * x, void *params, gsl_vector * f);
int assign_threeconstraints_fun(const gsl_vector * x, void *params, gsl_vector * f)
{
    double rhop = ((struct rparams *) params)->rhop;
    double del = ((struct rparams *) params)->del;
    double rhog = ((struct rparams *) params)->rhog;

    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);

    rhop = rhop*(1.-del)/2.;

    struct threeconstraints_fun functs;
    functs = eval_threeconstraints_fun(x0, del, x1, rhop, rhog);

    const double y0 = functs.f_stability;
    const double y1 = functs.f_press;

    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);

    return GSL_SUCCESS;
}

void print_state(gsl_multiroot_fsolver * s, double rhob_, struct rparams params_);
void print_state(gsl_multiroot_fsolver * s, double rhob_, struct rparams params_)
{
    double aa_eq;
    double rho0_eq;
    double press_equat;
    double zz;
    double rhop;

    aa_eq = gsl_vector_get(s->x,0);
    rho0_eq = gsl_vector_get(s->x, 1);
    press_equat = gsl_vector_get(s->f, 1);
    zz = aa_eq*(1.-params_.del)/2.;
    rhop = rhob_*(1.-params_.del)/2.;

    printf ("%g %g %g %g %g %g %g %g\n", rhob_,  aa_eq, params_.del, zz, rhop/rhob_, rho0_eq, params_.rhog/rhob_, press_equat);
}

// ==================== MAIN ====================

int main(void)
{
    double rhob;
    rhob = 1.e-3;

    double rhog;
    for(rhog = rhob/10.; rhog < rhob; rhog += rhob/100.)
    {

        const gsl_multiroot_fsolver_type *T;
        gsl_multiroot_fsolver *s;

        int status;
        size_t iter = 0;

        struct rparams p;
        p.del = 0.4;
        p.rhop = rhob; // modification in assign_threeconstraints_fun (DIRTY!)
        p.rhog = rhog;

        double x_init[2] = {82., 0.12};
        const size_t n = 2;
        gsl_vector *x = gsl_vector_alloc(n);
        gsl_vector_set(x, 0, x_init[0]);
        gsl_vector_set(x, 1, x_init[1]);

        gsl_multiroot_function f = {&assign_threeconstraints_fun, n, &p};

        T = gsl_multiroot_fsolver_dnewton;
        s = gsl_multiroot_fsolver_alloc(T,2);
        gsl_multiroot_fsolver_set(s, &f, x);

        do
        {
            /* print_state(s, rhob, p); */

            iter++;

            status = gsl_multiroot_fsolver_iterate(s);

            if (status)
                break;

            status = gsl_multiroot_test_residual (s->f, 9e-9);
        }

        while (status == GSL_CONTINUE && iter < 5000);

        print_state(s, rhob, p);

        gsl_multiroot_fsolver_free(s);
        gsl_vector_free(x);
    }

    return 0;
}
