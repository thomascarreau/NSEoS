#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "../nseos/nuclear_matter.h"
#include "../nseos/coulomb_en.h"

//====================================================================================
//====================================================================================

#define ASSIGN_PARAM (assign_param_sly4)

static const int taylor_exp_order = 2;

//====================================================================================
//====================================================================================

// ==================== FUNCTIONS ====================

double calc_core_fun(double del_, double rhob_);
struct rparams
{
    double rhob;
};
int assign_core_fun(const gsl_vector * x, void *params, gsl_vector * f);
double calc_epsws(double del_, double rhob_);
void print_state_core(gsl_multiroot_fsolver * s, double rhob_);

// ==================== MAIN ====================

int main(void)
{
    int irhob;
    double rhob;
    double del_new_guess;

    del_new_guess = 0.7;

    for(irhob = 40; irhob <= 500; irhob ++)
    {
        const gsl_multiroot_fsolver_type *T;
        gsl_multiroot_fsolver *s;

        int status;
        size_t iter = 0;

        double del_old, del_new, dstep;

        rhob = irhob/1000.;
        struct rparams p;
        p.rhob = rhob;

        double del_init = del_new_guess;
        const size_t n = 1;
        gsl_vector *x = gsl_vector_alloc(n);
        gsl_vector_set(x, 0, del_init);

        gsl_multiroot_function f = {&assign_core_fun, n, &p};

        T = gsl_multiroot_fsolver_broyden;
        s = gsl_multiroot_fsolver_alloc(T,1);
        gsl_multiroot_fsolver_set(s, &f, x);

        do
        {
            /* print_state_core(s, rhob); */

            del_old = gsl_vector_get (s->x, 0);

            iter++;

            status = gsl_multiroot_fsolver_iterate(s);

            del_new = gsl_vector_get (s->x, 0);
            dstep = del_new - del_old;

            if (status)
                break;

            // dirty backstepping
            while (gsl_vector_get (s->x, 0) < 0.1 || gsl_vector_get (s->x, 0) > 1.0) {
                dstep = dstep/4.;
                del_new = del_old + dstep;
                gsl_vector_set (x, 0, del_new);
                gsl_multiroot_fsolver_set (s, &f, x);
                del_new = gsl_vector_get (s->x, 0);
            }

            status = gsl_multiroot_test_residual (s->f, 9e-9);
        }

        while (status == GSL_CONTINUE && iter < 1000);

        if (gsl_vector_get(s->x, 0) != gsl_vector_get(s->x, 0))
            break;

        print_state_core(s, rhob);

        del_new_guess = gsl_vector_get(s->x, 0);

        gsl_multiroot_fsolver_free(s);
        gsl_vector_free(x);
    }

    return 0;
}

double calc_core_fun(double del_, double rhob_)
{
    struct parameters satdata;
    double epsd;
    struct hnm meta_p;
    struct hnm meta_m;
    double dehnmddel;
    double rhop;
    double mueltot;
    double beta_eq_fun;

    // energy per particle of infinite matter
    satdata = ASSIGN_PARAM(satdata);
    epsd = del_/1000.;
    meta_p = calc_meta_model_nuclear_matter(satdata, taylor_exp_order, rhob_, del_+epsd);
    meta_m = calc_meta_model_nuclear_matter(satdata, taylor_exp_order, rhob_, del_-epsd);
    dehnmddel = (meta_p.enpernuc - meta_m.enpernuc)/2./epsd;

    // electron gas chemical potential (including rest mass)
    rhop = rhob_*(1.-del_)/2.;
    mueltot = calc_egas_chemical_potential(rhop); 

    beta_eq_fun = 2.*dehnmddel + RMN - RMP - mueltot;

    return beta_eq_fun;
}

int assign_core_fun(const gsl_vector * x, void *params, gsl_vector * f)
{
    double rhob = ((struct rparams *) params)->rhob;

    const double x0 = gsl_vector_get (x, 0);

    double func = calc_core_fun(x0, rhob);

    const double y0 = func;

    gsl_vector_set(f, 0, y0);

    return GSL_SUCCESS;
}

double calc_epsws(double del_, double rhob_)
{
    struct parameters satdata;
    struct hnm meta;
    double rhop;
    double epseltot;
    double epsws;

    // energy per particle of infinite matter
    satdata = ASSIGN_PARAM(satdata);
    meta = calc_meta_model_nuclear_matter(satdata, taylor_exp_order, rhob_, del_);

    // electron gas energy density (including rest mass)
    rhop = rhob_*(1.-del_)/2.;
    epseltot = calc_egas_energy_density(rhop);

    epsws = rhob_*meta.enpernuc + epseltot + rhop*(RMP-RMN) + rhob_*(RMN-AMU);

    return epsws;
}

double calc_pressure(double del_, double rhob_)
{
    struct parameters satdata;
    struct hnm meta;
    double rhop;
    double epseltot;
    double mueltot;
    double pressure;

    // pressure of infinite matter
    satdata = ASSIGN_PARAM(satdata);
    meta = calc_meta_model_nuclear_matter(satdata, taylor_exp_order, rhob_, del_);

    // electron gas energy density and chemical potential (including rest mass)
    rhop = rhob_*(1.-del_)/2.;
    epseltot = calc_egas_energy_density(rhop);
    mueltot = calc_egas_chemical_potential(rhop);

    pressure = meta.p + rhop*mueltot - epseltot;

    return pressure;
}

void print_state_core(gsl_multiroot_fsolver * s, double rhob_)
{
    double del_eq;
    double yp_eq;
    double epsws;
    double pressure;

    del_eq = gsl_vector_get(s->x, 0);
    yp_eq = (1.-del_eq)/2.;
    epsws = calc_epsws(del_eq, rhob_);
    pressure = calc_pressure(del_eq, rhob_);

    printf("%g %g %g %g %g\n", rhob_, del_eq, yp_eq, epsws, pressure);
    /* printf("%g %g\n", rhob_, pressure); */
}
