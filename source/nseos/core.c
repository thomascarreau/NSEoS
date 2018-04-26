#include <gsl/gsl_multiroots.h>

#include "nuclear_matter.h"
#include "coulomb.h"
#include "modeling.h"
#include "core.h"

struct core_fun calc_core_fun(struct parameters satdata, double del_, double nu_, double nb_)
{
    struct core_fun result;
    double epsd;
    struct hnm meta_p;
    struct hnm meta_m;
    double dehnmddel;
    double np;
    double mueltot;
    double muutot;

    // energy per particle of infinite matter
    epsd = del_/1000.;
    meta_p = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, nb_, del_+epsd);
    meta_m = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, nb_, del_-epsd);
    dehnmddel = (meta_p.enpernuc - meta_m.enpernuc)/2./epsd;

    // electron(muon) chemical potential (including rest mass)
    np = nb_*(1.-del_)/2.;
    mueltot = calc_egas_chemical_potential(np - nu_); 
    muutot = calc_ugas_chemical_potential(nu_); 

    result.f_beta = 2.*dehnmddel + RMN - RMP - mueltot;
    result.f_mueq = mueltot - muutot;

    return result;
}

int assign_npecore_fun(const gsl_vector * x, void *params, gsl_vector * f)
{
    double nb = ((struct rparams_core *) params)->nb;
    struct parameters satdata = ((struct rparams_core *) params)->satdata;

    const double x0 = gsl_vector_get (x, 0);

    struct core_fun functs;
    functs = calc_core_fun(satdata, x0, 0., nb);

    const double y0 = functs.f_beta;

    gsl_vector_set(f, 0, y0);

    return GSL_SUCCESS;
}

struct core_compo calc_npecore_composition(double nb_, double *guess, struct parameters satdata)
{
    struct core_compo eq;

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t iter = 0;

    double del_old, del_new, dstep;

    struct rparams_core p;
    p.nb = nb_;
    p.satdata = satdata;

    double del_init = *guess;
    const size_t n = 1;
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector_set(x, 0, del_init);

    gsl_multiroot_function f = {&assign_npecore_fun, n, &p};

    T = gsl_multiroot_fsolver_broyden;
    s = gsl_multiroot_fsolver_alloc(T,1);
    gsl_multiroot_fsolver_set(s, &f, x);

    do
    {
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

    eq.del = gsl_vector_get(s->x, 0);
    eq.nu = 0.;

    *guess = eq.del;

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return eq;
}

int assign_npeucore_fun(const gsl_vector * x, void *params, gsl_vector * f)
{
    double nb = ((struct rparams_core *) params)->nb;
    struct parameters satdata = ((struct rparams_core *) params)->satdata;

    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);

    struct core_fun functs;
    functs = calc_core_fun(satdata, x0, x1, nb);

    const double y0 = functs.f_beta;
    const double y1 = functs.f_mueq;

    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);

    return GSL_SUCCESS;
}

struct core_compo calc_npeucore_composition(double nb_, double *guess, struct parameters satdata)
{
    struct core_compo eq;

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t iter = 0;

    double del_old, del_new, dstep;
    double nu_old, nu_new, ustep;

    struct rparams_core p;
    p.nb = nb_;
    p.satdata = satdata;

    const size_t n = 2;
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector_set(x, 0, guess[0]);
    gsl_vector_set(x, 1, guess[1]);

    gsl_multiroot_function f = {&assign_npeucore_fun, n, &p};

    T = gsl_multiroot_fsolver_broyden;
    s = gsl_multiroot_fsolver_alloc(T,2);
    gsl_multiroot_fsolver_set(s, &f, x);

    do
    {
        del_old = gsl_vector_get (s->x, 0);
        nu_old = gsl_vector_get (s->x, 1);

        iter++;

        status = gsl_multiroot_fsolver_iterate(s);

        del_new = gsl_vector_get (s->x, 0);
        nu_new = gsl_vector_get (s->x, 1);
        dstep = del_new - del_old;
        ustep = nu_new - nu_old;

        if (status)
            break;

        // dirty backstepping
        while (gsl_vector_get (s->x, 0) < 0.1 || gsl_vector_get (s->x, 0) > 1.0) {
            dstep = dstep/4.;
            del_new = del_old + dstep;
            gsl_vector_set (x, 0, del_new);
            gsl_vector_set (x, 1, nu_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            del_new = gsl_vector_get (s->x, 0);
        }
        while (gsl_vector_get (s->x, 1) < 0. || gsl_vector_get (s->x, 1) > nb_*(1.-gsl_vector_get(s->x, 0))/2.) {
            ustep = ustep/4.;
            nu_new = nu_old + ustep;
            gsl_vector_set (x, 0, del_new);
            gsl_vector_set (x, 1, nu_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            nu_new = gsl_vector_get (s->x, 1);
        }

        status = gsl_multiroot_test_residual (s->f, 9e-9);
    }
    while (status == GSL_CONTINUE && iter < 1000);

    eq.del = gsl_vector_get(s->x, 0);
    eq.nu = gsl_vector_get(s->x, 1);

    guess[0] = eq.del;
    guess[1] = eq.nu;

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return eq;
}

double calc_core_ws_cell_energy_density(struct parameters satdata, struct core_compo eq,
        double nb_)
{
    double np;
    double epseltot;
    double epsutot;
    struct hnm meta;
    double epsws;

    np = nb_*(1.-eq.del)/2.;
    epseltot = calc_egas_energy_density(np-eq.nu);
    epsutot = calc_ugas_energy_density(eq.nu);
    meta = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, nb_, eq.del);
    epsws = nb_*meta.enpernuc + epseltot + epsutot + np*(RMP-RMN) + nb_*RMN;

    return epsws;
}

double calc_core_ws_cell_pressure(struct parameters satdata, struct core_compo eq, 
        double nb_)
{
    double np;
    double egas_pressure;
    double ugas_pressure;
    struct hnm meta;
    double ws_cell_pressure;

    np = nb_*(1.-eq.del)/2.;
    egas_pressure = calc_egas_pressure(np-eq.nu);
    ugas_pressure = calc_ugas_pressure(eq.nu);

    meta = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, nb_, eq.del);
    ws_cell_pressure = meta.p + egas_pressure + ugas_pressure; 

    return ws_cell_pressure;
}

void print_state_core(struct parameters satdata, struct core_compo eq,
        double nb_, FILE *core, FILE *eos)
{
    double rhob;
    double np;
    double xp, xe, xu;
    double pressws;

    np = nb_*(1.-eq.del)/2.;
    xp = np/nb_;
    xe = (np - eq.nu)/nb_;
    xu = eq.nu/nb_;

    rhob = calc_core_ws_cell_energy_density(satdata, eq, nb_)*(ELEMC/1.e-19)/pow(SPEEDOFL/1.e8,2.)*1.e13;
    pressws = calc_core_ws_cell_pressure(satdata, eq, nb_);

    fprintf(core, "%g %g %g %g\n", nb_, xp, xe, xu);
    fprintf(eos, "%g %g %g\n", nb_, rhob, pressws);
}
