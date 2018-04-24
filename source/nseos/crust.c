#include <stdio.h>
#include <gsl/gsl_multiroots.h>

#include "nuclear_en.h"
#include "modeling.h"
#include "crust.h"

struct crust_fun_4d calc_crust_fun_4d(struct parameters satdata, struct sf_params sparams, 
        double aa_, double del_, double rho0_, double rhop_, double rhog_)
{
    struct crust_fun_4d result;
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

    enuc = CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, aa_, del_, rho0_, rhop_);
    epsa = 0.001;
    epsb = 0.0001;
    epsr = 0.0001;
    enuc_ap = CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, aa_+epsa, del_, rho0_, rhop_);
    enuc_am = CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, aa_-epsa, del_, rho0_, rhop_);
    enuc_bp = CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, aa_, del_+epsb, rho0_, rhop_);
    enuc_bm = CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, aa_, del_-epsb, rho0_, rhop_);
    enuc_rp = CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, aa_, del_, rho0_+epsr, rhop_);
    enuc_rm = CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, aa_, del_, rho0_-epsr, rhop_);
    denucdaa = (enuc_ap - enuc_am)/2./epsa; // 2 points derivatives
    denucddel = (enuc_bp - enuc_bm)/2./epsb;
    denucdrho0 = (enuc_rp - enuc_rm)/2./epsr;

    muel = calc_egas_chemical_potential(rhop_);
    dmu = calc_screening_derivative(satdata, aa_, del_, rho0_, rhop_);

    ngas = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, rhog_, 1.);

    result.f_stability = denucdaa - enuc/aa_;
    result.f_beta = denucddel*2./aa_ - muel - dmu - RMP + RMN;
    result.f_muneq = enuc/aa_ - (ngas.mun)*(1.-rhog_/rho0_) 
        + (1.-del_)/2.*(muel + dmu + RMP - RMN) - rhog_*(ngas.enpernuc)/rho0_ ;
    result.f_presseq = rho0_*rho0_*denucdrho0/aa_ - rhog_*ngas.mun + rhog_*ngas.enpernuc;

    return result;
}

int assign_ocrust_fun_3d(const gsl_vector * x, void *params, gsl_vector * f)
{
    double rhop = ((struct rparams_crust *) params)->rhop;
    struct parameters satdata = ((struct rparams_crust *) params)->satdata;
    struct sf_params sparams = ((struct rparams_crust *) params)->sparams;

    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);
    const double x2 = gsl_vector_get (x, 2);

    rhop = rhop*(1.-x1)/2.;

    struct crust_fun_4d functs;
    functs = calc_crust_fun_4d(satdata, sparams, x0, x1, x2, rhop, 0.);

    const double y0 = functs.f_stability;
    const double y1 = functs.f_beta;
    const double y2 = functs.f_presseq;

    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);
    gsl_vector_set(f, 2, y2);

    return GSL_SUCCESS;
}

struct compo calc_ocrust3d_composition(double rhob_, double *guess,
        struct parameters satdata, struct sf_params sparams)
{
    struct compo eq;

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t iter = 0;

    double aa_old, aa_new, astep;
    double basym_old, basym_new, bstep;
    double rho0_old, rho0_new, rstep;

    struct rparams_crust p;
    p.rhop = rhob_; // modification in assign_ocrust_fun_3d (DIRTY!)
    p.satdata = satdata;
    p.sparams = sparams;

    const size_t n = 3;
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector_set(x, 0, guess[0]);
    gsl_vector_set(x, 1, guess[1]);
    gsl_vector_set(x, 2, guess[2]);

    gsl_multiroot_function f = {&assign_ocrust_fun_3d, n, &p};

    T = gsl_multiroot_fsolver_broyden;
    s = gsl_multiroot_fsolver_alloc(T,3);
    gsl_multiroot_fsolver_set(s, &f, x);

    do
    {
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
        while (gsl_vector_get (s->x, 0) < 0.) {
            astep = astep/4.;
            aa_new = aa_old + astep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, rho0_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            aa_new = gsl_vector_get (s->x, 0.);
        }
        while (gsl_vector_get (s->x, 1) < 0.0 || gsl_vector_get (s->x, 1) > 0.5) {
            bstep = bstep/4.;
            basym_new = basym_old + bstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, rho0_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            basym_new = gsl_vector_get (s->x, 1);
        }
        while (gsl_vector_get (s->x, 2) < 0.) {
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

    eq.aa = gsl_vector_get(s->x, 0);
    eq.del = gsl_vector_get(s->x, 1);
    eq.rho0 = gsl_vector_get(s->x, 2);
    eq.rhog = 0.;

    guess[0] = eq.aa;
    guess[1] = eq.del;
    guess[2] = eq.rho0;

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return eq;
}

int assign_icrust_fun_4d(const gsl_vector * x, void *params, gsl_vector * f)
{
    double rhop = ((struct rparams_crust *) params)->rhop;
    struct parameters satdata = ((struct rparams_crust *) params)->satdata;
    struct sf_params sparams = ((struct rparams_crust *) params)->sparams;

    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);
    const double x2 = gsl_vector_get (x, 2);
    const double x3 = gsl_vector_get (x, 3);

    rhop = (rhop-x3)*(1.-x1)/2./(1.-x3/x2);

    struct crust_fun_4d functs;
    functs = calc_crust_fun_4d(satdata, sparams, x0, x1, x2, rhop, x3);

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

struct compo calc_icrust4d_composition(double rhob_, double *guess,
        struct parameters satdata, struct sf_params sparams)
{
    struct compo eq;

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t iter = 0;

    double aa_old, aa_new, astep;
    double basym_old, basym_new, bstep;
    double rho0_old, rho0_new, rstep;
    double rhog_old, rhog_new, gstep;

    struct rparams_crust p;
    p.rhop = rhob_; // modification in assign_icrust_fun_4d (DIRTY!)
    p.satdata = satdata;
    p.sparams = sparams;

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

    eq.aa = gsl_vector_get(s->x, 0);
    eq.del = gsl_vector_get(s->x, 1);
    eq.rho0 = gsl_vector_get(s->x, 2);
    eq.rhog = gsl_vector_get(s->x, 3);

    guess[0] = eq.aa;
    guess[1] = eq.del;
    guess[2] = eq.rho0;
    guess[3] = eq.rhog;

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return eq;
}

double calc_muncl(struct parameters satdata, struct sf_params sparams, struct compo eq, double rhob_)
{
    double rhop;
    double enuc;
    double epsb;
    double enuc_bp, enuc_bm;
    double denucddel;

    rhop = rhob_*(1.-eq.del)/2.;

    enuc = CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, eq.aa, eq.del, eq.rho0, rhop);
    epsb = 0.0001;
    enuc_bp = CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, eq.aa, eq.del+epsb, eq.rho0, rhop);
    enuc_bm = CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, eq.aa, eq.del-epsb, eq.rho0, rhop);
    denucddel = (enuc_bp - enuc_bm)/2./epsb;

    return enuc/eq.aa + (1.-eq.del)/eq.aa * denucddel;
}

double calc_crust_ws_cell_energy_density(struct parameters satdata, struct sf_params sparams, struct compo eq,
        double rhob_)
{
    double rhop;
    double vws;
    double epseltot;
    double enuc;
    struct hnm ngas;
    double epsg;
    double epsws; 

    vws = eq.aa*(1. - eq.rhog/eq.rho0)/(rhob_ - eq.rhog);
    rhop = eq.aa*(1.-eq.del)/2./vws;

    epseltot = calc_egas_energy_density(rhop);

    enuc = CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, eq.aa, eq.del, eq.rho0, rhop);

    ngas = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, eq.rhog, 1.);
    epsg = eq.rhog*ngas.enpernuc;

    epsws = enuc/vws + epseltot + epsg*(1.-eq.aa/eq.rho0/vws)
        + rhop*(RMP-RMN) + rhob_*(RMN-AMU);

    return epsws;
}

double calc_ngas_pressure(struct parameters satdata, double rhog_)
{
    struct hnm ngas;
    ngas = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, rhog_, 1.);

    return rhog_*ngas.mun - rhog_*ngas.enpernuc;
}

double calc_crust_ws_cell_pressure(struct parameters satdata, struct compo eq, double rhob_)
{
    double rhop;
    double egas_pressure;
    double lattice_pressure;
    double ngas_pressure;
    double ws_cell_pressure;


    rhop = (rhob_-eq.rhog)*(1.-eq.del)/2./(1.-eq.rhog/eq.rho0);
    egas_pressure = calc_egas_pressure(rhop);
    lattice_pressure = calc_lattice_pressure(satdata, eq.aa, eq.del, eq.rho0, rhop);
    ngas_pressure = calc_ngas_pressure(satdata, eq.rhog);
    ws_cell_pressure = egas_pressure + lattice_pressure + ngas_pressure;

    return ws_cell_pressure;
}

void print_state_crust(struct parameters satdata, struct sf_params sparams, struct compo eq, 
        double rhob_, FILE *compo, FILE *eos)
{
    double vws, rws;
    double epsws;
    double pressws;

    vws = eq.aa/(rhob_-eq.rhog)*(1.-eq.rhog/eq.rho0);
    rws = pow(3.*vws/4./PI,1./3.);

    epsws = calc_crust_ws_cell_energy_density(satdata, sparams, eq, rhob_);
    pressws = calc_crust_ws_cell_pressure(satdata, eq, rhob_);

    fprintf (compo, "%g %g %g %g %g %g %g %g\n", rhob_, eq.aa, eq.del, eq.aa*(1.-eq.del)/2., eq.rho0, eq.rhog, rws, epsws);
    fprintf(eos, "%g %g\n", rhob_, pressws);
}
