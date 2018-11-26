#include <stdio.h>
#include <gsl/gsl_multiroots.h>

#include "lepton.h"
#include "nuclear_en.h"
#include "modeling.h"
#include "crust.h"

double calc_zp_en(struct parameters satdata, struct sf_params sparams, 
        double aa_, double ii_, double n0_, double np_)
{
    double hbaromega_p;
    double zz;
    double mi;
    double vws;
    double u1;

    zz = aa_*(1.-ii_)/2.;
    mi = zz*RMP + (aa_-zz)*RMN 
        + CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, aa_, ii_, n0_);
    vws = zz/np_;

    hbaromega_p = sqrt(pow(HBARC,2.)*4.*PI*pow(zz,2.)*ALPHAFS*HBARC
            /mi/vws);

    u1 = 0.5113875; // bcc lattice; Table 2.4 of Haensel book

    return 3./2.*hbaromega_p*u1;
}

double calc_ion_en(struct parameters satdata, struct sf_params sparams,
        double aa_, double del_, double n0_, double np_)
{
    return CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, aa_, del_, n0_)
        + calc_zp_en(satdata, sparams, aa_, del_, n0_, np_)
        + calc_lattice_en(satdata, aa_, del_, n0_, np_);
}

struct crust_fun_4d calc_crust_fun_4d(struct parameters satdata, 
        struct sf_params sparams, 
        double aa_, double del_, double n0_, double np_, double ng_)
{
    struct crust_fun_4d result;
    double eion;
    double epsa;
    double epsb;
    double epsr;
    double epsp;
    double eion_ap, eion_am;
    double eion_bp, eion_bm;
    double eion_rp, eion_rm;
    double eion_pp, eion_pm;
    double deiondaa;
    double deionddel;
    double deiondn0;
    double deiondnp;
    double muel;
    struct hnm ngas;

    eion = calc_ion_en(satdata, sparams, aa_, del_, n0_, np_);
    epsa = 0.001;
    epsb = 0.0001;
    epsr = 0.0001;
    epsp = np_/1000.;
    eion_ap = calc_ion_en(satdata, sparams, aa_+epsa, del_, n0_, np_);
    eion_am = calc_ion_en(satdata, sparams, aa_-epsa, del_, n0_, np_);
    eion_bp = calc_ion_en(satdata, sparams, aa_, del_+epsb, n0_, np_);
    eion_bm = calc_ion_en(satdata, sparams, aa_, del_-epsb, n0_, np_);
    eion_rp = calc_ion_en(satdata, sparams, aa_, del_, n0_+epsr, np_);
    eion_rm = calc_ion_en(satdata, sparams, aa_, del_, n0_-epsr, np_);
    eion_pp = calc_ion_en(satdata, sparams, aa_, del_, n0_, np_+epsp);
    eion_pm = calc_ion_en(satdata, sparams, aa_, del_, n0_, np_-epsp);
    deiondaa = (eion_ap - eion_am)/2./epsa; // 2 points derivatives
    deionddel = (eion_bp - eion_bm)/2./epsb;
    deiondn0 = (eion_rp - eion_rm)/2./epsr;
    deiondnp = (eion_pp - eion_pm)/2./epsp;

    muel = calc_egas_chemical_potential(np_);

    ngas = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, ng_, 1.);

    result.f_stability = deiondaa - eion/aa_;
    result.f_beta = (deionddel + np_/(1.-del_)*deiondnp)*2./aa_ - muel 
        - RMP + RMN;
    result.f_muneq = eion/aa_ + (1.-del_)/aa_*deionddel 
        - ngas.mun + ng_/n0_*(ngas.mun - ngas.enpernuc);
    result.f_presseq = n0_*n0_*deiondn0/aa_ - ng_*ngas.mun + ng_*ngas.enpernuc;

    return result;
}

int assign_ocrust_fun_3d(const gsl_vector * x, void *params, gsl_vector * f)
{
    double np = ((struct rparams_crust *) params)->np;
    struct parameters satdata = ((struct rparams_crust *) params)->satdata;
    struct sf_params sparams = ((struct rparams_crust *) params)->sparams;

    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);
    const double x2 = gsl_vector_get (x, 2);

    np = np*(1.-x1)/2.;

    struct crust_fun_4d functs;
    functs = calc_crust_fun_4d(satdata, sparams, x0, x1, x2, np, 0.);

    const double y0 = functs.f_stability;
    const double y1 = functs.f_beta;
    const double y2 = functs.f_presseq;

    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);
    gsl_vector_set(f, 2, y2);

    return GSL_SUCCESS;
}

struct compo calc_ocrust3d_composition(double nb_, double *guess,
        struct parameters satdata, struct sf_params sparams)
{
    struct compo eq;

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t iter = 0;

    double aa_old, aa_new, astep;
    double basym_old, basym_new, bstep;
    double n0_old, n0_new, rstep;

    struct rparams_crust p;
    p.np = nb_; // modification in assign_ocrust_fun_3d (DIRTY!)
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
        n0_old = gsl_vector_get (s->x, 2);

        iter++;

        status = gsl_multiroot_fsolver_iterate(s);

        aa_new = gsl_vector_get (s->x, 0);
        basym_new = gsl_vector_get (s->x, 1);
        n0_new = gsl_vector_get (s->x, 2);
        astep = aa_new - aa_old;
        bstep = basym_new - basym_old;
        rstep = n0_new - n0_old;

        if (status)
            break;

        // dirty backstepping
        while (gsl_vector_get (s->x, 0) < 0.) {
            astep = astep/4.;
            aa_new = aa_old + astep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, n0_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            aa_new = gsl_vector_get (s->x, 0.);
        }
        while (gsl_vector_get (s->x, 1) < 0.01 
                || gsl_vector_get (s->x, 1) > 1.0) {
            bstep = bstep/4.;
            basym_new = basym_old + bstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, n0_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            basym_new = gsl_vector_get (s->x, 1);
        }
        while (gsl_vector_get (s->x, 2) < 0.) {
            rstep = rstep/4.;
            n0_new = n0_old + rstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, n0_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            n0_new = gsl_vector_get (s->x, 2);
        }

        status = gsl_multiroot_test_residual (s->f, 9e-9);
    }

    while (status == GSL_CONTINUE && iter < 1000);

    if (iter == 1000)
    {
        eq.aa = NAN;
        eq.del = NAN;
        eq.n0 = NAN;
        eq.ng = 0.0;
    } 
    else 
    {
        eq.aa = gsl_vector_get(s->x, 0);
        eq.del = gsl_vector_get(s->x, 1);
        eq.n0 = gsl_vector_get(s->x, 2);
        eq.ng = 0.0;
    }

    guess[0] = eq.aa;
    guess[1] = eq.del;
    guess[2] = eq.n0;

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return eq;
}

int assign_icrust_fun_4d(const gsl_vector * x, void *params, gsl_vector * f)
{
    double np = ((struct rparams_crust *) params)->np;
    struct parameters satdata = ((struct rparams_crust *) params)->satdata;
    struct sf_params sparams = ((struct rparams_crust *) params)->sparams;

    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);
    const double x2 = gsl_vector_get (x, 2);
    const double x3 = gsl_vector_get (x, 3);

    np = (np-x3)*(1.-x1)/2./(1.-x3/x2);

    struct crust_fun_4d functs;
    functs = calc_crust_fun_4d(satdata, sparams, x0, x1, x2, np, x3);

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

struct compo calc_icrust4d_composition(double nb_, double *guess,
        struct parameters satdata, struct sf_params sparams)
{
    struct compo eq;

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t iter = 0;

    double aa_old, aa_new, astep;
    double basym_old, basym_new, bstep;
    double n0_old, n0_new, rstep;
    double ng_old, ng_new, gstep;

    struct rparams_crust p;
    p.np = nb_; // modification in assign_icrust_fun_4d (DIRTY!)
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
        n0_old = gsl_vector_get (s->x, 2);
        ng_old = gsl_vector_get (s->x, 3);

        iter++;

        status = gsl_multiroot_fsolver_iterate(s);

        aa_new = gsl_vector_get (s->x, 0);
        basym_new = gsl_vector_get (s->x, 1);
        n0_new = gsl_vector_get (s->x, 2);
        ng_new = gsl_vector_get (s->x, 3);
        astep = aa_new - aa_old;
        bstep = basym_new - basym_old;
        rstep = n0_new - n0_old;
        gstep = ng_new - ng_old;

        if (status)
            break;

        // dirty backstepping
        while (gsl_vector_get (s->x, 0) < 0.) {
            astep = astep/4.;
            aa_new = aa_old + astep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, n0_new);
            gsl_vector_set (x, 3, ng_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            aa_new = gsl_vector_get (s->x, 0.);
        }
        while (gsl_vector_get (s->x, 1) < -1.0 
                || gsl_vector_get (s->x, 1) > 1.0) {
            bstep = bstep/4.;
            basym_new = basym_old + bstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, n0_new);
            gsl_vector_set (x, 3, ng_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            basym_new = gsl_vector_get (s->x, 1);
        }
        while (gsl_vector_get (s->x, 2) < 0. 
                || gsl_vector_get(s->x, 2) > 1.0) {
            rstep = rstep/4.;
            n0_new = n0_old + rstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, n0_new);
            gsl_vector_set (x, 3, ng_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            n0_new = gsl_vector_get (s->x, 2);
        }
        while (gsl_vector_get (s->x, 3) < 1.e-10 
                || gsl_vector_get (s->x, 3) > nb_) {
            gstep = gstep/4.;
            ng_new = ng_old + gstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, basym_new);
            gsl_vector_set (x, 2, n0_new);
            gsl_vector_set (x, 3, ng_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            ng_new = gsl_vector_get (s->x, 3);
        }
        
        if  (gsl_vector_get(s->x,1) > 0.99)
            iter = 1000;

        status = gsl_multiroot_test_residual (s->f, 9e-9);
    }

    while (status == GSL_CONTINUE && iter < 1000);

    if (iter == 1000)
    {
        eq.aa = NAN;
        eq.del = NAN;
        eq.n0 = NAN;
        eq.ng = NAN;
    } 
    else 
    {
        eq.aa = gsl_vector_get(s->x, 0);
        eq.del = gsl_vector_get(s->x, 1);
        eq.n0 = gsl_vector_get(s->x, 2);
        eq.ng = gsl_vector_get(s->x, 3);
    }

    guess[0] = eq.aa;
    guess[1] = eq.del;
    guess[2] = eq.n0;
    guess[3] = eq.ng;

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return eq;
}

double calc_muncl(struct parameters satdata, struct sf_params sparams, 
        struct compo eq, double nb_)
{
    double np;
    double eion;
    double epsb;
    double eion_bp, eion_bm;
    double deionddel;

    np = nb_*(1.-eq.del)/2.;

    eion = calc_ion_en(satdata, sparams, eq.aa, eq.del, eq.n0, np);
    epsb = 0.0001;
    eion_bp = calc_ion_en(satdata, sparams, eq.aa, eq.del+epsb, eq.n0, np);
    eion_bm = calc_ion_en(satdata, sparams, eq.aa, eq.del-epsb, eq.n0, np);
    deionddel = (eion_bp - eion_bm)/2./epsb;

    return eion/eq.aa + (1.-eq.del)/eq.aa * deionddel;
}

double calc_crust_ws_cell_energy_density(struct parameters satdata, 
        struct sf_params sparams, 
        struct compo eq, double nb_)
{
    double np;
    double vws;
    double epseltot;
    double eion;
    struct hnm ngas;
    double epsg;
    double epsws; 

    vws = eq.aa*(1. - eq.ng/eq.n0)/(nb_ - eq.ng);
    np = eq.aa*(1.-eq.del)/2./vws;

    epseltot = calc_egas_energy_density(np);

    eion = calc_ion_en(satdata, sparams, eq.aa, eq.del, eq.n0, np);

    ngas = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, 
            eq.ng, 1.);
    epsg = eq.ng*ngas.enpernuc;

    epsws = eion/vws + epseltot + epsg*(1.-eq.aa/eq.n0/vws)
        + np*(RMP-RMN) + nb_*RMN;

    return epsws;
}

double calc_ngas_pressure(struct parameters satdata, double ng_)
{
    struct hnm ngas;
    ngas = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, ng_, 1.);

    return ng_*ngas.mun - ng_*ngas.enpernuc;
}

double calc_ion_pressure(struct parameters satdata, struct sf_params sparams, 
        double aa_, double del_, double n0_, double np_)
{
    double epsp;
    double eion_pp, eion_pm;
    double deiondnp;

    epsp = np_/1000.;
    eion_pp = calc_ion_en(satdata, sparams, aa_, del_, n0_, np_+epsp);
    eion_pm = calc_ion_en(satdata, sparams, aa_, del_, n0_, np_-epsp);
    deiondnp = (eion_pp - eion_pm)/2./epsp;

    return 2.*np_*np_/aa_/(1.-del_)*deiondnp;
}

double calc_crust_ws_cell_pressure(struct parameters satdata, 
        struct sf_params sparams, 
        struct compo eq, double nb_)
{
    double np;
    double egas_pressure;
    double ion_pressure;
    double ngas_pressure;
    double ws_cell_pressure;

    np = (nb_-eq.ng)*(1.-eq.del)/2./(1.-eq.ng/eq.n0);
    egas_pressure = calc_egas_pressure(np);
    ion_pressure = calc_ion_pressure(satdata, sparams, 
            eq.aa, eq.del, eq.n0, np);
    ngas_pressure = calc_ngas_pressure(satdata, eq.ng);
    ws_cell_pressure = egas_pressure + ion_pressure + ngas_pressure;

    return ws_cell_pressure;
}

void print_state_crust(struct parameters satdata, struct sf_params sparams, 
        struct compo eq, double nb_, 
        FILE *compo, FILE *eos)
{
    double vws, rws;
    double rhob;
    double pressws;

    vws = eq.aa/(nb_-eq.ng)*(1.-eq.ng/eq.n0);
    rws = pow(3.*vws/4./PI,1./3.);

    rhob = calc_crust_ws_cell_energy_density(satdata, sparams, eq, nb_)
        *(ELEMC/1.e-19)/pow(SPEEDOFL/1.e8,2.)*1.e13;
    pressws = calc_crust_ws_cell_pressure(satdata, sparams, eq, nb_);

    fprintf (compo, "%g %g %g %g %g %g %g\n", nb_, 
            eq.aa, eq.del, eq.aa*(1.-eq.del)/2., eq.n0, eq.ng, rws);
    fprintf(eos, "%g %g\n", rhob, pressws);
}
