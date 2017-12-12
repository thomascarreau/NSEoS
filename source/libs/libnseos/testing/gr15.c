#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "../nseos/hnm.h"

/* // MATHEMATICAL CONSTANTS */
/* static const float pi = 3.14159265359; */
/* static const float pi2 = 9.86960440109; */

/* // PHYSICAL CONSTANTS */
/* static const float alphafs = 7.2973525664e-3; */
/* static const float hbarc = 197.3269788; // MeV.fm */
/* static const float rmn = 939.5654133; // MeV */
/* static const float rmp = 938.2720813; // MeV */
/* static const float mel = 0.5109989461; // MeV */

// EMPIRICAL PARAMETERS (SLy4)
static const float nsat0 = 0.160; // /fm^3
static const float lasat0 = -15.97; // MeV
static const float ksat0 = 229.91; // MeV
static const float jsym0 = 32.00; // MeV
static const float lsym0 = 45.94; // MeV
static const float ksym0 = -119.73; // MeV
static const float effm0 = 0.69;

// MODEL
struct nuclear_matter
{
    double enpernuc;
    double mung;
};
struct gas
{
    double energy_density;
    double chemical_potential;
};
double calc_asymmetry_factor(double m_, double ii_);
struct nuclear_matter calc_nuclear_matter(double nn_, double ii_);
double calc_asymmetric_saturation_density(double ii_);
double calc_nuclear_energy(double aa_, double ii_, double np_); // including rest mass
struct gas calc_egas(double np_); // including rest mass
struct gas calc_ngas(double ngas_); // including rest mass

// NEWTON METHOD
struct icrust_fun_3d
{
    double f_stability;
    double f_beta;
    double f_muneq;
};
struct icrust_fun_3d calc_icrust_fun_3d(double aa_, double ii_, double ngas_, double np_);
struct rparams
{
    double np;
};
int assign_icrust_fun_3d(const gsl_vector * x, void *params, gsl_vector * f);
void print_icrust(gsl_multiroot_fsolver * s, double nb_);

int main(void)
{
    double inb;
    double nb;

    double aa_new_guess;
    double ii_new_guess;
    double ng_new_guess;

    aa_new_guess = 75.;
    ii_new_guess = 0.25;
    ng_new_guess = 1.e-6;

    for(inb = 3; inb <= 921; inb += 1) // last step: 1300
    {
        const gsl_multiroot_fsolver_type *T;
        gsl_multiroot_fsolver *s;

        int status;
        size_t iter = 0;

        double aa_old, aa_new, astep;
        double ii_old, ii_new, istep;
        double ng_old, ng_new, gstep;

        nb = inb/10000.;
        struct rparams p;
        p.np = nb; // modification in assign_icrust_fun_4d (DIRTY!)

        double x_init[3] = {aa_new_guess, ii_new_guess, ng_new_guess};
        const size_t n = 3;
        gsl_vector *x = gsl_vector_alloc(n);
        gsl_vector_set(x, 0, x_init[0]);
        gsl_vector_set(x, 1, x_init[1]);
        gsl_vector_set(x, 2, x_init[2]);

        gsl_multiroot_function f = {&assign_icrust_fun_3d, n, &p};

        T = gsl_multiroot_fsolver_broyden;
        s = gsl_multiroot_fsolver_alloc(T,3);
        gsl_multiroot_fsolver_set(s, &f, x);

        do
        {
            /* print_icrust(s, nb); */

            aa_old = gsl_vector_get (s->x, 0);
            ii_old = gsl_vector_get (s->x, 1);
            ng_old = gsl_vector_get (s->x, 2);

            iter++;

            status = gsl_multiroot_fsolver_iterate(s);

            aa_new = gsl_vector_get (s->x, 0);
            ii_new = gsl_vector_get (s->x, 1);
            ng_new = gsl_vector_get (s->x, 2);
            astep = aa_new - aa_old;
            istep = ii_new - ii_old;
            gstep = ng_new - ng_old;

            if (status)
                break;

            // dirty backstepping
            while (gsl_vector_get (s->x, 0) < 50. || gsl_vector_get (s->x, 0) > 100000.) {
                astep = astep/4.;
                aa_new = aa_old + astep;
                gsl_vector_set (x, 0, aa_new);
                gsl_vector_set (x, 1, ii_new);
                gsl_vector_set (x, 2, ng_new);
                gsl_multiroot_fsolver_set (s, &f, x);
                aa_new = gsl_vector_get (s->x, 0.);
            }
            while (gsl_vector_get (s->x, 1) < 0.1 || gsl_vector_get (s->x, 1) > 0.9) {
                istep = istep/4.;
                ii_new = ii_old + istep;
                gsl_vector_set (x, 0, aa_new);
                gsl_vector_set (x, 1, ii_new);
                gsl_vector_set (x, 2, ng_new);
                gsl_multiroot_fsolver_set (s, &f, x);
                ii_new = gsl_vector_get (s->x, 1);
            }
            while (gsl_vector_get (s->x, 2) < nb/1000. || gsl_vector_get (s->x, 2) > nb) { // SLy4 value
                gstep = gstep/4.;
                ng_new = ng_old + gstep;
                gsl_vector_set (x, 0, aa_new);
                gsl_vector_set (x, 1, ii_new);
                gsl_vector_set (x, 2, ng_new);
                gsl_multiroot_fsolver_set (s, &f, x);
                ng_new = gsl_vector_get (s->x, 2);
            }

            status = gsl_multiroot_test_residual (s->f, 9e-5);
        }

        while (status == GSL_CONTINUE && iter < 5000);

        print_icrust(s, nb);

        aa_new_guess = gsl_vector_get(s->x, 0);
        ii_new_guess = gsl_vector_get(s->x, 1);
        ng_new_guess = gsl_vector_get(s->x, 2);

        gsl_multiroot_fsolver_free(s);
        gsl_vector_free(x);
    }

    return 0;
}

double calc_asymmetry_factor(double m_, double ii_)
{
    return 0.5*(pow(1.+ii_,m_) + pow(1.-ii_,m_));
}

struct nuclear_matter calc_nuclear_matter(double nn_, double ii_)
{ 
    struct nuclear_matter nucl_matter;
    float t0, t1, t2, t3; 
    float sigma;
    float x0, x1, x2, x3; 
    double f53, f2, f83;

    // Skyrme interaction parameters (SLy4)
    t0 = -2488.91; 
    t1 = 486.82;
    t2 = -546.39;
    t3 = 13777.0;
    sigma = 1./6.;
    x0 = 0.834;
    x1 = -0.344;
    x2 = -1.0;
    x3 = 1.354;

    f53 = calc_asymmetry_factor(5./3., ii_);
    f2 = calc_asymmetry_factor(2., ii_);
    f83 = calc_asymmetry_factor(8./3., ii_);

    nucl_matter.enpernuc = 3./5.*hbarc*hbarc/2./rmn*pow(1.5*pi2,2./3.)*pow(nn_,2./3.)*f53
        + 1./8.*t0*nn_*(2.*(x0 + 2.) - (2.*x0 + 1.)*f2)
        + 1./48.*t3*pow(nn_,sigma+1.)*(2.*(x3 + 2.) - (2.*x3 + 1.)*f2)
        + 3./40.*pow(1.5*pi2,2./3.)*pow(nn_,5./3.)*((t1*(x1 + 2.) + t2*(x2 + 2.))*f53
                + 0.5*(t2*(2.*x2 + 1.) - t1*(2.*x1 + 1.))*f83);

    nucl_matter.mung = hbarc*hbarc/2./rmn*pow(1.5*pi2,2./3.)*pow(nn_,2./3.)*pow(2.,2./3.)
        + 1./4.*t0*nn_*(2.*(x0 + 2.) - (2.*x0 + 1.)*2.)
        + 1./48.*t3*(sigma+2.)*pow(nn_,sigma+1.)*(2.*(x3 + 2.) - (2.*x3 + 1.)*2.)
        + 3./40.*8./3.*pow(1.5*pi2,2./3.)*pow(nn_,5./3.)*((t1*(x1 + 2.) + t2*(x2 + 2.))*pow(2.,2./3.)
                + 0.5*(t2*(2.*x2 + 1.) - t1*(2.*x1 + 1.))*pow(2.,5./3.));

    return nucl_matter;
}

double calc_asymmetric_saturation_density(double ii_)
{
    return nsat0*(1. - 3.*lsym0*ii_*ii_/(ksat0 + ksym0*ii_*ii_));
}

double calc_nuclear_energy(double aa_, double ii_, double np_)
{
    double n0;
    double sym_en;
    struct nuclear_matter nucl_matter;
    double bulk_en;

    double surf_en;

    double rsat;
    double ac;
    double rpt;
    double fws;
    double coul_en;

    double nucl_en;

    // Sly4 based
    n0 = calc_asymmetric_saturation_density(ii_);
    /* sym_en = jsym0/(1.+jsym0/16.60/pow(aa_,1./3.)); */
    /* bulk_en = (lasat0 + sym_en*ii_*ii_)*aa_; */ 
    /* nucl_matter = calc_nuclear_matter(n0, ii_); */
    /* bulk_en = nucl_matter.enpernuc*aa_; */

    struct hnm numat;
    struct parameters satdata;
    satdata = assign_param(satdata);
    numat = calc_hnm(satdata, n0, ii_);
    bulk_en = numat.enpernuc*aa_;

    surf_en = 18.24*pow(aa_,2./3.);

    rsat = pow(3./4./pi/nsat0,1./3.);
    ac = 3./5.*alphafs*hbarc/rsat;
    rpt = 2.*np_/n0/(1.-ii_);
    fws = -1.5*pow(rpt,1./3.) + 0.5*rpt;
    coul_en = ac*(1.-ii_)*(1.-ii_)/4.*pow(aa_,5./3.)*(1.+fws);

    nucl_en = bulk_en + surf_en + coul_en
        + aa_*(1.+ii_)/2.*rmn + aa_*(1.-ii_)/2.*rmp;

    return nucl_en;
}

struct gas calc_egas(double np_)
{
    struct gas egas;
    double t;

    t = hbarc*pow(3.0*pi2*np_,1./3.)/mel;
    egas.energy_density = pow(mel,4.)/8./pi2/pow(hbarc,3.)*((2.*t*t+1.)*t*pow(t*t+1.,1./2.) - log(t + pow(t*t+1.,1./2.)));
    egas.chemical_potential = pow(mel,3.)/8./pow(3.*np_*pi2,2./3.)/pow(hbarc,2.)
        *(pow(t*t+1.,1./2.)*(1.+6.*t*t) + t*t*(2.*t*t+1.)/pow(t*t+1.,1./2.)
                - 1./pow(t*t+1.,1./2.));

    return egas;
}

struct gas calc_ngas(double ngas_)
{
    /* struct nuclear_matter nucl_matter; */
    struct hnm nucl_matter;
    struct parameters satdata;
    satdata = assign_param(satdata);
    struct gas ngas;

    nucl_matter = calc_hnm(satdata, ngas_ ,1.);
    ngas.energy_density = ngas_*(nucl_matter.enpernuc + rmn);
    ngas.chemical_potential = nucl_matter.mug + rmn;

    return ngas;
}

struct icrust_fun_3d calc_icrust_fun_3d(double aa_, double ii_, double ngas_, double np_)
{
    struct icrust_fun_3d functions;

    double enuc;
    double epsa;
    double epsi;
    double enuc_ap, enuc_am;
    double enuc_ip, enuc_im;
    double denucdaa;
    double denucdii;
    enuc = calc_nuclear_energy(aa_, ii_, np_);
    epsa = 0.001;
    epsi = 0.0001;
    enuc_ap = calc_nuclear_energy(aa_+epsa, ii_, np_);
    enuc_am = calc_nuclear_energy(aa_-epsa, ii_, np_);
    enuc_ip = calc_nuclear_energy(aa_, ii_+epsi, np_);
    enuc_im = calc_nuclear_energy(aa_, ii_-epsi, np_);
    denucdaa = (enuc_ap - enuc_am)/2./epsa;
    denucdii = (enuc_ip - enuc_im)/2./epsi;

    struct gas egas;
    struct gas ngas;
    double mueltot;
    double lambda;
    double pgas;
    egas = calc_egas(np_);
    mueltot = egas.chemical_potential;
    ngas = calc_ngas(ngas_);
    lambda = ngas.chemical_potential;
    pgas = ngas_*lambda - ngas.energy_density;

    double n0, dn0dii;
    n0 = calc_asymmetric_saturation_density(ii_);
    dn0dii = -6.*nsat0*lsym0*ksat0*ii_/pow((ksat0 + ksym0*ii_*ii_),2.);

    double rsat;
    double rpt;
    double ac;
    double delta;
    rpt = 2.*np_/n0/(1.-ii_);
    rsat = pow(3./4./pi/nsat0,1./3.);
    ac = 3./5.*alphafs*hbarc/rsat;
    delta = 0.5*ac*pow(aa_,2./3.)*(1.-ii_)*(1.-ii_)/4.*rpt*(pow(rpt,-2./3.) - 1.);

    functions.f_stability = denucdaa/aa_ - enuc/aa_/aa_;
    functions.f_beta = 2./aa_*denucdii - mueltot + 2./(1.-ii_)*delta - 2.*dn0dii*pgas/n0/n0;
    functions.f_muneq = enuc/aa_ + (1.-ii_)/aa_*denucdii - lambda + pgas/n0 - dn0dii/n0/n0*(1.-ii_)*pgas;

    return functions;
}

int assign_icrust_fun_3d(const gsl_vector * x, void *params, gsl_vector * f)
{
    double np = ((struct rparams *) params)->np;

    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);
    const double x2 = gsl_vector_get (x, 2);

    double n0 = calc_asymmetric_saturation_density(x1);

    np = (np-x2)*(1.-x1)/2./(1.-x2/n0);

    struct icrust_fun_3d functs;
    functs = calc_icrust_fun_3d(x0, x1, x2, np);

    const double y0 = functs.f_stability;
    const double y1 = functs.f_beta;
    const double y2 = functs.f_muneq;

    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);
    gsl_vector_set(f, 2, y2);

    return GSL_SUCCESS;
}

void print_icrust(gsl_multiroot_fsolver * s, double nb_)
{
    double aa, ii, ng;
    double zz, n0;
    double f0, f1, f2;

    aa = gsl_vector_get(s->x, 0);
    ii = gsl_vector_get(s->x, 1);
    ng = gsl_vector_get(s->x, 2);
    zz = aa*(1.-ii)/2.;
    n0 = calc_asymmetric_saturation_density(ii);
    f0 = gsl_vector_get(s->f, 0);
    f1 = gsl_vector_get(s->f, 1);
    f2 = gsl_vector_get(s->f, 2);

    printf("%g %g %g %g %g %g %g %g %g\n", nb_, aa, zz, ii, n0, ng, f0, f1, f2);
}
