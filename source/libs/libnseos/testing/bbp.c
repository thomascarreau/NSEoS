#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

// MATHEMATICAL CONSTANTS
static const float pi = 3.14159265359;
static const float pi2 = 9.86960440109;

// PHYSICAL CONSTANTS
static const float alphafs = 7.2973525664e-3;
static const float hbarc = 197.3269788; // MeV.fm
static const float rmn = 939.5654133; // MeV
static const float rmp = 938.2720813; // MeV
static const float mel = 0.5109989461; // MeV

// EMPIRICAL PARAMETERS
static const float k0 = 1.43; // /fm
static const float w0 = 16.5; // MeV
static const float k = 143.; // MeV
static const float s = 33.; // MeV

// MODEL
struct gas
{
    double energy_density;
    double chemical_potential;
};
double calc_nuclear_matter_energy(double kk_, double proton_fraction_);
double calc_nuclear_surface_energy(double aa_, double ii_, double n0_, double ngas_);
double calc_coulomb_energy(double aa_, double ii_, double n0_, double np_);
double calc_nuclear_energy(double aa_, double ii_, double n0_, double ngas_,
        double np_);
struct gas calc_egas(double np_); // for completely relativistic electrons
struct gas calc_ngas(double ngas_);

// NEWTON METHOD
struct icrust_fun_4d
{
    double f_stability;
    double f_beta;
    double f_muneq;
    double f_presseq;
};
struct icrust_fun_4d calc_icrust_fun_4d(double aa_, double ii_, double n0_, double ngas_, double np_);
struct rparams
{
    double np;
};
int assign_icrust_fun_4d(const gsl_vector * x, void *params, gsl_vector * f);
void print_icrust(gsl_multiroot_fsolver * s, double nb_);

int main(void)
{

    // ======================================= PLOT NUCLEAR MATTER ENERGY =======================================

    /* int ikk; */
    /* double kk; */
    /* int ixx; */
    /* double xx; */
    /* double nuc_mat_en; */
    /* for(ikk = 0; ikk < 16; ikk++) */
    /* { */
    /*     kk = ikk/10.; */
    /*     for(ixx = 0; ixx < 51; ixx += 5) */
    /*     { */
    /*         xx = ixx/100.; */
    /*         nuc_mat_en = calc_nuclear_matter_energy(kk, xx); */
    /*         printf("%g %g %g\n", nuc_mat_en); */
    /*     } */
    /* } */


    // ======================================= PLOT NUCLEAR SURFACE ENERGY =======================================

    /* double aa; */
    /* int iii; */
    /* double ii; */
    /* double zz; */
    /* double n0; */
    /* double ng; */
    /* double esurf; */

    /* aa = 200.; */
    /* n0 = pow(k0,3.)/1.5/pi2; */

    /* for(iii = 30; iii <= 90; iii++) */
    /* { */
    /*     ii = iii/100.; */
    /*     zz = aa*(1.-ii)/2.; */
    /*     for(ng = 0.; ng < n0; ng+= n0/1000.) */
    /*     { */
    /*         esurf = calc_nuclear_surface_energy(aa, zz, n0, ng); */ 
    /*         printf("%g %g %g\n", ii, ng, esurf/pow(aa,2./3.)); */
    /*     } */
    /* } */


    // ======================================= PLOT NUCLEAR ENERGY =======================================

    /* double aa; */
    /* double iii, ii; */
    /* double n0; */
    /* double evac0, evac1, evac2, evac3, evac4; */

    /* aa = 200.; */
    /* n0 = pow(k0,3.)/1.5/pi2; */
    /* double ng[5] = {0., n0/10., n0/4., n0/2., n0}; */

    /* for (iii = 0; iii <= 90; iii++) */
    /* { */
    /*     ii = iii/100.; */
    /*     evac0 = calc_nuclear_energy(aa, ii, n0, ng[0], 0.); */
    /*     evac1 = calc_nuclear_energy(aa, ii, n0, ng[1], 0.); */
    /*     evac2 = calc_nuclear_energy(aa, ii, n0, ng[2], 0.); */
    /*     evac3 = calc_nuclear_energy(aa, ii, n0, ng[3], 0.); */
    /*     evac4 = calc_nuclear_energy(aa, ii, n0, ng[4], 0.); */
    /*     printf("%g %g %g %g %g %g\n", ii, evac0/aa, evac1/aa, */ 
    /*             evac2/aa, evac3/aa, evac4/aa); */
    /* } */


    // ======================================= NEWTON METHOD =======================================

    double nb;
    nb = 5.e-4;
    struct rparams p;
    p.np = nb; // modification in assign_icrust_fun_4d

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t iter = 0;

    double aa_old, aa_new, astep;
    double ii_old, ii_new, istep;
    double n0_old, n0_new, nstep;
    double ng_old, ng_new, gstep;

    double x_init[4] = {300., 0.56, 0.15, 1.e-4};
    const size_t n = 4;
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector_set(x, 0, x_init[0]);
    gsl_vector_set(x, 1, x_init[1]);
    gsl_vector_set(x, 2, x_init[2]);
    gsl_vector_set(x, 3, x_init[3]);

    gsl_multiroot_function f = {&assign_icrust_fun_4d, n, &p};

    T = gsl_multiroot_fsolver_dnewton;
    s = gsl_multiroot_fsolver_alloc(T,4);
    gsl_multiroot_fsolver_set(s, &f, x);

    do
    {
        print_icrust(s, nb);

        aa_old = gsl_vector_get (s->x, 0);
        ii_old = gsl_vector_get (s->x, 1);
        n0_old = gsl_vector_get (s->x, 2);
        ng_old = gsl_vector_get (s->x, 3);

        iter++;

        status = gsl_multiroot_fsolver_iterate(s);

        aa_new = gsl_vector_get (s->x, 0);
        ii_new = gsl_vector_get (s->x, 1);
        n0_new = gsl_vector_get (s->x, 2);
        ng_new = gsl_vector_get (s->x, 3);
        astep = aa_new - aa_old;
        istep = ii_new - ii_old;
        nstep = n0_new - n0_old;
        gstep = ng_new - ng_old;

        if (status)
            break;

        // dirty backstepping
        while (gsl_vector_get (s->x, 0) < 100. || gsl_vector_get (s->x, 0) > 1000.) {
            astep = astep/4.;
            aa_new = aa_old + astep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, ii_new);
            gsl_vector_set (x, 2, n0_new);
            gsl_vector_set (x, 3, ng_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            aa_new = gsl_vector_get (s->x, 0.);
        }
        while (gsl_vector_get (s->x, 1) < 0.3 || gsl_vector_get (s->x, 1) > 0.9) {
            istep = istep/4.;
            ii_new = ii_old + istep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, ii_new);
            gsl_vector_set (x, 2, n0_new);
            gsl_vector_set (x, 3, ng_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            ii_new = gsl_vector_get (s->x, 1);
        }
        while (gsl_vector_get (s->x, 2) < 0.1 || gsl_vector_get (s->x, 2) > 0.2) {
            nstep = nstep/4.;
            n0_new = n0_old + nstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, ii_new);
            gsl_vector_set (x, 2, n0_new);
            gsl_vector_set (x, 3, ng_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            n0_new = gsl_vector_get (s->x, 2);
        }
        while (gsl_vector_get (s->x, 3) < 1.e-6 || gsl_vector_get (s->x, 3) > nb) {
            gstep = gstep/4.;
            ng_new = ng_old + gstep;
            gsl_vector_set (x, 0, aa_new);
            gsl_vector_set (x, 1, ii_new);
            gsl_vector_set (x, 2, n0_new);
            gsl_vector_set (x, 3, ng_new);
            gsl_multiroot_fsolver_set (s, &f, x);
            ng_new = gsl_vector_get (s->x, 3);
        }

        status = gsl_multiroot_test_residual (s->f, 9e-9);
    }

    while (status == GSL_CONTINUE && iter < 10000);

    print_icrust(s, nb);

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return 0;
}

double calc_nuclear_matter_energy(double kk_, double proton_fraction_)
{
    // construction of the interpolation
    double bulk_sym_en; // W(k,1/2)
    double pnm_en; // W(k,0)
    double kin_en; // W_kin(k,x)
    double mup0;
    double dpnm_endk; // dW(k,0)/dk
    double mun0;
    double t;
    double ii;
    double nuc_matter_en;

    bulk_sym_en = 3.*hbarc*hbarc*kk_*kk_/10./rmn*pow(1.-kk_/k0,3.)
        - w0*pow(kk_/k0,3.)*(1.+(1.-kk_/k0)*(9.-6.*kk_/k0))
        + 0.5*k*pow(1.-kk_/k0,2.)*pow(kk_/k0,3.);

    pnm_en = 19.74*kk_*kk_
        - kk_*kk_*kk_*(40.4 - 1.088*kk_*kk_*kk_)/(1.+2.545*kk_);

    kin_en = 3.*pow(pow(2.,1./3.)*hbarc*kk_,2.)/10./rmn
        *(pow(proton_fraction_,5./3.) + pow(1.-proton_fraction_,5./3.));

    mup0 = -pow(kk_,3.)*(218.+277.*kk_)/(1+8.57*kk_*kk_);
    /* dpnm_endk = 39.48*kk_ - (-2.545*pow(kk_,3.)*(40.4-1.088*pow(kk_,3.)) */
    /*         + (1.+2.545*kk_)*(121.2*kk_*kk_-6.528*pow(kk_,5.)))/pow(1.+2.545*kk_,2.); */
    dpnm_endk = kk_*(2.13752*pow(kk_,5.) + 1.00787*pow(kk_,4.) + 7.73147*kk_*kk_
            + 12.3132*kk_ + 6.09539)/pow(kk_+0.392927,2.);
    mun0 = pnm_en + kk_/3.*dpnm_endk;

    t = hbarc*hbarc*kk_*kk_/2./rmn;
    ii = 1.-2.*proton_fraction_;
    nuc_matter_en = (bulk_sym_en - 3./5.*t)*(1. - 3.*pow(ii,4.) + 2.*pow(ii,6.))
        + (s*pow(kk_/k0,2.) - t/3.)*ii*ii*pow(1.-ii*ii,2.)
        + (pnm_en - 3./5.*t*pow(2.,2./3.))*(3.*pow(ii,4.) - 2.*pow(ii,6.))
        + 0.25*(mup0 - mun0 + t*pow(2.,2./3.))*(pow(ii,4.) - pow(ii,6.))
        + kin_en;

    return nuc_matter_en;
}

double calc_nuclear_surface_energy(double aa_, double ii_, double n0_, double ngas_)
{
    float sigma;
    double kk_in;
    double kk_out;
    double nuc_matter_en_in;
    double nuc_matter_en_out;
    double nuc_surf_en;

    sigma = 21.;

    kk_in = pow(3./2.*pi2*n0_,1./3.);
    kk_out = pow(3./2.*pi2*ngas_,1./3.);
    nuc_matter_en_in = calc_nuclear_matter_energy(kk_in, (1.-ii_)/2.);
    nuc_matter_en_out = calc_nuclear_matter_energy(kk_out, 0.);

    nuc_surf_en = sigma*(nuc_matter_en_out - nuc_matter_en_in)/w0*pow((1.-ngas_/n0_)*aa_,2./3.);

    return nuc_surf_en;
}

double calc_coulomb_energy(double aa_, double ii_, double n0_, double np_)
{
    double r0;
    double ac;
    double rpt;
    double fws;
    double coul_en;

    r0 = pow(3./4./pi/n0_,1./3.);
    ac = 3./5.*alphafs*hbarc/r0;
    rpt = 2.*np_/n0_/(1.-ii_);
    fws = -1.5*pow(rpt,1./3.) + 0.5*rpt;
    coul_en = ac*(1.-ii_)*(1.-ii_)/4.*pow(aa_,5./3.)*(1.+fws);

    return coul_en;
}

double calc_nuclear_energy(double aa_, double ii_, double n0_, double ngas_, double np_)
{
    double k0_;
    double zz_;
    double bulk_en;
    double surf_en;
    double coul_en;
    double d;
    double wthick_times_aa;
    double wexch_times_aa;
    double nucl_en;

    k0_ = pow(1.5*pi2*n0_,1./3.);
    zz_ = aa_*(1.-ii_)/2.;
    bulk_en = calc_nuclear_matter_energy(k0_, (1.-ii_)/2.);
    bulk_en = bulk_en*aa_;
    surf_en = calc_nuclear_surface_energy(aa_, ii_, n0_, ngas_);
    coul_en = calc_coulomb_energy(aa_, ii_, n0_, np_);
    d = 0.74/1.5/pi2/(n0_-ngas_);
    wthick_times_aa = -4./9.*pi*zz_*alphafs*hbarc*d*d*k0_*k0_*k0_*zz_/aa_;
    wexch_times_aa = -3./4./pi*zz_*alphafs*hbarc*pow(2.*zz_/aa_,1./3.)*k0_;
    nucl_en = bulk_en + surf_en + coul_en 
        + aa_*((1.-zz_/aa_)*rmn + zz_/aa_*rmp);
        /* + wthick_times_aa */ 
        /* + wexch_times_aa; */

    return nucl_en;
}

struct gas calc_egas(double np_)
{
    struct gas egas;

    egas.energy_density = 3./4.*hbarc*pow(1.5*pi2*np_,1./3.)*np_;
    egas.chemical_potential = hbarc*pow(1.5*pi2*np_,1./3.);

    return egas;
}

struct gas calc_ngas(double ngas_)
{
    double kn;
    double pnm_en; // W(k,0)
    double dpnm_endk; // dW(k,0)/dk
    struct gas ngas;

    kn = pow(1.5*pi2*ngas_,1./3.);
    ngas.energy_density = ngas_*(calc_nuclear_matter_energy(kn, 0.) + rmn);

    pnm_en = 19.74*kn*kn
        - kn*kn*kn*(40.4 - 1.088*kn*kn*kn)/(1.+2.545*kn);
    /* dpnm_endk = 39.48*kn - (-2.545*pow(kn,3.)*(40.4-1.088*pow(kn,3.)) */
    /*         + (1.+2.545*kn)*(121.2*kn*kn-6.528*pow(kn,5.)))/pow(1.+2.545*kn,2.); */
    dpnm_endk = kn*(2.13752*pow(kn,5.) + 1.00787*pow(kn,4.) + 7.73147*kn*kn
            + 12.3132*kn + 6.09539)/pow(kn+0.392927,2.);
    ngas.chemical_potential = pnm_en + kn/3.*dpnm_endk; // for infinite matter

    return ngas;
}

struct icrust_fun_4d calc_icrust_fun_4d(double aa_, double ii_, double n0_, double ngas_, double np_)
{
    struct icrust_fun_4d functions;
    double zz;
    double enuc;
    double epsa;
    double epsi;
    double epsn;
    double enuc_ap, enuc_am;
    double enuc_ip, enuc_im;
    double enuc_np, enuc_nm;
    double denucdaa;
    double denucdii;
    double denucdn0;
    double munn;
    struct gas egas;
    struct gas ngas;
    double epsg;
    double es_gp, es_gm;
    double desdngas;
    double mung;

    zz = aa_*(1.-ii_)/2.;

    enuc = calc_nuclear_energy(aa_, ii_, n0_, ngas_, np_);
    epsa = 0.001;
    epsi = 0.0001;
    epsn = 0.0001;
    enuc_ap = calc_nuclear_energy(aa_+epsa, ii_, n0_, ngas_, np_);
    enuc_am = calc_nuclear_energy(aa_-epsa, ii_, n0_, ngas_, np_);
    enuc_ip = calc_nuclear_energy(aa_, ii_+epsi, n0_, ngas_, np_);
    enuc_im = calc_nuclear_energy(aa_, ii_-epsi, n0_, ngas_, np_);
    enuc_np = calc_nuclear_energy(aa_, ii_, n0_+epsn, ngas_, np_);
    enuc_nm = calc_nuclear_energy(aa_, ii_, n0_-epsn, ngas_, np_);
    denucdaa = (enuc_ap - enuc_am)/2./epsa;
    denucdii = (enuc_ip - enuc_im)/2./epsi;
    denucdn0 = (enuc_np - enuc_nm)/2./epsn;

    double es, ec;
    es = calc_nuclear_surface_energy(aa_, ii_, n0_, ngas_);
    ec = calc_coulomb_energy(aa_, ii_, n0_, np_);
    functions.f_stability = es - 2.*ec;

    egas = calc_egas(np_);
    ngas = calc_ngas(ngas_);
    epsg = ngas_/100.;
    es_gp = calc_nuclear_surface_energy(aa_, ii_, n0_, ngas_+epsg);
    es_gm = calc_nuclear_surface_energy(aa_, ii_, n0_, ngas_-epsg);
    desdngas = (es_gp - es_gm)/2./epsg;
    desdngas = 0.;
    mung = ngas.chemical_potential + rmn;
        /* + 1./(zz/np_ - aa_/n0_)*desdngas; */ 
    munn = (enuc/aa_ + (1.-ii_)/aa_*denucdii - ngas.energy_density/n0_)/(1.-ngas_/n0_);

    double rpt;
    double r0;
    double ac;
    double delta;
    rpt = 2.*np_/n0_/(1.-ii_);
    r0 = pow(3./4./pi/n0_,1./3.);
    ac = 3./5.*alphafs*hbarc/r0;
    delta = 0.5*ac*pow(aa_,2./3.)*(1.-ii_)*(1.-ii_)/4.*rpt*(pow(rpt,-2./3.) - 1.);

    /* functions.f_stability = denucdaa/aa_ - enuc/aa_/aa_; */
    functions.f_beta = denucdii*2./aa_ - egas.chemical_potential - delta;
    functions.f_muneq = munn - mung;
    functions.f_presseq = n0_*n0_*denucdn0/aa_ - ngas_*mung + ngas.energy_density - ngas_*rmn;

    return functions;
}

int assign_icrust_fun_4d(const gsl_vector * x, void *params, gsl_vector * f)
{
    double np = ((struct rparams *) params)->np;

    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);
    const double x2 = gsl_vector_get (x, 2);
    const double x3 = gsl_vector_get (x, 3);

    np = (np-x3)*(1.-x1)/2./(1.-x3/x2);

    struct icrust_fun_4d functs;
    functs = calc_icrust_fun_4d(x0, x1, x2, x3, np);

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

void print_icrust(gsl_multiroot_fsolver * s, double nb_)
{
    double aa, ii, n0, ng;
    double zz;
    double f0, f1, f2, f3;

    aa = gsl_vector_get(s->x, 0);
    ii = gsl_vector_get(s->x, 1);
    n0 = gsl_vector_get(s->x, 2);
    ng = gsl_vector_get(s->x, 3);
    zz = aa*(1.-ii)/2.;
    f0 = gsl_vector_get(s->f, 0);
    f1 = gsl_vector_get(s->f, 1);
    f2 = gsl_vector_get(s->f, 2);
    f3 = gsl_vector_get(s->f, 3);

    printf("%g %g %g %g %g %g %g %g %g %g\n", nb_, aa, zz, ii, n0, ng, f0, f1, f2, f3);
}
