#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

// some constants
static const float hbarc = 197.3269788; // [MeV fm]
static const float alphafs = 7.2973525664e-3; // fine-structure constant
static const float pi = 3.14159265359; // pi
static const float pi2 = 9.86960440109; // pi squared
static const float rmn = 939.5654133; // [MeV]

// empirical parameters reference values (see: arXiv:1708.06894v3)
static const float nsat = 0.154; // [/fm^3]
static const float lasat = -16.04; // [MeV]
static const float ksat = 255.91; // [MeV]
static const float qsat = 0.; // [MeV]
static const float zsat = 0.; // [MeV]
static const float jsym = 33.43; // [MeV]
static const float lsym = 77.92; // [MeV]
static const float ksym = -2.19; // [MeV]
static const float qsym = 0.; // [MeV]
static const float zsym = 0.; // [MeV]
static const float barm = 1./0.7-1.; // 1/effm - 1
static const float cfin = 59.; // [MeV fm^5]

struct diffop
{
    double function;
    double gradient;
    double laplacian;
};
double calc_basym(double aa_, double zz_);
double calc_n0(double aa_, double zz_);
double calc_n0p(double aa_, double zz_);
struct diffop calc_n(double r_, double aa_, double zz_);
struct diffop calc_np(double r_, double aa_, double zz_);
double calc_k(double r_, double aa_, double zz_);
double calc_meta_model_low_density_correction(int order, double xx_);
double calc_hmeta(double r_, double aa_, double zz_);
double calc_hfin(double r_, double aa_, double zz_); // w/o dfin term
double calc_h(double r_, double aa_, double zz_);
struct f_params
{
    double aa;
    double zz;
};
double f(double x, void * params);
double calc_energy(double aa_, double zz_);
double calc_coulomb_energy(double aa_, double zz_);
double calc_binding_energy_per_nucleon(double aa_, double zz_);
double calc_bulk_energy_per_nucleon(double aa_, double zz_);
double calc_surface_energy_per_nucleon(double aa_, double zz_);
double calc_evacpernucleon_etf_para(double aa_, double zz_);

int main(void)
{

    //========================== TABLE ============================
    double aa;
    double zz;
    int iaa;
    int izz;
    double surface_energy_per_nucleon;
    for(iaa = 40; iaa < 209; iaa += 4)
    {
        for(izz = iaa/4; izz < iaa/2+1; izz += 1)
        {
            aa = iaa;
            zz = izz;
            surface_energy_per_nucleon = calc_surface_energy_per_nucleon(aa,zz);
            printf("%g %g %g\n", aa, zz, surface_energy_per_nucleon);
        }
    }

    return 0;
}

double calc_basym(double aa_, double zz_)
{
    double ii;
    double rsat;
    double ac;
    double qq;

    ii = 1.-2.*zz_/aa_; 
    rsat = pow(3./4./pi/nsat,1./3.);
    ac = 3./5.*alphafs*hbarc/rsat;
    qq = 144.5*jsym/(lsym + 55.5); 

    return (ii + 3.*ac*zz_*zz_/8./qq/pow(aa_,5./3.))/(1.+9.*jsym/4./qq/pow(aa_,1./3.));
    /* return ii; */ // to explore symmetric nuclei according to the 1st paper
}

double calc_n0(double aa_, double zz_)
{
    double basym;

    basym = calc_basym(aa_,zz_);

    return nsat*(1.-3.*lsym*basym*basym/(ksat+ksym*basym*basym));
}

double calc_n0p(double aa_, double zz_)
{
    double basym;
    double n0;

    basym = calc_basym(aa_,zz_);
    n0 = calc_n0(aa_,zz_);

    return n0*(1.-basym)/2.;
}

struct diffop calc_n(double r_, double aa_, double zz_)
{
    struct diffop mass_density_profile;
    double n0;
    double r0;
    double rhs;
    double basym;
    double a; // see: Papakonstantinou et al. (2013)
    double R;

    n0 = calc_n0(aa_,zz_);
    r0 = pow(3./4./pi/n0,1./3.);
    rhs = pow(aa_,1./3.)*r0;
    basym = calc_basym(aa_,zz_);
    a = 0.53+1.06*basym*basym; // [fm]
    R = rhs*(1.-pi2/3.*pow(a/rhs,2.));

    mass_density_profile.function = n0/(1.+exp((r_-R)/a));
    mass_density_profile.gradient = -n0/2./a/(1.+cosh((R-r_)/a));
    mass_density_profile.laplacian = mass_density_profile.gradient/a*sinh((R-r_)/a)/(1.+cosh((R-r_)/a))
        + 2./r_*mass_density_profile.gradient;

    return mass_density_profile;
}

struct diffop calc_np(double r_, double aa_, double zz_)
{
    struct diffop proton_density_profile;
    double n0p;
    double r0p;
    double rhsp;
    double basym;
    double ap; // see: Papakonstantinou et al. (2013)
    double Rp;

    n0p = calc_n0p(aa_,zz_);
    r0p = pow(3./4./pi/n0p,1./3.);
    rhsp = pow(zz_,1./3.)*r0p;
    basym = calc_basym(aa_,zz_);
    ap = 0.53 + 0.35*basym*basym; // [fm]
    Rp = rhsp*(1.-pi2/3.*pow(ap/rhsp,2.));

    proton_density_profile.function = n0p/(1.+exp((r_-Rp)/ap));
    proton_density_profile.gradient = -n0p/2./ap/(1.+cosh((Rp-r_)/ap));
    proton_density_profile.laplacian = proton_density_profile.gradient/ap*sinh((Rp-r_)/ap)/(1.+cosh((Rp-r_)/ap))
        + 2./r_*proton_density_profile.gradient;

    return proton_density_profile;
}

double calc_k(double r_, double aa_, double zz_)
{
    // calculation of k0
    double k0;
    struct diffop mass_density;
    struct diffop proton_density;
    float t0fac, t0fg;
    double tmp, xx;
    double delta;

    mass_density = calc_n(r_,aa_,zz_);
    proton_density = calc_np(r_,aa_,zz_);
    delta = 1.-2.*proton_density.function/mass_density.function;
    t0fac = 3.*pi2/2.*nsat;
    t0fg = 3./10./rmn*(pow(t0fac,2./3.))*(pow(hbarc,2.));
    tmp = log(mass_density.function) - log(nsat) - log(3.);
    xx = exp(tmp) - 1./3.;
    k0 = mass_density.function/2.*t0fg*(1.+barm*(1.+3.*xx))
        *(pow(1.+delta,5./3.) + pow(1.-delta,5./3.))*pow(mass_density.function/nsat,2./3.);

    // calculation of k2
    double k2;
    double k2n;
    double k2p;
    double tau2n_l, tau2n_nl, tau2n;
    double tau2p_l, tau2p_nl, tau2p;
    struct diffop f_eff;

    tau2n_l = 1./36.*pow(mass_density.gradient-proton_density.gradient,2.)/(mass_density.function-proton_density.function)
        + (mass_density.laplacian-proton_density.laplacian)/3.;
    f_eff.function = 1.+barm*(1.+3.*xx);
    f_eff.gradient = barm/nsat*mass_density.gradient;
    f_eff.laplacian = barm/nsat*mass_density.laplacian;
    tau2n_nl = 1./6.*(mass_density.gradient-proton_density.gradient)*f_eff.gradient/f_eff.function
        + 1./6.*(mass_density.function-proton_density.function)*f_eff.laplacian/f_eff.function
        - 1./12.*(mass_density.function-proton_density.function)*pow(f_eff.gradient/f_eff.function,2.);
    tau2n = tau2n_l + tau2n_nl;
    k2n = hbarc*hbarc/2./rmn*f_eff.function*tau2n;

    tau2p_l = 1./36.*pow(proton_density.gradient,2.)/(proton_density.function)
        + (proton_density.laplacian)/3.;
    tau2p_nl = 1./6.*(proton_density.gradient)*f_eff.gradient/f_eff.function
        + 1./6.*(proton_density.function)*f_eff.laplacian/f_eff.function
        - 1./12.*(proton_density.function)*pow(f_eff.gradient/f_eff.function,2.);
    tau2p = tau2p_l + tau2p_nl;
    k2p = hbarc*hbarc/2./rmn*f_eff.function*tau2p;
    k2 = k2n + k2p;

    return k0+k2;
}

double calc_meta_model_low_density_correction(int order, double xx_)
{
    double bb;
    double bexp;
    double corr;
    bb = 10.*log(2);
    bexp = exp(-bb*(3.*xx_+1.));
    corr = 1. - pow(-3.*xx_,2+1-order)*bexp;
    return corr;
}

double calc_hmeta(double r_, double aa_, double zz_)
{
    struct diffop mass_density;
    struct diffop proton_density;
    double tmp, xx;
    double delta;
    float t0fac, t0fg;
    double a00, a10, a20, a30, a40;
    double a02, a12, a22, a32, a42;
    double u0, u1, u2, u3, u4;

    t0fac = 3.*pi2/2.*nsat;
    t0fg = 3./10./rmn*(pow(t0fac,2./3.))*(pow(hbarc,2.));

    a00 = lasat - t0fg*(1. + barm);
    a10 = - t0fg*(2. + 5.*barm);
    a20 = ksat - 2.*t0fg*(5.*barm - 1.);
    a30 = qsat - 2.*t0fg*(4.-5.*barm);
    a40 = zsat - 8.*t0fg*(-7.+5.*barm);
    a02 = jsym - 5./9.*t0fg*(1. + barm);
    a12 = lsym - 5./9.*t0fg*(2. + 5.*barm);
    a22 = ksym - 10./9.*t0fg*(-1. + 5.*barm);
    a32 = qsym - 10./9.*t0fg*(4.-5.*barm);
    a42 = zsym - 40./9.*t0fg*(-7.+5.*barm);

    mass_density = calc_n(r_,aa_,zz_);
    proton_density = calc_np(r_,aa_,zz_);
    tmp = log(mass_density.function) - log(nsat) - log(3.);
    xx = exp(tmp) - 1./3.;
    delta = 1.-2.*proton_density.function/mass_density.function;

    u0 = calc_meta_model_low_density_correction(0, xx);
    u1 = calc_meta_model_low_density_correction(1, xx);
    u2 = calc_meta_model_low_density_correction(2, xx);
    /* u3 = calc_meta_model_low_density_correction(3, xx); */
    /* u4 = calc_meta_model_low_density_correction(4, xx); */
    u3 = 0.; // up to N=2
    u4 = 0.;

    return mass_density.function*(a00*u0 + a10*xx*u1 + 0.5*a20*xx*xx*u2
            + 1./6.*a30*xx*xx*xx*u3 + 1./24.*a40*xx*xx*xx*xx*u4
            + delta*delta*(a02*u0 + a12*xx*u1 + 0.5*a22*xx*xx*u2
                + 1./6.*a32*xx*xx*xx*u3 + 1./24.*a42*xx*xx*xx*xx*u4));
}

double calc_hfin(double r_, double aa_, double zz_)
{
    struct diffop mass_density;

    mass_density = calc_n(r_,aa_,zz_);

    return cfin*pow(mass_density.gradient,2.);
}

double calc_h(double r_, double aa_, double zz_)
{
    double k, hmeta, hfin;

    k = calc_k(r_,aa_,zz_);
    hmeta = calc_hmeta(r_,aa_,zz_);
    hfin = calc_hfin(r_,aa_,zz_);

    return k + hmeta + hfin;
}

double f(double x, void * params_ptr)
{
    struct f_params * params = (struct f_params *) params_ptr;
    double aa = (params->aa);
    double zz = (params->zz);
    double h = calc_h(x,aa,zz);
    double my_integrand = 4.*pi*x*x*h;

    return my_integrand;
}

double calc_energy(double aa_, double zz_)
{
    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

    double lower_limit = 0.;      /* start integral from 0 */
    double upper_limit = 100.;      /* end integral at 100 */
    double abs_error = 1.0e-8;    /* to avoid round-off problems */
    double rel_error = 1.0e-8;    /* the result will usually be much better */
    double result;                /* the result from the integration */
    double error;                 /* the estimated error from the integration */
    struct f_params params = {aa_, zz_};

    gsl_function My_function;
    My_function.function = &f;
    My_function.params = &params;
    gsl_integration_qag (&My_function, lower_limit, upper_limit,
            abs_error, rel_error, 1000, 6, work_ptr, &result,
            &error);

    return result;
}

double calc_coulomb_energy(double aa_, double zz_)
{
    double rsat;
    double ac;
    double coulomb_energy;

    rsat = pow(3./4./pi/nsat,1./3.);
    ac = 3./5.*alphafs*hbarc/rsat;
    coulomb_energy = ac*zz_*zz_/pow(aa_,1./3.);

    return coulomb_energy;
}

double calc_binding_energy_per_nucleon(double aa_, double zz_)
{
    double eb_plus_es;
    double ec;

    eb_plus_es = calc_energy(aa_,zz_);
    ec = calc_coulomb_energy(aa_,zz_);
    
    return (eb_plus_es + ec)/aa_;
}

double calc_bulk_energy_per_nucleon(double aa_, double zz_)
{
    double k0_sat;
    float t0fac, t0fg;
    double tmp, xx;
    double delta;
    double saturation_density;

    delta = calc_basym(aa_,zz_);
    saturation_density = calc_n0(aa_,zz_);
    t0fac = 3.*pi2/2.*nsat;
    t0fg = 3./10./rmn*(pow(t0fac,2./3.))*(pow(hbarc,2.));
    tmp = log(saturation_density) - log(nsat) - log(3.);
    xx = exp(tmp) - 1./3.;
    k0_sat = saturation_density/2.*t0fg*(1.+barm*(1.+3.*xx))
        *(pow(1.+delta,5./3.) + pow(1.-delta,5./3.))*pow(saturation_density/nsat,2./3.);

    double a00, a10, a20, a30, a40;
    double a02, a12, a22, a32, a42;
    double u0, u1, u2, u3, u4;
    double hmeta_sat;

    t0fac = 3.*pi2/2.*nsat;
    t0fg = 3./10./rmn*(pow(t0fac,2./3.))*(pow(hbarc,2.));

    a00 = lasat - t0fg*(1. + barm);
    a10 = - t0fg*(2. + 5.*barm);
    a20 = ksat - 2.*t0fg*(5.*barm - 1.);
    a30 = qsat - 2.*t0fg*(4.-5.*barm);
    a40 = zsat - 8.*t0fg*(-7.+5.*barm);
    a02 = jsym - 5./9.*t0fg*(1. + barm);
    a12 = lsym - 5./9.*t0fg*(2. + 5.*barm);
    a22 = ksym - 10./9.*t0fg*(-1. + 5.*barm);
    a32 = qsym - 10./9.*t0fg*(4.-5.*barm);
    a42 = zsym - 40./9.*t0fg*(-7.+5.*barm);

    u0 = calc_meta_model_low_density_correction(0, xx);
    u1 = calc_meta_model_low_density_correction(1, xx);
    u2 = calc_meta_model_low_density_correction(2, xx);
    /* u3 = calc_meta_model_low_density_correction(3, xx); */
    /* u4 = calc_meta_model_low_density_correction(4, xx); */
    u3 = 0.; // up to N=2
    u4 = 0.;

    hmeta_sat = saturation_density*(a00*u0 + a10*xx*u1 + 0.5*a20*xx*xx*u2
            + 1./6.*a30*xx*xx*xx*u3 + 1./24.*a40*xx*xx*xx*xx*u4
            + delta*delta*(a02*u0 + a12*xx*u1 + 0.5*a22*xx*xx*u2
                + 1./6.*a32*xx*xx*xx*u3 + 1./24.*a42*xx*xx*xx*xx*u4));

    return (k0_sat + hmeta_sat)/saturation_density;
}

double calc_surface_energy_per_nucleon(double aa_, double zz_)
{
    double energy;
    double energy_per_nucleon;
    double bulk_energy_per_nucleon;
    
    energy = calc_energy(aa_,zz_);
    energy_per_nucleon = energy/aa_;
    bulk_energy_per_nucleon = calc_bulk_energy_per_nucleon(aa_,zz_);

    return energy_per_nucleon - bulk_energy_per_nucleon;
}

double calc_evacpernucleon_etf_para(double aa_, double zz_)
{
    double basym;
    double bulk_energy_per_nucleon;
    double sigmas, ss;
    double n0, r0;
    double surface_energy_per_nucleon;
    double coulomb_energy;
    double coulomb_energy_per_nucleon;
    double evacpernucleon_etf_para;

    basym = calc_basym(aa_,zz_);
    bulk_energy_per_nucleon = calc_bulk_energy_per_nucleon(aa_,zz_);
    sigmas = 1.0723; // [MeV/fm^2]
    ss = 47.7291; // [MeV]
    n0 = calc_n0(aa_,zz_);
    r0 = pow(3./4./pi/n0,1./3.);
    surface_energy_per_nucleon = (4.*pi*r0*r0*sigmas + ss*basym*basym)*pow(aa_,-1./3.);
    coulomb_energy = calc_coulomb_energy(aa_,zz_);
    coulomb_energy_per_nucleon = coulomb_energy/aa_;
    evacpernucleon_etf_para = bulk_energy_per_nucleon + surface_energy_per_nucleon + coulomb_energy_per_nucleon;

    return evacpernucleon_etf_para;
}
