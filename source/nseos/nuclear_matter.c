#include <math.h>
#include <complex.h>

#include "nuclear_matter.h"

double calc_meta_model_low_density_correction(struct parameters satdata, int max_order, int order, double xx_)
{
    double bexp;
    bexp = exp(-satdata.b*(3.*xx_+1.));
    return 1. - pow(-3.*xx_,max_order+1-order)*bexp;
}

double calc_meta_model_low_density_correction_derivative(struct parameters satdata, int max_order, int order, double xx_)
{
    double bexp;
    bexp = exp(-satdata.b*(3.*xx_+1.));
    return bexp*pow(3.,-order+max_order+1)*pow(-xx_,max_order-order)
        *(-order - 3.*satdata.b*xx_ + max_order + 1.);
}

double calc_meta_model_low_density_correction_second_derivative(struct parameters satdata, int max_order, int order, double xx_)
{
    double bexp;
    bexp = exp(-satdata.b*(3.*xx_+1.));
    return bexp*(-pow(-3.,-order+max_order+1))*pow(-xx_,-order+max_order-1)
        *(9.*satdata.b*satdata.b*xx_*xx_ + (-order + max_order + 1)*(-order - 6.*satdata.b*xx_ + max_order));
}

struct hnm calc_meta_model_nuclear_matter(struct parameters satdata, int max_order, double nn_, double ii_)
{
    double tmp;
    double xx;
    float t0fac, t0fg;
    double a00, a10, a20, a30, a40;
    double barfac;
    double a02, a12, a22, a32, a42;
    double rmns, rmps;
    double u0, u1, u2, u3, u4;
    double u0p, u1p, u2p, u3p, u4p;
    double u0pp, u1pp, u2pp, u3pp, u4pp;
    double epotpernuc;
    double dekinpernucdx;
    double depotpernucdx;
    double denpernucdx;
    double dekinpernucdi;
    double depotpernucdi;
    double denpernucdi;
    struct hnm result;

    t0fac = 3.*PI2/2.*satdata.rhosat0;
    t0fg = 3./10./RMN*(pow(t0fac,2./3.))*(pow(HBARC,2.));

    // coeff
    a00 = satdata.lasat0 - t0fg*(1. + satdata.barm);
    a10 = - t0fg*(2. + 5.*satdata.barm);
    a20 = satdata.ksat0 - 2.*t0fg*(5.*satdata.barm - 1.);
    barfac = (satdata.barm + 3.*satdata.bardel);
    a02 = satdata.jsym0 - 5./9.*t0fg*(1. + barfac);
    a12 = satdata.lsym0 - 5./9.*t0fg*(2. + 5.*barfac);
    a22 = satdata.ksym0 - 10./9.*t0fg*(-1. + 5.*barfac);

    // x_bulk
    tmp = log(nn_) - log(satdata.rhosat0) - log(3.);
    xx = exp(tmp) - 1./3.;

    // correction at nn=0
    u0 = calc_meta_model_low_density_correction(satdata, max_order, 0, xx);
    u1 = calc_meta_model_low_density_correction(satdata, max_order, 1, xx);
    u2 = calc_meta_model_low_density_correction(satdata, max_order, 2, xx);
    u0p = calc_meta_model_low_density_correction_derivative(satdata, max_order, 0, xx);
    u1p = calc_meta_model_low_density_correction_derivative(satdata, max_order, 1, xx);
    u2p = calc_meta_model_low_density_correction_derivative(satdata, max_order, 2, xx);
    u0pp = calc_meta_model_low_density_correction_second_derivative(satdata, max_order, 0, xx);
    u1pp = calc_meta_model_low_density_correction_second_derivative(satdata, max_order, 1, xx);
    u2pp = calc_meta_model_low_density_correction_second_derivative(satdata, max_order, 2, xx);

    rmns = RMN/(1.+ (satdata.barm + ii_*satdata.bardel)*(1.+3.*xx));
    rmps = RMN/(1.+ (satdata.barm - ii_*satdata.bardel)*(1.+3.*xx));

    dekinpernucdx = t0fg*cpow((1.+3.*xx),-1./3.)
        *( cpow((1.+ii_),5./3.)*(1.+(satdata.barm+satdata.bardel*ii_)*(1.+3.*xx))
                +  cpow((1.-ii_),5./3.)*(1.+(satdata.barm-satdata.bardel*ii_)*(1.+3.*xx)) )
        + 0.5*t0fg*cpow((1.+3.*xx),2./3.)
        *( cpow((1.+ii_),5./3.)*3.*(satdata.barm+satdata.bardel*ii_)
                + cpow((1.-ii_),5./3.)*3.*(satdata.barm-satdata.bardel*ii_));

    dekinpernucdi = 5./6.*t0fg*cpow((1.+3.*xx),2./3.)
        *(cpow((1.+ii_),2./3.)*RMN/rmns - cpow((1.-ii_),2./3.)*RMN/rmps);

    if(max_order == 3)
    {
        a30 = satdata.qsat0 - 2.*t0fg*(4.-5.*satdata.barm);
        a32 = satdata.qsym0 - 10./9.*t0fg*(4.-5.*barfac);
        u3 = calc_meta_model_low_density_correction(satdata, max_order, 3, xx);
        u3p = calc_meta_model_low_density_correction_derivative(satdata, max_order, 3, xx);
        u3pp = calc_meta_model_low_density_correction_second_derivative(satdata, max_order, 3, xx);

        epotpernuc = a00*u0 + a10*xx*u1 + 0.5*a20*xx*xx*u2 + 1./6.*a30*xx*xx*xx*u3
            + (a02*u0 + a12*xx*u1 + 0.5*a22*xx*xx*u2 + 1./6.*a32*xx*xx*xx*u3)*ii_*ii_;
        depotpernucdx = a00*u0p + a10*u1 + a10*xx*u1p
            +  a20*u2*xx + 0.5*a20*xx*xx*u2p 
            + a30/6.*(3.*xx*xx*u3 + xx*xx*xx*u3p)
            +  ii_*ii_*(a02*u0p + a12*u1 + a12*xx*u1p
                    + a22*xx*u2 + 0.5*a22*xx*xx*u2p
                    + a32/6.*(3.*xx*xx*u3 + xx*xx*xx*u3p));
        depotpernucdi = 2.*ii_*(a02*u0 + a12*xx*u1 + 0.5*a22*xx*xx*u2 + 1./6.*a32*xx*xx*xx*u3);

        result.jsym = 5./9.*t0fg*pow(1.+3.*xx,2./3.)*(1.+satdata.barm*(1.+3.*xx))
            + a02*u0 + a12*xx*u1 + 0.5*a22*xx*xx*u2 + 1./6*a32*xx*xx*xx*u3;

        result.lsym = 5./9.*t0fg*(1.+3.*xx)*(2.*pow(1.+3.*xx,-1./3.)*(1.+satdata.barm*(1.+3.*xx))
                + 3.*satdata.barm*pow(1.+3.*xx,2./3.))
            + (1.+3.*xx)*(a02*u0p + a12*(u1 + xx*u1p) + a22/2.*(2.*xx*u2 + xx*xx*u2p)
                    + a32/6.*(3.*xx*xx*u3 + xx*xx*xx*u3p));

        result.ksym = 10./9.*t0fg*pow(1.+3.*xx,2.)*(6.*satdata.barm*pow(1.+3.*xx,-1./3.)
                - pow(1.+3.*xx,-4./3.)*(1.+satdata.barm*(1.+3.*xx)))
            + pow(1.+3.*xx,2.)*(a02*u0pp + a12*(2.*u1p + xx*u1pp) + a22/2.*(2.*u2 + 4.*xx*u2p + xx*xx*u2pp)
                    + a32/6.*(6.*xx*u3 + 6.*xx*xx*u3p + xx*xx*xx*u3pp));
    }
    else if(max_order == 4)
    {
        a30 = satdata.qsat0 - 2.*t0fg*(4.-5.*satdata.barm);
        a40 = satdata.zsat0 - 8.*t0fg*(-7.+5.*satdata.barm);
        a32 = satdata.qsym0 - 10./9.*t0fg*(4.-5.*barfac);
        a42 = satdata.zsym0 - 40./9.*t0fg*(-7.+5.*barfac);
        u3 = calc_meta_model_low_density_correction(satdata, max_order, 3, xx);
        u4 = calc_meta_model_low_density_correction(satdata, max_order, 4, xx);
        u3p = calc_meta_model_low_density_correction_derivative(satdata, max_order, 3, xx);
        u4p = calc_meta_model_low_density_correction_derivative(satdata, max_order, 4, xx);
        u3pp = calc_meta_model_low_density_correction_second_derivative(satdata, max_order, 3, xx);
        u4pp = calc_meta_model_low_density_correction_second_derivative(satdata, max_order, 4, xx);

        epotpernuc = a00*u0 + a10*xx*u1 + 0.5*a20*xx*xx*u2 + 1./6.*a30*xx*xx*xx*u3 + 1./24.*a40*xx*xx*xx*xx*u4
            + (a02*u0 + a12*xx*u1 + 0.5*a22*xx*xx*u2 + 1./6.*a32*xx*xx*xx*u3 + 1./24.*a42*xx*xx*xx*xx*u4)*ii_*ii_;
        depotpernucdx = a00*u0p + a10*u1 + a10*xx*u1p
            +  a20*u2*xx + 0.5*a20*xx*xx*u2p 
            + a30/6.*(3.*xx*xx*u3 + xx*xx*xx*u3p)
            + a40/24.*(4.*xx*xx*xx*u4+xx*xx*xx*xx*u4p)
            +  ii_*ii_*(a02*u0p + a12*u1 + a12*xx*u1p
                    + a22*xx*u2 + 0.5*a22*xx*xx*u2p
                    + a32/6.*(3.*xx*xx*u3 + xx*xx*xx*u3p)
                    + a42/24.*(4.*xx*xx*xx*u4+xx*xx*xx*xx*u4p));
        depotpernucdi = 2.*ii_*(a02*u0 + a12*xx*u1 + 0.5*a22*xx*xx*u2 + 1./6.*a32*xx*xx*xx*u3 + 1./24.*a42*xx*xx*xx*xx*u4);

        result.jsym = 5./9.*t0fg*pow(1.+3.*xx,2./3.)*(1.+satdata.barm*(1.+3.*xx))
            + a02*u0 + a12*xx*u1 + 0.5*a22*xx*xx*u2 + 1./6.*a32*xx*xx*xx*u3 + 1./24.*a42*xx*xx*xx*xx*u4;

        result.lsym = 5./9.*t0fg*(1.+3.*xx)*(2.*pow(1.+3.*xx,-1./3.)*(1.+satdata.barm*(1.+3.*xx))
                + 3.*satdata.barm*pow(1.+3.*xx,2./3.))
            + (1.+3.*xx)*(a02*u0p + a12*(u1 + xx*u1p) + a22/2.*(2.*xx*u2 + xx*xx*u2p)
                    + a32/6.*(3.*xx*xx*u3 + xx*xx*xx*u3p) + a42/24.*(4.*xx*xx*xx*u4 + xx*xx*xx*xx*u4p));

        result.ksym = 10./9.*t0fg*pow(1.+3.*xx,2.)*(6.*satdata.barm*pow(1.+3.*xx,-1./3.)
                - pow(1.+3.*xx,-4./3.)*(1.+satdata.barm*(1.+3.*xx)))
            + pow(1.+3.*xx,2.)*(a02*u0pp + a12*(2.*u1p + xx*u1pp) + a22/2.*(2.*u2 + 4.*xx*u2p + xx*xx*u2pp)
                    + a32/6.*(6.*xx*u3 + 6.*xx*xx*u3p + xx*xx*xx*u3pp)
                    + a42/24.*(12.*xx*xx*u4 + 8.*xx*xx*xx*u4p + xx*xx*xx*xx*u4pp));
    }
    else
    {
        epotpernuc = a00*u0 + a10*xx*u1 + 0.5*a20*xx*xx*u2
            + (a02*u0 + a12*xx*u1 + 0.5*a22*xx*xx*u2)*ii_*ii_;
        depotpernucdx = a00*u0p + a10*u1 + a10*xx*u1p
            +  a20*u2*xx + 0.5*a20*xx*xx*u2p
            +  ii_*ii_*(a02*u0p + a12*u1 + a12*xx*u1p
                    + a22*xx*u2 + 0.5*a22*xx*xx*u2p);
        depotpernucdi = 2.*ii_*(a02*u0 + a12*xx*u1 + 0.5*a22*xx*xx*u2);

        result.jsym = 5./9.*t0fg*pow(1.+3.*xx,2./3.)*(1.+satdata.barm*(1.+3.*xx))
            + a02*u0 + a12*xx*u1 + 0.5*a22*xx*xx*u2;

        result.lsym = 5./9.*t0fg*(1.+3.*xx)*(2.*pow(1.+3.*xx,-1./3.)*(1.+satdata.barm*(1.+3.*xx))
                + 3.*satdata.barm*pow(1.+3.*xx,2./3.))
            + (1.+3.*xx)*(a02*u0p + a12*(u1 + xx*u1p) + a22/2.*(2.*xx*u2 + xx*xx*u2p));

        result.ksym = 10./9.*t0fg*pow(1.+3.*xx,2.)*(6.*satdata.barm*pow(1.+3.*xx,-1./3.)
                - pow(1.+3.*xx,-4./3.)*(1.+satdata.barm*(1.+3.*xx)))
            + pow(1.+3.*xx,2.)*(a02*u0pp + a12*(2.*u1p + xx*u1pp) + a22/2.*(2.*u2 + 4.*xx*u2p + xx*xx*u2pp));
    }

    result.enpernuc = 0.5*t0fg*cpow((1.+3.*xx),2./3.)
        *(cpow((1.+ii_),5./3.)*RMN/rmns + cpow((1.-ii_),5./3.)*RMN/rmps)
        + epotpernuc;

    denpernucdx = dekinpernucdx + depotpernucdx;
    denpernucdi = dekinpernucdi + depotpernucdi;

    result.p = 1./3.*satdata.rhosat0*pow(1.+3.*xx,2.)*denpernucdx;

    if (nn_ == 0.)
        result.mun = 0.;
    else
        result.mun = result.enpernuc + nn_*(1./3./satdata.rhosat0*denpernucdx
                + (1.-ii_)/2.*denpernucdi);

    return result;
}

double calc_asymmetry_factor(double m_, double ii_)
{
    return 0.5*(pow(1.+ii_,m_) + pow(1.-ii_,m_));
}

double calc_asymmetry_factor_derivative(double m_, double ii_)
{
    return 0.5*m_*(pow(1.+ii_,m_-1.) - pow(1.-ii_,m_-1.));
}

struct hnm calc_skyrme_nuclear_matter(struct skyrme_parameters coeff, double nn_, double ii_)
{
    double f53, f2, f83;
    double f53p, f2p, f83p;
    double denpernucdn;
    double denpernucdi;
    struct hnm result;

    f53 = calc_asymmetry_factor(5./3., ii_);
    f2 = calc_asymmetry_factor(2., ii_);
    f83 = calc_asymmetry_factor(8./3., ii_);
    f53p = calc_asymmetry_factor_derivative(5./3., ii_);
    f2p = calc_asymmetry_factor_derivative(2., ii_);
    f83p = calc_asymmetry_factor_derivative(8./3., ii_);

    result.enpernuc = 3./5.*HBARC*HBARC/2./RMN*pow(1.5*PI2,2./3.)*pow(nn_,2./3.)*f53
        + 1./8.*coeff.t0*nn_*(2.*(coeff.x0 + 2.) - (2.*coeff.x0 + 1.)*f2)
        + 1./48.*coeff.t3*pow(nn_,coeff.sigma+1.)*(2.*(coeff.x3 + 2.) - (2.*coeff.x3 + 1.)*f2)
        + 3./40.*pow(1.5*PI2,2./3.)*pow(nn_,5./3.)*((coeff.t1*(coeff.x1 + 2.) + coeff.t2*(coeff.x2 + 2.))*f53
                + 0.5*(coeff.t2*(2.*coeff.x2 + 1.) - coeff.t1*(2.*coeff.x1 + 1.))*f83);

    denpernucdn = 2./5.*HBARC*HBARC/2./RMN*pow(1.5*PI2,2./3.)*pow(nn_,-1./3.)*f53
        + 1./8.*coeff.t0*(2.*(coeff.x0 + 2.) - (2.*coeff.x0 + 1.)*f2)
        + 1./48.*(coeff.sigma+1.)*coeff.t3*pow(nn_,coeff.sigma)*(2.*(coeff.x3 + 2.) - (2.*coeff.x3 + 1.)*f2)
        + 1./8.*pow(1.5*PI2,2./3.)*pow(nn_,2./3.)*((coeff.t1*(coeff.x1 + 2.) + coeff.t2*(coeff.x2 + 2.))*f53
                + 0.5*(coeff.t2*(2.*coeff.x2 + 1.) - coeff.t1*(2.*coeff.x1 + 1.))*f83);

    denpernucdi = 3./5.*HBARC*HBARC/2./RMN*pow(1.5*PI2,2./3.)*pow(nn_,2./3.)*f53p
        + 1./8.*coeff.t0*nn_*(2.*(coeff.x0 + 2.) - (2.*coeff.x0 + 1.)*f2p)
        + 1./48.*coeff.t3*pow(nn_,coeff.sigma+1.)*(2.*(coeff.x3 + 2.) - (2.*coeff.x3 + 1.)*f2p)
        + 3./40.*pow(1.5*PI2,2./3.)*pow(nn_,5./3.)*((coeff.t1*(coeff.x1 + 2.) + coeff.t2*(coeff.x2 + 2.))*f53p
                + 0.5*(coeff.t2*(2.*coeff.x2 + 1.) - coeff.t1*(2.*coeff.x1 + 1.))*f83p);

    result.mun = result.enpernuc + nn_*(denpernucdn + (1.-ii_)/nn_*denpernucdi);

    result.p = nn_*nn_*denpernucdn;

    result.jsym = HBARC*HBARC/2./RMN/3.*pow(1.5*PI2,2./3.)*pow(nn_,2./3.) + 1./8.*coeff.t0*nn_*(1.-coeff.x0)
        + 1./48.*coeff.t3*pow(nn_,coeff.sigma+1.)*(1.-coeff.x3) + 1./24.*pow(1.5*PI2,2./3.)*pow(nn_,5./3.)
        *(coeff.t2*(5.*coeff.x2 + 4.) - 3.*coeff.t1*coeff.x1);

    result.lsym = HBARC*HBARC/RMN/3.*pow(1.5*PI2,2./3.)*pow(nn_,2./3.) + 3./8.*coeff.t0*nn_*(1.-coeff.x0)
        + 1./16.*(coeff.sigma + 1.)*coeff.t3*pow(nn_,coeff.sigma+1.)*(1.-coeff.x3) + 5./24.*pow(1.5*PI2,2./3.)*pow(nn_,5./3.)
        *(coeff.t2*(5.*coeff.x2 + 4.) - 3.*coeff.t1*coeff.x1);

    result.ksym = -HBARC*HBARC/RMN/3.*pow(1.5*PI2,2./3.)*pow(nn_,2./3.)
        + 3./16.*coeff.sigma*(coeff.sigma + 1.)*coeff.t3*pow(nn_,coeff.sigma+1.)*(1.-coeff.x3) 
        + 5./12.*pow(1.5*PI2,2./3.)*pow(nn_,5./3.)*(coeff.t2*(5.*coeff.x2 + 4.) - 3.*coeff.t1*coeff.x1);

    return result;
}
