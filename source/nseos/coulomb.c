#include <math.h>

#include "nuclear_en.h"
#include "modeling.h"
#include "mathconst.h"
#include "phyconst.h"
#include "coulomb.h"

#define CM (0.895929255682) // bcc lattice; Table 2.4 of Haensel book
#define U1 (0.5113875) // bcc lattice; Table 2.4 of Haensel book

double calc_coulomb_en(struct parameters satdata, double aa_, double ii_)
{
    double rsat;
    double ac;
    double zz;

    rsat = pow(3./4./PI/satdata.rhosat0,1./3.);
    ac = 3./5.*ALPHAFS*HBARC/rsat;
    zz = aa_*(1.-ii_)/2.;

    return ac*zz*zz*pow(aa_,-1./3.);
}

double calc_lattice_en(struct parameters satdata, 
        double aa_, double ii_, double n0_, double np_)
{
    double rsat;
    double rpt;
    double zz;

    rsat = pow(3./4./PI/satdata.rhosat0,1./3.);
    rpt = 2.*np_/(1.-ii_)/n0_;
    zz = aa_*(1.-ii_)/2.;

    return ALPHAFS*HBARC/rsat*zz*zz*pow(aa_,-1./3.)
        *(-CM*pow(rpt,1./3.) + 3./10.*rpt);
}

double calc_zp_en(struct parameters satdata, struct sf_params sparams, 
        double aa_, double ii_, double n0_, double np_)
{
    double hbaromega_p;
    double zz;
    double mi;
    double vws;

    zz = aa_*(1.-ii_)/2.;
    mi = CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, aa_, ii_, n0_)
        + zz*RMP + (aa_-zz)*RMN;
    vws = zz/np_;

    hbaromega_p = sqrt(pow(HBARC,2.)*4.*PI*pow(zz,2.)*ALPHAFS*HBARC
            /mi/vws);

    return 1.5*hbaromega_p*U1;
}

double calc_harmonic_contrib(
        struct parameters satdata, struct sf_params sparams, 
        double aa_, double del_, double n0_, double np_, 
        double tt_)
{
    // see: Baiko et al. (2001) for details
    double alpha1 = 0.932446;
    double alpha2 = 0.334547;
    double alpha3 = 0.265764;
    double alpha6 = 4.757014e-3;
    double alpha8 = 4.7770935e-3;

    double a[9] = {1., 0.1839, 0.593586, 
        5.4814e-3, 5.01813e-4, 0., 
        3.9247e-7, 0., 5.8356e-11};
    double b[8] = {261.66, 0., 7.07997, 0., 
        0.0409484, 3.97355e-4, 
        5.11148e-5, 2.19749e-6};

    double theta = calc_zp_en(satdata, sparams, aa_, del_, n0_, np_)
        /1.5/U1/tt_;

    double sum_aa = 0.;
    double sum_bb = alpha6*a[6]*pow(theta,9.) + alpha8*a[8]*pow(theta,11.);

    for(int n = 0; n < 9; n++)
    {
        sum_aa += a[n]*pow(theta,n);
        if (n < 8)
        {
            sum_bb += b[n]*pow(theta,n);
        }
    }

    double sum_ln = log(1.-exp(-alpha1*theta))
        + log(1.-exp(-alpha2*theta))
        + log(1.-exp(-alpha3*theta));

    return tt_*(sum_ln - sum_aa/sum_bb);
}

double calc_translational_free_en(
        struct parameters satdata, struct sf_params sparams, 
        double aa_, double del_, double n0_, double np_, 
        double tt_)
{
    // see: eq. (2.71) of "Neutron Stars 1: Equation of State and Structure"
    double zz = aa_*(1.-del_)/2.;
    double vws = zz/np_;
    double mi = 
        CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, aa_, del_, n0_)
        + zz*RMP + (aa_-zz)*RMN;
    double lambdai = pow(2.*PI*HBARC*HBARC/mi/tt_,0.5);
    double gi = 1.;

    return tt_*(log(pow(lambdai,3.)/gi/vws) - 1.);
}

double calc_total_coulomb_contrib(
        double zz_, double np_, double tt_)
{
    // see: Potekhin and Chabrier (2000) for details
    double a1 = -0.9070;
    double a2 = 0.62954;
    double a3 = -sqrt(3.)/2. - a1/sqrt(a2);
    double b1 = 4.56e-3;
    double b2 = 211.6;
    double b3 = -1.e-4;
    double b4 = 4.62e-3;

    double vws = zz_/np_;
    double a = pow(4./3.*PI/vws,-1./3.);
    double gamma = zz_*zz_*ALPHAFS*HBARC/tt_/a;

    return tt_*(a1*(pow(gamma*(a2+gamma),0.5) 
                - a2*log(pow(gamma/a2,0.5) + pow(1.+gamma/a2,0.5)))
            + 2.*a3*(pow(gamma,0.5) - atan(pow(gamma,0.5)))
            + b1*(gamma - b2*log(1. + gamma/b2))
            + b3/2.*log(1. + gamma*gamma/b4));
}
