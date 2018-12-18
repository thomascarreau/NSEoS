#include <math.h>
#include <gsl/gsl_spline.h>

#include "mathconst.h"
#include "phyconst.h"
#include "tov.h"

#define G_CGS (G*1000.)
#define SPEEDOFL_CGS (SPEEDOFL*100.)
#define MSUN_CGS (1.989e33)
#define P_FACTOR_NU_TO_CGS (1.6022e33) // MeV/fm^3 to dyn/cm^2

double calc_dm(double rho_, double r_, double dr_)
{
    return 4.*PI*r_*r_*rho_*dr_;
}

double calc_dp(double rho_, double p_, double r_, double dr_, double m_)
{
    return -G_CGS*m_*rho_/r_/r_*(1.+p_/rho_/SPEEDOFL_CGS/SPEEDOFL_CGS)
        *(1.+4.*PI*r_*r_*r_*p_/m_/SPEEDOFL_CGS/SPEEDOFL_CGS)
        /(1.-2.*G_CGS*m_/r_/SPEEDOFL_CGS/SPEEDOFL_CGS)*dr_;
}

double calc_dy(double rho_, double p_, double r_, double dr_, double m_, 
        double y_, double vs2_)
{
    double q1, q2, q;

    q1 = 4.*PI*G_CGS*((5.-y_)*rho_+(9.+y_)*p_/SPEEDOFL_CGS/SPEEDOFL_CGS
            +(p_/SPEEDOFL_CGS/SPEEDOFL_CGS+rho_)
            /(vs2_/SPEEDOFL_CGS/SPEEDOFL_CGS))
        /(SPEEDOFL_CGS*SPEEDOFL_CGS-2.*G_CGS*m_/r_);
    q2 = pow(2.*G_CGS*(m_+4.*PI*r_*r_*r_*p_/SPEEDOFL_CGS/SPEEDOFL_CGS)
            /r_/(r_*SPEEDOFL_CGS*SPEEDOFL_CGS-2*G_CGS*m_),2.);
    q=q1-q2;

    return (-y_*y_/r_ - (y_ - 6.)/(r_-2.*G_CGS*m_/SPEEDOFL_CGS/SPEEDOFL_CGS) 
            - r_*q)*dr_;
}

double calc_tidal_love_number(double r_, double m_, double y_)
{
    double beta;
    double denom;

    beta = G_CGS*m_/r_/SPEEDOFL_CGS/SPEEDOFL_CGS; // compactness
    denom = 6.*beta*(2. - y_ + beta*(5.*y_-8.)) // see: arXiv:1512.07820
        + 4.*pow(beta,3.)*(13.-11.*y_ + beta*(3.*y_-2.) + 2.*beta*beta*(1.+y_))
        + 3.*pow(1.-2.*beta,2.)*(2. - y_ + 2.*beta*(y_-1.))*log(1.-2.*beta);

    return 8./5.*pow(beta,5.)*pow(1.-2.*beta,2.)
        *(2. - y_ + 2.*beta*(y_-1.))/denom;
}

double calc_dimensionless_tidal_deformability(double r_, double m_, 
        double k2_)
{
    double beta;

    beta = G_CGS*m_/r_/SPEEDOFL_CGS/SPEEDOFL_CGS; // compactness

    return 2.*k2_/3.*pow(beta,-5.);
}

// see: PRC82,025810(2010)
double calc_normalized_moment_of_inertia(double r_, double m_)
{
    double rs = 2.*G_CGS*m_;
    double i_over_mr2 = 0.21/(1.-rs/r_/SPEEDOFL_CGS/SPEEDOFL_CGS);
    return i_over_mr2;
}

double calc_normalized_crustal_moment_of_inertia(double r_, double m_, 
        double i_over_mr2, double epst, double pt, double rcore)
{
    double rs = 2.*G_CGS*m_;
    double icrust_over_mr2 = 16.*PI/3.*pow(rcore,6.)*pt*P_FACTOR_NU_TO_CGS/rs
        *(1. - rs/r_/SPEEDOFL_CGS/SPEEDOFL_CGS*i_over_mr2)
        *(1. + 48./5.*(rcore/rs*SPEEDOFL_CGS*SPEEDOFL_CGS - 1.)*(pt/epst))
        /m_/r_/r_;
    return icrust_over_mr2;
}

double get_observable_for_a_given_mass(double m, double mm, double mp, 
        double om, double op)
{
    return (op-om)/(mp-mm)*(m-mm) + om; // stupid linear interpolation
}

double solve_tov_equation(int lines, double pt, double epst, FILE *eos, 
        struct tov_solution *tovs_m, double fixed_m, FILE *tov)
{
    double Rho[lines], P[lines];
    double Rho_tmp, P_tmp;

    int j = 0;
    int l = lines;

    // to avoid 'insufficient number of points for interpolation type' error
    if (l < 10) 
        return 0.;

    // reading the EoS table
    for(int i = 0; i < lines; i++)
    {
        if (j < l)
        {
            fscanf(eos, "%lf %lf", &Rho_tmp, &P_tmp);
            P_tmp *= P_FACTOR_NU_TO_CGS;
        }

        if (j == 0)
        {
            Rho[j] = Rho_tmp;
            P[j] = P_tmp;
            j += 1;
        }
        else if (j > 0 
                && Rho_tmp > Rho[j-1] && P_tmp > P[j-1])
            // to avoid interpolation issues
        {
            Rho[j] = Rho_tmp;
            P[j] = P_tmp;
            j += 1;
        }
        else
            l -= 1;
    }

    gsl_interp_accel *acc
        = gsl_interp_accel_alloc ();
    gsl_spline *spline
        = gsl_spline_alloc (gsl_interp_linear, l);

    double rhosat = 2.3e14; // saturation density in g/cm^3
    double dr = 100.; // radius step in cm
    double rhoc, pc;
    double rho, p, vs2;
    double m, r;
    double y;
    double rho_sav, p_sav, m_sav, r_sav;
    double mcore, rcore;
    double i_over_mr2, icrust_over_mr2;
    double k2, lambda_dimless;

    double rhoc_last = 0.;
    double pc_last = 0.;
    double r_last = 0.;
    double m_last = 0.;
    double i_over_mr2_last = 0.;
    double rcore_last = 0.;
    double mcore_last = 0.;
    double icrust_over_mr2_last = 0.;
    double k2_last = 0.;
    double lambda_dimless_last = 0.;
    tovs_m->mmax = 0.;
    tovs_m->rhoc = 0.;
    tovs_m->pc = 0.;
    tovs_m->r = 0.;
    tovs_m->i_over_mr2 = 0;
    tovs_m->rcore = 0;
    tovs_m->mcore = 0;
    tovs_m->icrust_over_mr2 = 0;
    tovs_m->k2 = 0;
    tovs_m->lambda_dimless = 0;

    int N = 100; // number of points

    if (Rho[l-1] < rhosat) // to avoid the 'interpolation error' error
        return 0.;

    for(int j = 0; j < N; j++)
    {
        rhoc = rhosat + ((double)j / (double)N)*(Rho[l-1] - rhosat);

        gsl_spline_init (spline, Rho, P, l);
        pc = gsl_spline_eval (spline, rhoc, acc);

        rho = rhoc;
        rho_sav = rho;
        p = pc;
        p_sav = p;

        m = 0.;
        r = 10.;
        y = 2.;
        m_sav = m;
        r_sav = r;

        gsl_spline_init (spline, P, Rho, l);

        while(rho > 2.e5)
        {
            r += dr;
            m += calc_dm(rho, r, dr);
            p += calc_dp(rho, p, r, dr, m);

            if (pt*P_FACTOR_NU_TO_CGS != P[0])
            {
                if(p_sav > pt * P_FACTOR_NU_TO_CGS 
                        && p < pt * P_FACTOR_NU_TO_CGS)
                {
                    mcore = (m_sav + m)/2.;
                    rcore = (r_sav + r)/2.;
                }
            }
            else
            {
                mcore = m;
                rcore = r;
            }

            if(p > P[0])
                rho = gsl_spline_eval(spline, p, acc);
            else
                break;

            vs2 = (p_sav - p)/(rho_sav-rho);
            y += calc_dy(rho_sav, p, r, dr, m, y, vs2);

            rho_sav = rho;
            p_sav = p;
            r_sav = r;
            m_sav = m;
        }

        i_over_mr2 = calc_normalized_moment_of_inertia(r, m);
        if (pt*P_FACTOR_NU_TO_CGS != P[0])
            icrust_over_mr2 = calc_normalized_crustal_moment_of_inertia(r, m, 
                    i_over_mr2, epst, pt, rcore);
        else
            icrust_over_mr2 = 0.;

        k2 = calc_tidal_love_number(r, m, y);
        lambda_dimless = calc_dimensionless_tidal_deformability(r, m, k2);

        r /= 100000.;
        m /= MSUN_CGS;
        rcore /= 100000.;
        mcore /= MSUN_CGS;

        if (m > tovs_m->mmax)
            tovs_m->mmax = m;

        if (m > fixed_m && m_last < fixed_m)
        {
            tovs_m->rhoc = get_observable_for_a_given_mass(fixed_m, 
                    m_last, m, rhoc_last, rhoc);
            tovs_m->pc = get_observable_for_a_given_mass(fixed_m, 
                    m_last, m, pc_last, pc);
            tovs_m->r = get_observable_for_a_given_mass(fixed_m, 
                    m_last, m, r_last, r);
            tovs_m->i_over_mr2 = get_observable_for_a_given_mass(fixed_m, 
                    m_last, m, i_over_mr2_last, i_over_mr2);
            tovs_m->rcore = get_observable_for_a_given_mass(fixed_m, 
                    m_last, m, rcore_last, rcore);
            tovs_m->mcore = get_observable_for_a_given_mass(fixed_m, 
                    m_last, m, mcore_last, mcore);
            tovs_m->icrust_over_mr2 = get_observable_for_a_given_mass(fixed_m, 
                    m_last, m, icrust_over_mr2_last, icrust_over_mr2);
            tovs_m->k2 = get_observable_for_a_given_mass(fixed_m, 
                    m_last, m, k2_last, k2);
            tovs_m->lambda_dimless = get_observable_for_a_given_mass(fixed_m, 
                    m_last, m, lambda_dimless_last, lambda_dimless);
        }

        rhoc_last = rhoc;
        pc_last = pc;
        r_last = r;
        m_last = m;
        i_over_mr2_last = i_over_mr2;
        rcore_last = rcore;
        mcore_last = mcore;
        icrust_over_mr2_last = icrust_over_mr2;
        k2_last = k2;
        lambda_dimless_last = lambda_dimless;

        fprintf(tov, "%g %g %g %g %g %g %g %g %g %g\n", rhoc, pc, 
                r, m,
                rcore, mcore,
                i_over_mr2, icrust_over_mr2/i_over_mr2,
                k2, lambda_dimless);
    }

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    return tovs_m->mmax;
}
