#include <math.h>

#include <gsl/gsl_spline.h>

#include "mathconst.h"
#include "phyconst.h"
#include "tov.h"

#define G_CGS (G * 1000.)
#define SPEEDOFL_CGS (SPEEDOFL * 100.)
#define MSUN_CGS (1.989e33)
#define P_FACTOR_NU_TO_CGS (1.6022e33) // MeV/fm^3 to dyn/cm^2

double calc_dm(double rho_, double r_, double dr_) {
  return 4. * PI * r_ * r_ * rho_ * dr_;
}

double calc_dp(double rho_, double p_, double r_, double dr_, double m_) {
  return -G_CGS * m_ * rho_ / r_ / r_ *
         (1. + p_ / rho_ / SPEEDOFL_CGS / SPEEDOFL_CGS) *
         (1. + 4. * PI * r_ * r_ * r_ * p_ / m_ / SPEEDOFL_CGS / SPEEDOFL_CGS) /
         (1. - 2. * G_CGS * m_ / r_ / SPEEDOFL_CGS / SPEEDOFL_CGS) * dr_;
}

double calc_dw(
    double rho_, double p_, double r_, double dr_, double m_, double w_) {
  return (4. * PI * G_CGS / SPEEDOFL_CGS / SPEEDOFL_CGS *
                 (rho_ * SPEEDOFL_CGS * SPEEDOFL_CGS + p_) * (4. + w_) * r_ /
                 (SPEEDOFL_CGS * SPEEDOFL_CGS - 2. * G_CGS * m_ / r_) -
             w_ / r_ * (3. + w_)) *
         dr_;
}

double calc_dy(double rho_, double p_, double r_, double dr_, double m_,
    double y_, double vs2_) {
  double q1, q2, q;

  q1 = 4. * PI * G_CGS *
       ((5. - y_) * rho_ + (9. + y_) * p_ / SPEEDOFL_CGS / SPEEDOFL_CGS +
           (p_ / SPEEDOFL_CGS / SPEEDOFL_CGS + rho_) /
               (vs2_ / SPEEDOFL_CGS / SPEEDOFL_CGS)) /
       (SPEEDOFL_CGS * SPEEDOFL_CGS - 2. * G_CGS * m_ / r_);
  q2 =
      pow(2. * G_CGS *
              (m_ + 4. * PI * r_ * r_ * r_ * p_ / SPEEDOFL_CGS / SPEEDOFL_CGS) /
              r_ / (r_ * SPEEDOFL_CGS * SPEEDOFL_CGS - 2 * G_CGS * m_),
          2.);
  q = q1 - q2;

  return (-y_ * y_ / r_ -
             (y_ - 6.) / (r_ - 2. * G_CGS * m_ / SPEEDOFL_CGS / SPEEDOFL_CGS) -
             r_ * q) *
         dr_;
}

double calc_moment_of_inertia(double r_, double w_) {
  return SPEEDOFL_CGS * SPEEDOFL_CGS / G_CGS * w_ * pow(r_, 3.) /
         (6. + 2. * w_); // see: Phys. Rep. 621, 2016, 127
}

double calc_normalized_moment_of_inertia_approx(double r_, double m_) {
  double beta = G_CGS * m_ / r_ / SPEEDOFL_CGS / SPEEDOFL_CGS;
  double i_over_mr2 =
      0.237 *
      (1. + 2.844 * beta + 18.91 * pow(beta, 4.)); // Phys. Rep. 621, 2016, 127

  return i_over_mr2;
}

double calc_normalized_crustal_moment_of_inertia_approx(double r_, double m_,
    double i_over_mr2, double epst, double pt, double rcore) {
  double rs = 2. * G_CGS * m_;
  double icrust_over_mr2 =
      16. * PI / 3. * pow(rcore, 6.) * pt * P_FACTOR_NU_TO_CGS / rs *
      (1. - rs / r_ / SPEEDOFL_CGS / SPEEDOFL_CGS * i_over_mr2) *
      (1. +
          48. / 5. * (rcore / rs * SPEEDOFL_CGS * SPEEDOFL_CGS - 1.) *
              (pt / epst)) /
      m_ / r_ / r_;

  return icrust_over_mr2;
}

double calc_tidal_love_number(double r_, double m_, double y_) {
  double beta;
  double denom;

  beta  = G_CGS * m_ / r_ / SPEEDOFL_CGS / SPEEDOFL_CGS; // compactness
  denom = 6. * beta * (2. - y_ + beta * (5. * y_ - 8.)) // see: arXiv:1512.07820
          +
          4. * pow(beta, 3.) * (13. - 11. * y_ + beta * (3. * y_ - 2.) +
                                   2. * beta * beta * (1. + y_)) +
          3. * pow(1. - 2. * beta, 2.) * (2. - y_ + 2. * beta * (y_ - 1.)) *
              log(1. - 2. * beta);

  return 8. / 5. * pow(beta, 5.) * pow(1. - 2. * beta, 2.) *
         (2. - y_ + 2. * beta * (y_ - 1.)) / denom;
}

double calc_dimensionless_tidal_deformability(
    double r_, double m_, double k2_) {
  double beta = G_CGS * m_ / r_ / SPEEDOFL_CGS / SPEEDOFL_CGS; // compactness

  return 2. * k2_ / 3. * pow(beta, -5.);
}

double get_observable_for_a_given_mass(
    double m, double mm, double mp, double om, double op) {
  return (op - om) / (mp - mm) * (m - mm) + om; // stupid linear interpolation
}

double solve_tov_equation(int lines, double pt, FILE *eos,
    struct tov_solution *tovs_m, double fixed_m, FILE *tov) {
  double Rho[lines], P[lines];
  double Rho_tmp, P_tmp;

  int j = 0;
  int l = lines;

  // to avoid 'insufficient number of points for interpolation type' error
  if (l < 10)
    return 0.;

  // reading the EoS table
  for (int i = 0; i < lines; i++) {
    if (j < l) {
      fscanf(eos, "%lf %lf", &Rho_tmp, &P_tmp);
      P_tmp *= P_FACTOR_NU_TO_CGS;
    }

    if (j == 0) {
      Rho[j] = Rho_tmp;
      P[j]   = P_tmp;
      j += 1;
    } else if (j > 0 && Rho_tmp > Rho[j - 1] && P_tmp > P[j - 1])
    // to avoid interpolation issues
    {
      Rho[j] = Rho_tmp;
      P[j]   = P_tmp;
      j += 1;
    } else
      l -= 1;
  }

  gsl_interp_accel *acc    = gsl_interp_accel_alloc();
  gsl_spline *      spline = gsl_spline_alloc(gsl_interp_linear, l);

  double rhosat = 2.3e14; // saturation density in g/cm^3
  double dr     = 100.;   // radius step in cm
  double rhoc, pc;
  double rho, p, vs2;
  double m, r;
  double w;
  double y;
  double rho_sav, p_sav, m_sav, r_sav, w_sav;
  double mcore, rcore, wcore;
  double i_over_mr2, icrust_over_mr2;
  double k2, lambda_dimless;

  double rhoc_last            = 0.;
  double pc_last              = 0.;
  double r_last               = 0.;
  double m_last               = 0.;
  double i_over_mr2_last      = 0.;
  double rcore_last           = 0.;
  double mcore_last           = 0.;
  double icrust_over_mr2_last = 0.;
  double k2_last              = 0.;
  double lambda_dimless_last  = 0.;
  tovs_m->mmax                = 0.;
  tovs_m->rhoc                = 0.;
  tovs_m->pc                  = 0.;
  tovs_m->r                   = 0.;
  tovs_m->i_over_mr2          = 0;
  tovs_m->rcore               = 0;
  tovs_m->mcore               = 0;
  tovs_m->icrust_over_mr2     = 0;
  tovs_m->k2                  = 0;
  tovs_m->lambda_dimless      = 0;

  int N = 100; // number of points

  if (Rho[l - 1] < rhosat) // to avoid the 'interpolation error' error
    return 0.;

  for (int j = 0; j < N; j++) {
    rhoc = rhosat + ((double)j / (double)N) * (Rho[l - 1] - rhosat);

    gsl_spline_init(spline, Rho, P, l);
    pc = gsl_spline_eval(spline, rhoc, acc);

    rho     = rhoc;
    rho_sav = rho;
    p       = pc;
    p_sav   = p;

    m     = 0.;
    r     = 10.;
    w     = 0.;
    y     = 2.;
    m_sav = m;
    r_sav = r;
    w_sav = w;

    gsl_spline_init(spline, P, Rho, l);

    while (rho > 2.e5) {
      r += dr;
      m += calc_dm(rho, r, dr);
      p += calc_dp(rho, p, r, dr, m);
      w += calc_dw(rho, p, r, dr, m, w);

      if (pt * P_FACTOR_NU_TO_CGS != P[0]) {
        if (p_sav > pt * P_FACTOR_NU_TO_CGS && p < pt * P_FACTOR_NU_TO_CGS) {
          mcore = (m_sav + m) / 2.;
          rcore = (r_sav + r) / 2.;
          wcore = (w_sav + w) / 2.;
        }
      } else {
        mcore = m;
        rcore = r;
        wcore = w;
      }

      if (p > P[0])
        rho = gsl_spline_eval(spline, p, acc);
      else
        break;

      vs2 = (p_sav - p) / (rho_sav - rho);
      y += calc_dy(rho_sav, p, r, dr, m, y, vs2);

      rho_sav = rho;
      p_sav   = p;
      r_sav   = r;
      m_sav   = m;
      w_sav   = w;
    }

    i_over_mr2 = calc_moment_of_inertia(r, w) / m / r / r;
    if (pt * P_FACTOR_NU_TO_CGS != P[0]) {
      icrust_over_mr2 =
          i_over_mr2 - calc_moment_of_inertia(rcore, wcore) / m / r / r;
    } else
      icrust_over_mr2 = 0.;

    k2             = calc_tidal_love_number(r, m, y);
    lambda_dimless = calc_dimensionless_tidal_deformability(r, m, k2);

    r /= 100000.;
    m /= MSUN_CGS;
    rcore /= 100000.;
    mcore /= MSUN_CGS;

    if (m > tovs_m->mmax)
      tovs_m->mmax = m;

    if (m > fixed_m && m_last < fixed_m) {
      tovs_m->rhoc =
          get_observable_for_a_given_mass(fixed_m, m_last, m, rhoc_last, rhoc);
      tovs_m->pc =
          get_observable_for_a_given_mass(fixed_m, m_last, m, pc_last, pc);
      tovs_m->r =
          get_observable_for_a_given_mass(fixed_m, m_last, m, r_last, r);
      tovs_m->i_over_mr2 = get_observable_for_a_given_mass(
          fixed_m, m_last, m, i_over_mr2_last, i_over_mr2);
      tovs_m->rcore = get_observable_for_a_given_mass(
          fixed_m, m_last, m, rcore_last, rcore);
      tovs_m->mcore = get_observable_for_a_given_mass(
          fixed_m, m_last, m, mcore_last, mcore);
      tovs_m->icrust_over_mr2 = get_observable_for_a_given_mass(
          fixed_m, m_last, m, icrust_over_mr2_last, icrust_over_mr2);
      tovs_m->k2 =
          get_observable_for_a_given_mass(fixed_m, m_last, m, k2_last, k2);
      tovs_m->lambda_dimless = get_observable_for_a_given_mass(
          fixed_m, m_last, m, lambda_dimless_last, lambda_dimless);
    }

    rhoc_last            = rhoc;
    pc_last              = pc;
    r_last               = r;
    m_last               = m;
    i_over_mr2_last      = i_over_mr2;
    rcore_last           = rcore;
    mcore_last           = mcore;
    icrust_over_mr2_last = icrust_over_mr2;
    k2_last              = k2;
    lambda_dimless_last  = lambda_dimless;

    fprintf(tov, "%g %g %g %g %g %g %g %g %g %g\n", rhoc, pc, r, m, rcore,
        mcore, i_over_mr2, icrust_over_mr2 / i_over_mr2, k2, lambda_dimless);
  }

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  return tovs_m->mmax;
}

int eval_observables_from_glitch_activity(int lines, double pt, double epst,
    FILE *eos, double g, double sigma_g, FILE *f_pulsar) {
  double Rho[lines], P[lines];
  double Rho_tmp, P_tmp;

  int j = 0;
  int l = lines;

  // to avoid 'insufficient number of points for interpolation type' error
  if (l < 10)
    return 0;

  // reading the EoS table
  for (int i = 0; i < lines; i++) {
    if (j < l) {
      fscanf(eos, "%lf %lf", &Rho_tmp, &P_tmp);
      P_tmp *= P_FACTOR_NU_TO_CGS;
    }

    if (j == 0) {
      Rho[j] = Rho_tmp;
      P[j]   = P_tmp;
      j += 1;
    } else if (j > 0 && Rho_tmp > Rho[j - 1] && P_tmp > P[j - 1])
    // to avoid interpolation issues
    {
      Rho[j] = Rho_tmp;
      P[j]   = P_tmp;
      j += 1;
    } else
      l -= 1;
  }

  gsl_interp_accel *acc    = gsl_interp_accel_alloc();
  gsl_spline *      spline = gsl_spline_alloc(gsl_interp_linear, l);

  double rhosat = 2.3e14; // saturation density in g/cm^3
  double dr     = 100.;   // radius step in cm
  double rhoc, pc;
  double rho, p, vs2;
  double m, r;
  double y;
  double rho_sav, p_sav, m_sav, r_sav;
  double mcore, rcore;
  double i_over_mr2, icrust_over_mr2;
  double k2, lambda_dimless;

  double gp         = g + sigma_g;
  double gm         = g - sigma_g;
  double gp_checker = 0;
  double gm_checker = 0;

  int N = 100; // number of points

  if (Rho[l - 1] < rhosat) // to avoid the 'interpolation error' error
    return 0.;

  for (int j = 0; j < N; j++) {
    rhoc = rhosat + ((double)j / (double)N) * (Rho[l - 1] - rhosat);

    gsl_spline_init(spline, Rho, P, l);
    pc = gsl_spline_eval(spline, rhoc, acc);

    rho     = rhoc;
    rho_sav = rho;
    p       = pc;
    p_sav   = p;

    m     = 0.;
    r     = 10.;
    y     = 2.;
    m_sav = m;
    r_sav = r;

    gsl_spline_init(spline, P, Rho, l);

    while (rho > 2.e5) {
      r += dr;
      m += calc_dm(rho, r, dr);
      p += calc_dp(rho, p, r, dr, m);

      if (pt * P_FACTOR_NU_TO_CGS != P[0]) {
        if (p_sav > pt * P_FACTOR_NU_TO_CGS && p < pt * P_FACTOR_NU_TO_CGS) {
          mcore = (m_sav + m) / 2.;
          rcore = (r_sav + r) / 2.;
        }
      } else {
        mcore = m;
        rcore = r;
      }

      if (p > P[0])
        rho = gsl_spline_eval(spline, p, acc);
      else
        break;

      vs2 = (p_sav - p) / (rho_sav - rho);
      y += calc_dy(rho_sav, p, r, dr, m, y, vs2);

      rho_sav = rho;
      p_sav   = p;
      r_sav   = r;
      m_sav   = m;
    }

    i_over_mr2 = calc_normalized_moment_of_inertia_approx(r, m);
    if (pt * P_FACTOR_NU_TO_CGS != P[0]) {
      icrust_over_mr2 = calc_normalized_crustal_moment_of_inertia_approx(
          r, m, i_over_mr2, epst, pt, rcore);
    } else
      icrust_over_mr2 = 0.;

    k2             = calc_tidal_love_number(r, m, y);
    lambda_dimless = calc_dimensionless_tidal_deformability(r, m, k2);

    r /= 100000.;
    m /= MSUN_CGS;
    rcore /= 100000.;
    mcore /= MSUN_CGS;

    if (gp_checker == 0 && icrust_over_mr2 / i_over_mr2 < gp && m > 1.3) {
      fprintf(f_pulsar, "%g %g %g %g %g %g %g %g %g %g\n", rhoc, pc, r, m,
          rcore, mcore, i_over_mr2, icrust_over_mr2 / i_over_mr2, k2,
          lambda_dimless);
      gp_checker = 1;
    }

    if (gm_checker == 0 && icrust_over_mr2 / i_over_mr2 < gm && m > 1.3) {
      fprintf(f_pulsar, "%g %g %g %g %g %g %g %g %g %g\n", rhoc, pc, r, m,
          rcore, mcore, i_over_mr2, icrust_over_mr2 / i_over_mr2, k2,
          lambda_dimless);
      gm_checker = 1;
    }

    if (gp_checker + gm_checker == 2) {
      return gp_checker + gm_checker;
    }
  }

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  return gp_checker + gm_checker;
}

int equation_for_m2(const gsl_vector *x, void *params, gsl_vector *f) {
  double m1     = ((struct rparams_tov *)params)->m1;
  double mchirp = ((struct rparams_tov *)params)->mchirp;

  const double x0 = gsl_vector_get(x, 0);
  const double y0 = mchirp - pow(m1 * x0, 3. / 5.) / pow(m1 + x0, 1. / 5.);

  gsl_vector_set(f, 0, y0);

  return GSL_SUCCESS;
}

double m2_for_m1_mchirp(const double mchirp, const double m1, double *guess) {
  double m2;

  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *           s;

  int    status;
  size_t iter = 0;

  double m2_old, m2_new, mstep;

  struct rparams_tov p;
  p.m1     = m1;
  p.mchirp = mchirp;

  double       m2_init = *guess;
  const size_t n       = 1;
  gsl_vector * x       = gsl_vector_alloc(n);
  gsl_vector_set(x, 0, m2_init);

  gsl_multiroot_function f = {&equation_for_m2, n, &p};

  T = gsl_multiroot_fsolver_dnewton;
  s = gsl_multiroot_fsolver_alloc(T, 1);
  gsl_multiroot_fsolver_set(s, &f, x);

  do {
    m2_old = gsl_vector_get(s->x, 0);

    iter++;

    status = gsl_multiroot_fsolver_iterate(s);

    m2_new = gsl_vector_get(s->x, 0);
    mstep  = m2_new - m2_old;

    if (status)
      break;

    // dirty backstepping
    int count = 0;
    while (gsl_vector_get(s->x, 0) < 0.) {
      mstep  = mstep / 4.;
      m2_new = m2_old + mstep;
      gsl_vector_set(x, 0, m2_new);
      gsl_multiroot_fsolver_set(s, &f, x);
      m2_new = gsl_vector_get(s->x, 0);
      count += 1;
      if (count > 100) {
        m2     = NAN;
        *guess = m2;
        gsl_multiroot_fsolver_free(s);
        gsl_vector_free(x);
        return m2;
      }
    }

    status = gsl_multiroot_test_residual(s->f, 9e-9);
  }

  while (status == GSL_CONTINUE && iter < 100);

  m2 = gsl_vector_get(s->x, 0);

  *guess = m2;

  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);

  return m2;
}

void dimensionless_lambda1_lambda2_for_m1_mchirp(const int lines,
    const char *eos, const double mchirp, const double m1, double *m2,
    struct Lambda *Lambdai) {
  double guess = m1 - 0.1;

  *m2 = m2_for_m1_mchirp(mchirp, m1, &guess);

  double pt = -1; // we do not care about crust properties here
  struct tov_solution TOVSolutionM1, TOVSolutionM2;

  FILE *tov   = fopen("tov.out", "w+");
  FILE *myeos = fopen(eos, "r");
  solve_tov_equation(lines, pt, myeos, &TOVSolutionM1, m1, tov);
  fclose(myeos);
  myeos = fopen(eos, "r");
  solve_tov_equation(lines, pt, myeos, &TOVSolutionM2, *m2, tov);
  fclose(myeos);
  fclose(tov);

  Lambdai->s1 = TOVSolutionM1.lambda_dimless;
  Lambdai->s2 = TOVSolutionM2.lambda_dimless;
}

void dimensionless_lambda1_lambda2_relation(const int lines, const char *eos,
    const double mchirp, const char *outfile) {
  double pt = -1; // we do not care about crust properties here
  struct tov_solution tovs;

  FILE * tov   = fopen("tov.out", "w+");
  FILE * myeos = fopen(eos, "r");
  double mmax  = solve_tov_equation(lines, pt, myeos, &tovs, 1.4, tov);
  fclose(myeos);
  fclose(tov);

  double m1 = mmax;
  double m2;

  struct Lambda Lambdai;

  FILE *output_file = fopen(outfile, "w+");

  while (m1 > 1.0) {
    dimensionless_lambda1_lambda2_for_m1_mchirp(
        lines, eos, mchirp, m1, &m2, &Lambdai);

    fprintf(output_file, "%6f %6f %6f %6f\n", m1, m2, Lambdai.s1, Lambdai.s2);

    m1 -= 0.02;
  }

  fclose(output_file);
}
