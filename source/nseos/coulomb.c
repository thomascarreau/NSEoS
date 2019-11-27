#include <math.h>

#include "coulomb.h"
#include "mathconst.h"
#include "modeling.h"
#include "nuclear_en.h"
#include "phyconst.h"

#define CM (0.895929255682) // bcc lattice; Table 2.4 of Haensel book
#define U1 (0.5113875)      // bcc lattice; Table 2.4 of Haensel book

double calc_coulomb_en(double aa_, double ii_, double n0_) {
  double rsat;
  double ac;
  double zz;

  rsat = pow(3. / 4. / PI / n0_, 1. / 3.);
  ac   = 3. / 5. * ALPHAFS * HBARC / rsat;
  zz   = aa_ * (1. - ii_) / 2.;

  return ac * zz * zz * pow(aa_, -1. / 3.);
}

double calc_fd(double u_, int d_) {
  if (d_ == 2)
    return 0.25 * (log(1. / u_) - 1. + u_);
  else
    return 1. / (d_ + 2.) *
           (2. / (d_ - 2.) * (1. - d_ * pow(u_, 1. - 2. / d_) / 2.) + u_);
}

double calc_lattice_en(double aa_, double ii_, double n0_, double np_) {
  double rsat;
  double rpt;
  double zz;

  rsat = pow(3. / 4. / PI / n0_, 1. / 3.);
  rpt  = 2. * np_ / (1. - ii_) / n0_;
  zz   = aa_ * (1. - ii_) / 2.;

  return -CM * ALPHAFS * HBARC / rsat * zz * zz * pow(aa_, -1. / 3.) *
         pow(rpt, 1. / 3.);
}

double calc_lattice_en_for_tm(double zz_, double np_) {
  double vws;
  double an;

  vws = zz_ / np_;
  an  = pow(4. * PI / 3. / vws, -1. / 3.);

  return -CM * zz_ * zz_ * ALPHAFS * HBARC / an;
}

double calc_finite_size_contrib(
    double aa_, double ii_, double n0_, double np_) {
  double rsat;
  double rpt;
  double zz;

  rsat = pow(3. / 4. / PI / n0_, 1. / 3.);
  rpt  = 2. * np_ / (1. - ii_) / n0_;
  zz   = aa_ * (1. - ii_) / 2.;

  return ALPHAFS * HBARC / rsat * zz * zz * pow(aa_, -1. / 3.) * 3. / 10. * rpt;
}

double calc_zp_en(double zz_, double np_, double mi_) {
  double hbaromega_p;
  double vws;

  vws         = zz_ / np_;
  hbaromega_p = sqrt(
      pow(HBARC, 2.) * 4. * PI * pow(zz_, 2.) * ALPHAFS * HBARC / mi_ / vws);

  return 1.5 * hbaromega_p * U1;
}

double calc_harmonic_contrib(double zz_, double np_, double mi_, double tt_) {
  // see: Baiko et al. (2001) for details
  if (tt_ == 0) {
    return 0.;
  } else {
    double alpha1 = 0.932446;
    double alpha2 = 0.334547;
    double alpha3 = 0.265764;
    double alpha6 = 4.757014e-3;
    double alpha8 = 4.7770935e-3;

    double a[9] = {1., 0.1839, 0.593586, 5.4814e-3, 5.01813e-4, 0., 3.9247e-7,
        0., 5.8356e-11};
    double b[8] = {
        261.66, 0., 7.07997, 0., 0.0409484, 3.97355e-4, 5.11148e-5, 2.19749e-6};

    double theta = calc_zp_en(zz_, np_, mi_) / 1.5 / U1 / tt_;

    double sum_aa = 0.;
    double sum_bb =
        alpha6 * a[6] * pow(theta, 9.) + alpha8 * a[8] * pow(theta, 11.);

    for (int n = 0; n < 9; n++) {
      sum_aa += a[n] * pow(theta, n);
      if (n < 8) {
        sum_bb += b[n] * pow(theta, n);
      }
    }

    double sum_ln = log(1. - exp(-alpha1 * theta)) +
                    log(1. - exp(-alpha2 * theta)) +
                    log(1. - exp(-alpha3 * theta));

    return tt_ * (sum_ln - sum_aa / sum_bb);
  }
}

double calc_anharmonic_contrib(double zz_, double np_, double tt_) {
  // see: eq. (2.117) of "Neutron Stars 1: Equation of State and Structure"
  double a[3] = {10.9, 247.0, 1.765e5};

  double sum = 0.;

  double vws   = zz_ / np_;
  double an    = pow(4. / 3. * PI / vws, -1. / 3.);
  double gamma = zz_ * zz_ * ALPHAFS * HBARC / tt_ / an;

  for (int p = 1; p < 4; p++) {
    sum += a[p - 1] / p / pow(gamma, p);
  }

  return -tt_ * sum;
}

double calc_translational_free_en(
    double zz_, double np_, double mi_, double tt_) {
  // see: eq. (2.71) of "Neutron Stars 1: Equation of State and Structure"
  double vws     = zz_ / np_;
  double lambdai = pow(2. * PI * HBARC * HBARC / mi_ / tt_, 0.5);
  double gi      = 1.;

  return tt_ * (log(pow(lambdai, 3.) / gi / vws) - 1.);
}

double calc_total_coulomb_contrib(double zz_, double np_, double tt_) {
  // see: Potekhin and Chabrier (2000) for details
  double a1 = -0.9070;
  double a2 = 0.62954;
  double a3 = -sqrt(3.) / 2. - a1 / sqrt(a2);
  double b1 = 4.56e-3;
  double b2 = 211.6;
  double b3 = -1.e-4;
  double b4 = 4.62e-3;

  double vws   = zz_ / np_;
  double an    = pow(4. / 3. * PI / vws, -1. / 3.);
  double gamma = zz_ * zz_ * ALPHAFS * HBARC / tt_ / an;

  return tt_ *
         (a1 * (pow(gamma * (a2 + gamma), 0.5) -
                   a2 * log(pow(gamma / a2, 0.5) + pow(1. + gamma / a2, 0.5))) +
             2. * a3 * (pow(gamma, 0.5) - atan(pow(gamma, 0.5))) +
             b1 * (gamma - b2 * log(1. + gamma / b2)) +
             b3 / 2. * log(1. + gamma * gamma / b4));
}
