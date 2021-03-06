#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include "nuclear_matter.h"

double calc_meta_model_low_density_correction(
    struct parameters satdata, int max_order, int order, double xx_) {
  double bexp;
  bexp = exp(-satdata.b * (3. * xx_ + 1.));
  return 1. - pow(-3. * xx_, max_order + 1 - order) * bexp;
}

double calc_meta_model_low_density_correction_derivative(
    struct parameters satdata, int max_order, int order, double xx_) {
  double bexp;
  bexp = exp(-satdata.b * (3. * xx_ + 1.));
  return bexp * pow(3., -order + max_order + 1) * pow(-xx_, max_order - order) *
         (-order - 3. * satdata.b * xx_ + max_order + 1.);
}

double calc_meta_model_low_density_correction_second_derivative(
    struct parameters satdata, int max_order, int order, double xx_) {
  double bexp;
  bexp = exp(-satdata.b * (3. * xx_ + 1.));
  return -pow(3., -order + max_order + 1) * bexp *
         pow(-xx_, -order + max_order - 1) *
         (order * order + order * (6. * satdata.b * xx_ - 2. * max_order - 1.) +
             9. * satdata.b * satdata.b * xx_ * xx_ +
             (max_order + 1.) * (max_order - 6. * satdata.b * xx_));
}

struct hnm calc_meta_model_nuclear_matter(struct parameters satdata,
    int max_order, double nn_, double ii_, double tt_) {
  double     tmp;
  double     xx;
  float      t0fac, t0fg;
  double     a00, a10, a20, a30, a40;
  double     barfac;
  double     a02, a12, a22, a32, a42;
  double     rmns, rmps;
  double     u0, u1, u2, u3, u4;
  double     u0p, u1p, u2p, u3p, u4p;
  double     u0pp, u1pp, u2pp, u3pp, u4pp;
  double     epotpernuc;
  double     depotpernucdx;
  double     depotpernucdi;
  double     nn_n, nn_p;
  double     kfn, kfp;
  double     efn, efp;
  double     diidnn_n, diidnn_p;
  double     un, up; // local mean field potential
  struct hnm result;

  t0fac = 3. * PI2 / 2. * satdata.rhosat0;
  t0fg  = 3. / 10. / RMN * (pow(t0fac, 2. / 3.)) * (pow(HBARC, 2.));

  // coeff
  a00    = satdata.lasat0 - t0fg * (1. + satdata.barm);
  a10    = -t0fg * (2. + 5. * satdata.barm);
  a20    = satdata.ksat0 - 2. * t0fg * (5. * satdata.barm - 1.);
  barfac = (satdata.barm + 3. * satdata.bardel);
  a02    = satdata.jsym0 - 5. / 9. * t0fg * (1. + barfac);
  a12    = satdata.lsym0 - 5. / 9. * t0fg * (2. + 5. * barfac);
  a22    = satdata.ksym0 - 10. / 9. * t0fg * (-1. + 5. * barfac);

  // x_bulk
  tmp = log(nn_) - log(satdata.rhosat0) - log(3.);
  xx  = exp(tmp) - 1. / 3.;

  // correction at nn=0
  u0  = calc_meta_model_low_density_correction(satdata, max_order, 0, xx);
  u1  = calc_meta_model_low_density_correction(satdata, max_order, 1, xx);
  u2  = calc_meta_model_low_density_correction(satdata, max_order, 2, xx);
  u0p = calc_meta_model_low_density_correction_derivative(
      satdata, max_order, 0, xx);
  u1p = calc_meta_model_low_density_correction_derivative(
      satdata, max_order, 1, xx);
  u2p = calc_meta_model_low_density_correction_derivative(
      satdata, max_order, 2, xx);
  u0pp = calc_meta_model_low_density_correction_second_derivative(
      satdata, max_order, 0, xx);
  u1pp = calc_meta_model_low_density_correction_second_derivative(
      satdata, max_order, 1, xx);
  u2pp = calc_meta_model_low_density_correction_second_derivative(
      satdata, max_order, 2, xx);

  rmns = RMN / (1. + (satdata.barm + ii_ * satdata.bardel) * (1. + 3. * xx));
  rmps = RMN / (1. + (satdata.barm - ii_ * satdata.bardel) * (1. + 3. * xx));

  if (max_order == 3) {
    a30 = satdata.qsat0 - 2. * t0fg * (4. - 5. * satdata.barm);
    a32 = satdata.qsym0 - 10. / 9. * t0fg * (4. - 5. * barfac);
    u3  = calc_meta_model_low_density_correction(satdata, max_order, 3, xx);
    u3p = calc_meta_model_low_density_correction_derivative(
        satdata, max_order, 3, xx);
    u3pp = calc_meta_model_low_density_correction_second_derivative(
        satdata, max_order, 3, xx);

    epotpernuc = a00 * u0 + a10 * xx * u1 + 0.5 * a20 * xx * xx * u2 +
                 1. / 6. * a30 * xx * xx * xx * u3 +
                 (a02 * u0 + a12 * xx * u1 + 0.5 * a22 * xx * xx * u2 +
                     1. / 6. * a32 * xx * xx * xx * u3) *
                     ii_ * ii_;
    depotpernucdx =
        a00 * u0p + a10 * u1 + a10 * xx * u1p + a20 * u2 * xx +
        0.5 * a20 * xx * xx * u2p +
        a30 / 6. * (3. * xx * xx * u3 + xx * xx * xx * u3p) +
        ii_ * ii_ * (a02 * u0p + a12 * u1 + a12 * xx * u1p + a22 * xx * u2 +
                        0.5 * a22 * xx * xx * u2p +
                        a32 / 6. * (3. * xx * xx * u3 + xx * xx * xx * u3p));
    depotpernucdi =
        2. * ii_ * (a02 * u0 + a12 * xx * u1 + 0.5 * a22 * xx * xx * u2 +
                       1. / 6. * a32 * xx * xx * xx * u3);

    result.jsym =
        5. / 9. * t0fg * pow(1. + 3. * xx, 2. / 3.) *
            (1. + (satdata.barm + 3. * satdata.bardel) * (1. + 3. * xx)) +
        a02 * u0 + a12 * xx * u1 + 0.5 * a22 * xx * xx * u2 +
        1. / 6 * a32 * xx * xx * xx * u3;

    result.lsym = 5. / 9. * t0fg * (1. + 3. * xx) *
                      (2. * pow(1. + 3. * xx, -1. / 3.) *
                              (1. + satdata.barm * (1. + 3. * xx)) +
                          3. * satdata.barm * pow(1. + 3. * xx, 2. / 3.)) +
                  (1. + 3. * xx) *
                      (a02 * u0p + a12 * (u1 + xx * u1p) +
                          a22 / 2. * (2. * xx * u2 + xx * xx * u2p) +
                          a32 / 6. * (3. * xx * xx * u3 + xx * xx * xx * u3p));

    result.ksym =
        10. / 9. * t0fg * pow(1. + 3. * xx, 2.) *
            (6. * satdata.barm * pow(1. + 3. * xx, -1. / 3.) -
                pow(1. + 3. * xx, -4. / 3.) *
                    (1. + satdata.barm * (1. + 3. * xx))) +
        pow(1. + 3. * xx, 2.) *
            (a02 * u0pp + a12 * (2. * u1p + xx * u1pp) +
                a22 / 2. * (2. * u2 + 4. * xx * u2p + xx * xx * u2pp) +
                a32 / 6. *
                    (6. * xx * u3 + 6. * xx * xx * u3p + xx * xx * xx * u3pp));
  } else if (max_order == 4) {
    a30 = satdata.qsat0 - 2. * t0fg * (4. - 5. * satdata.barm);
    a40 = satdata.zsat0 - 8. * t0fg * (-7. + 5. * satdata.barm);
    a32 = satdata.qsym0 - 10. / 9. * t0fg * (4. - 5. * barfac);
    a42 = satdata.zsym0 - 40. / 9. * t0fg * (-7. + 5. * barfac);
    u3  = calc_meta_model_low_density_correction(satdata, max_order, 3, xx);
    u4  = calc_meta_model_low_density_correction(satdata, max_order, 4, xx);
    u3p = calc_meta_model_low_density_correction_derivative(
        satdata, max_order, 3, xx);
    u4p = calc_meta_model_low_density_correction_derivative(
        satdata, max_order, 4, xx);
    u3pp = calc_meta_model_low_density_correction_second_derivative(
        satdata, max_order, 3, xx);
    u4pp = calc_meta_model_low_density_correction_second_derivative(
        satdata, max_order, 4, xx);

    epotpernuc = a00 * u0 + a10 * xx * u1 + 0.5 * a20 * xx * xx * u2 +
                 1. / 6. * a30 * xx * xx * xx * u3 +
                 1. / 24. * a40 * xx * xx * xx * xx * u4 +
                 (a02 * u0 + a12 * xx * u1 + 0.5 * a22 * xx * xx * u2 +
                     1. / 6. * a32 * xx * xx * xx * u3 +
                     1. / 24. * a42 * xx * xx * xx * xx * u4) *
                     ii_ * ii_;
    depotpernucdx =
        a00 * u0p + a10 * u1 + a10 * xx * u1p + a20 * u2 * xx +
        0.5 * a20 * xx * xx * u2p +
        a30 / 6. * (3. * xx * xx * u3 + xx * xx * xx * u3p) +
        a40 / 24. * (4. * xx * xx * xx * u4 + xx * xx * xx * xx * u4p) +
        ii_ * ii_ *
            (a02 * u0p + a12 * u1 + a12 * xx * u1p + a22 * xx * u2 +
                0.5 * a22 * xx * xx * u2p +
                a32 / 6. * (3. * xx * xx * u3 + xx * xx * xx * u3p) +
                a42 / 24. * (4. * xx * xx * xx * u4 + xx * xx * xx * xx * u4p));
    depotpernucdi =
        2. * ii_ * (a02 * u0 + a12 * xx * u1 + 0.5 * a22 * xx * xx * u2 +
                       1. / 6. * a32 * xx * xx * xx * u3 +
                       1. / 24. * a42 * xx * xx * xx * xx * u4);

    result.jsym =
        5. / 9. * t0fg * pow(1. + 3. * xx, 2. / 3.) *
            (1. + (satdata.barm + 3. * satdata.bardel) * (1. + 3. * xx)) +
        a02 * u0 + a12 * xx * u1 + 0.5 * a22 * xx * xx * u2 +
        1. / 6. * a32 * xx * xx * xx * u3 +
        1. / 24. * a42 * xx * xx * xx * xx * u4;

    result.lsym =
        5. / 9. * t0fg * (1. + 3. * xx) *
            (2. * pow(1. + 3. * xx, -1. / 3.) *
                    (1. + satdata.barm * (1. + 3. * xx)) +
                3. * satdata.barm * pow(1. + 3. * xx, 2. / 3.)) +
        (1. + 3. * xx) *
            (a02 * u0p + a12 * (u1 + xx * u1p) +
                a22 / 2. * (2. * xx * u2 + xx * xx * u2p) +
                a32 / 6. * (3. * xx * xx * u3 + xx * xx * xx * u3p) +
                a42 / 24. * (4. * xx * xx * xx * u4 + xx * xx * xx * xx * u4p));

    result.ksym =
        10. / 9. * t0fg * pow(1. + 3. * xx, 2.) *
            (6. * satdata.barm * pow(1. + 3. * xx, -1. / 3.) -
                pow(1. + 3. * xx, -4. / 3.) *
                    (1. + satdata.barm * (1. + 3. * xx))) +
        pow(1. + 3. * xx, 2.) *
            (a02 * u0pp + a12 * (2. * u1p + xx * u1pp) +
                a22 / 2. * (2. * u2 + 4. * xx * u2p + xx * xx * u2pp) +
                a32 / 6. *
                    (6. * xx * u3 + 6. * xx * xx * u3p + xx * xx * xx * u3pp) +
                a42 / 24. * (12. * xx * xx * u4 + 8. * xx * xx * xx * u4p +
                                xx * xx * xx * xx * u4pp));
  } else {
    epotpernuc =
        a00 * u0 + a10 * xx * u1 + 0.5 * a20 * xx * xx * u2 +
        (a02 * u0 + a12 * xx * u1 + 0.5 * a22 * xx * xx * u2) * ii_ * ii_;
    depotpernucdx = a00 * u0p + a10 * u1 + a10 * xx * u1p + a20 * u2 * xx +
                    0.5 * a20 * xx * xx * u2p +
                    ii_ * ii_ * (a02 * u0p + a12 * u1 + a12 * xx * u1p +
                                    a22 * xx * u2 + 0.5 * a22 * xx * xx * u2p);
    depotpernucdi =
        2. * ii_ * (a02 * u0 + a12 * xx * u1 + 0.5 * a22 * xx * xx * u2);

    result.jsym =
        5. / 9. * t0fg * pow(1. + 3. * xx, 2. / 3.) *
            (1. + (satdata.barm + 3. * satdata.bardel) * (1. + 3. * xx)) +
        a02 * u0 + a12 * xx * u1 + 0.5 * a22 * xx * xx * u2;

    result.lsym =
        5. / 9. * t0fg * (1. + 3. * xx) *
            (2. * pow(1. + 3. * xx, -1. / 3.) *
                    (1. + satdata.barm * (1. + 3. * xx)) +
                3. * satdata.barm * pow(1. + 3. * xx, 2. / 3.)) +
        (1. + 3. * xx) * (a02 * u0p + a12 * (u1 + xx * u1p) +
                             a22 / 2. * (2. * xx * u2 + xx * xx * u2p));

    result.ksym =
        10. / 9. * t0fg * pow(1. + 3. * xx, 2.) *
            (6. * satdata.barm * pow(1. + 3. * xx, -1. / 3.) -
                pow(1. + 3. * xx, -4. / 3.) *
                    (1. + satdata.barm * (1. + 3. * xx))) +
        pow(1. + 3. * xx, 2.) *
            (a02 * u0pp + a12 * (2. * u1p + xx * u1pp) +
                a22 / 2. * (2. * u2 + 4. * xx * u2p + xx * xx * u2pp));
  }

  if (nn_ == 0.) // to avoid singularities
  {
    result.enpernuc = 0.;
    result.spernuc  = 0.;
    result.mun      = 0.;
    result.mup      = 0.;
    result.fpernuc  = 0.;
    result.p        = 0.;

    return result;
  }

  nn_n = nn_ * (1. + ii_) / 2.;
  nn_p = nn_ * (1. - ii_) / 2.;

  kfn = cpow(3. * PI2 * nn_n, 1. / 3.);
  kfp = cpow(3. * PI2 * nn_p, 1. / 3.);

  efn = HBARC * HBARC * kfn * kfn / 2. / rmns;
  efp = HBARC * HBARC * kfp * kfp / 2. / rmps;

  diidnn_n = (1. - ii_) / nn_;
  diidnn_p = -(1. + ii_) / nn_;

  if (ii_ == 1.) // Pure Neutron Matter (PNM)
  {
    result.enpernuc = 1. / nn_ * 3. / 5. * nn_n * efn *
                          (1. + 5. * PI2 / 12. * cpow(tt_ / efn, 2.)) +
                      epotpernuc;

    result.spernuc = 1. / nn_ * PI2 / 2. * tt_ * nn_n / efn;

    un = epotpernuc +
         nn_ * (1. / 3. / satdata.rhosat0 * depotpernucdx +
                   diidnn_n * depotpernucdi) +
         3. / 5. * nn_n * (cpow(HBARC * kfn, 2.) / 2. / RMN -
                              5. * PI2 / 12. * tt_ * tt_ / efn * rmns / RMN) *
             ((satdata.barm + ii_ * satdata.bardel) / satdata.rhosat0 +
                 satdata.bardel * (1. + 3. * xx) * diidnn_n);

    result.mun = efn * (1. - PI2 / 12. * cpow(tt_ / efn, 2.)) + un;

    result.mup = 0.;

    result.fpernuc = result.enpernuc - tt_ * result.spernuc;

    result.p = nn_n * result.mun - nn_ * result.fpernuc;

    return result;
  } else // use of cpow for ii_ > 1 or ii_ < -1...
  {
    result.enpernuc =
        1. / nn_ * 3. / 5. *
            (nn_n * efn * (1. + 5. * PI2 / 12. * cpow(tt_ / efn, 2.)) +
                nn_p * efp * (1. + 5. * PI2 / 12. * cpow(tt_ / efp, 2.))) +
        epotpernuc;

    result.spernuc = 1. / nn_ * PI2 / 2. * tt_ * (nn_n / efn + nn_p / efp);

    un = epotpernuc +
         nn_ * (1. / 3. / satdata.rhosat0 * depotpernucdx +
                   diidnn_n * depotpernucdi) +
         3. / 5. * nn_n * efn * (1. - 5 * PI2 / 12. * cpow(tt_ / efn, 2.)) *
             rmns / RMN *
             ((satdata.barm + ii_ * satdata.bardel) / satdata.rhosat0 +
                 satdata.bardel * (1. + 3. * xx) * diidnn_n) +
         3. / 5. * nn_p * efp * (1. - 5 * PI2 / 12. * cpow(tt_ / efp, 2.)) *
             rmps / RMN *
             ((satdata.barm - ii_ * satdata.bardel) / satdata.rhosat0 -
                 satdata.bardel * (1. + 3. * xx) * diidnn_n);
    up = epotpernuc +
         nn_ * (1. / 3. / satdata.rhosat0 * depotpernucdx +
                   diidnn_p * depotpernucdi) +
         3. / 5. * nn_n * efn * (1. - 5 * PI2 / 12. * cpow(tt_ / efn, 2.)) *
             rmns / RMN *
             ((satdata.barm + ii_ * satdata.bardel) / satdata.rhosat0 +
                 satdata.bardel * (1. + 3. * xx) * diidnn_p) +
         3. / 5. * nn_p * efp * (1. - 5 * PI2 / 12. * cpow(tt_ / efp, 2.)) *
             rmps / RMN *
             ((satdata.barm - ii_ * satdata.bardel) / satdata.rhosat0 -
                 satdata.bardel * (1. + 3. * xx) * diidnn_p);

    result.mun = efn * (1. - PI2 / 12. * cpow(tt_ / efn, 2.)) + un;

    result.mup = efp * (1. - PI2 / 12. * cpow(tt_ / efp, 2.)) + up;

    result.fpernuc = result.enpernuc - tt_ * result.spernuc;

    result.p = nn_n * result.mun + nn_p * result.mup - nn_ * result.fpernuc;

    return result;
  }
}

double calc_squared_nucleon_sound_velocity(struct parameters satdata,
    int max_order, double nn_, double ii_, double tt_) {
  double dnn = nn_ / 1000.;
  double en_p =
      calc_meta_model_nuclear_matter(satdata, max_order, nn_ + dnn, ii_, tt_)
          .enpernuc;
  struct hnm t =
      calc_meta_model_nuclear_matter(satdata, max_order, nn_, ii_, tt_);
  double en_m =
      calc_meta_model_nuclear_matter(satdata, max_order, nn_ - dnn, ii_, tt_)
          .enpernuc;
  double d2enpernucdnn2 = (en_p - 2. * t.enpernuc + en_m) / pow(dnn, 2.);
  double kis            = 9. * nn_ * nn_ * d2enpernucdnn2 + 18. * t.p / nn_;
  return kis / 9. / (RMN + t.enpernuc + t.p / nn_);
}

double calc_asymmetry_factor(double m_, double ii_) {
  return 0.5 * (cpow(1. + ii_, m_) + cpow(1. - ii_, m_));
}

struct hnm calc_skyrme_nuclear_matter(
    struct skyrme_parameters coeff, double nn_, double ii_) {
  struct hnm result;

  double f23, f53, f83;
  double denpernucdi;
  double kf;
  double endens;

  f23 = calc_asymmetry_factor(2. / 3., ii_);
  f53 = calc_asymmetry_factor(5. / 3., ii_);
  f83 = calc_asymmetry_factor(8. / 3., ii_);

  kf = pow(3. * PI2 * nn_ / 2., 1. / 3.);
  endens =
      3. * HBARC * HBARC / 20. *
          (cpow(1. + ii_, 5. / 3.) / RMN + cpow(1 - ii_, 5. / 3.) / RMP) * nn_ *
          kf * kf +
      1. / 8. * coeff.t0 * (3. - (1. + 2. * coeff.x0) * ii_ * ii_) * nn_ * nn_ +
      3. / 40. * coeff.t1 * ((2. + coeff.x1) * f53 - (0.5 + coeff.x1) * f83) *
          nn_ * nn_ * kf * kf +
      3. / 40. * coeff.t2 * ((2. + coeff.x2) * f53 + (0.5 + coeff.x2) * f83) *
          nn_ * nn_ * kf * kf +
      1. / 48. * coeff.t3 * (3. - (1. + 2. * coeff.x3) * ii_ * ii_) *
          pow(nn_, coeff.alpha + 2.) +
      3. / 40. * coeff.t4 * ((2. + coeff.x4) * f53 - (0.5 + coeff.x4) * f83) *
          pow(nn_, coeff.beta + 2.) * kf * kf +
      3. / 40. * coeff.t5 * ((2. + coeff.x5) * f53 + (0.5 + coeff.x5) * f83) *
          pow(nn_, coeff.gamma + 2.) * kf * kf;

  result.enpernuc = endens / nn_;

  result.spernuc = 0.;

  result.fpernuc = result.enpernuc; // T=0

  result.p =
      HBARC * HBARC / 10. *
          (cpow(1. + ii_, 5. / 3.) / RMN + cpow(1. - ii_, 5. / 3.) / RMP) *
          nn_ * kf * kf +
      1. / 8. * coeff.t0 * (3. - (1. + 2. * coeff.x0) * ii_ * ii_) * nn_ * nn_ +
      1. / 8. * coeff.t1 * ((2. + coeff.x1) * f53 - (0.5 + coeff.x1) * f83) *
          nn_ * nn_ * kf * kf +
      1. / 8. * coeff.t2 * ((2. + coeff.x2) * f53 + (0.5 + coeff.x2) * f83) *
          nn_ * nn_ * kf * kf +
      (coeff.alpha + 1.) / 48. * coeff.t3 *
          (3. - (1. + 2. * coeff.x3) * ii_ * ii_) * pow(nn_, coeff.alpha + 2.) +
      (3. * coeff.beta + 5.) / 40. * coeff.t4 *
          ((2. + coeff.x4) * f53 - (0.5 + coeff.x4) * f83) *
          pow(nn_, coeff.beta + 2.) * kf * kf +
      (3. * coeff.gamma + 5.) / 40. * coeff.t5 *
          ((2. + coeff.x5) * f53 + (0.5 + coeff.x5) * f83) *
          pow(nn_, coeff.gamma + 2.) * kf * kf;

  denpernucdi =
      HBARC * HBARC / 4.0 *
          (cpow(1. + ii_, 2. / 3.) / RMN - cpow(1. - ii_, 2. / 3.) / RMP) * kf *
          kf -
      1. / 4. * coeff.t0 * (1. + 2. * coeff.x0) * ii_ * nn_ +
      1. / 40. * coeff.t1 *
          (-(10. + 5 * coeff.x1) * f23 + (14. + 13. * coeff.x1) * f53 -
              (4. + 8. * coeff.x1) * f83) *
          nn_ * kf * kf / ii_ +
      1. / 40. * coeff.t2 *
          (-(10. + 5 * coeff.x2) * f23 + (6. - 3. * coeff.x2) * f53 +
              (4. + 8. * coeff.x2) * f83) *
          nn_ * kf * kf / ii_ -
      1. / 24. * coeff.t3 * (1. + 2. * coeff.x3) * ii_ *
          pow(nn_, coeff.alpha + 1.) +
      1. / 40. * coeff.t4 *
          (-(10. + 5 * coeff.x4) * f23 + (14. + 13. * coeff.x4) * f53 -
              (4. + 8. * coeff.x4) * f83) *
          pow(nn_, coeff.beta + 1.) * kf * kf / ii_ +
      1. / 40. * coeff.t5 *
          (-(10. + 5 * coeff.x5) * f23 + (6. - 3. * coeff.x5) * f53 +
              (4. + 8. * coeff.x5) * f83) *
          pow(nn_, coeff.gamma + 1.) * kf * kf / ii_;

  result.mun = result.enpernuc + 1. / nn_ * result.p + (1. - ii_) * denpernucdi;

  result.mup = result.enpernuc + 1. / nn_ * result.p - (1. + ii_) * denpernucdi;

  result.jsym = 0.;

  result.lsym = 0.;

  result.ksym = 0.;

  return result;
}
