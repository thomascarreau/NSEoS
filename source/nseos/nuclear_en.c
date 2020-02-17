#include <math.h>
#include <string.h>

#include <gsl/gsl_roots.h>

#include "nuclear_en.h"

double calc_ldm_meta_model_nuclear_en(struct parameters satdata, int max_order,
    double aa_, double ii_, double n0_) {
  struct hnm meta;
  double     bulk_en;
  double     surf_en;
  double     coul_en;

  meta    = calc_meta_model_nuclear_matter(satdata, max_order, n0_, ii_, 0.);
  bulk_en = meta.enpernuc * aa_;
  surf_en = calc_ldm_surface_en(satdata, aa_);
  coul_en = calc_coulomb_en(aa_, ii_, n0_);

  return bulk_en + surf_en + coul_en;
}

double calc_etf_ana_meta_model_nuclear_en(struct parameters satdata,
    int max_order, double aa_, double ii_, double n0_) {
  struct hnm meta;
  double     bulk_en;
  double     surf_en;
  double     coul_en;

  meta    = calc_meta_model_nuclear_matter(satdata, max_order, n0_, ii_, 0.);
  bulk_en = meta.enpernuc * aa_;
  surf_en = calc_etf_ana_surface_en(satdata, aa_, ii_, n0_);
  coul_en = calc_coulomb_en(aa_, ii_, n0_);

  return bulk_en + surf_en + coul_en;
}

double calc_ls_meta_model_nuclear_free_en(struct parameters satdata,
    struct sf_params sparams, int max_order, double aa_, double ii_, double n0_,
    double tt_) {
  struct hnm meta;
  double     bulk_free_en;
  double     surf_free_en;
  double     coul_en;

  meta = calc_meta_model_nuclear_matter(satdata, max_order, n0_, ii_, tt_);
  bulk_free_en = meta.fpernuc * aa_;
  surf_free_en = calc_ls_surface_free_en(satdata, sparams, aa_, ii_, n0_, tt_);
  coul_en      = calc_coulomb_en(aa_, ii_, n0_);

  return bulk_free_en + surf_free_en + coul_en;
}

double get_shell_en_from_myers_table(int aa_, int zz_) {
  int   aa, zz;
  float shell_en;
  int   checker;
  FILE *myers_table = fopen("../../input/shell_corr/myers.table", "r");

  aa      = 0;
  zz      = 0;
  checker = 0;

  while (fscanf(myers_table, "%d %d %f", &zz, &aa, &shell_en) == 3)
    if (aa == aa_ && zz == zz_) {
      checker = 1;
      break;
    }

  if (checker == 0) {
    fprintf(stderr, "ERROR: out of Myers table! ; ESH = 0\n");
    shell_en = 0.;
  }

  fclose(myers_table);

  return shell_en;
}

double calc_pairing_en(int aa_, int zz_) {
  double ap = 12.;
  int    nn = aa_ - zz_;

  if (nn % 2 == 1 && zz_ % 2 == 1)
    return ap / pow(aa_, 1. / 2.);
  else if (aa_ % 2 == 1)
    return 0.;
  else if (nn % 2 == 0 && zz_ % 2 == 0)
    return -ap / pow(aa_, 1. / 2.);
  else {
    fprintf(stderr, "ERROR: check nuclear pairing function!");
    return 0.;
  }
}

double calc_ls_meta_model_nuclear_en_micro(struct parameters satdata,
    struct sf_params sparams, int max_order, double aa_, double ii_,
    double n0_) {
  struct hnm meta;
  double     bulk_en;
  double     surf_en;
  double     coul_en;
  double     shell_en;
  double     pairing_en;

  meta    = calc_meta_model_nuclear_matter(satdata, max_order, n0_, ii_, 0.);
  bulk_en = meta.enpernuc * aa_;
  surf_en = calc_ls_surface_free_en(satdata, sparams, aa_, ii_, n0_, 0.);
  coul_en = calc_coulomb_en(aa_, ii_, n0_);
  shell_en =
      get_shell_en_from_myers_table((int)aa_, (int)aa_ * (1. - ii_) / 2.);
  pairing_en = calc_pairing_en((int)aa_, (int)aa_ * (1. - ii_) / 2.);

  return bulk_en + surf_en + coul_en + shell_en + pairing_en;
}

double calc_electron_binding_energy(int zz_) {
  return 1.44381e-5 * pow(zz_, 2.39) + 1.55468e-12 * pow(zz_, 5.35);
}

double calc_nuclear_mass_from_mass_excess(int aa_, int zz_, double deps_) {
  return deps_ + aa_ * AMU - zz_ * MEL + calc_electron_binding_energy(zz_);
}

double function_rn(double x, void *params) {
  struct pars *p = (struct pars *)params;

  double a = p->a;
  double b = p->b;
  double c = p->c;

  return a * x * x * x * x + b * x + c;
}

double evaluate_rn(double nb_, double ii_, double n0_, double ng_, int d_,
    struct sf_params sparams, char phase[]) {
  int                          status;
  int                          iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *           s;
  double                       r    = 0;
  double                       x_lo = 0.1, x_hi = 1000.0;
  gsl_function                 F;

  double u;
  double ypnuc;
  double sigmas;
  double sigmac;

  if (strcmp(phase, "nuclei") == 0)
    u = (nb_ - ng_) / (n0_ - ng_);
  else if (strcmp(phase, "bubbles") == 0)
    u = (n0_ - nb_) / (n0_ - ng_);
  else {
    fprintf(stderr, "ERROR: phase must be either 'nuclei' or 'bubbles'!\n");
    exit(1);
  }

  ypnuc  = (1. - ii_) / 2.;
  sigmas = sparams.sigma0 * (pow(2., sparams.p + 1.) + sparams.b) /
           (pow(ypnuc, -sparams.p) + sparams.b + pow(1. - ypnuc, -sparams.p));
  sigmac = sigmas * sparams.sigma0c / sparams.sigma0 * sparams.alpha *
           (sparams.beta - ypnuc);

  double a = 4. * PI * ALPHAFS * HBARC * pow(ypnuc * n0_, 2.) * calc_fd(u, d_);
  double b = -d_ * sigmas;
  double c = -2 * d_ * (d_ - 1.) * sigmac;

  struct pars params = {a, b, c};

  F.function = &function_rn;
  F.params   = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);

  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    r      = gsl_root_fsolver_root(s);
    x_lo   = gsl_root_fsolver_x_lower(s);
    x_hi   = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);
  } while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free(s);

  return r;
}

double calc_surface_plus_coulomb_energy_density(double nb_, double ii_,
    double n0_, double ng_, int d_, struct sf_params sparams, char phase[]) {
  double ypnuc;
  double sigmas;
  double sigmac;
  double u;
  double rn;
  double surface_endens;
  double curvature_endens;
  double coulomb_endens;

  ypnuc  = (1. - ii_) / 2.;
  sigmas = sparams.sigma0 * (pow(2., sparams.p + 1.) + sparams.b) /
           (pow(ypnuc, -sparams.p) + sparams.b + pow(1. - ypnuc, -sparams.p));
  sigmac = sigmas * sparams.sigma0c / sparams.sigma0 * sparams.alpha *
           (sparams.beta - ypnuc);

  if (strcmp(phase, "nuclei") == 0)
    u = (nb_ - ng_) / (n0_ - ng_);
  else if (strcmp(phase, "bubbles") == 0)
    u = (n0_ - nb_) / (n0_ - ng_);
  else {
    fprintf(stderr, "ERROR: phase must be either 'nuclei' or 'bubbles'!\n");
    exit(1);
  }

  rn = evaluate_rn(nb_, ii_, n0_, ng_, d_, sparams, phase);

  surface_endens   = u * d_ * sigmas / rn;
  curvature_endens = u * d_ * (d_ - 1.) * sigmac / rn / rn;

  coulomb_endens = u * 2. * PI * ALPHAFS * HBARC * pow(ypnuc * n0_ * rn, 2.) *
                   calc_fd(u, d_);

  return surface_endens + curvature_endens + coulomb_endens;
}
