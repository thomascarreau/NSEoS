#include <gsl/gsl_multiroots.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "crust.h"
#include "lepton.h"
#include "modeling.h"
#include "nuclear_en.h"

double get_shell_energy_per_nucleon(float nb_, float zz_) {
  FILE *bsk_table = fopen("../../input/shell_corr/bsk24_nb_zz_esh.data", "r");

  float esh_nb;

  if (nb_ < 2.7e-4) {
    fprintf(stderr, "WARNING: out of table! nb = %g \t zz = %g ; eSH = 0\n",
        nb_, zz_);
    esh_nb = 0.;
    return esh_nb;
  }

  float nb, nb_sav, nb_m, nb_p;
  float zz;
  float esh, esh_sav, esh_m;
  int   checker;

  checker = 0;

  while (fscanf(bsk_table, "%f %f %f", &nb, &zz, &esh) == 3) {
    if (nb > nb_ && checker == 0) {
      nb_m    = nb_sav;
      nb_p    = nb;
      esh_m   = esh_sav;
      checker = 1;
    }

    if (zz == zz_ && nb == nb_p && checker == 1) {
      esh_nb = (esh - esh_m) / (nb - nb_m) * (nb_ - nb) + esh;
      break;
    }

    if (checker == 1 && nb != nb_p) {
      fprintf(stderr, "WARNING: out of table! nb = %g \t zz = %g ; eSH = 0\n",
          nb_, zz_);
      esh_nb = 0.;
      break;
    }

    if (zz == zz_) {
      nb_sav  = nb;
      esh_sav = esh;
    }
  }

  if (checker == 0) {
    fprintf(stderr, "WARNING: out of table! ; eSH = 0\n");
    esh_nb = 0.;
  }

  fclose(bsk_table);

  return esh_nb;
}

double calc_ion_free_en(struct parameters satdata, struct sf_params sparams,
    double aa_, double del_, double n0_, double np_, double ng_, double tt_,
    char phase[]) {
  double zz = aa_ * (1. - del_) / 2.;
  double mi =
      CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, aa_, del_, n0_, tt_) +
      zz * RMP + (aa_ * (1. - ng_ / n0_) - zz) * RMN -
      aa_ / n0_ * ng_ *
          calc_meta_model_nuclear_matter(
              satdata, TAYLOR_EXP_ORDER, ng_, 1., tt_)
              .fpernuc;

  if (strcmp(phase, "sol") == 0) {
    if (tt_ == 0.) {
      return CALC_NUCLEAR_EN(
                 satdata, sparams, TAYLOR_EXP_ORDER, aa_, del_, n0_, 0.) +
             calc_lattice_en(aa_, del_, n0_, np_) +
             calc_finite_size_contrib(aa_, del_, n0_, np_);
    } else {
      return CALC_NUCLEAR_EN(
                 satdata, sparams, TAYLOR_EXP_ORDER, aa_, del_, n0_, tt_) +
             calc_lattice_en_for_tm(zz, np_) +
             calc_finite_size_contrib(aa_, del_, n0_, np_) +
             calc_zp_en(zz, np_, mi) + calc_harmonic_contrib(zz, np_, mi, tt_) +
             calc_anharmonic_contrib(zz, np_, tt_);
    }
  } else if (strcmp(phase, "liq") == 0) {
    return CALC_NUCLEAR_EN(
               satdata, sparams, TAYLOR_EXP_ORDER, aa_, del_, n0_, tt_) +
           calc_finite_size_contrib(aa_, del_, n0_, np_) +
           calc_translational_free_en(zz, np_, mi, tt_) +
           calc_total_coulomb_contrib(zz, np_, tt_);
  } else {
    fprintf(stderr, "ERROR: phase must be either 'sol' or 'liq'!\n");
    exit(EXIT_FAILURE);
  }
}

struct crust_fun_4d calc_crust_fun_4d(struct parameters satdata,
    struct sf_params sparams, double aa_, double del_, double n0_, double np_,
    double ng_, double tt_, char phase[]) {
  struct crust_fun_4d result;
  double              fion;
  double              epsa;
  double              epsb;
  double              epsr;
  double              epsp;
  double              epsg;
  double              fion_ap, fion_am;
  double              fion_bp, fion_bm;
  double              fion_rp, fion_rm;
  double              fion_pp, fion_pm;
  double              fion_gp, fion_gm;
  double              dfiondaa;
  double              dfionddel;
  double              dfiondn0;
  double              dfiondnp;
  double              dfiondng;
  double              muel;
  struct hnm          ngas;
  double              mu;

  epsa = 0.001;
  epsb = 0.0001;
  epsr = 0.0001;
  epsp = np_ / 1000.;
  epsg = ng_ / 1000.;

  fion =
      calc_ion_free_en(satdata, sparams, aa_, del_, n0_, np_, ng_, tt_, phase);
  fion_ap = calc_ion_free_en(
      satdata, sparams, aa_ + epsa, del_, n0_, np_, ng_, tt_, phase);
  fion_am = calc_ion_free_en(
      satdata, sparams, aa_ - epsa, del_, n0_, np_, ng_, tt_, phase);
  fion_bp = calc_ion_free_en(
      satdata, sparams, aa_, del_ + epsb, n0_, np_, ng_, tt_, phase);
  fion_bm = calc_ion_free_en(
      satdata, sparams, aa_, del_ - epsb, n0_, np_, ng_, tt_, phase);
  fion_rp = calc_ion_free_en(
      satdata, sparams, aa_, del_, n0_ + epsr, np_, ng_, tt_, phase);
  fion_rm = calc_ion_free_en(
      satdata, sparams, aa_, del_, n0_ - epsr, np_, ng_, tt_, phase);
  fion_pp = calc_ion_free_en(
      satdata, sparams, aa_, del_, n0_, np_ + epsp, ng_, tt_, phase);
  fion_pm = calc_ion_free_en(
      satdata, sparams, aa_, del_, n0_, np_ - epsp, ng_, tt_, phase);
  fion_gp = calc_ion_free_en(
      satdata, sparams, aa_, del_, n0_, np_, ng_ + epsg, tt_, phase);
  fion_gm = calc_ion_free_en(
      satdata, sparams, aa_, del_, n0_, np_, ng_ - epsg, tt_, phase);

  dfiondaa  = (fion_ap - fion_am) / 2. / epsa; // 2 points derivatives
  dfionddel = (fion_bp - fion_bm) / 2. / epsb;
  dfiondn0  = (fion_rp - fion_rm) / 2. / epsr;
  dfiondnp  = (fion_pp - fion_pm) / 2. / epsp;
  dfiondng  = (fion_gp - fion_gm) / 2. / epsg;

  muel = calc_egas_chemical_potential(np_, tt_);

  ngas =
      calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, ng_, 1., tt_);

  if (ng_ > 0.) {
    mu = RMN + ngas.mun +
         2. * np_ / aa_ / (1. - del_) / (1. - 2. * np_ / n0_ / (1. - del_)) *
             dfiondng;
  } else {
    mu = RMN;
  }

  result.f_stability = dfiondaa - fion / aa_;
  result.f_beta =
      (dfionddel - np_ / (1. - del_) * dfiondnp) * 2. / aa_ - muel - RMP + RMN;
  result.f_muneq = fion / aa_ + (1. - del_) / aa_ * dfionddel - (mu - RMN) +
                   ng_ / n0_ * (mu - RMN - ngas.fpernuc);
  result.f_presseq =
      n0_ * n0_ * dfiondn0 / aa_ - ng_ * (mu - RMN) + ng_ * ngas.fpernuc;

  return result;
}

struct crust_fun_4d calc_crust_fun_zz_fixed(struct parameters satdata,
    struct sf_params sparams, double aa_, double del_, double n0_, double np_,
    double ng_, double tt_, char phase[]) {
  struct crust_fun_4d result;
  double              fion;
  double              epsa;
  double              epsb;
  double              epsr;
  double              epsp;
  double              epsg;
  double              fion_ap, fion_am;
  double              fion_bp, fion_bm;
  double              fion_rp, fion_rm;
  double              fion_pp, fion_pm;
  double              fion_gp, fion_gm;
  double              dfiondaa;
  double              dfionddel;
  double              dfiondn0;
  double              dfiondnp;
  double              dfiondng;
  double              muel;
  struct hnm          ngas;
  double              mu;

  epsa = 0.001;
  epsb = 0.0001;
  epsr = 0.0001;
  epsp = np_ / 1000.;
  epsg = ng_ / 1000.;

  fion =
      calc_ion_free_en(satdata, sparams, aa_, del_, n0_, np_, ng_, tt_, phase);
  fion_ap = calc_ion_free_en(
      satdata, sparams, aa_ + epsa, del_, n0_, np_, ng_, tt_, phase);
  fion_am = calc_ion_free_en(
      satdata, sparams, aa_ - epsa, del_, n0_, np_, ng_, tt_, phase);
  fion_bp = calc_ion_free_en(
      satdata, sparams, aa_, del_ + epsb, n0_, np_, ng_, tt_, phase);
  fion_bm = calc_ion_free_en(
      satdata, sparams, aa_, del_ - epsb, n0_, np_, ng_, tt_, phase);
  fion_rp = calc_ion_free_en(
      satdata, sparams, aa_, del_, n0_ + epsr, np_, ng_, tt_, phase);
  fion_rm = calc_ion_free_en(
      satdata, sparams, aa_, del_, n0_ - epsr, np_, ng_, tt_, phase);
  fion_pp = calc_ion_free_en(
      satdata, sparams, aa_, del_, n0_, np_ + epsp, ng_, tt_, phase);
  fion_pm = calc_ion_free_en(
      satdata, sparams, aa_, del_, n0_, np_ - epsp, ng_, tt_, phase);
  fion_gp = calc_ion_free_en(
      satdata, sparams, aa_, del_, n0_, np_, ng_ + epsg, tt_, phase);
  fion_gm = calc_ion_free_en(
      satdata, sparams, aa_, del_, n0_, np_, ng_ - epsg, tt_, phase);

  dfiondaa  = (fion_ap - fion_am) / 2. / epsa; // 2 points derivatives
  dfionddel = (fion_bp - fion_bm) / 2. / epsb;
  dfiondn0  = (fion_rp - fion_rm) / 2. / epsr;
  dfiondnp  = (fion_pp - fion_pm) / 2. / epsp;
  dfiondng  = (fion_gp - fion_gm) / 2. / epsg;

  muel = calc_egas_chemical_potential(np_, tt_);

  ngas =
      calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, ng_, 1., tt_);

  if (ng_ > 0.) {
    mu = RMN + ngas.mun +
         2. * np_ / aa_ / (1. - del_) / (1. - 2. * np_ / n0_ / (1. - del_)) *
             dfiondng;
  } else {
    mu = RMN;
  }

  result.f_beta = 0.;

  result.f_stability = dfiondaa - fion / aa_ -
                       (1. - del_) / 2. *
                           (muel + 2 * np_ / aa_ / (1. - del_) * dfiondnp -
                               2. / aa_ * dfionddel + RMP - RMN);
  result.f_muneq = dfiondaa + (1. - del_) / aa_ * dfionddel - (mu - RMN) +
                   ng_ / n0_ * (mu - RMN - ngas.fpernuc);
  result.f_presseq =
      n0_ * n0_ * dfiondn0 / aa_ - ng_ * (mu - RMN) + ng_ * ngas.fpernuc;

  return result;
}

int assign_ocrust_fun_3d(const gsl_vector *x, void *params, gsl_vector *f) {
  double np = ((struct rparams_crust *)params)->np;
  double tt = ((struct rparams_crust *)params)->tt;
  char   phase[3];
  strcpy(phase, ((struct rparams_crust *)params)->phase);
  struct parameters satdata = ((struct rparams_crust *)params)->satdata;
  struct sf_params  sparams = ((struct rparams_crust *)params)->sparams;

  const double x0 = gsl_vector_get(x, 0);
  const double x1 = gsl_vector_get(x, 1);
  const double x2 = gsl_vector_get(x, 2);

  np = np * (1. - x1) / 2.;

  struct crust_fun_4d functs;
  functs = calc_crust_fun_4d(satdata, sparams, x0, x1, x2, np, 0., tt, phase);

  const double y0 = functs.f_stability;
  const double y1 = functs.f_beta;
  const double y2 = functs.f_presseq;

  gsl_vector_set(f, 0, y0);
  gsl_vector_set(f, 1, y1);
  gsl_vector_set(f, 2, y2);

  return GSL_SUCCESS;
}

struct compo calc_ocrust3d_composition(double nb_, double tt_, char phase[],
    double *guess, struct parameters satdata, struct sf_params sparams) {
  struct compo eq;

  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *           s;

  int    status;
  size_t iter = 0;

  double aa_old, aa_new, astep;
  double basym_old, basym_new, bstep;
  double n0_old, n0_new, rstep;

  struct rparams_crust p;
  p.np = nb_; // modification in assign_ocrust_fun_3d (DIRTY!)
  p.tt = tt_;
  strcpy(p.phase, phase);
  p.satdata = satdata;
  p.sparams = sparams;

  const size_t n = 3;
  gsl_vector * x = gsl_vector_alloc(n);
  gsl_vector_set(x, 0, guess[0]);
  gsl_vector_set(x, 1, guess[1]);
  gsl_vector_set(x, 2, guess[2]);

  gsl_multiroot_function f = {&assign_ocrust_fun_3d, n, &p};

  T = gsl_multiroot_fsolver_broyden;
  s = gsl_multiroot_fsolver_alloc(T, 3);
  gsl_multiroot_fsolver_set(s, &f, x);

  do {
    aa_old    = gsl_vector_get(s->x, 0);
    basym_old = gsl_vector_get(s->x, 1);
    n0_old    = gsl_vector_get(s->x, 2);

    iter++;

    status = gsl_multiroot_fsolver_iterate(s);

    /* if (gsl_vector_get (s->f, 0) != gsl_vector_get (s->f, 0) */
    /*         && gsl_vector_get (s->f, 1) != gsl_vector_get (s->f, 1) */
    /*         && gsl_vector_get (s->f, 2) != gsl_vector_get (s->f, 2)) */
    /*     iter = 1000; // to avoid 'matrix is singular' error */

    aa_new    = gsl_vector_get(s->x, 0);
    basym_new = gsl_vector_get(s->x, 1);
    n0_new    = gsl_vector_get(s->x, 2);
    astep     = aa_new - aa_old;
    bstep     = basym_new - basym_old;
    rstep     = n0_new - n0_old;

    if (status)
      break;

    // dirty backstepping
    while (gsl_vector_get(s->x, 0) < 0.) {
      astep  = astep / 4.;
      aa_new = aa_old + astep;
      gsl_vector_set(x, 0, aa_new);
      gsl_vector_set(x, 1, basym_new);
      gsl_vector_set(x, 2, n0_new);
      gsl_multiroot_fsolver_set(s, &f, x);
      aa_new = gsl_vector_get(s->x, 0.);
    }
    while (gsl_vector_get(s->x, 1) < 0.01 || gsl_vector_get(s->x, 1) > 1.0) {
      bstep     = bstep / 4.;
      basym_new = basym_old + bstep;
      gsl_vector_set(x, 0, aa_new);
      gsl_vector_set(x, 1, basym_new);
      gsl_vector_set(x, 2, n0_new);
      gsl_multiroot_fsolver_set(s, &f, x);
      basym_new = gsl_vector_get(s->x, 1);
    }
    while (gsl_vector_get(s->x, 2) < 0.) {
      rstep  = rstep / 4.;
      n0_new = n0_old + rstep;
      gsl_vector_set(x, 0, aa_new);
      gsl_vector_set(x, 1, basym_new);
      gsl_vector_set(x, 2, n0_new);
      gsl_multiroot_fsolver_set(s, &f, x);
      n0_new = gsl_vector_get(s->x, 2);
    }

    status = gsl_multiroot_test_residual(s->f, 9e-9);
  }

  while (status == GSL_CONTINUE && iter < 1000);

  if (iter == 1000) {
    eq.aa  = NAN;
    eq.del = NAN;
    eq.n0  = NAN;
    eq.ng  = 0.0;
  } else {
    eq.aa  = gsl_vector_get(s->x, 0);
    eq.del = gsl_vector_get(s->x, 1);
    eq.n0  = gsl_vector_get(s->x, 2);
    eq.ng  = 0.0;
  }

  guess[0] = eq.aa;
  guess[1] = eq.del;
  guess[2] = eq.n0;

  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);

  return eq;
}

struct compo calc_ocrust_composition_with_mass_table(
    double nb_, double tt_, char phase[], char mass_table[]) {
  char path_of_mass_table[128] = "../../input/mass_tables/";
  strcat(path_of_mass_table, mass_table);

  FILE *table = fopen(path_of_mass_table, "r");

  struct compo eq;

  eq.n0 = -1.;
  eq.ng = 0.;

  int   nn, zz;
  float deps;

  int    aa;
  double ii;
  double np;
  double vws;

  double mi;
  double fdenswsmin = 1.e99;
  double fdensws;

  while (fscanf(table, "%d %d %f", &zz, &nn, &deps) == 3) {
    aa  = nn + zz;
    ii  = 1. - 2. * (double)zz / (double)aa;
    vws = aa / nb_;
    np  = zz / vws;

    mi = calc_nuclear_mass_from_mass_excess(aa, zz, deps);

    if (strcmp(phase, "sol") == 0 || tt_ == 0.) {
      fdensws = mi / vws + calc_lattice_en_for_tm((double)zz, np) / vws +
                calc_zp_en((double)zz, np, mi) / vws +
                calc_harmonic_contrib((double)zz, np, mi, tt_) / vws +
                calc_egas_free_energy_density(np, tt_);
    } else if (strcmp(phase, "liq") == 0) {
      fdensws = mi / vws +
                calc_translational_free_en((double)zz, np, mi, tt_) / vws +
                calc_total_coulomb_contrib((double)zz, np, tt_) / vws +
                calc_egas_free_energy_density(np, tt_);
    } else {
      fprintf(stderr, "ERROR: phase must be either 'sol' or 'liq'!\n");
      exit(EXIT_FAILURE);
    }

    if (fdensws < fdenswsmin) {
      eq.aa      = aa;
      eq.del     = ii;
      fdenswsmin = fdensws;
    }
  }

  fclose(table);

  return eq;
}

int assign_icrust_fun_4d(const gsl_vector *x, void *params, gsl_vector *f) {
  double np = ((struct rparams_crust *)params)->np;
  double tt = ((struct rparams_crust *)params)->tt;
  char   phase[3];
  strcpy(phase, ((struct rparams_crust *)params)->phase);
  struct parameters satdata = ((struct rparams_crust *)params)->satdata;
  struct sf_params  sparams = ((struct rparams_crust *)params)->sparams;

  const double x0 = gsl_vector_get(x, 0);
  const double x1 = gsl_vector_get(x, 1);
  const double x2 = gsl_vector_get(x, 2);
  const double x3 = gsl_vector_get(x, 3);

  np = (np - x3) * (1. - x1) / 2. / (1. - x3 / x2);

  struct crust_fun_4d functs;
  functs = calc_crust_fun_4d(satdata, sparams, x0, x1, x2, np, x3, tt, phase);

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

struct compo calc_icrust4d_composition(double nb_, double tt_, char phase[],
    double *guess, struct parameters satdata, struct sf_params sparams) {
  struct compo eq;

  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *           s;

  int    status;
  size_t iter = 0;

  double aa_old, aa_new, astep;
  double basym_old, basym_new, bstep;
  double n0_old, n0_new, rstep;
  double ng_old, ng_new, gstep;

  struct rparams_crust p;
  p.np = nb_; // modification in assign_icrust_fun_4d (DIRTY!)
  p.tt = tt_;
  strcpy(p.phase, phase);
  p.satdata = satdata;
  p.sparams = sparams;

  const size_t n = 4;
  gsl_vector * x = gsl_vector_alloc(n);
  gsl_vector_set(x, 0, guess[0]);
  gsl_vector_set(x, 1, guess[1]);
  gsl_vector_set(x, 2, guess[2]);
  gsl_vector_set(x, 3, guess[3]);

  gsl_multiroot_function f = {&assign_icrust_fun_4d, n, &p};

  T = gsl_multiroot_fsolver_broyden;
  s = gsl_multiroot_fsolver_alloc(T, 4);
  gsl_multiroot_fsolver_set(s, &f, x);

  do {
    aa_old    = gsl_vector_get(s->x, 0);
    basym_old = gsl_vector_get(s->x, 1);
    n0_old    = gsl_vector_get(s->x, 2);
    ng_old    = gsl_vector_get(s->x, 3);

    iter++;

    status = gsl_multiroot_fsolver_iterate(s);

    /* if (gsl_vector_get (s->f, 0) != gsl_vector_get (s->f, 0) */
    /*         && gsl_vector_get (s->f, 1) != gsl_vector_get (s->f, 1) */
    /*         && gsl_vector_get (s->f, 2) != gsl_vector_get (s->f, 2) */
    /*         && gsl_vector_get (s->f, 3) != gsl_vector_get (s->f, 3)) */
    /*     iter = 1000; // to avoid 'matrix is singular' error */

    aa_new    = gsl_vector_get(s->x, 0);
    basym_new = gsl_vector_get(s->x, 1);
    n0_new    = gsl_vector_get(s->x, 2);
    ng_new    = gsl_vector_get(s->x, 3);
    astep     = aa_new - aa_old;
    bstep     = basym_new - basym_old;
    rstep     = n0_new - n0_old;
    gstep     = ng_new - ng_old;

    if (status)
      break;

    // dirty backstepping
    while (gsl_vector_get(s->x, 0) < 0.) {
      astep  = astep / 4.;
      aa_new = aa_old + astep;
      gsl_vector_set(x, 0, aa_new);
      gsl_vector_set(x, 1, basym_new);
      gsl_vector_set(x, 2, n0_new);
      gsl_vector_set(x, 3, ng_new);
      gsl_multiroot_fsolver_set(s, &f, x);
      aa_new = gsl_vector_get(s->x, 0.);
    }
    while (gsl_vector_get(s->x, 1) < -1.0 || gsl_vector_get(s->x, 1) > 1.0) {
      bstep     = bstep / 4.;
      basym_new = basym_old + bstep;
      gsl_vector_set(x, 0, aa_new);
      gsl_vector_set(x, 1, basym_new);
      gsl_vector_set(x, 2, n0_new);
      gsl_vector_set(x, 3, ng_new);
      gsl_multiroot_fsolver_set(s, &f, x);
      basym_new = gsl_vector_get(s->x, 1);
    }
    while (gsl_vector_get(s->x, 2) < 0. || gsl_vector_get(s->x, 2) > 1.0) {
      rstep  = rstep / 4.;
      n0_new = n0_old + rstep;
      gsl_vector_set(x, 0, aa_new);
      gsl_vector_set(x, 1, basym_new);
      gsl_vector_set(x, 2, n0_new);
      gsl_vector_set(x, 3, ng_new);
      gsl_multiroot_fsolver_set(s, &f, x);
      n0_new = gsl_vector_get(s->x, 2);
    }
    while (gsl_vector_get(s->x, 3) < 1.e-10 || gsl_vector_get(s->x, 3) > nb_) {
      /* if (nb_ < 3.e-4) // to avoid 'matrix is singular' error */
      /* { */
      /*     iter = 1000; */
      /*     break; */
      /* } */
      gstep  = gstep / 4.;
      ng_new = ng_old + gstep;
      gsl_vector_set(x, 0, aa_new);
      gsl_vector_set(x, 1, basym_new);
      gsl_vector_set(x, 2, n0_new);
      gsl_vector_set(x, 3, ng_new);
      gsl_multiroot_fsolver_set(s, &f, x);
      ng_new = gsl_vector_get(s->x, 3);
    }

    if (gsl_vector_get(s->x, 1) > 0.99)
      iter = 1000;

    status = gsl_multiroot_test_residual(s->f, 9e-9);
  }

  while (status == GSL_CONTINUE && iter < 1000);

  if (iter == 1000) {
    eq.aa  = NAN;
    eq.del = NAN;
    eq.n0  = NAN;
    eq.ng  = NAN;
  } else {
    eq.aa  = gsl_vector_get(s->x, 0);
    eq.del = gsl_vector_get(s->x, 1);
    eq.n0  = gsl_vector_get(s->x, 2);
    eq.ng  = gsl_vector_get(s->x, 3);
  }

  guess[0] = eq.aa;
  guess[1] = eq.del;
  guess[2] = eq.n0;
  guess[3] = eq.ng;

  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);

  return eq;
}

int assign_icrust_fun_zz_fixed(
    const gsl_vector *x, void *params, gsl_vector *f) {
  double np = ((struct rparams_crust *)params)->np;
  double tt = ((struct rparams_crust *)params)->tt;
  int    zz = ((struct rparams_crust *)params)->zz;
  char   phase[3];
  strcpy(phase, ((struct rparams_crust *)params)->phase);
  struct parameters satdata = ((struct rparams_crust *)params)->satdata;
  struct sf_params  sparams = ((struct rparams_crust *)params)->sparams;

  const double x0 = gsl_vector_get(x, 0);
  const double x1 = gsl_vector_get(x, 1);
  const double x2 = gsl_vector_get(x, 2);

  np = (np - x2) * (float)zz / x0 / (1. - x2 / x1);

  struct crust_fun_4d functs;
  functs = calc_crust_fun_zz_fixed(
      satdata, sparams, x0, 1. - 2. * (float)zz / x0, x1, np, x2, tt, phase);

  const double y0 = functs.f_stability;
  const double y1 = functs.f_muneq;
  const double y2 = functs.f_presseq;

  gsl_vector_set(f, 0, y0);
  gsl_vector_set(f, 1, y1);
  gsl_vector_set(f, 2, y2);

  return GSL_SUCCESS;
}
struct compo calc_icrust_composition_zz_fixed(double nb_, double tt_, int zz_,
    char phase[], double *guess, struct parameters satdata,
    struct sf_params sparams) {
  struct compo eq;

  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *           s;

  int    status;
  size_t iter = 0;

  double aa_old, aa_new, astep;
  double n0_old, n0_new, rstep;
  double ng_old, ng_new, gstep;

  struct rparams_crust p;
  p.np = nb_; // modification in assign_icrust_fun_4d (DIRTY!)
  p.tt = tt_;
  p.zz = zz_;
  strcpy(p.phase, phase);
  p.satdata = satdata;
  p.sparams = sparams;

  const size_t n = 3;
  gsl_vector * x = gsl_vector_alloc(n);
  gsl_vector_set(x, 0, guess[0]);
  gsl_vector_set(x, 1, guess[1]);
  gsl_vector_set(x, 2, guess[2]);

  gsl_multiroot_function f = {&assign_icrust_fun_zz_fixed, n, &p};

  T = gsl_multiroot_fsolver_broyden;
  s = gsl_multiroot_fsolver_alloc(T, 3);
  gsl_multiroot_fsolver_set(s, &f, x);

  do {
    aa_old = gsl_vector_get(s->x, 0);
    n0_old = gsl_vector_get(s->x, 1);
    ng_old = gsl_vector_get(s->x, 2);

    iter++;

    status = gsl_multiroot_fsolver_iterate(s);

    aa_new = gsl_vector_get(s->x, 0);
    n0_new = gsl_vector_get(s->x, 1);
    ng_new = gsl_vector_get(s->x, 2);
    astep  = aa_new - aa_old;
    rstep  = n0_new - n0_old;
    gstep  = ng_new - ng_old;

    if (status)
      break;

    // dirty backstepping
    while (gsl_vector_get(s->x, 0) < zz_) {
      astep  = astep / 4.;
      aa_new = aa_old + astep;
      gsl_vector_set(x, 0, aa_new);
      gsl_vector_set(x, 1, n0_new);
      gsl_vector_set(x, 2, ng_new);
      gsl_multiroot_fsolver_set(s, &f, x);
      aa_new = gsl_vector_get(s->x, 0);
    }
    while (gsl_vector_get(s->x, 1) < 0. || gsl_vector_get(s->x, 1) > 1.0) {
      rstep  = rstep / 4.;
      n0_new = n0_old + rstep;
      gsl_vector_set(x, 0, aa_new);
      gsl_vector_set(x, 1, n0_new);
      gsl_vector_set(x, 2, ng_new);
      gsl_multiroot_fsolver_set(s, &f, x);
      n0_new = gsl_vector_get(s->x, 1);
    }
    while (gsl_vector_get(s->x, 2) < 1.e-10 || gsl_vector_get(s->x, 2) > nb_) {
      gstep  = gstep / 4.;
      ng_new = ng_old + gstep;
      gsl_vector_set(x, 0, aa_new);
      gsl_vector_set(x, 1, n0_new);
      gsl_vector_set(x, 2, ng_new);
      gsl_multiroot_fsolver_set(s, &f, x);
      ng_new = gsl_vector_get(s->x, 2);
    }

    status = gsl_multiroot_test_residual(s->f, 9e-9);
  }

  while (status == GSL_CONTINUE && iter < 1000);

  if (iter == 1000) {
    eq.aa  = NAN;
    eq.del = NAN;
    eq.n0  = NAN;
    eq.ng  = NAN;
  } else {
    eq.aa  = gsl_vector_get(s->x, 0);
    eq.del = 1. - 2. * (float)zz_ / eq.aa;
    eq.n0  = gsl_vector_get(s->x, 1);
    eq.ng  = gsl_vector_get(s->x, 2);
  }

  guess[0] = eq.aa;
  guess[1] = eq.n0;
  guess[2] = eq.ng;

  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);

  return eq;
}

struct compo calc_icrust_composition_w_shl(double nb_, double tt_, char phase[],
    double *guess_zz20, struct parameters satdata, struct sf_params sparams) {
  struct compo eq, eq_opt;
  double       fdensws, fdensws_min;

  fdensws_min = 1e99;

  double guess[3];
  guess[0] = guess_zz20[0];
  guess[1] = guess_zz20[1];
  guess[2] = guess_zz20[2];

  double x    = 1.;
  double tt0  = 1.;
  double beta = 1. / tt0 * tan(PI / 2. * (1. - x));

  for (int zz = 20; zz < 61; zz += 2) {
    eq = calc_icrust_composition_zz_fixed(
        nb_, tt_, zz, phase, guess, satdata, sparams);

    if (zz == 20) {
      guess_zz20[0] = eq.aa;
      guess_zz20[1] = eq.n0;
      guess_zz20[2] = eq.ng;
    }

    fdensws = calc_crust_ws_cell_free_energy_density(
                  satdata, sparams, eq, nb_, tt_, phase) +
              get_shell_energy_per_nucleon(nb_, zz) * nb_ *
                  (1. - 2. / PI * atan(beta * tt_));

    if (fdensws < fdensws_min) {
      eq_opt      = eq;
      fdensws_min = fdensws;
    }

    /* if (guess[0] != guess[0]) // old */
    /* { */
    /*     guess[0] = guess_zz20[0]; */
    /*     guess[1] = guess_zz20[1]; */
    /*     guess[2] = guess_zz20[2]; */
    /* } */
  }

  return eq_opt;
}

double calc_muncl(struct parameters satdata, struct sf_params sparams,
    struct compo eq, double nb_, double tt_, char phase[]) {
  double np;
  double epsa, epsb;
  double fion_ap, fion_am;
  double fion_bp, fion_bm;
  double dfiondaa;
  double dfionddel;

  np = nb_ * (1. - eq.del) / 2.;

  epsa    = 0.001;
  epsb    = 0.0001;
  fion_ap = calc_ion_free_en(
      satdata, sparams, eq.aa + epsa, eq.del, eq.n0, np, eq.ng, tt_, phase);
  fion_am = calc_ion_free_en(
      satdata, sparams, eq.aa - epsa, eq.del, eq.n0, np, eq.ng, tt_, phase);
  fion_bp = calc_ion_free_en(
      satdata, sparams, eq.aa, eq.del + epsb, eq.n0, np, eq.ng, tt_, phase);
  fion_bm = calc_ion_free_en(
      satdata, sparams, eq.aa, eq.del - epsb, eq.n0, np, eq.ng, tt_, phase);
  dfiondaa  = (fion_ap - fion_am) / 2. / epsa;
  dfionddel = (fion_bp - fion_bm) / 2. / epsb;

  return dfiondaa + (1. - eq.del) / eq.aa * dfionddel;
}

double calc_mupcl(struct parameters satdata, struct sf_params sparams,
    struct compo eq, double nb_, double tt_, char phase[]) {
  double np;
  double epsa, epsb;
  double fion_ap, fion_am;
  double fion_bp, fion_bm;
  double dfiondaa;
  double dfionddel;

  np = nb_ * (1. - eq.del) / 2.;

  epsa    = 0.001;
  epsb    = 0.0001;
  fion_ap = calc_ion_free_en(
      satdata, sparams, eq.aa + epsa, eq.del, eq.n0, np, eq.ng, tt_, phase);
  fion_am = calc_ion_free_en(
      satdata, sparams, eq.aa - epsa, eq.del, eq.n0, np, eq.ng, tt_, phase);
  fion_bp = calc_ion_free_en(
      satdata, sparams, eq.aa, eq.del + epsb, eq.n0, np, eq.ng, tt_, phase);
  fion_bm = calc_ion_free_en(
      satdata, sparams, eq.aa, eq.del - epsb, eq.n0, np, eq.ng, tt_, phase);
  dfiondaa  = (fion_ap - fion_am) / 2. / epsa;
  dfionddel = (fion_bp - fion_bm) / 2. / epsb;

  return dfiondaa - (1. + eq.del) / eq.aa * dfionddel;
}

double calc_crust_ws_cell_free_energy_density(struct parameters satdata,
    struct sf_params sparams, struct compo eq, double nb_, double tt_,
    char phase[]) {
  double     np;
  double     vws;
  double     feldenstot;
  double     fion;
  struct hnm ngas;
  double     fgdens;
  double     epsws;

  vws = eq.aa * (1. - eq.ng / eq.n0) / (nb_ - eq.ng);
  np  = eq.aa * (1. - eq.del) / 2. / vws;

  feldenstot = calc_egas_free_energy_density(np, tt_);

  fion = calc_ion_free_en(
      satdata, sparams, eq.aa, eq.del, eq.n0, np, eq.ng, tt_, phase);

  ngas =
      calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, eq.ng, 1., tt_);
  fgdens = eq.ng * ngas.fpernuc;

  epsws = fion / vws + feldenstot + fgdens * (1. - eq.aa / eq.n0 / vws) +
          np * (RMP - RMN) + nb_ * RMN;

  return epsws;
}

double calc_ion_pressure(struct parameters satdata, struct sf_params sparams,
    double aa_, double del_, double n0_, double np_, double ng_, double tt_,
    char phase[]) {
  double epsp;
  double fion_pp, fion_pm;
  double dfiondnp;

  epsp = np_ / 1000.;

  fion_pp = calc_ion_free_en(
      satdata, sparams, aa_, del_, n0_, np_ + epsp, ng_, tt_, phase);
  fion_pm = calc_ion_free_en(
      satdata, sparams, aa_, del_, n0_, np_ - epsp, ng_, tt_, phase);

  dfiondnp = (fion_pp - fion_pm) / 2. / epsp;

  return 2. * np_ * np_ / aa_ / (1. - del_) * dfiondnp;
}

double calc_crust_ws_cell_pressure(struct parameters satdata,
    struct sf_params sparams, struct compo eq, double nb_, double tt_,
    char phase[]) {
  double np;
  double egas_pressure;
  double ion_pressure;
  double ngas_pressure;
  double ws_cell_pressure;

  np            = (nb_ - eq.ng) * (1. - eq.del) / 2. / (1. - eq.ng / eq.n0);
  egas_pressure = calc_egas_pressure(np, tt_);
  ion_pressure  = calc_ion_pressure(
      satdata, sparams, eq.aa, eq.del, eq.n0, np, eq.ng, tt_, phase);
  ngas_pressure =
      calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, eq.ng, 1., tt_)
          .p;
  ws_cell_pressure = egas_pressure + ion_pressure + ngas_pressure;

  return ws_cell_pressure;
}

double approximate_melting_temperature(struct compo comp, double nb_) {
  // see eq. (2.28) of "Neutron Stars 1: Equation of State and Structure"
  double gamma_m = 175.; // ion coupling parameter
  double zz      = comp.aa * (1. - comp.del) / 2.;
  double vws     = comp.aa / (nb_ - comp.ng) * (1. - comp.ng / comp.n0);

  return zz * zz * ALPHAFS * HBARC / gamma_m * pow(4. * PI / 3. / vws, 1. / 3.);
}

double eval_melting_temperature(struct parameters satdata,
    struct sf_params sparams, double nb_, struct compo *eq, int nd_checker) {
  struct compo comp;
  double       guess_oc[3] = {eq->aa, eq->del, eq->n0};
  double       guess_ic[4] = {eq->aa, eq->del, eq->n0, eq->ng};

  /* double guess_zz20[3] = {80., eq->n0, eq->ng}; */

  if (nd_checker == 0) {
    comp =
        calc_ocrust3d_composition(nb_, 0., "sol", guess_oc, satdata, sparams);
  } else if (nd_checker == 1) {
    comp =
        calc_icrust4d_composition(nb_, 0., "sol", guess_ic, satdata, sparams);
    /* comp = calc_icrust_composition_w_shl(nb_, 0., */
    /*         "sol", guess_zz20, satdata, sparams); */
  }

  double tt = approximate_melting_temperature(comp, nb_) * 1.2;

  double fws_sol, fws_liq;

  do {
    if (nd_checker == 0) {
      comp =
          calc_ocrust3d_composition(nb_, tt, "liq", guess_oc, satdata, sparams);
    } else {
      comp =
          calc_icrust4d_composition(nb_, tt, "liq", guess_ic, satdata, sparams);
      /* comp = calc_icrust_composition_w_shl(nb_, tt, */
      /*         "liq", guess_zz20, satdata, sparams); */
    }

    fws_sol = calc_crust_ws_cell_free_energy_density(
        satdata, sparams, comp, nb_, tt, "sol");
    fws_liq = calc_crust_ws_cell_free_energy_density(
        satdata, sparams, comp, nb_, tt, "liq");

    tt -= 0.0002;

    if (tt <= 0.) {
      fprintf(stderr, "WARNING: sign of T_m cannot be negative! "
                      "Tm set to -1!\n");
      return -1;
    }
  } while (fws_liq < fws_sol);

  eq->aa  = comp.aa;
  eq->del = comp.del;
  eq->n0  = comp.n0;
  eq->ng  = comp.ng;

  return tt;
}

void print_state_crust(struct parameters satdata, struct sf_params sparams,
    struct compo eq, double nb_, double tt_, char phase[], FILE *compo,
    FILE *eos) {
  double vws, rws;
  double rhob;
  double pressws;

  vws = eq.aa / (nb_ - eq.ng) * (1. - eq.ng / eq.n0);
  rws = pow(3. * vws / 4. / PI, 1. / 3.);

  rhob = calc_crust_ws_cell_free_energy_density(
             satdata, sparams, eq, nb_, tt_, phase) *
         (ELEMC / 1.e-19) / pow(SPEEDOFL / 1.e8, 2.) * 1.e13;
  pressws = calc_crust_ws_cell_pressure(satdata, sparams, eq, nb_, tt_, phase);

  if (tt_ != 0.) {
    fprintf(compo, "%g %g %g %g %g %g %g %g\n", nb_, tt_, eq.aa, eq.del,
        eq.aa * (1. - eq.del) / 2., eq.n0, eq.ng, rws);
  } else {
    fprintf(compo, "%g %g %g %g %g %g %g\n", nb_, eq.aa, eq.del,
        eq.aa * (1. - eq.del) / 2., eq.n0, eq.ng, rws);
  }
  fprintf(eos, "%g %g\n", rhob, pressws);
}
