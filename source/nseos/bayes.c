#include <math.h>

#include <gsl/gsl_rng.h>

#include "bayes.h"
#include "eos.h"
#include "modeling.h"
#include "nuclear_matter.h"
#include "nuclear_surface_en.h"
#include "tov.h"

void get_low_density_posterior(FILE *prior, FILE *posterior) {
  struct parameters satdata;
  float             m, dm;

  double nn[10] = {0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20};
  double e_sm_min[10] = {-4.4921, -7.6542, -10.1682, -12.1990, -13.8226,
      -15.0815, -16.0017, -16.6000, -16.8879, -16.8733};
  double e_sm_max[10] = {-3.8713, -6.9232, -9.2607, -11.0292, -12.3161,
      -13.0155, -13.2269, -13.0000, -12.3471, -11.2774};
  double e_nm_min[10] = {4.2123, 6.0424, 7.4706, 8.8545, 10.3097, 11.8712,
      13.5397, 15.2000, 16.8821, 18.5931};
  double e_nm_max[10] = {4.3001, 6.2685, 7.9385, 9.6124, 11.3935, 13.3080,
      15.3503, 17.5000, 19.7299, 22.0096};
  struct hnm test_hnm_sm;
  struct hnm test_hnm_nm;
  int        ld_check;

  while (read_table_of_sets(prior, &satdata, &m, &dm) == 0) {
    ld_check = 0;

    for (int i = 0; i < 10; i++) {
      test_hnm_sm = calc_meta_model_nuclear_matter(
          satdata, TAYLOR_EXP_ORDER, nn[i], 0., 0.);
      test_hnm_nm = calc_meta_model_nuclear_matter(
          satdata, TAYLOR_EXP_ORDER, nn[i], 1., 0.);

      if (test_hnm_sm.enpernuc < 1.05 * e_sm_min[i] ||
          test_hnm_sm.enpernuc > 0.95 * e_sm_max[i] ||
          test_hnm_nm.enpernuc < 0.95 * e_nm_min[i] ||
          test_hnm_nm.enpernuc > 1.05 * e_nm_max[i]) {
        ld_check = 1;
        break;
      }
    }

    if (ld_check == 0)
      fprintf(posterior, "%g %g %g %g %g %g %g %g %g %g %g %g %g\n",
          satdata.rhosat0, satdata.lasat0, satdata.ksat0, satdata.qsat0,
          satdata.zsat0, satdata.jsym0, satdata.lsym0, satdata.ksym0,
          satdata.qsym0, satdata.zsym0, m, dm, satdata.b);
  }
}

void get_high_density_posterior(FILE *prior, FILE *posterior_par,
    FILE *posterior_spar, FILE *posterior_obs, const double mmax_obs,
    size_t posterior_size) {
  struct parameters satdata;
  float             m, dm;

  size_t set_no = 0;
  size_t count  = 0;
  char   eos[128];
  char   tov[128];

  double           p = 3.0;
  struct sf_params sparams;

  fprintf(posterior_par,
      "nsat,esat,ksat,qsat,zsat,esym,lsym,ksym,qsym,zsym,m,dm,b\n");
  fprintf(posterior_spar, "sigma0,b,sigma0c,beta,chi2\n");
  fprintf(posterior_obs,
      "mmax,rhoc14,pc14,r14,nt,pt,rcore14,mcore14,inorm14,ifrac14,lambda14\n");

  while (read_table_of_sets(prior, &satdata, &m, &dm) == 0 &&
         count <= posterior_size) {
    set_no++;
    fprintf(stderr, "- Set no %zu -\n", set_no);

    // CAUSALITY, STABILITY, SYM EN ===================================
    sparams = fit_sf_params(satdata, p, TABLE_FOR_SFPAR);
    struct transition_qtt tqtt;
    tqtt.nt          = 0.0005;
    int   hd_checker = 0;
    FILE *mycrust    = fopen("crust.out", "w+");
    FILE *mycore     = fopen("core.out", "w+");
    snprintf(eos, sizeof(eos), "eos/eos%zu.out", count);
    FILE *myeos = fopen(eos, "w+");
    int   lines = calc_zero_temperature_equation_of_state(
        satdata, sparams, &tqtt, &hd_checker, mycrust, mycore, myeos);
    fclose(mycrust);
    fclose(mycore);
    fclose(myeos);

    // TEST MAX MASS, OBSERVABLES =====================================
    if (hd_checker == 0 && tqtt.nt > 0.01) {
      myeos = fopen(eos, "r");
      struct tov_solution tovs14;
      snprintf(tov, sizeof(tov), "tov/tov%zu.out", count);
      FILE * mytov = fopen(tov, "w+");
      double mmax =
          solve_tov_equation(lines, tqtt.pt, myeos, &tovs14, 1.4, mytov);
      fclose(myeos);
      fclose(mytov);
      if (mmax < mmax_obs) {
        hd_checker = 1;
        fprintf(stderr, "Mmax = %g Msun < Mmax(obs)\n", mmax);
      } else {
        fprintf(posterior_par, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
            satdata.rhosat0, satdata.lasat0, satdata.ksat0, satdata.qsat0,
            satdata.zsat0, satdata.jsym0, satdata.lsym0, satdata.ksym0,
            satdata.qsym0, satdata.zsym0, m, dm, satdata.b);
        fprintf(posterior_spar, "%g,%g,%g,%g,%g\n", sparams.sigma0, sparams.b,
            sparams.sigma0c, sparams.beta, sparams.chi2);
        fprintf(posterior_obs, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n", mmax,
            tovs14.rhoc, tovs14.pc, tovs14.r, tqtt.nt, tqtt.pt, tovs14.rcore,
            tovs14.mcore, tovs14.i_over_mr2,
            tovs14.icrust_over_mr2 / tovs14.i_over_mr2, tovs14.lambda_dimless);

        count++;
        if (count == posterior_size) {
          return;
        }
      }
    }
  }
}

void calc_observables(
    FILE *posterior, int p_switch, FILE *observables, FILE *new_posterior) {
  struct parameters satdata;
  float             m, dm;

  struct sf_params prms;

  FILE * ftov = NULL;
  double mmax;

  if (p_switch == 0) // then p=3
  {
    while (read_table_of_sets(posterior, &satdata, &m, &dm) == 0) {
      prms = fit_sf_params(satdata, 3.0, TABLE_FOR_SFPAR);
      struct transition_qtt tqtt;
      int                   hd_checker = 0;

      FILE *fcrust = fopen("crust.out", "w+");
      FILE *fcore  = fopen("core.out", "w+");
      FILE *feos   = fopen("eos.out", "w+");
      int   lines  = calc_zero_temperature_equation_of_state(
          satdata, prms, &tqtt, &hd_checker, fcrust, fcore, feos);
      fclose(fcrust);
      fclose(fcore);
      fclose(feos);

      feos = fopen("eos.out", "r");
      struct tov_solution tovs14;
      ftov = fopen("tov.out", "w+");
      mmax = solve_tov_equation(lines, tqtt.pt, feos, &tovs14, 1.4, ftov);
      fclose(feos);
      fclose(ftov);

      if (mmax > 1.4 && tqtt.nt > 0.0005 && tqtt.pt > 0.) {
        fprintf(observables, "%g %g %g %g %g %g %g %g %g %g %g\n", mmax,
            tovs14.rhoc, tovs14.pc, tovs14.r, tqtt.nt, tqtt.pt, tovs14.rcore,
            tovs14.mcore, tovs14.i_over_mr2,
            tovs14.icrust_over_mr2 / tovs14.i_over_mr2, tovs14.lambda_dimless);
        fprintf(new_posterior,
            "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
            satdata.rhosat0, satdata.lasat0, satdata.ksat0, satdata.qsat0,
            satdata.zsat0, satdata.jsym0, satdata.lsym0, satdata.ksym0,
            satdata.qsym0, satdata.zsym0, m, dm, satdata.b, prms.p, prms.sigma0,
            prms.b, prms.chi2);
      }
    }
  } else if (p_switch == 1) // then p={2.5, 3.0, 3.5}
  {
    double p[3] = {2.5, 3., 3.5};

    while (read_table_of_sets(posterior, &satdata, &m, &dm) == 0) {
      for (int i = 0; i < 3; i++) {
        prms = fit_sf_params(satdata, p[i], TABLE_FOR_SFPAR);
        struct transition_qtt tqtt;
        int                   hd_checker = 0;

        FILE *fcrust = fopen("crust.out", "w+");
        FILE *fcore  = fopen("core.out", "w+");
        FILE *feos   = fopen("eos.out", "w+");
        int   lines  = calc_zero_temperature_equation_of_state(
            satdata, prms, &tqtt, &hd_checker, fcrust, fcore, feos);
        fclose(fcrust);
        fclose(fcore);
        fclose(feos);

        feos = fopen("eos.out", "r");
        struct tov_solution tovs14;
        ftov = fopen("tov.out", "w+");
        mmax = solve_tov_equation(lines, tqtt.pt, feos, &tovs14, 1.4, ftov);
        fclose(feos);
        fclose(ftov);

        if (mmax > 1.4 && tqtt.nt > 0.0005 && tqtt.pt > 0.) {
          fprintf(observables, "%g %g %g %g %g %g %g %g %g %g %g\n", mmax,
              tovs14.rhoc, tovs14.pc, tovs14.r, tqtt.nt, tqtt.pt, tovs14.rcore,
              tovs14.mcore, tovs14.i_over_mr2,
              tovs14.icrust_over_mr2 / tovs14.i_over_mr2,
              tovs14.lambda_dimless);
          fprintf(new_posterior,
              "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
              satdata.rhosat0, satdata.lasat0, satdata.ksat0, satdata.qsat0,
              satdata.zsat0, satdata.jsym0, satdata.lsym0, satdata.ksym0,
              satdata.qsym0, satdata.zsym0, m, dm, satdata.b, prms.p,
              prms.sigma0, prms.b, prms.chi2);
        }
      }
    }
  } else {
    fprintf(stderr, "ERROR: p_switch should be either 0 or 1");
    return;
  }
}
