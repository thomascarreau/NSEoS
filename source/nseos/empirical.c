#include <math.h>
#include <stdio.h>
#include <string.h>

#include "empirical.h"
#include "modeling.h"

// (effm,isosplit) -> (barm,bardel)
void eval_barm_and_bardel(
    float effm, float isosplit, float *barm, float *bardel) {
  float kv;

  *barm = 1. / effm - 1.; // = ks

  if (isosplit == 0.0) {
    kv = (*barm) - 0.5 * isosplit * (1. + (*barm)) * (1. + (*barm));
    // eq (8) of arXiv:1708:06894
  } else {
    kv = (sqrt(((*barm) + 1.) * ((*barm) + 1.) * isosplit * isosplit + 1.) +
             (*barm) * isosplit - 1.) /
         isosplit;
    // see: ref [113] of arXiv:1708:06894
  }

  *bardel = *barm - kv;
}

struct parameters assign_param(char set[], float b) {
  struct parameters satdata;
  FILE *            fin;
  float             effm;
  float             isosplit;

  char path_of_set[50] = "../../input/satdata/";
  strcat(path_of_set, set);

  fin = fopen(path_of_set, "r");

  fscanf(fin, "%f %f %f %f %f", &satdata.rhosat0, &satdata.lasat0,
      &satdata.ksat0, &satdata.qsat0, &satdata.zsat0);
  fscanf(fin, "%f %f %f %f %f", &satdata.jsym0, &satdata.lsym0, &satdata.ksym0,
      &satdata.qsym0, &satdata.zsym0);
  fscanf(fin, "%f %f", &effm, &isosplit);

  eval_barm_and_bardel(effm, isosplit, &satdata.barm, &satdata.bardel);

  satdata.b = b;

  fclose(fin);

  return satdata;
}

void print_parameters(struct parameters satdata) {
  fprintf(stderr, "%g %g %g %g %g\n", satdata.rhosat0, satdata.lasat0,
      satdata.ksat0, satdata.qsat0, satdata.zsat0);
  fprintf(stderr, "%g %g %g %g %g\n", satdata.jsym0, satdata.lsym0,
      satdata.ksym0, satdata.qsym0, satdata.zsym0);
  fprintf(stderr, "%g %g\n", satdata.barm, satdata.bardel);
  fprintf(stderr, "%g\n", satdata.b);
}

struct skyrme_parameters assign_skyrme_param(char skyrme_params[]) {
  struct skyrme_parameters skypars;
  FILE *                   fin;

  char path_of_skyrme_params[50] = "../../input/skyrme_params/";
  strcat(path_of_skyrme_params, skyrme_params);

  fin = fopen(path_of_skyrme_params, "r");

  fscanf(fin, "%f %f %f %f %f %f", &skypars.t0, &skypars.t1, &skypars.t2,
      &skypars.t3, &skypars.t4, &skypars.t5);
  fscanf(fin, "%f %f %f %f %f %f", &skypars.x0, &skypars.x1, &skypars.x2,
      &skypars.x3, &skypars.x4, &skypars.x5);
  fscanf(fin, "%f %f %f", &skypars.alpha, &skypars.beta, &skypars.gamma);
  fclose(fin);

  return skypars;
}

void print_skyrme_parameters(struct skyrme_parameters skypars) {
  fprintf(stderr, "t0    = %12f\t MeV.fm^3\n", skypars.t0);
  fprintf(stderr, "t1    = %12f\t MeV.fm^5\n", skypars.t1);
  fprintf(stderr, "t2    = %12f\t MeV.fm^5\n", skypars.t2);
  fprintf(stderr, "t3    = %12f\t MeV.fm^(3+3*alpha)\n", skypars.t3);
  fprintf(stderr, "t4    = %12f\t MeV.fm^(5+3*beta)\n", skypars.t4);
  fprintf(stderr, "t5    = %12f\t MeV.fm^(5+3*gamma)\n", skypars.t5);
  fprintf(stderr, "x0    = %12f\n", skypars.x0);
  fprintf(stderr, "x1    = %12f\n", skypars.x1);
  fprintf(stderr, "x2    = %12f\n", skypars.x2);
  fprintf(stderr, "x3    = %12f\n", skypars.x3);
  fprintf(stderr, "x4    = %12f\n", skypars.x4);
  fprintf(stderr, "x5    = %12f\n", skypars.x5);
  fprintf(stderr, "alpha = %12f\n", skypars.alpha);
  fprintf(stderr, "beta  = %12f\n", skypars.beta);
  fprintf(stderr, "gamma = %12f\n", skypars.gamma);
}

int read_table_of_sets(
    FILE *sets, struct parameters *satdata, float *m, float *dm) {
  float effm, isosplit;

  int buffer = fscanf(sets, "%f %f %f %f %f %f %f %f %f %f %f %f %f",
      &satdata->rhosat0, &satdata->lasat0, &satdata->ksat0, &satdata->qsat0,
      &satdata->zsat0, &satdata->jsym0, &satdata->lsym0, &satdata->ksym0,
      &satdata->qsym0, &satdata->zsym0, &effm, &isosplit, &satdata->b);

  if (buffer == 13) {
    eval_barm_and_bardel(effm, isosplit, &satdata->barm, &satdata->bardel);

    *m  = effm;
    *dm = isosplit;

    return 0;
  } else
    return 1;
}
