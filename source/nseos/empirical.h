#ifndef H_EMPIRICAL
#define H_EMPIRICAL

#include <stdio.h>

struct parameters {
  float rhosat0;
  float lasat0;
  float ksat0;
  float qsat0;
  float zsat0;
  float jsym0;
  float lsym0;
  float ksym0;
  float qsym0;
  float zsym0;
  float barm;
  float bardel;
  float b;
};

struct skyrme_parameters {
  float t0;
  float t1;
  float t2;
  float t3;
  float t4;
  float t5;
  float x0;
  float x1;
  float x2;
  float x3;
  float x4;
  float x5;
  float alpha;
  float beta;
  float gamma;
};

void eval_barm_and_bardel(float effm, float isosplit, float *ks, float *kv);

struct parameters assign_param(char set[], float b);
void print_parameters(struct parameters satdata);

struct skyrme_parameters assign_skyrme_param(char skyrme_params[]);
void print_skyrme_parameters(struct skyrme_parameters skypars);

int read_table_of_sets(
    FILE *sets, struct parameters *satdata, float *m, float *dm);

#endif // H_EMPIRICAL
