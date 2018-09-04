#ifndef H_EMPIRICAL
#define H_EMPIRICAL

#include <stdio.h>

struct parameters
{
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

struct skyrme_parameters
{
    float t0;
    float t1;
    float t2;
    float t3;
    float sigma;
    float x0;
    float x1;
    float x2;
    float x3;
};

void eval_barm_and_bardel(float effm, float isosplit,
        float *ks, float *kv);
struct parameters assign_param(char set[], float b);
void print_parameters(struct parameters satdata);
struct skyrme_parameters assign_skyrme_param(struct skyrme_parameters coeff);
int read_table_of_sets(FILE *, struct parameters *satdata, float *m, float *dm);

#endif // H_EMPIRICAL
