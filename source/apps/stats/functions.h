#ifndef H_FUNCTIONS
#define H_FUNCTIONS

#define N (325) // wc -l jm_sets.data
#define N_PARAMS (17) // number of parameters in the matrix

#include <stdio.h>

#include "../../nseos/empirical.h"

struct parameters read_table_of_sets(FILE *, float *, float *);

struct transtion_qtt
{
    double nt;
    double pt;
};
struct transtion_qtt eval_transition_qtt(struct parameters);

void calc_weights_for_masses_filter(double chi2[N], double w[N]);

struct stats
{
    double average[N_PARAMS];
    double variance[N_PARAMS];
    double deviation[N_PARAMS];
    double correlation[N_PARAMS][N_PARAMS];
};
struct stats calc_stats(double data[N_PARAMS][N], double w[N]);

#endif // H_FUNCTIONS
