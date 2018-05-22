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

void calc_weights_for_masses_filter(float chi2[N], float w[N]);

struct stats
{
    float average[N_PARAMS];
    float variance[N_PARAMS];
    float deviation[N_PARAMS];
    float correlation[N_PARAMS][N_PARAMS];
};
struct stats calc_stats(float data[N_PARAMS][N], float w[N]);

#endif // H_FUNCTIONS
