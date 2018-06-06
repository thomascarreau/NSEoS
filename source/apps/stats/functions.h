#ifndef H_FUNCTIONS
#define H_FUNCTIONS

#include <stdio.h>

#include "../../nseos/empirical.h"

int read_table_of_sets(FILE *, struct parameters *satdata, float *m, float *dm);

struct transition_qtt
{
    double nt;
    double pt;
};
void eval_transition_qtt(struct parameters, double p,
        struct transition_qtt *tqtt);

#endif // H_FUNCTIONS
