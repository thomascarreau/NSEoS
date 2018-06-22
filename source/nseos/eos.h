#ifndef H_EOS
#define H_EOS

#include <stdio.h>

int read_table_of_sets(FILE *, struct parameters *satdata, float *m, float *dm);

void calc_equation_of_state(struct parameters satdata, double p, char *argv[]);

struct transition_qtt
{
    double nt;
    double pt;
};
void eval_transition_qtt(struct parameters satdata, double p,
        struct transition_qtt *tqtt);

#endif // H_EOS
