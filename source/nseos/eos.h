#ifndef H_EOS
#define H_EOS

#include "empirical.h"

struct transition_qtt
{
    double nt;
    double pt;
};

int calc_zero_temperature_equation_of_state(
        struct parameters satdata, double p, 
        struct transition_qtt *tqtt, double *epst, int *hd_checker, 
        FILE *crust, FILE *core, FILE *eos);

void eval_transition_qtt(struct parameters satdata, double p,
        struct transition_qtt *tqtt, double *epst);

#endif // H_EOS
