#ifndef H_EOS
#define H_EOS

#include "nuclear_surface_en.h"
#include "empirical.h"

struct transition_qtt {
  double nt;
  double pt;
};

int calc_zero_temperature_equation_of_state(struct parameters satdata,
    struct sf_params sparams, struct transition_qtt *tqtt, int *hd_checker,
    FILE *crust, FILE *core, FILE *eos);

void eval_transition_qtt(struct parameters satdata, struct sf_params sparams,
    struct transition_qtt *tqtt);

#endif // H_EOS
