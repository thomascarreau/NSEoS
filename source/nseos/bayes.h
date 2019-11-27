#ifndef H_BAYES
#define H_BAYES

#include <stdio.h>

#include "empirical.h"

// ====================
//  Low-density filter:
//  Drischler el al. ab-initio calculations
//  PRC93,054314(2016)
// ====================
void get_low_density_posterior(FILE *prior, FILE *posterior);

// ====================
//  High-density filter:
//  * causality
//  * stability
//  * positive symmetry energy
//  * maximum mass
// ====================
void get_high_density_posterior(FILE *prior, FILE *posterior_par,
    FILE *posterior_spar, FILE *posterior_obs, size_t posterior_size);

void calc_observables(
    FILE *posterior, int p_switch, FILE *observables, FILE *new_posterior);

#endif // H_BAYES
