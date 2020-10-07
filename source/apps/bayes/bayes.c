#include <math.h>

#include "../../nseos/bayes.h"

int main(void) {

  // LD posterior
  FILE *prior = fopen("prior.in", "r");
  FILE *ld_posterior = fopen("ld_posterior.out", "w+");
  get_low_density_posterior(prior, ld_posterior);
  fclose(prior);
  fclose(ld_posterior);

  // LD+HD posterior
  ld_posterior = fopen("ld_posterior.out", "r");
  FILE *pars = fopen("ldhdp3_pars.csv", "w+");
  FILE *spars = fopen("ldhdp3_spars.csv", "w+");
  FILE *obs = fopen("ldhdp3_obs.csv", "w+");
  get_high_density_posterior(ld_posterior, pars, spars, obs, 1.97, 1000);
  fclose(ld_posterior);
  fclose(pars);
  fclose(spars);
  fclose(obs);

  return 0;
}
