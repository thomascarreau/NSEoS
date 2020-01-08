#include <math.h>

#include "../../nseos/coulomb.h"
#include "../../nseos/crust.h"
#include "../../nseos/mathconst.h"
#include "../../nseos/modeling.h"
#include "../../nseos/nuclear_en.h"
#include "../../nseos/phyconst.h"

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "ERROR: syntax is './crust_crystallization "
                    "set.in crust.out eos.out'\n");
    return 1;
  }

  // nuclear modeling
  float             b       = 10. * log(2.);
  struct parameters satdata = assign_param(argv[1], b);
  double            p       = 3.0;
  struct sf_params  sparams = fit_sf_params(satdata, p, TABLE_FOR_SFPAR);

  double       nb = 1.0e-8; // initial baryon density
  double       tm;
  struct compo comp;

  // initial guess
  comp.aa  = 60.;
  comp.del = 0.15;
  comp.n0  = 0.1595;
  comp.ng  = 1.e-4;

  FILE *crust = fopen(argv[2], "w+");
  FILE *eos   = fopen(argv[3], "w+");

  int nd_checker = 0;

  while (nd_checker == 0) // up to neutron drip density
  {
    tm = eval_melting_temperature(satdata, sparams, nb, &comp, 0);

    print_state_crust(satdata, sparams, comp, nb, tm, "sol", crust, eos);

    nb += nb / 20.;

    if (calc_crust_ws_cell_neutron_chemical_potential(
            satdata, sparams, comp, nb, tm, "sol") > RMN) {
      nd_checker = 1;
    }
  }

  fprintf(stderr, "INFO: neutron drip reached!\n");

  if (nb < 2.7e-4) {
    // to avoid the 'out of table' issue when sh. corr. are added
    nb = 2.8e-4;
  }

  while (nb < 0.04) {
    tm = eval_melting_temperature(satdata, sparams, nb, &comp, 1);

    if (tm == -1 || tm != tm) {
      fprintf(stderr, "ERROR: sign of T_m is negative or T_m = NaN!\n");
      break;
    }

    print_state_crust(satdata, sparams, comp, nb, tm, "sol", crust, eos);

    nb += nb / 10.;
  }

  fclose(crust);
  fclose(eos);

  return 0;
}
