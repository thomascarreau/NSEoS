#include <math.h>

#include "../../nseos/eos.h"
#include "../../nseos/modeling.h"
#include "../../nseos/tov.h"

int main(int argc, char *argv[]) {
  if (argc != 6) {
    fprintf(stderr, "ERROR: syntax is './nseos set.in "
                    "crust.out core.out eos.out tov.out'\n");
    return 1;
  }

  // =================== EOS ===================

  float                 b       = 10. * log(2.);
  struct parameters     satdata = assign_param(argv[1], b);
  double                p       = 3.0;
  struct transition_qtt tqtt;
  int                   hd_checker = 0;

  FILE *mycrust = fopen(argv[2], "w+");
  FILE *mycore  = fopen(argv[3], "w+");
  FILE *myeos   = fopen(argv[4], "w+");

  int lines = calc_zero_temperature_equation_of_state(
      satdata, p, &tqtt, &hd_checker, mycrust, mycore, myeos);

  fclose(mycrust);
  fclose(mycore);
  fclose(myeos);

  // =================== TOV ===================

  struct tov_solution tovs14;
  double              fixed_m = 1.4;

  myeos       = fopen(argv[4], "r");
  FILE *mytov = fopen(argv[5], "w+");

  double mmax =
      solve_tov_equation(lines, tqtt.pt, myeos, &tovs14, fixed_m, mytov);

  fclose(myeos);
  fclose(mytov);

  if (mmax > fixed_m) {
    fprintf(stderr, "==============================================\n\n");
    fprintf(stderr, "Mmax        = %g Msun\n", mmax);
    fprintf(stderr, "rhoc_14     = %g g/cm^3\n", tovs14.rhoc);
    fprintf(stderr, "pc_14       = %g dyn/cm^2\n", tovs14.pc);
    fprintf(stderr, "R_14        = %g km\n", tovs14.r);
    fprintf(stderr, "Rcore_14    = %g km\n", tovs14.rcore);
    fprintf(stderr, "Mcore_14    = %g Msun\n", tovs14.mcore);
    fprintf(stderr, "I/MR^2_14   = %g\n", tovs14.i_over_mr2);
    fprintf(stderr, "Icrust/I_14 = %g\n",
        tovs14.icrust_over_mr2 / tovs14.i_over_mr2);
    fprintf(stderr, "k2_14       = %g\n", tovs14.k2);
    fprintf(stderr, "lambda_14   = %g\n", tovs14.lambda_dimless);
  } else {
    fprintf(stderr, "==============================================\n\n");
    fprintf(stderr, "Mmax = %g Msun\n", mmax);
  }

  return 0;
}
