#include <math.h>

#include "../../nseos/empirical.h"
#include "../../nseos/eos.h"
#include "../../nseos/tov.h"
#include "../../nseos/modeling.h"

int main(int argc, char* argv[])
{
    if (argc != 6)
    {
        fprintf(stderr, "ERROR: syntax is './nseos set.in crust.out core.out eos.out tov.out'\n");
        return 1;
    }

    // =================== EOS ===================

    float b = 10.*log(2.);
    struct parameters satdata = assign_param(argv[1], b);
    double p = 3.;
    struct transition_qtt tqtt;
    double epst;

    FILE *mycrust = fopen(argv[2], "w+");
    FILE *mycore = fopen(argv[3], "w+");
    FILE *myeos = fopen(argv[4], "w+");

    int hd_checker = 0;
    int lines = calc_equation_of_state(satdata, p, &tqtt, &epst, &hd_checker, mycrust, mycore, myeos);

    fclose(mycrust);
    fclose(mycore);
    fclose(myeos);

    // =================== TOV ===================

    myeos = fopen(argv[4], "r");
    struct tov_solution tovs14;
    FILE *mytov = fopen(argv[5], "w+");

    solve_tov_equation(lines, tqtt.pt, epst, myeos, &tovs14, mytov);

    fprintf(stderr, "==============================================\n\n");
    fprintf(stderr, "Mmax        = %g Msun\n", tovs14.mmax);
    fprintf(stderr, "rhoc_14     = %g g/cm^3\n", tovs14.rhoc);
    fprintf(stderr, "pc_14       = %g dyn/cm^2\n", tovs14.pc);
    fprintf(stderr, "R_14        = %g km\n", tovs14.r);
    fprintf(stderr, "Rcore_14    = %g km\n", tovs14.rcore);
    fprintf(stderr, "Mcore_14    = %g Msun\n", tovs14.mcore);
    fprintf(stderr, "I/MR^2_14   = %g\n", tovs14.i_over_mr2);
    fprintf(stderr, "Icrust/I_14 = %g\n", tovs14.icrust_over_mr2/tovs14.i_over_mr2);

    fclose(myeos);
    fclose(mytov);

    return 0;
}
