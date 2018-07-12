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

    /* float b = 10.*log(2.); */
    float b = 1.1;
    struct parameters satdata = assign_param(argv[1], b);
    double p = 3.;
    struct transition_qtt tqtt;
    double epst;

    FILE *mycrust = fopen(argv[2], "w+");
    FILE *mycore = fopen(argv[3], "w+");
    FILE *myeos = fopen(argv[4], "w+");

    int lines = calc_equation_of_state(satdata, p, &tqtt, &epst, mycrust, mycore, myeos);

    fclose(mycrust);
    fclose(mycore);
    fclose(myeos);

    // =================== TOV ===================

    myeos = fopen(argv[4], "r");
    FILE *mytov = fopen(argv[5], "w+");

    solve_tov_equation(lines, tqtt.pt, epst, myeos, mytov);

    fclose(myeos);
    fclose(mytov);

    return 0;
}
