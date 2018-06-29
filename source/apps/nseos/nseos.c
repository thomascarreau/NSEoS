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

    struct parameters satdata = assign_param(argv[1]);
    double p = 3.;

    int lines = calc_equation_of_state(satdata, p, argv);

    solve_tov_equation(lines, argv);

    return 0;
}
