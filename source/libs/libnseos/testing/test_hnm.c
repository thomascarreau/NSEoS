#include <stdio.h>

#include "../nseos/hnm.h"

int main(int argc, char* argv[])
{
    FILE *fout_plot;
    struct parameters satdata0;
    struct hnm my_hnm;
    double rho;
    double del;

    del = 0.;
    satdata0 = assign_param(satdata0);

    fout_plot = fopen(argv[1], "a+");

    for(rho = 1.e-4; rho < 3.*satdata0.rhosat0; rho += 0.01) {
        my_hnm = calc_hnm(satdata0,rho,del); 
        fprintf(fout_plot, "%g %g %g\n", rho, del, my_hnm.enpernuc);
    }

    fclose(fout_plot);

    return 0;
}
