#include <math.h>

#include "../../nseos/nuclear_surface_en.h"

#include "functions.h"

int main(void)
{
    FILE *sets = NULL;

    sets = fopen("../../input/jm_sets.data", "r");
    if(sets == NULL)
    {
        fprintf(stderr, "ERROR: file issue\n");
        return 1;
    }

    FILE *posterior = NULL;
    struct parameters satdata;
    double p[3] = {2., 3., 4.};
    float m, dm;
    struct sf_params sparams;
    struct transtion_qtt tqtt;

    posterior = fopen("posterior.out", "w+"); 

    for(int i = 0; i < N; i++)
    {
        fprintf(stderr, "Set %d:\n", i+1);

        satdata = read_table_of_sets(sets, &m, &dm);

        print_parameters(satdata); // test
        fprintf(stderr, "\n==============================================\n\n");

        for(int j = 0; j < 3; j++)
        {
            sparams = fit_sf_params(satdata, p[j]);

            tqtt = eval_transition_qtt(satdata, p[j]);

            if (tqtt.nt > 0.001) // =if transition occurs
                fprintf(posterior, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                        tqtt.nt, tqtt.pt, satdata.rhosat0, satdata.lasat0, satdata.ksat0,
                        satdata.qsat0, satdata.zsat0, satdata.jsym0, satdata.lsym0, 
                        satdata.ksym0, satdata.qsym0, satdata.zsym0, m, dm, satdata.b, 
                        sparams.p, sparams.sigma0, sparams.b, sparams.chi2);
        }
    }

    fclose(sets);
    fclose(posterior);

    fprintf(stderr, "\\o/\n");

    return 0;
}
