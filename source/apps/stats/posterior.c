#include <math.h>

#include "../../nseos/nuclear_surface_en.h"

#include "functions.h"

int main(int argc, char* argv[])
{
    if(argc != 3)
    {
        fprintf(stderr, "ERROR: syntax is './posterior list_of_sets.data posterior.out\n");
    }

    FILE *sets = NULL;

    sets = fopen(argv[1], "r");

    if(sets == NULL)
    {
        fprintf(stderr, "ERROR: file issue\n");
        return 1;
    }

    struct parameters satdata;
    double p[3] = {2.5, 3., 3.5};
    float m, dm;
    struct sf_params sparams;
    struct transition_qtt tqtt;
    FILE *posterior = NULL;

    posterior = fopen(argv[2], "w+"); 

    int line = 1;
    while(read_table_of_sets(sets, &satdata, &m, &dm) == 0)
    {
        print_parameters(satdata);
        fprintf(stderr, "\n==============================================\n\n");

        for(int j = 0; j < 3; j++)
        {
            fprintf(stderr, "Set %d (p=%g):\n", line, p[j]);

            sparams = fit_sf_params(satdata, p[j]);

            tqtt.nt = 0.0005;

            eval_transition_qtt(satdata, sparams.p, &tqtt);

            if (tqtt.nt > 0.0005)
                fprintf(posterior, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                        tqtt.nt, tqtt.pt, satdata.rhosat0, satdata.lasat0, satdata.ksat0,
                        satdata.qsat0, satdata.zsat0, satdata.jsym0, satdata.lsym0, 
                        satdata.ksym0, satdata.qsym0, satdata.zsym0, m, dm, satdata.b, 
                        sparams.p, sparams.sigma0, sparams.b, sparams.chi2);
        }

        line += 1;
    }

    fclose(sets);
    fclose(posterior);

    fprintf(stderr, "\\o/\n");

    return 0;
}
