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
    FILE *statistics = NULL;
    FILE *matrix = NULL;
    struct parameters satdata;
    float m, dm;
    struct sf_params sparams;
    struct transtion_qtt tqtt;
    double p[N_PARAMS][N];
    double chi2[N];
    struct stats st;

    posterior = fopen("posterior.out", "w+"); 

    for(int i = 0; i < N; i++)
    {
        fprintf(stderr, "Set %d:\n", i+1);

        satdata = read_table_of_sets(sets, &m, &dm);

        print_parameters(satdata); // test
        fprintf(stderr, "\n==============================================\n\n");

        sparams = fit_sf_params(satdata);

        tqtt = eval_transition_qtt(satdata);
        fprintf(posterior, "%g %g %g %g\n", tqtt.nt, tqtt.pt, sparams.sigma0, sparams.b);

        p[0][i] = tqtt.nt;
        p[1][i] = tqtt.pt;
        p[2][i] = satdata.rhosat0;
        p[3][i] = satdata.lasat0;
        p[4][i] = satdata.ksat0;
        p[5][i] = satdata.qsat0;
        p[6][i] = satdata.zsat0;
        p[7][i] = satdata.jsym0;
        p[8][i] = satdata.lsym0;
        p[9][i] = satdata.ksym0;
        p[10][i] = satdata.qsym0;
        p[11][i] = satdata.zsym0;
        p[12][i] = m;
        p[13][i] = dm;
        p[14][i] = satdata.b;
        p[15][i] = sparams.sigma0;
        p[16][i] = sparams.b;

        chi2[i] = sparams.chi2;
    }

    double w[N];

    /* // w/o any filter */
    /* for(int i = 0; i < N; i++) */
    /*     w[i] = 1.0; */

    // w/ masses filter
    calc_weights_for_masses_filter(chi2, w);

    st = calc_stats(p, w);

    statistics = fopen("statistics.out", "w+"); 
    matrix = fopen("matrix.out", "w+"); 

    for(int i = 0; i < N_PARAMS; i++)
        fprintf(statistics, "%g %g %g\n", st.average[i], st.variance[i], st.deviation[i]);

    for(int i = 0; i < N_PARAMS; i++)
        fprintf(matrix, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", st.correlation[i][0],
                st.correlation[i][1], st.correlation[i][2], st.correlation[i][3], st.correlation[i][4], 
                st.correlation[i][5], st.correlation[i][6], st.correlation[i][7], st.correlation[i][8], 
                st.correlation[i][9], st.correlation[i][10], st.correlation[i][11], st.correlation[i][12], 
                st.correlation[i][13], st.correlation[i][14], st.correlation[i][15], st.correlation[i][16]);

    fclose(sets);
    fclose(posterior);
    fclose(statistics);
    fclose(matrix);

    fprintf(stderr, "\\o/\n");

    return 0;
}
