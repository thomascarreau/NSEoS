#include <math.h>

#include "../../nseos/nuclear_surface_en.h"

#include "functions.h"

#define MASSES (1) // 0 -> w/o masses filter

int main(void)
{
    FILE *posterior = NULL;

    posterior = fopen("posterior.out", "r");
    if(posterior == NULL)
    {
        fprintf(stderr, "ERROR: file issue\n");
        return 1;
    }

    FILE *statistics = NULL;
    FILE *matrix = NULL;
    float p[N_PARAMS][N];
    float chi2[N];
    float w[N];
    struct stats st;

    /*
        ___________________
       /
       | 0  -> nt
       | 1  -> pt
       | 2  -> rhosat0
       | 3  -> lasat0
       | 4  -> ksat0
       | 5  -> qsat0
       | 6  -> zsat0
       | 7  -> jsym0
       | 8  -> lsym0
       | 9  -> ksym0
       | 10 -> qsym0
       | 11 -> zsym0
       | 12 -> m
       | 13 -> dm
       | 14 -> b
       | 15 -> sigma0
       | 16 -> bs
       \___________________

       */


    for(int i = 0; i < N; i++)
    {
        fscanf(posterior, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
                &p[0][i], &p[1][i], &p[2][i], &p[3][i], &p[4][i], &p[5][i], &p[6][i],
                &p[7][i], &p[8][i], &p[9][i], &p[10][i], &p[11][i], &p[12][i], &p[13][i],
                &p[14][i], &p[15][i], &p[16][i], &chi2[i]);
    }

    if (MASSES == 1)
        calc_weights_for_masses_filter(chi2, w);
    else
        for(int i = 0; i < N; i++)
            w[i] = 1.0;

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

    fclose(posterior);
    fclose(statistics);
    fclose(matrix);

    fprintf(stderr, "\\o/\n");

    return 0;
}
