#include <math.h>

#include "../../nseos/nuclear_en.h"
#include "../../nseos/crust.h"
#include "../../nseos/core.h"
#include "../../nseos/modeling.h"

#define N (211) // using wc -l
#define N_PARAMS (14) // number of parameters in the matrix

struct parameters read_table_of_sets(FILE *, float *, float *, float *);

struct stats
{
    double average[N_PARAMS];
    double variance[N_PARAMS];
    double deviation[N_PARAMS];
    double correlation[N_PARAMS][N_PARAMS];
};
struct stats calc_stats(double data[N_PARAMS][N], double w[N]);

int main(void)
{
    FILE *sets = NULL;

    sets = fopen("../../input/mw_sets.data", "r");
    if(sets == NULL)
    {
        fprintf(stderr, "ERROR: file issue\n");
        return 1;
    }

    FILE *statistics = NULL;
    FILE *matrix = NULL;
    struct parameters satdata;
    float m, dm;
    float mmax;
    double p[N_PARAMS][N];
    struct stats st;

    for(int i = 0; i < N; i++)
    {
        fprintf(stderr, "Set %d:\n", i+1);

        satdata = read_table_of_sets(sets, &m, &dm, &mmax);

        print_parameters(satdata); // test
        fprintf(stderr, "\n==============================================\n\n");

        p[0][i] = satdata.rhosat0;
        p[1][i] = satdata.lasat0;
        p[2][i] = satdata.ksat0;
        p[3][i] = satdata.qsat0;
        p[4][i] = satdata.zsat0;
        p[5][i] = satdata.jsym0;
        p[6][i] = satdata.lsym0;
        p[7][i] = satdata.ksym0;
        p[8][i] = satdata.qsym0;
        p[9][i] = satdata.zsym0;
        p[10][i] = m;
        p[11][i] = dm;
        p[12][i] = satdata.b;
        p[13][i] = mmax;
    }

    double w[N];

    for(int i = 0; i < N; i++)
        w[i] = 1.;

    st = calc_stats(p, w);

    statistics = fopen("statistics_tov.out", "w+"); 
    matrix = fopen("matrix_tov.out", "w+"); 

    for(int i = 0; i < N_PARAMS; i++)
        fprintf(statistics, "%g %g %g\n", st.average[i], st.variance[i], st.deviation[i]);

    for(int i = 0; i < N_PARAMS; i++)
        fprintf(matrix, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", st.correlation[i][0],
                st.correlation[i][1], st.correlation[i][2], st.correlation[i][3], st.correlation[i][4], 
                st.correlation[i][5], st.correlation[i][6], st.correlation[i][7], st.correlation[i][8], 
                st.correlation[i][9], st.correlation[i][10], st.correlation[i][11], st.correlation[i][12], 
                st.correlation[i][13]);

    fclose(sets);
    fclose(statistics);
    fclose(matrix);

    fprintf(stderr, "\\o/\n");

    return 0;
}

struct parameters read_table_of_sets(FILE *sets, float *m, float *dm, float *mmax)
{
    struct parameters satdata;
    float effm, isosplit, kv;

    // reading the parameters
    fscanf(sets, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f", 
            &satdata.lasat0, &satdata.rhosat0, &satdata.ksat0,
            &satdata.qsat0, &satdata.zsat0, &satdata.jsym0, 
            &satdata.lsym0, &satdata.ksym0, &satdata.qsym0, 
            &satdata.zsym0, &effm, &isosplit, &satdata.b, mmax);

    // (effm,isosplit) <-> (barm,bardel)
    satdata.barm = 1./effm - 1.;

    if (isosplit == 0.0)
        kv = satdata.barm - 0.5*isosplit*(1.+satdata.barm)*(1.+satdata.barm); // eq (8) of arXiv:1708:06894
    else
        kv = (sqrt((satdata.barm + 1.)*(satdata.barm + 1.)*isosplit*isosplit + 1.)
                + satdata.barm * isosplit - 1.)/isosplit; // see: ref [103] of arXiv:1708:06894

    satdata.bardel = satdata.barm - kv;

    *m = effm;
    *dm = isosplit;

    return satdata;
}

struct stats calc_stats(double data[N_PARAMS][N], double w[N])
{
    struct stats result;

    // initializing the sum
    double sumi[N_PARAMS];
    for(int i = 0; i < N_PARAMS; i++)
        sumi[i] = 0.;

    // normalizing constant (=N if w[i]=1 for each set)
    double norm_const = 0.;
    for(int i = 0; i < N; i++)
        norm_const += w[i];

    // average
    for(int i = 0; i < N_PARAMS; i++)
        for(int j = 0; j < N; j++)
            sumi[i] += data[i][j]*w[j];
    for(int i = 0; i < N_PARAMS; i++)
        result.average[i] = sumi[i] / norm_const;

    // initializing the sum
    for(int i = 0; i < N_PARAMS; i++)
        sumi[i] = 0.;

    // variance and deviation
    for(int i = 0; i < N_PARAMS; i++)
        for(int j = 0; j < N; j++)
            sumi[i] += pow(data[i][j] - result.average[i], 2.)*w[j];
    for(int i = 0; i < N_PARAMS; i++)
    {
        result.variance[i] = sumi[i] / norm_const;
        result.deviation[i] = sqrt(result.variance[i]);
    }

    // initializing the sum
    double sumij[N_PARAMS][N_PARAMS];
    for(int i = 0; i < N_PARAMS; i++)
        for(int j = 0; j < N_PARAMS; j++)
            sumij[i][j] = 0.;

    // correlation matrix: cor(x,y) = cov(x,y)/(sigx*sigy)
    for(int i = 0; i < N_PARAMS; i++)
        for(int j = 0; j < N_PARAMS; j++)
            for(int k = 0; k < N; k++)
                sumij[i][j] += (data[i][k] - result.average[i])
                    *(data[j][k] - result.average[j])*w[k];
    for(int i = 0; i < N_PARAMS; i++)
        for(int j = 0; j < N_PARAMS; j++)
            result.correlation[i][j] = fabs(sumij[i][j] / norm_const)
                /(result.deviation[i]*result.deviation[j]);

    return result;
}
