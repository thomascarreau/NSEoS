#include <math.h>

#include "../../nseos/nuclear_en.h"
#include "../../nseos/crust.h"
#include "../../nseos/core.h"
#include "../../nseos/modeling.h"

#define N (325) // wc -l list_of_good_sets.data
#define N_PARAMS (17) // number of parameters in the matrix

struct parameters read_table_of_good_sets(FILE *, float *, float *);

struct transtion_qtt
{
    double nt;
    double pt;
};
struct transtion_qtt eval_transition_qtt(struct parameters);

void calc_weights_for_masses_filter(double chi2[N], double w[N]);

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
    FILE *good_sets = NULL;

    good_sets = fopen("../../input/list_of_good_sets.data", "r");
    if(good_sets == NULL)
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
        satdata = read_table_of_good_sets(good_sets, &m, &dm);

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

    fclose(good_sets);
    fclose(posterior);
    fclose(statistics);
    fclose(matrix);

    return 0;
}

struct parameters read_table_of_good_sets(FILE *good_sets, float *m, float *dm)
{
    struct parameters satdata;
    float effm, isosplit, kv;

    // reading the parameters
    fscanf(good_sets, "%f %f %f %f %f %f %f %f %f %f %f %f %f", 
            &satdata.lasat0, &satdata.rhosat0, &satdata.ksat0,
            &satdata.qsat0, &satdata.zsat0, &satdata.jsym0, 
            &satdata.lsym0, &satdata.ksym0, &satdata.qsym0, 
            &satdata.zsym0, &effm, &isosplit, &satdata.b);

    // deleting factors
    satdata.lasat0 = -satdata.lasat0/10.;
    satdata.rhosat0 /= 1000.;
    effm /= 100.;
    isosplit /= 100.;
    satdata.b /= 10.;

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

struct transtion_qtt eval_transition_qtt(struct parameters satdata)
{
    struct transtion_qtt tqtt;
    struct sf_params sparams = fit_sf_params(satdata);
    struct compo comp;
    double epsws_ic;
    double epsws_core;
    int transition = 0;
    double guess_ic[4] = {100., 0.4, 0.14, 1.e-4}; // initial guess for the inner crust
    struct core_compo ccomp;
    double guess_npecore = 0.7; // initial guess for the core
    double nb = 0.001;

    while(1)
    {
        comp = calc_icrust4d_composition(nb, guess_ic, satdata, sparams);
        if (guess_ic[0] != guess_ic[0]) // break if nan
            break;

        if (nb > 0.01)
        {
            // calculation of the energy density in the cell in the inner crust
            epsws_ic = calc_crust_ws_cell_energy_density(satdata, sparams, comp, nb);

            ccomp = calc_npecore_composition(nb, &guess_npecore, satdata);
            if (guess_npecore != guess_npecore) // break if nan
                break;

            // calculation of the energy density in the cell in the core
            epsws_core = calc_core_ws_cell_energy_density(satdata, ccomp, nb);

            if (epsws_core < epsws_ic) // crust-core transition
            {
                fprintf(stderr, "e_core < e_crust\n");
                transition = 1; // transition occurs
                break;
            }
        }

        nb += 0.0001;
    }

    if (transition == 0)
        fprintf(stderr, "e_core - e_crust = %g MeV/fm^3\n", epsws_core - epsws_ic);

    fprintf(stderr, "\n==============================================\n");
    fprintf(stderr, "==============================================\n\n");

    tqtt.nt = nb; 
    tqtt.pt = calc_core_ws_cell_pressure(satdata, ccomp, nb);

    return tqtt;
}

void calc_weights_for_masses_filter(double chi2[N], double w[N]) // where chi2 is associated to the masses filter
{
    for(int i = 0; i < N; i++)
        w[i] = exp(-chi2[i]/2.);
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
