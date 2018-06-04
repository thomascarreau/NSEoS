#include <math.h>

#include "../../nseos/crust.h"
#include "../../nseos/core.h"
#include "../../nseos/modeling.h"

#include "functions.h"

int read_table_of_sets(FILE *sets, struct parameters *satdata, float *m, float *dm)
{
    float effm, isosplit, kv;

    int buffer = 
        fscanf(sets, "%f %f %f %f %f %f %f %f %f %f %f %f %f", 
                &satdata->lasat0, &satdata->rhosat0, &satdata->ksat0,
                &satdata->qsat0, &satdata->zsat0, &satdata->jsym0, 
                &satdata->lsym0, &satdata->ksym0, &satdata->qsym0, 
                &satdata->zsym0, &effm, &isosplit, &satdata->b);

    if (buffer == 13)
    {
        // (effm,isosplit) <-> (barm,bardel)
        satdata->barm = 1./effm - 1.;

        if (isosplit == 0.0)
            kv = satdata->barm - 0.5*isosplit*(1.+satdata->barm)*(1.+satdata->barm); // eq (8) of arXiv:1708:06894
        else
            kv = (sqrt((satdata->barm + 1.)*(satdata->barm + 1.)*isosplit*isosplit + 1.)
                    + satdata->barm * isosplit - 1.)/isosplit; // see: ref [103] of arXiv:1708:06894

        satdata->bardel = satdata->barm - kv;

        *m = effm;
        *dm = isosplit;

        return 0;
    }
    else
        return 1;
}

void eval_transition_qtt(struct parameters satdata, double p,
        struct transition_qtt *tqtt)
{
    struct sf_params sparams = fit_sf_params(satdata, p);
    struct compo comp;
    double guess_oc[3] = {60., 0.15, 0.9*satdata.rhosat0}; // initial guess for the outer crust
    double muncl = -1.; // sign of muncl is negative is the outer crust
    double nb = 1.e-6;

    while(1)
    {
        comp = calc_ocrust3d_composition(nb, guess_oc, satdata, sparams);
        if (guess_oc[0] != guess_oc[0]) // exit if nan
            return;

        muncl = calc_muncl(satdata, sparams, comp, nb);
        if (muncl > 0.) // neutron drip -> transtion to inner crust
            break;

        nb += nb/50.;
    }

    double guess_ic[4] = {guess_oc[0], guess_oc[1], guess_oc[2], 1.e-4}; // initial guess for the inner crust
    struct core_compo ccomp;
    double guess_npecore = 0.7; // initial guess for the core
    double epsws_ic;
    double epsws_core;
    int transition = 0;

    while(1)
    {
        comp = calc_icrust4d_composition(nb, guess_ic, satdata, sparams);
        if (guess_ic[0] != guess_ic[0]) // break if nan
            return;

        // calculation of the energy density in the cell in the inner crust
        epsws_ic = calc_crust_ws_cell_energy_density(satdata, sparams, comp, nb);

        ccomp = calc_npecore_composition(nb, &guess_npecore, satdata);
        if (guess_npecore != guess_npecore) // break if nan
            return;

        // calculation of the energy density in the cell in the core
        epsws_core = calc_core_ws_cell_energy_density(satdata, ccomp, nb);

        if (epsws_core < epsws_ic) // crust-core transition
        {
            fprintf(stderr, "e_core < e_crust\n");
            transition = 1; // transition occurs
            break;
        }

        nb += 0.0001;
    }

    if (transition == 0)
        fprintf(stderr, "e_core - e_crust = %g MeV/fm^3\n", epsws_core - epsws_ic);

    fprintf(stderr, "\n==============================================\n");
    fprintf(stderr, "==============================================\n\n");

    tqtt->nt = nb; 
    tqtt->pt = calc_core_ws_cell_pressure(satdata, ccomp, nb);
}

void calc_weights_for_masses_filter(float chi2[N], float w[N]) // where chi2 is associated to the masses filter
{
    for(int i = 0; i < N; i++)
        w[i] = exp(-chi2[i]/2.);
}

struct stats calc_stats(float data[N_PARAMS][N], float w[N])
{
    struct stats result;

    // initializing the sum
    float sumi[N_PARAMS];
    for(int i = 0; i < N_PARAMS; i++)
        sumi[i] = 0.;

    // normalizing constant (=N if w[i]=1 for each set)
    float norm_const = 0.;
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
    float sumij[N_PARAMS][N_PARAMS];
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
