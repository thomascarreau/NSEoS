#include <stdio.h>
#include <string.h>
#include <math.h>

#include "modeling.h"
#include "empirical.h"

/* (effm,isosplit) <-> (barm,bardel) */
void eval_barm_and_bardel(float effm, float isosplit,
        float *barm, float *bardel)
{
    float kv;

    *barm = 1./effm - 1.; // = ks

    if (isosplit == 0.0)
        kv = (*barm) - 0.5*isosplit*(1.+(*barm))*(1.+(*barm)); // eq (8) of arXiv:1708:06894
    else
        kv = (sqrt(((*barm) + 1.)*((*barm) + 1.)*isosplit*isosplit + 1.)
                + (*barm) * isosplit - 1.)/isosplit; // see: ref [103] of arXiv:1708:06894

    *bardel = *barm - kv;
}

struct parameters assign_param(char set[], float b)
{
    struct parameters satdata;
    FILE *fin;
    float effm;
    float isosplit;

    char path_of_set[50] = "../../input/satdata/";
    strcat(path_of_set, set);

    fin = fopen(path_of_set, "r");

    fscanf(fin, "%f %f %f %f %f", &satdata.rhosat0, &satdata.lasat0, &satdata.ksat0, &satdata.qsat0, &satdata.zsat0);
    fscanf(fin, "%f %f %f %f %f", &satdata.jsym0, &satdata.lsym0, &satdata.ksym0, &satdata.qsym0, &satdata.zsym0);
    fscanf(fin, "%f %f", &effm, &isosplit);

    eval_barm_and_bardel(effm, isosplit, &satdata.barm, &satdata.bardel);

    satdata.b = b;

    fclose(fin);

    return satdata;
}

void print_parameters(struct parameters satdata)
{
    fprintf(stderr, "%g %g %g %g %g\n", satdata.rhosat0, satdata.lasat0, satdata.ksat0, satdata.qsat0, satdata.zsat0);
    fprintf(stderr, "%g %g %g %g %g\n", satdata.jsym0, satdata.lsym0, satdata.ksym0, satdata.qsym0, satdata.zsym0);
    fprintf(stderr, "%g %g\n", satdata.barm, satdata.bardel);
    fprintf(stderr, "%g\n", satdata.b);
}

struct skyrme_parameters assign_skyrme_param(struct skyrme_parameters coeff)
{
    FILE *fin;
    float sigmadenom;

    fin = fopen("../../input/coeff.in", "r");
    fscanf(fin, "%f %f %f %f", &coeff.t0, &coeff.t1, &coeff.t2, &coeff.t3);
    fscanf(fin, "%f", &sigmadenom);
    fscanf(fin, "%f %f %f %f", &coeff.x0, &coeff.x1, &coeff.x2, &coeff.x3);
    coeff.sigma = 1/sigmadenom;
    fclose(fin);

    return coeff;
}

int read_table_of_sets(FILE *sets, struct parameters *satdata, float *m, float *dm)
{
    float effm, isosplit;

    int buffer = 
        fscanf(sets, "%f %f %f %f %f %f %f %f %f %f %f %f %f", 
                &satdata->rhosat0, &satdata->lasat0, &satdata->ksat0,
                &satdata->qsat0, &satdata->zsat0, &satdata->jsym0, 
                &satdata->lsym0, &satdata->ksym0, &satdata->qsym0, 
                &satdata->zsym0, &effm, &isosplit, &satdata->b);

    if (buffer == 13)
    {
        eval_barm_and_bardel(effm, isosplit, &satdata->barm, &satdata->bardel);

        *m = effm;
        *dm = isosplit;

        return 0;
    }
    else
        return 1;
}
