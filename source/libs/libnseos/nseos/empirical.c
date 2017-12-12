#include <stdio.h>

#include "empirical.h"

struct parameters assign_param(struct parameters satdata)
{
    FILE *fin;
    float effm;

    fin = fopen("../testing/satdata.in", "r");
    fscanf(fin, "%f %f %f %f %f", &satdata.rhosat0, &satdata.lasat0, &satdata.ksat0, &satdata.qsat0, &satdata.zsat0);
    fscanf(fin, "%f %f %f %f %f", &satdata.jsym0, &satdata.lsym0, &satdata.ksym0, &satdata.qsym0, &satdata.zsym0);
    fscanf(fin, "%f %f", &effm, &satdata.bardel);
    satdata.barm = 1./effm - 1.;
    fclose(fin);

    return satdata;
}

struct skyrme_parameters assign_skyrme_param(struct skyrme_parameters coeff)
{
    FILE *fin;
    float sigmadenom;

    fin = fopen("../testing/coeff.in", "r");
    fscanf(fin, "%f %f %f %f", &coeff.t0, &coeff.t1, &coeff.t2, &coeff.t3);
    fscanf(fin, "%f", &sigmadenom);
    fscanf(fin, "%f %f %f %f", &coeff.x0, &coeff.x1, &coeff.x2, &coeff.x3);
    coeff.sigma = 1/sigmadenom;
    fclose(fin);

    return coeff;
}
