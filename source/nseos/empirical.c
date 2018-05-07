#include <stdio.h>
#include <string.h>
#include <math.h>

#include "modeling.h"
#include "empirical.h"

struct parameters assign_param(char set[])
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
    satdata.barm = 1./effm - 1.;
    satdata.bardel = satdata.barm - isosplit;
    satdata.b = 10.*log(2.);

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
