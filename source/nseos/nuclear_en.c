#include <math.h>

#include "nuclear_en.h"

double calc_ldm_meta_model_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_, double np_)
{
    struct hnm meta;
    double bulk_en;
    double surf_en;
    double coul_en;

    meta = calc_meta_model_nuclear_matter(satdata, max_order, n0_, ii_);
    bulk_en = meta.enpernuc*aa_;
    surf_en = calc_ldm_surface_en(satdata, aa_);
    coul_en = calc_coulomb_en(satdata, aa_, ii_, n0_, np_);

    return bulk_en + surf_en + coul_en;
}

double calc_etf_ana_meta_model_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_, double np_)
{
    struct hnm meta;
    double bulk_en;
    double surf_en;
    double coul_en;

    meta = calc_meta_model_nuclear_matter(satdata, max_order, n0_, ii_);
    bulk_en = meta.enpernuc*aa_;
    surf_en = calc_etf_ana_surface_en(satdata, aa_, ii_, n0_);
    coul_en = calc_coulomb_en(satdata, aa_, ii_, n0_, np_);

    return bulk_en + surf_en + coul_en;
}

double calc_ls_meta_model_nuclear_en(struct parameters satdata, struct sf_params sparams, 
        int max_order, double aa_, double ii_, double n0_, double np_)
{
    struct hnm meta;
    double bulk_en;
    double surf_en;
    double coul_en;

    meta = calc_meta_model_nuclear_matter(satdata, max_order, n0_, ii_);
    bulk_en = meta.enpernuc*aa_;
    surf_en = calc_ls_surface_en(sparams, aa_, ii_, n0_);
    coul_en = calc_coulomb_en(satdata, aa_, ii_, n0_, np_);

    return bulk_en + surf_en + coul_en;
}

double get_shell_en_from_myers_table(int aa_, int zz_)
{
    int aa, zz;
    float shell_en;
    int checker;
    FILE *myers_table = fopen("../../input/myers.table", "r");

    aa = 0;
    zz = 0;
    checker = 0;

    while(fscanf(myers_table, "%d %d %f", &zz, &aa, &shell_en) == 3)
        if(aa == aa_ && zz == zz_)
        {
            checker = 1;
            break;
        }

    if (checker == 0)
    {
        fprintf(stderr, "ERROR: out of Myers table!");
        shell_en = 0.;
    }

    fclose(myers_table);

    return shell_en;
}

double calc_pairing_en(int aa_, int zz_)
{
    double ap = 12.;
    int nn = aa_ - zz_;

    if(nn % 2 == 1 && zz_ % 2 == 1)
        return ap/pow(aa_,1./2.);
    else if(aa_ % 2 == 1)
        return 0.;
    else if(nn % 2 == 0 && zz_ % 2 == 0)
        return -ap/pow(aa_,1./2.);
    else
    {
        fprintf(stderr, "ERROR: check nuclear pairing function!");
        return 0.;
    }
}

double calc_ls_meta_model_nuclear_en_micro(struct parameters satdata, struct sf_params sparams, 
        int max_order, double aa_, double ii_, double n0_, double np_)
{
    struct hnm meta;
    double bulk_en;
    double surf_en;
    double coul_en;
    double shell_en;
    double pairing_en;

    meta = calc_meta_model_nuclear_matter(satdata, max_order, n0_, ii_);
    bulk_en = meta.enpernuc*aa_;
    surf_en = calc_ls_surface_en(sparams, aa_, ii_, n0_);
    coul_en = calc_coulomb_en(satdata, aa_, ii_, n0_, np_);
    shell_en = get_shell_en_from_myers_table((int)aa_, (int)aa_*(1.-ii_)/2.);
    pairing_en = calc_pairing_en((int)aa_, (int)aa_*(1.-ii_)/2.);

    return bulk_en + surf_en + coul_en + shell_en + pairing_en;
}
