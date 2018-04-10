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
