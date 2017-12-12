#include "nuclear_en.h"

double calc_cldm_skyrme_based_nuclear_en(struct parameters satdata, struct skyrme_parameters coeff, 
        double aa_, double ii_, double n0_, double np_)
{
    struct hnm skyrme;
    double bulk_en;
    double surf_en;
    double coul_en;

    skyrme = calc_skyrme_nuclear_matter(coeff, n0_, ii_);
    bulk_en = skyrme.enpernuc*aa_;
    surf_en = calc_ldm_surf_en(aa_);
    coul_en = calc_coulomb_en(satdata, aa_, ii_, n0_, np_);

    return bulk_en + surf_en + coul_en;
}

double calc_cldm_meta_model_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_, double np_)
{
    struct hnm meta;
    double bulk_en;
    double surf_en;
    double coul_en;

    meta = calc_meta_model_nuclear_matter(satdata, max_order, n0_, ii_);
    bulk_en = meta.enpernuc*aa_;
    surf_en = calc_ldm_surf_en(aa_);
    coul_en = calc_coulomb_en(satdata, aa_, ii_, n0_, np_);

    return bulk_en + surf_en + coul_en;
}
