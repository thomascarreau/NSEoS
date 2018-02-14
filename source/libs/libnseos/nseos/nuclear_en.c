#include <math.h>

#include "nuclear_en.h"

double calc_sly4_ldm_meta_model_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_, double np_)
{
    struct hnm meta;
    double bulk_en;
    double surf_en;
    double coul_en;

    meta = calc_meta_model_nuclear_matter(satdata, max_order, n0_, ii_);
    bulk_en = meta.enpernuc*aa_;
    surf_en = calc_sly4_ldm_surface_en(aa_);
    coul_en = calc_coulomb_en(satdata, aa_, ii_, n0_, np_);

    return bulk_en + surf_en + coul_en;
}

double calc_etf_meta_model_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_, double np_)
{
    struct hnm meta;
    double bulk_en;
    double surf_en;
    double coul_en;

    meta = calc_meta_model_nuclear_matter(satdata, max_order, n0_, ii_);
    bulk_en = meta.enpernuc*aa_;
    surf_en = calc_etf_surface_en(satdata, aa_, ii_, n0_);
    coul_en = calc_coulomb_en(satdata, aa_, ii_, n0_, np_);

    return bulk_en + surf_en + coul_en;
}

double calc_ls_meta_model_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_, double np_)
{
    struct hnm meta;
    double bulk_en;
    double surf_en;
    double coul_en;

    meta = calc_meta_model_nuclear_matter(satdata, max_order, n0_, ii_);
    bulk_en = meta.enpernuc*aa_;
    surf_en = calc_ls_surface_en(satdata, aa_, ii_);
    coul_en = calc_coulomb_en(satdata, aa_, ii_, n0_, np_);

    return bulk_en + surf_en + coul_en;
}

double calc_ls_etf_meta_model_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_, double np_)
{
    double basym;
    struct hnm meta;
    double bulk_en;
    double surf_en;
    double coul_en;

    basym = calc_bulk_asymmetry(satdata, aa_, ii_);
    meta = calc_meta_model_nuclear_matter(satdata, max_order, n0_, basym);
    bulk_en = meta.enpernuc*aa_;
    surf_en = calc_ls_etf_surface_en(satdata, aa_, ii_,n0_);
    coul_en = calc_coulomb_en(satdata, aa_, ii_, n0_, np_);

    return bulk_en + surf_en + coul_en;
}

double calc_dl_etf_meta_model_nuclear_en(struct parameters satdata, int max_order,
        double aa_, double ii_, double n0_, double np_)
{
    struct hnm meta0;
    double eb0, jsym;
    double r0;
    double sigmas, ss;
    double coul_en;

    meta0 = calc_meta_model_nuclear_matter(satdata, max_order, n0_, 0.);
    eb0 = meta0.enpernuc;
    jsym = meta0.jsym;

    r0 = pow(3./4./PI/n0_,1./3.);       // !!!!!!!!!!!!!!!!!!!!!!!!!!
    sigmas = 1.07763; // in MeV/fm^2    // !!!! Reference values !!!!
    ss = 11.5514; // in MeV             // !!!!!!!!!!!!!!!!!!!!!!!!!!

    coul_en = calc_coulomb_en(satdata, aa_, ii_, n0_, np_);

    return eb0*aa_ + 4.*PI*r0*r0*sigmas*pow(aa_,2./3.) + jsym/(1.+jsym/ss/pow(aa_,1./3.))*ii_*ii_*aa_ + coul_en; 
}
