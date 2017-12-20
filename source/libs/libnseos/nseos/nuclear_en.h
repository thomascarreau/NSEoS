#ifndef H_NUCLEAR_EN
#define H_NUCLEAR_EN

#include "nuclear_matter.h"
#include "nuclear_surface_en.h"
#include "coulomb_en.h"

// CLDM meta (ELFc)
double calc_cldm_meta_model_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_, double np_);
double calc_cldm_meta_model_wbbpcorr_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_, double np_, double ng_);

// ETF meta (ELFc)
double calc_etf_meta_model_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_, double np_);
double calc_etf_meta_model_wbbpcorr_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_, double np_, double ng_);

// LS meta (ELFc; SkI' values)
double calc_ls_meta_model_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_, double np_);
double calc_ls_meta_model_wbbpcorr_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_, double np_, double ng_);

#endif // H_NUCLEAR_EN
