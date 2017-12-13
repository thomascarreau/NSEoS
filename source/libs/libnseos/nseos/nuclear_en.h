#ifndef H_NUCLEAR_EN
#define H_NUCLEAR_EN

#include "nuclear_matter.h"
#include "nuclear_surface_en.h"
#include "coulomb_en.h"

double calc_cldm_skyrme_based_nuclear_en(struct parameters satdata, struct skyrme_parameters coeff, 
        double aa_, double ii_, double n0_, double np_);
double calc_cldm_meta_model_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_, double np_);
double calc_cldm_skyrme_based_wbbpcorr_nuclear_en(struct parameters satdata, struct skyrme_parameters coeff, 
        double aa_, double ii_, double n0_, double np_, double ng_);

#endif // H_NUCLEAR_EN
