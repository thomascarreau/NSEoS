#ifndef H_NUCLEAR_EN
#define H_NUCLEAR_EN

#include "nuclear_matter.h"
#include "nuclear_surface_en.h"
#include "coulomb.h"

//=============================
//      LDM(ETF,ref)-ELFc(N)     
//=============================
// bulk: ELFc(N) meta-modeling
// surface: LDM (SLy4 value)
double calc_ldm_meta_model_nuclear_en(struct parameters satdata, int max_order, 
        double aa_, double ii_, double n0_);

//=============================
//      ETF-ELFc(N)         
//=============================
// bulk: ELFc(N) meta-modeling
// surface: ETF analytical formula w/ Cfin reference value
double calc_etf_ana_meta_model_nuclear_en(struct parameters satdata, 
        int max_order, 
        double aa_, double ii_, double n0_);

//=============================
//      LS-ELFc(N)     
//=============================
// bulk: ELFc(N) meta-modeling
// surface: LS corrected surface energy
double calc_ls_meta_model_nuclear_free_en(struct parameters satdata, 
        struct sf_params sparams, int max_order, 
        double aa_, double ii_, double n0_, double tt_);

// micro components
double get_shell_en_from_myers_table(int aa_, int zz_);
double calc_pairing_en(int aa_, int zz_);

//=============================
//      LS-MM+micro
//=============================
// bulk: MM
// surface: LS corrected surface energy
// micro: Myers shell energy + LDM pairing term
double calc_ls_meta_model_nuclear_en_micro(struct parameters satdata, 
        struct sf_params sparams, int max_order, 
        double aa_, double ii_, double n0_);

//=============================
//      mass table
//=============================
double calc_electron_binding_energy(int zz_);
double calc_nuclear_mass_from_mass_excess(int aa_, int zz_, double deps_);

//=============================
//      pasta
//=============================
struct pars
{
    double a, b, c;
};
double function_rn(double x, void *params);
double evaluate_rn(double nb_, double ii_, double n0_, double ng_,
        int d_, struct sf_params sparams, char phase[]);
double calc_surface_plus_coulomb_energy_density(double nb_, 
        double ii_, double n0_, double ng_, 
        int d_, struct sf_params sparams, char phase[]);


#endif // H_NUCLEAR_EN
