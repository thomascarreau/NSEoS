#include <stdio.h>

#include "../nseos/ws_cell.h"

struct functions
{
    double f_stability;
    double f_beta;
    double f_press;
    double esurf, ecoul;
    double dmunp, muel;
};

struct functions calc_functions(double aa_, double del_, double rho0_, double rhop_);
struct functions calc_functions(double aa_, double del_, double rho0_, double rhop_)
{
    struct functions my_functions;
    struct parameters satdata;
    double enuc;
    double epsa;
    double epsd;
    double epsr;
    double enuc_ap, enuc_am;
    double enuc_dp, enuc_dm;
    double enuc_rp, enuc_rm;
    double denucdaa;
    double denucddel;
    double denucdrho0;
    struct gas egas;
    struct coulomb_energy_shift coul_shift;
    double dmu;
    double mueltot;

    satdata = assign_param(satdata);
    enuc = calc_enuc(satdata, aa_, del_, rho0_, rhop_);
    epsa = 0.001;
    epsd = 0.0001;
    epsr = 0.0001;
    enuc_ap = calc_enuc(satdata, aa_+epsa, del_, rho0_, rhop_);
    enuc_am = calc_enuc(satdata, aa_-epsa, del_, rho0_, rhop_);
    enuc_dp = calc_enuc(satdata, aa_, del_+epsd, rho0_, rhop_);
    enuc_dm = calc_enuc(satdata, aa_, del_-epsd, rho0_, rhop_);
    enuc_rp = calc_enuc(satdata, aa_, del_, rho0_+epsr, rhop_);
    enuc_rm = calc_enuc(satdata, aa_, del_, rho0_-epsr, rhop_);
    denucdaa = (enuc_ap - enuc_am)/2./epsa;
    denucddel = (enuc_dp - enuc_dm)/2./epsd;
    denucdrho0 = (enuc_rp - enuc_rm)/2./epsr;
    egas = calc_egas(rhop_);
    coul_shift = calc_coulomb_energy_shift(aa_, del_, rho0_, rhop_);
    dmu = coul_shift.derivative;
    mueltot = egas.mutot;

    /* my_functions.f_stability = denucdaa/aa_ - enuc/aa_/aa_; */
    double eca, ecb;
    struct coulomb_energy_shift screen;
    screen = calc_coulomb_energy_shift(aa_, del_, rho0_, rhop_);
    my_functions.esurf = calc_surface_energy(satdata, aa_, del_, rho0_);
    eca = calc_coulomb_energy(aa_, del_, rho0_);
    ecb = screen.energy;
    my_functions.ecoul = eca + ecb;
    my_functions.f_stability = my_functions.esurf - 2.*my_functions.ecoul; 

    my_functions.dmunp = denucddel*2./aa_ - rmp + rmn; 
    my_functions.muel = mueltot + dmu;
    my_functions.f_beta = denucddel*2./aa_ - mueltot - dmu - rmp + rmn;

    my_functions.f_press = rho0_*rho0_*denucdrho0/aa_;

    return my_functions;
}

int main(int argc, char* argv[])
{
    FILE *fout_test, *fout_plot;

    /* some declarations */
    /* int irho0; */
    int iaa;
    int idel;
    double rhob;
    double aa, del, rho0, rhop;
    struct functions my_functions;

    /* initialization */
    rhob = 1.2e-4;
    aa = 68.5151;
    /* del = 0.427596; */
    rho0 = 0.0931774;

    fout_test = fopen(argv[1], "w");
    fout_plot = fopen(argv[2], "w");

    /* loop on rho0 (saturation density of anm) */
    for(idel = 10; idel <= 50; idel++)
    {
        /* rho0 = irho0/1000.; */
        /* aa = iaa; */
        del = idel/100.;

        rhop = (1.-del)/2.*rhob;

        my_functions = calc_functions(aa, del, rho0, rhop);
        fprintf(fout_test, "%g %g %g %g\n", aa, del, rho0, rhop);
        /* fprintf(fout_plot, "%g %g %g %g\n", aa, my_functions.f_stability, */ 
        /*         my_functions.f_beta, my_functions.f_press); */
        fprintf(fout_plot, "%g %g %g %g %g %g %g %g %g\n", aa, del, rho0, rhop,
                my_functions.esurf, 2.*my_functions.ecoul,
                my_functions.dmunp, my_functions.muel,
                my_functions.f_press);
    }

    fclose(fout_plot);
    fclose(fout_test);

    return 0;
}
