// ****************************************************************************
//
//     source/apps/mcp_icrust/mcp_icrust.C
//
// Author:         Thomas Carreau
// Last update:    07/10/2020
// Description:
//  Multicomponent Coulomb plasma calculation in the free neutron regime with
//  perturbative implementation of the nuclear statistical equilibrium
//
// Build:
//  make
//  Run (see list of sets in source/input/satdata):
//  ./mcp_icrust set.in > file.csv
//
// ****************************************************************************

#include <cmath>
#include <iostream>
#include <vector>

#include <gsl/gsl_multiroots.h>

extern "C" {
#include "../../nseos/core.h"
#include "../../nseos/crust.h"
#include "../../nseos/lepton.h"
#include "../../nseos/modeling.h"
#include "../../nseos/nuclear_en.h"
}

// input variables to calculate the equilibrium nuclear density
struct rparams {
  double            np_ocp;
  double            tt;
  double            aa;
  double            del;
  double            ng;
  double            mun_ocp;
  struct parameters satdata;
  struct sf_params  sparams;
};

// function to define the pressure equilibrium equation
int PressureEquilibriumEquation(
    const gsl_vector *x, void *params, gsl_vector *f);

// function to evaluate the equilibrium nuclear density
double NuclearDensity(struct parameters satdata, struct sf_params sparams,
    double np_ocp, double tt, double aa, double del, double ng, double mun_ocp,
    double *guess);

// function to evaluate the rearrangment term
double RearrangementTerm(struct parameters satdata, struct sf_params sparams,
    double nb, double tt, struct compo ocp, double zz);

// function to define the single ion grand canonical potential
double SingleIonGrandCanonicalPotential(struct parameters satdata,
    struct sf_params sparams, double nb, double tt, struct compo ocp,
    double mun_ocp, double mup_ocp, double aa, double zz);

// function to calculate probabilities pj at a given baryon density nb
// with perturbative implementation of nuclear statistical equilibrium
double Probabilities(struct parameters satdata, struct sf_params sparams,
    double nb, double tt, struct compo ocp, double pres_ocp, double np_ocp,
    double mun_ocp, double mup_ocp);

// function to print probabilities pj
void PrintProbabilities(struct parameters satdata, struct sf_params sparams,
    double nb, double tt, char phase[]);

// perform mcp calcuation in free neutron regime at given temperature tt
void InnerCrustMCPCalculation(struct parameters satdata,
    struct sf_params sparams, double tt, char phase[]);

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "ERROR: Syntax is './mcp_icrust set.in'\n";
    return 1;
  }

  // modeling
  struct parameters satdata = assign_param(argv[1], 10.0 * log(2.0));
  struct sf_params  sparams =
      fit_sf_params(satdata, 3.0, (char *)TABLE_FOR_SFPAR);

  // calculation vs baryon density
  double tt      = 1.0 * 1e10 / 1.16e10;
  char   phase[] = "liq";
  std::cout << "pres_ocp,nb,tt,np,ng,mun,mup,vws_avg,aa_avg,aacell_avg,zz_avg,"
            << "qimp,aa_mp,aacell_mp,zz_mp,aa_ocp,aacell_ocp,zz_ocp\n";
  InnerCrustMCPCalculation(satdata, sparams, tt, phase);

  return 0;
}

int PressureEquilibriumEquation(
    const gsl_vector *x, void *params, gsl_vector *f) {
  double            np_ocp  = ((struct rparams *)params)->np_ocp;
  double            tt      = ((struct rparams *)params)->tt;
  double            aa      = ((struct rparams *)params)->aa;
  double            del     = ((struct rparams *)params)->del;
  double            ng      = ((struct rparams *)params)->ng;
  double            mun_ocp = ((struct rparams *)params)->mun_ocp;
  struct parameters satdata = ((struct rparams *)params)->satdata;
  struct sf_params  sparams = ((struct rparams *)params)->sparams;

  const double x0 = gsl_vector_get(x, 0);

  double epsr    = 0.0001;
  double fion_rp = calc_ion_free_en(
      satdata, sparams, aa, del, x0 + epsr, np_ocp, ng, tt, (char *)"liq");
  double fion_rm = calc_ion_free_en(
      satdata, sparams, aa, del, x0 - epsr, np_ocp, ng, tt, (char *)"liq");
  double dfiondn0 = (fion_rp - fion_rm) / 2. / epsr;

  struct hnm ngas =
      calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, ng, 1., tt);
  double presg = ng * (mun_ocp - RMN) - ng * ngas.fpernuc;

  const double y0 = x0 * x0 / aa * dfiondn0 - presg;

  gsl_vector_set(f, 0, y0);

  return GSL_SUCCESS;
}

double NuclearDensity(struct parameters satdata, struct sf_params sparams,
    double np_ocp, double tt, double aa, double del, double ng, double mun_ocp,
    double *guess) {
  double n0;

  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *           s;

  int    status;
  size_t iter = 0;

  struct rparams p;
  p.np_ocp  = np_ocp;
  p.tt      = tt;
  p.aa      = aa;
  p.del     = del;
  p.ng      = ng;
  p.mun_ocp = mun_ocp;
  p.satdata = satdata;
  p.sparams = sparams;

  double       n0_init = *guess;
  const size_t n       = 1;
  gsl_vector * x       = gsl_vector_alloc(n);
  gsl_vector_set(x, 0, n0_init);

  gsl_multiroot_function f = {&PressureEquilibriumEquation, n, &p};

  T = gsl_multiroot_fsolver_dnewton;
  s = gsl_multiroot_fsolver_alloc(T, 1);
  gsl_multiroot_fsolver_set(s, &f, x);

  do {
    iter++;

    status = gsl_multiroot_fsolver_iterate(s);

    if (status)
      break;

    status = gsl_multiroot_test_residual(s->f, 9e-9);
  }

  while (status == GSL_CONTINUE && iter < 100);

  n0 = gsl_vector_get(s->x, 0);

  *guess = n0;

  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);

  return n0;
}

double RearrangementTerm(struct parameters satdata, struct sf_params sparams,
    double nb, double tt, struct compo ocp, double zz) {
  double np_ocp  = (nb - ocp.ng) * (1. - ocp.del) / 2. / (1. - ocp.ng / ocp.n0);
  double zz_ocp  = ocp.aa * (1.0 - ocp.del) / 2.0;
  double vws_ocp = zz_ocp / np_ocp;

  // numerical derivative gives the same values
  double epsp    = np_ocp / 1000.;
  double fion_pp = calc_ion_free_en(satdata, sparams, ocp.aa, ocp.del, ocp.n0,
      np_ocp + epsp, ocp.ng, tt, (char *)"liq");
  double fion_pm = calc_ion_free_en(satdata, sparams, ocp.aa, ocp.del, ocp.n0,
      np_ocp - epsp, ocp.ng, tt, (char *)"liq");
  double dfintdnp_ocp = (fion_pp - fion_pm) / 2. / epsp - tt / np_ocp;

  return zz / vws_ocp * dfintdnp_ocp;
}

double SingleIonGrandCanonicalPotential(struct parameters satdata,
    struct sf_params sparams, double nb, double tt, struct compo ocp,
    double mun_ocp, double mup_ocp, double aa, double zz) {
  // useful variables
  double del    = 1. - 2. * zz / aa;
  double np_ocp = (nb - ocp.ng) * (1. - ocp.del) / 2. / (1. - ocp.ng / ocp.n0);
  double guess  = satdata.rhosat0;
  double n0     = NuclearDensity(
      satdata, sparams, np_ocp, tt, aa, del, ocp.ng, mun_ocp, &guess);

  // nuclear mass
  double mi =
      CALC_NUCLEAR_EN(satdata, sparams, TAYLOR_EXP_ORDER, aa, del, n0, tt) +
      zz * RMP + (aa * (1. - ocp.ng / n0) - zz) * RMN -
      aa / n0 * ocp.ng *
          calc_meta_model_nuclear_matter(
              satdata, TAYLOR_EXP_ORDER, ocp.ng, 1., tt)
              .fpernuc;

  double omega_tilde =
      mi + tt * log(pow(pow(2. * PI * HBARC * HBARC / mi / tt, 0.5), 3.)) +
      calc_finite_size_contrib(aa, del, n0, np_ocp) +
      calc_total_coulomb_contrib(zz, np_ocp, tt) +
      RearrangementTerm(satdata, sparams, nb, tt, ocp, zz) -
      mun_ocp * (aa * (1. - ocp.ng / n0) - zz) - mup_ocp * zz;

  return omega_tilde;
}

double Probabilities(struct parameters satdata, struct sf_params sparams,
    double nb, double tt, struct compo ocp, double pres_ocp, double np_ocp,
    double mun_ocp, double mup_ocp) {
  // initial chemical potential values are the ocp ones
  double mun = mun_ocp;
  double mup = mup_ocp;

  // probability
  double omega_tilde;
  double normalisation;
  double pj_sum;

  // most probable cluster
  double aa_mp;
  double aacell_mp;
  int    zz_mp;
  double pj_mp;

  // average and variance of several quantities
  double vws_avg;
  double aa_avg;
  double aacell_avg;
  double zz_avg;
  double qimp;

  // electron gas pressure
  double np = np_ocp;

  // aa cell ocp
  double zz_ocp     = ocp.aa * (1. - ocp.del) / 2.; // central point
  double aacell_ocp = ocp.aa * (1. - ocp.ng / ocp.n0) + ocp.ng * zz_ocp / np;

  // for equilibrium nuclear density calculation
  double del;
  double guess = satdata.rhosat0;
  double n0;

  std::vector<double> v_zz;
  std::vector<double> v_aa, v_aacell;
  std::vector<double> v_pj;

  // criterion for keeping nuclei
  double pj_non_norm;
  double omega_tilde_ocp = SingleIonGrandCanonicalPotential(
      satdata, sparams, nb, tt, ocp, mun, mup, ocp.aa, zz_ocp);
  double pj_non_norm_ocp = exp(-omega_tilde_ocp / tt);

  // loop for normalisation factor
  normalisation = 0.;

  for (int zz = (int)zz_ocp - 15; zz <= (int)zz_ocp + 15; zz++) {
    for (double aa = 2.5 * zz; aa <= 9.0 * zz; aa += 0.5) {
      omega_tilde = SingleIonGrandCanonicalPotential(
          satdata, sparams, nb, tt, ocp, mun, mup, aa, (double)zz);
      pj_non_norm = exp(-omega_tilde / tt);
      if (pj_non_norm > pj_non_norm_ocp / 1e9) {
        normalisation += pj_non_norm;
        v_zz.push_back((double)zz);
        v_aa.push_back(aa);
        del = 1. - 2. * (double)zz / aa;
        n0  = NuclearDensity(
            satdata, sparams, np, tt, aa, del, ocp.ng, mun, &guess);
        v_aacell.push_back(aa * (1. - ocp.ng / n0) + ocp.ng * (double)zz / np);
        v_pj.push_back(pj_non_norm);
      }
    }
  }

  // loop for pj and average quantities
  pj_sum     = 0.;
  pj_mp      = 0.;
  vws_avg    = 0.;
  aa_avg     = 0.;
  aacell_avg = 0.;
  zz_avg     = 0.;

  for (size_t i = 0; i < v_zz.size(); i++) {
    v_pj[i] /= normalisation;
    pj_sum += v_pj[i];
    if (v_pj[i] > pj_mp) {
      pj_mp     = v_pj[i];
      zz_mp     = v_zz[i];
      aa_mp     = v_aa[i];
      aacell_mp = v_aacell[i];
    }

    vws_avg += v_pj[i] * v_zz[i] / np;
    aa_avg += v_pj[i] * v_aa[i];
    aacell_avg += v_pj[i] * v_aacell[i];
    zz_avg += v_pj[i] * v_zz[i];
  }

  // loop for qimp
  qimp = 0.;

  for (size_t i = 0; i < v_zz.size(); i++) {
    qimp += v_pj[i] * pow(v_zz[i] - zz_avg, 2.); // impurity parameter
  }

  std::cerr << "normalisation factor = " << normalisation << "\n";
  std::cerr << "sum of pj            = " << pj_sum << "\n";
  std::cerr << "<V>       = " << vws_avg << " fm^3\n";
  std::cerr << "<A>       = " << aa_avg << "\n";
  std::cerr << "<Z>       = " << zz_avg << "\n";
  std::cerr << "Qimp      = " << qimp << "\n";
  std::cerr << "nb        = " << nb << " /fm^3\n";
  std::cerr << "np        = " << np << " /fm^3\n";
  std::cerr << "ng        = " << ocp.ng << " /fm^3\n";
  std::cerr << "mun       = " << mun << " MeV\n";
  std::cerr << "mup       = " << mup << " MeV\n\n";

  std::cout << pres_ocp << "," << nb << "," << tt << "," << np << "," << ocp.ng
            << "," << mun << "," << mup << "," << vws_avg << "," << aa_avg
            << "," << aacell_avg << "," << zz_avg << "," << qimp << "," << aa_mp
            << "," << aacell_mp << "," << zz_mp << "," << ocp.aa << ","
            << aacell_ocp << "," << zz_ocp << "\n";

  return mun;
}

void PrintProbabilities(struct parameters satdata, struct sf_params sparams,
    double nb, double tt, char phase[]) {
  double       guess[4] = {100., 0.35, satdata.rhosat0, nb / 5.};
  struct compo ocp;
  ocp.aa  = guess[0];
  ocp.del = guess[1];
  ocp.n0  = guess[2];
  ocp.ng  = guess[3];

  /* tt = eval_melting_temperature(satdata, sparams, nb, &ocp, 1); */
  ocp = calc_icrust4d_composition(nb, tt, phase, guess, satdata, sparams);

  // calculate ocp pressure and chemical potentials
  double fdensws = calc_crust_ws_cell_free_energy_density(
      satdata, sparams, ocp, nb, tt, phase);
  double pres_ocp =
      calc_crust_ws_cell_pressure(satdata, sparams, ocp, nb, tt, phase);
  double np_ocp = (nb - ocp.ng) * (1. - ocp.del) / 2. / (1. - ocp.ng / ocp.n0);
  double fdensel_ocp = calc_egas_free_energy_density(np_ocp, tt);
  double presel_ocp  = calc_egas_pressure(np_ocp, tt);
  double mun_ocp     = (fdensws + pres_ocp) / nb;
  double mup_ocp     = mun_ocp - (fdensel_ocp + presel_ocp) / np_ocp;

  Probabilities(
      satdata, sparams, nb, tt, ocp, pres_ocp, np_ocp, mun_ocp, mup_ocp);
}

void InnerCrustMCPCalculation(struct parameters satdata,
    struct sf_params sparams, double tt, char phase[]) {
  double nb = 3.e-4; // initial baryon density

  double       guess[4] = {100., 0.35, satdata.rhosat0, nb / 5.};
  struct compo ocp;
  ocp.aa  = guess[0];
  ocp.del = guess[1];
  ocp.n0  = guess[2];
  ocp.ng  = guess[3];

  double fdensws = 0.;
  double pres_ocp;
  double np_ocp;
  double fdensel_ocp;
  double presel_ocp;
  double mun_ocp;
  double mup_ocp;

  while (nb < 0.04) {
    // to uncomment if we want results at the crystallization temperature
    /* tt = eval_melting_temperature(satdata, sparams, nb, &ocp, 1); */

    // calculate ocp composition
    ocp = calc_icrust4d_composition(nb, tt, phase, guess, satdata, sparams);

    if (ocp.aa != ocp.aa) {
      std::cerr << "WARNING: Cannot find the OCP solution!\n";
      break;
    }

    fdensws = calc_crust_ws_cell_free_energy_density(
        satdata, sparams, ocp, nb, tt, phase);

    // calculate ocp pressure and chemical potentials
    pres_ocp =
        calc_crust_ws_cell_pressure(satdata, sparams, ocp, nb, tt, phase);
    np_ocp      = (nb - ocp.ng) * (1. - ocp.del) / 2. / (1. - ocp.ng / ocp.n0);
    fdensel_ocp = calc_egas_free_energy_density(np_ocp, tt);
    presel_ocp  = calc_egas_pressure(np_ocp, tt);
    mun_ocp     = (fdensws + pres_ocp) / nb;
    mup_ocp     = mun_ocp - (fdensel_ocp + presel_ocp) / np_ocp;

    // mcp calculation to get average quantities
    Probabilities(
        satdata, sparams, nb, tt, ocp, pres_ocp, np_ocp, mun_ocp, mup_ocp);

    nb += nb / 10.;
  }

  return;
}
