#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "../nseos/ws_cell.h"

// ==================== FUNCTIONS ====================

struct ocrust_fun_3d
{
	double f_stability;
	double f_beta;
    double f_press;
};

struct ocrust_fun_3d calc_ocrust_fun_3d(double aa_, double del_, double rho0_, double rhop_);
struct ocrust_fun_3d calc_ocrust_fun_3d(double aa_, double del_, double rho0_, double rhop_)
{
	struct ocrust_fun_3d result;
	struct parameters satdata;
	double enuc;
	double epsa;
	double epsb;
    double epsr;
	double enuc_ap, enuc_am;
	double enuc_bp, enuc_bm;
	double enuc_rp, enuc_rm;
	double denucdaa;
	double denucddel;
    double denucdrho0;
	struct gas egas;
	double mueltot;

	satdata = assign_param(satdata);

	enuc = calc_enuc(satdata, aa_, del_, rho0_, rhop_);
	epsa = 0.001;
	epsb = 0.0001;
	epsr = 0.0001;
	enuc_ap = calc_enuc(satdata, aa_+epsa, del_, rho0_, rhop_);
	enuc_am = calc_enuc(satdata, aa_-epsa, del_, rho0_, rhop_);
	enuc_bp = calc_enuc(satdata, aa_, del_+epsb, rho0_, rhop_);
	enuc_bm = calc_enuc(satdata, aa_, del_-epsb, rho0_, rhop_);
	enuc_rp = calc_enuc(satdata, aa_, del_, rho0_+epsr, rhop_);
	enuc_rm = calc_enuc(satdata, aa_, del_, rho0_-epsr, rhop_);
	denucdaa = (enuc_ap - enuc_am)/2./epsa; // 2 points derivatives
	denucddel = (enuc_bp - enuc_bm)/2./epsb;
	denucdrho0 = (enuc_rp - enuc_rm)/2./epsr;

	egas = calc_egas(rhop_);
	mueltot = egas.mutot;

	result.f_stability = denucdaa/aa_ - enuc/aa_/aa_;
	result.f_beta = denucddel*2./aa_ - mueltot - rmp + rmn;
	result.f_press = rho0_*rho0_*denucdrho0/aa_;

	return result;
}

struct rparams
{
	double rhop;
};

int assign_ocrust_fun_3d(const gsl_vector * x, void *params, gsl_vector * f);
int assign_ocrust_fun_3d(const gsl_vector * x, void *params, gsl_vector * f)
{
	double rhop = ((struct rparams *) params)->rhop;

	const double x0 = gsl_vector_get (x, 0);
	const double x1 = gsl_vector_get (x, 1);
	const double x2 = gsl_vector_get (x, 2);

	rhop = rhop*(1.-x1)/2.;

	struct ocrust_fun_3d functs;
	functs = calc_ocrust_fun_3d(x0, x1, x2, rhop);

	const double y0 = functs.f_stability;
	const double y1 = functs.f_beta;
	const double y2 = functs.f_press;

	gsl_vector_set(f, 0, y0);
	gsl_vector_set(f, 1, y1);
	gsl_vector_set(f, 2, y2);

	return GSL_SUCCESS;
}

void print_state_ocrust(gsl_multiroot_fsolver * s, double rhob_);
void print_state_ocrust(gsl_multiroot_fsolver * s, double rhob_)
{
	double aa_eq, del_eq, rhop_eq;
	double zz_eq;
	struct parameters satdata;
	double rho0_eq;
	double epsws;

	aa_eq = gsl_vector_get(s->x, 0);
	del_eq = gsl_vector_get(s->x, 1);
	rhop_eq = rhob_*(1.-del_eq)/2.;
	zz_eq = aa_eq*(1.-del_eq)/2.;
	satdata = assign_param(satdata);
	rho0_eq = gsl_vector_get(s->x, 2);
	epsws = calc_epsws(satdata, aa_eq, del_eq, rho0_eq, 0., rhop_eq, rhob_);
	printf ("%g %g %g %g %g %g %g %g\n", aa_eq, del_eq, zz_eq, rhop_eq/rhob_, rho0_eq, epsws,
            gsl_vector_get(s->f, 0), gsl_vector_get(s->f, 1));
}

// ==================== MAIN ====================

int main(void)
{
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;

	int status;
	size_t iter = 0;

	/* double aa_sav, del_sav; */

	double rhob;
	rhob = 1.e-6;

	struct rparams p;
	p.rhop = rhob; // modification in assign_ocrust_fun_3d (DIRTY!)

	double x_init[3] = {82., 0.15, 0.1540};
	const size_t n = 3;
	gsl_vector *x = gsl_vector_alloc(n);
	gsl_vector_set(x, 0, x_init[0]);
	gsl_vector_set(x, 1, x_init[1]);
	gsl_vector_set(x, 2, x_init[2]);

	gsl_multiroot_function f = {&assign_ocrust_fun_3d, n, &p};

	T = gsl_multiroot_fsolver_dnewton;
	s = gsl_multiroot_fsolver_alloc(T,3);
	gsl_multiroot_fsolver_set(s, &f, x);

	do
	{
		print_state_ocrust(s, rhob);

		iter++;

		status = gsl_multiroot_fsolver_iterate(s);

        // bracketting
		/* aa_sav = gsl_vector_get(s->x, 0); */
		/* del_sav = gsl_vector_get(s->x, 1); */

		/* while (gsl_vector_get(s->x, 0) < 40. || gsl_vector_get(s->x, 0) > 150. */
		/* 		|| gsl_vector_get(s->x, 1) < 0. || gsl_vector_get (s->x, 1) > 0.7) */
		/* { */
		/* 	if (gsl_vector_get (s->x, 0) < 40.) { */
		/* 		gsl_vector_set (x, 0, 80.-aa_sav); */
		/* 		gsl_vector_set (x, 1, del_sav); */
		/* 		gsl_multiroot_fsolver_set (s, &f, x); */
		/* 		aa_sav = gsl_vector_get (s->x, 0); */
		/* 	} */
		/* 	else if (gsl_vector_get (s->x, 0) > 150.) { */
		/* 		gsl_vector_set (x, 0, aa_sav-150.); */
		/* 		gsl_vector_set (x, 1, del_sav); */
		/* 		gsl_multiroot_fsolver_set (s, &f, x); */
		/* 		aa_sav = gsl_vector_get (s->x, 0); */
		/* 	} */
		/* 	if (gsl_vector_get (s->x, 1) < 0.) { */
		/* 		gsl_vector_set (x, 0, aa_sav); */
		/* 		gsl_vector_set (x, 1, -del_sav); */
		/* 		gsl_multiroot_fsolver_set (s, &f, x); */
		/* 		del_sav = gsl_vector_get (s->x, 1); */
		/* 	} */
		/* 	else if (gsl_vector_get (s->x, 1) > 0.7) { */
		/* 		gsl_vector_set (x, 0, aa_sav); */
		/* 		gsl_vector_set (x, 1, del_sav-0.7); */
		/* 		gsl_multiroot_fsolver_set (s, &f, x); */
		/* 		del_sav = gsl_vector_get (s->x, 1); */
		/* 	} */
		/* } */

		if (status)
			break;

		status = gsl_multiroot_test_residual (s->f, 9e-9);
	}

	while (status == GSL_CONTINUE && iter < 5000);

	print_state_ocrust(s, rhob);

	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);

	return 0;
}
