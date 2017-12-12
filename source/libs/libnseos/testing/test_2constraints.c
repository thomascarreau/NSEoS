#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "../nseos/ws_cell.h"

// ==================== FUNCTIONS ====================

struct twoconstraints_fun
{
    double f_stability;
    double f_press;
};

struct twoconstraints_fun eval_twoconstraints_fun(double aa_, double del_, double rho0_, double rhop_);
struct twoconstraints_fun eval_twoconstraints_fun(double aa_, double del_, double rho0_, double rhop_)
{
    struct twoconstraints_fun result;
    struct parameters satdata;
    double enuc;
    double epsa;
    double epsr;
    double enuc_ap, enuc_am;
    double enuc_rp, enuc_rm;
    double denucdaa;
    double denucdrho0;

    satdata = assign_param(satdata);

    enuc = calc_enuc(satdata, aa_, del_, rho0_, rhop_);
    epsa = 0.001;
    epsr = 0.0001;
    enuc_ap = calc_enuc(satdata, aa_+epsa, del_, rho0_, rhop_);
    enuc_am = calc_enuc(satdata, aa_-epsa, del_, rho0_, rhop_);
    enuc_rp = calc_enuc(satdata, aa_, del_, rho0_+epsr, rhop_);
    enuc_rm = calc_enuc(satdata, aa_, del_, rho0_-epsr, rhop_);
    denucdaa = (enuc_ap - enuc_am)/2./epsa; // 2 points derivatives
    denucdrho0 = (enuc_rp - enuc_rm)/2./epsr;

    result.f_stability = denucdaa/aa_ - enuc/aa_/aa_;
    result.f_press = rho0_*rho0_*denucdrho0/aa_;

    return result;
}

struct rparams
{
    double del;
    double rhop;
};

int assign_twoconstraints_fun(const gsl_vector * x, void *params, gsl_vector * f);
int assign_twoconstraints_fun(const gsl_vector * x, void *params, gsl_vector * f)
{
    double rhop = ((struct rparams *) params)->rhop;
    double del = ((struct rparams *) params)->del;

    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);

    rhop = rhop*(1.-del)/2.;

    struct twoconstraints_fun functs;
    functs = eval_twoconstraints_fun(x0, del, x1, rhop);

    const double y0 = functs.f_stability;
    const double y1 = functs.f_press;

    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);

    return GSL_SUCCESS;
}

void print_state(gsl_multiroot_fsolver * s, double rhob_, struct rparams params_);
void print_state(gsl_multiroot_fsolver * s, double rhob_, struct rparams params_)
{
    double aa_eq;
    double rho0_eq;
    double zz;
    double rhop;

    aa_eq = gsl_vector_get(s->x,0);
    rho0_eq = gsl_vector_get(s->x, 1);
    zz = aa_eq*(1.-params_.del)/2.;
    rhop = rhob_*(1.-params_.del)/2.;

    printf ("%g %g %g %g %g %g\n", rhob_,  aa_eq, params_.del, zz, rhop/rhob_, rho0_eq);
}

// ==================== MAIN ====================

int main(void)
{
    double rhob;
    rhob = 1.e-6;

    double del;
    for(del = 0.1; del < 0.41; del += 0.01)
    {

        const gsl_multiroot_fsolver_type *T;
        gsl_multiroot_fsolver *s;

        int status;
        size_t iter = 0;

        struct rparams p;
        p.del = del;
        p.rhop = rhob; // modification in assign_ocrust_fun_3d (DIRTY!)

        double x_init[2] = {82., 0.1540};
        const size_t n = 2;
        gsl_vector *x = gsl_vector_alloc(n);
        gsl_vector_set(x, 0, x_init[0]);
        gsl_vector_set(x, 1, x_init[1]);

        gsl_multiroot_function f = {&assign_twoconstraints_fun, n, &p};

        T = gsl_multiroot_fsolver_dnewton;
        s = gsl_multiroot_fsolver_alloc(T,2);
        gsl_multiroot_fsolver_set(s, &f, x);

        do
        {
            /* print_state(s, rhob, p); */

            iter++;

            status = gsl_multiroot_fsolver_iterate(s);

            if (status)
                break;

            status = gsl_multiroot_test_residual (s->f, 9e-9);
        }

        while (status == GSL_CONTINUE && iter < 5000);

        print_state(s, rhob, p);

        gsl_multiroot_fsolver_free(s);
        gsl_vector_free(x);
    }

    return 0;
}
