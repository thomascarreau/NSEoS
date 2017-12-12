#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "../nseos/ws_cell.h"

// ==================== FUNCTIONS ====================

double eval_press_eq_fun(double aa_, double del_, double rho0_, double rhop_);
double eval_press_eq_fun(double aa_, double del_, double rho0_, double rhop_)
{
    struct parameters satdata;
    double epsr;
    double enuc_rp, enuc_rm;
    double denucdrho0;
    double press_eq_fun;

    satdata = assign_param(satdata);
    epsr = 0.0001;
    enuc_rp = calc_enuc(satdata, aa_, del_, rho0_+epsr, rhop_);
    enuc_rm = calc_enuc(satdata, aa_, del_, rho0_-epsr, rhop_);
    denucdrho0 = (enuc_rp - enuc_rm)/2./epsr;
    press_eq_fun = rho0_*rho0_*denucdrho0/aa_;

    return press_eq_fun;
}

struct rparams
{
    double aa;
    double del;
    double rhop;
};

int assign_press_eq_fun(const gsl_vector * x, void *params, gsl_vector * f);
int assign_press_eq_fun(const gsl_vector * x, void *params, gsl_vector * f)
{
    double aa = ((struct rparams *) params)->aa;
    double del = ((struct rparams *) params)->del;
    double rhop = ((struct rparams *) params)->rhop;

    const double x0 = gsl_vector_get (x, 0);

    rhop = rhop*(1.-del)/2.;

    double press_eq_fun;
    press_eq_fun = eval_press_eq_fun(aa, del, x0, rhop);

    const double y0 = press_eq_fun;

    gsl_vector_set(f, 0, y0);

    return GSL_SUCCESS;
}

void print_state(gsl_multiroot_fsolver * s, double rhob_, struct rparams params_);
void print_state(gsl_multiroot_fsolver * s, double rhob_, struct rparams params_)
{
    double rho0_eq;
    double zz;
    double rhop;

    rho0_eq = gsl_vector_get(s->x, 0);
    zz = params_.aa*(1.-params_.del)/2.;
    rhop = rhob_*(1.-params_.del)/2.;

    printf ("%g %g %g %g %g %g\n", rhob_,  params_.aa, params_.del, zz, rhop/rhob_, rho0_eq);
}

// ==================== MAIN ====================

int main(void)
{
    double rhob;
    rhob = 1.e-6;
    double del;

    for (del = 0.; del <= 0.4; del += 0.01)
    {

        const gsl_multiroot_fsolver_type *T;
        gsl_multiroot_fsolver *s;

        int status;
        size_t iter = 0;

        struct rparams p;
        p.rhop = rhob; // modification in assign_ocrust_fun_2d (DIRTY!)
        p.aa = 100.;
        p.del = del;

        double x_init = 0.1540;
        const size_t n = 1;
        gsl_vector *x = gsl_vector_alloc(n);
        gsl_vector_set(x, 0, x_init);

        gsl_multiroot_function f = {&assign_press_eq_fun, n, &p};

        T = gsl_multiroot_fsolver_dnewton;
        s = gsl_multiroot_fsolver_alloc(T,1);
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
