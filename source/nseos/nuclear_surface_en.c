#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "nuclear_matter.h"
#include "nuclear_surface_en.h"
#include "coulomb.h"
#include "modeling.h"

double calc_ldm_surface_en(struct parameters satdata, double aa_)
{
    double rsat;
    float sigmas;

    rsat = pow(3./4./PI/satdata.rhosat0,1./3.); // SNM saturation density reference value
    sigmas = 1.02999; // in MeV/fm^2; SLy4 value

    return 4.*PI*rsat*rsat*sigmas*pow(aa_,2./3.);
}

double my_integrand(double x, void *params_ptr)
{
    struct f_params * params =
        (struct f_params *) params_ptr;
    int k = (params->k);
    double gamma = (params->gamma);
    double my_integrand =  pow(-1,k)*(( 1.+pow(-1,k)*exp(-gamma*x) )/( pow(1+exp(-x),gamma) ) - 1.)*pow(x,k);

    return my_integrand;
}

double eta_function(int a, double b)
{
    gsl_integration_workspace *work_ptr =
        gsl_integration_workspace_alloc (1000);

    double lower_limit = 0.;      /* start integral from 0 (to infinity) */
    double abs_error = 1.0e-8;    /* to avoid round-off problems */
    double rel_error = 1.0e-8;    /* the result will usually be much better */
    double result;                /* the result from the integration */
    double error;                 /* the estimated error from the integration */
    struct f_params params = {a, b};

    gsl_function My_function;
    My_function.function = &my_integrand;
    My_function.params = &params;
    gsl_integration_qagiu (&My_function, lower_limit,
            abs_error, rel_error, 1000, work_ptr, &result,
            &error);

    return result;
}

double calc_etf_ana_surface_en(struct parameters satdata, double aa_, double ii_, double n0_)
{
    double ckin;
    double alphanum, alphadenom, alpha;
    double c0num, c0denom, c0;
    double c3num, c3denom, c3;
    double xsat, rmsat, delmsat;
    struct hnm hnm12;
    double clsurf, clcurv, clind, clsurf0;
    int i, imax;
    double sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8;
    double cfin, vso;
    double cnlsurf, cnlcurv, cnlind, cnlsurf0;
    double zz;
    double rhs, rhosatp, rhsp, drhs;
    double aslab;
    double ais, agaus2, agaus;
    double r0;
    double adif, aratio;
    double drskin;
    double esurfissurf, esurfiscurv, esurfisind, esurfis;
    double ratioskin;
    double esurfivcoeff, esurfivsurf, esurfivind, esurfiv;
    double esurf;

    r0 = pow(3./4./PI/n0_,1./3.);

    ckin = 3.*HBARC*HBARC/10./RMN*pow((3.*PI2/2.),2./3.) ; // in MeV.fm^2

    alphanum = -(satdata.ksat0 + 9.*satdata.lasat0 - ckin*pow(satdata.rhosat0,2./3.)*(1. + 4.*satdata.barm));
    alphadenom =  9.*(satdata.lasat0 - ckin*pow(satdata.rhosat0,2./3.)*((1.+satdata.barm)/3. - satdata.barm));
    alpha = (alphanum/alphadenom);

    c0num = satdata.lasat0*satdata.ksat0 - satdata.ksat0*ckin*pow(satdata.rhosat0,2./3.)*(1.+satdata.barm)
        - satdata.lasat0*ckin*pow(satdata.rhosat0,2./3.)*(4.*(1.+satdata.barm) + 21.*satdata.barm) 
        + 9.*ckin*ckin*pow(satdata.rhosat0,4./3.)*satdata.barm;
    c0denom = satdata.rhosat0*(satdata.ksat0 + 9.*satdata.lasat0
            - ckin*pow(satdata.rhosat0,2./3.)*(1. + 4.*satdata.barm));
    c0 = c0num/c0denom;

    c3num = 9.*pow((satdata.lasat0 - ckin*pow(satdata.rhosat0,2./3.)*((1.+satdata.barm)/3. - satdata.barm)),2.);
    c3denom = pow(satdata.rhosat0,alpha)*c0denom;
    c3 = c3num/c3denom;

    xsat = (n0_ - satdata.rhosat0)/3./satdata.rhosat0;
    rmsat = RMN/(1. + satdata.barm*(1. + 3.*xsat));
    delmsat = (RMN - rmsat)/rmsat;

    hnm12 = calc_meta_model_nuclear_matter(satdata, 4, satdata.rhosat0/2., ii_);

    clsurf = 3.*ckin*pow(n0_,2./3.)*(eta_function(0,5./3.)*RMN/rmsat - 3./5.*delmsat) 
        - 3.*c0*n0_ + 3.*c3*pow(n0_,alpha+1)*eta_function(0,alpha+2) ; 

    clcurv = 6.*ckin*pow(n0_,2./3.)*((eta_function(1,5./3.) - PI2/6.)*RMN/rmsat 
            - 3./5.*eta_function(0,5./3.)*delmsat)
        +  6.*c3*pow(n0_,alpha+1)*(eta_function(1,alpha+2) - PI2/6.) ; 

    clind = 3.*ckin*pow(n0_,2./3.)*( (eta_function(2,5./3.) 
                - 2.*PI2/3.*eta_function(0,5./3.))*RMN/rmsat 
            - 2./5.*(3.*eta_function(1,5./3.) - PI2)*delmsat )
        + PI2*c0*n0_ + 3.*c3*pow(n0_,alpha+1)*(eta_function(2,alpha+2) 
                - 2.*PI2/3.*eta_function(0,alpha+2)) ;

    clsurf0 = 3.*ckin*pow(satdata.rhosat0,2./3.)*(eta_function(0,5./3.)*(1.+satdata.barm) - 3./5.*satdata.barm)
        - 3.*c0*satdata.rhosat0 + 3.*c3*pow(satdata.rhosat0,alpha+1)*eta_function(0,alpha+2) ;

    imax = 7;
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    sum6 = 0.0;
    sum7 = 0.0;
    sum8 = 0.0;
    for( i = 0; i <= imax; i = i + 1 ){
        sum1 = sum1 + pow(-1,i)*pow(delmsat,i+2)/(i+3)/(i+4);
        sum2 = sum2 + pow(-1,i)*pow(delmsat,i)/(i+3)/(i+4);
        sum3 = sum3 + pow(-1,i)*pow(delmsat,i+2)/(i+3)/(i+4)*(eta_function(0,i+2) + 1.);
        sum4 = sum4 + pow(-1,i)*pow(delmsat,i)/(i+3)/(i+4)*(eta_function(0,i+3) + 1.);
        sum5 = sum5 + pow(-1,i)*pow(delmsat,i+2)/(i+3)/(i+4)*(eta_function(1,i+2) + eta_function(0,i+2)-PI2/3.);
        sum6 = sum6 + pow(-1,i)*pow(delmsat,i)/(i+3)/(i+4)*(eta_function(1,i+3) + eta_function(0,i+3)-PI2/3.);
        sum7 = sum1 + pow(-1,i)*pow(satdata.barm,i+2)/(i+3)/(i+4);
        sum8 = sum2 + pow(-1,i)*pow(satdata.barm,i)/(i+3)/(i+4);
    }

    cfin = 59.; // reference value; see: arXiv:1709.00189
    vso = 0.; // from chi2 minimization

    cnlsurf =  HBARC*HBARC/4./RMN*(1./12. - 11./36.*delmsat - 1./2.*sum1)
        + 1./2.*cfin*n0_ + 3.*vso*n0_*n0_*sum2 ;
    cnlcurv = HBARC*HBARC/2./RMN*(1./12. - 1./2.*sum3)
        + 6.*vso*n0_*n0_*sum4 ;
    cnlind = HBARC*HBARC/2./RMN*(-PI2/12./6. + 11./36.*(1. + PI2/6.)*delmsat)
        - HBARC*HBARC/4./RMN*sum5 - (1.+PI2/6.)*cfin*n0_ 
        + 6.*vso*n0_*n0_*sum6 ;
    cnlsurf0 =  HBARC*HBARC/4./RMN*(1./12. - 11./36.*satdata.barm - 1./2.*sum7)
        + 1./2.*cfin*satdata.rhosat0 + 3.*vso*satdata.rhosat0*satdata.rhosat0*sum8 ;

    rhs = pow(aa_/(4./3.*PI*n0_),1./3.);
    rhosatp = n0_*(1.- ii_)/2.;
    zz = aa_*(1.-ii_)/2.;
    rhsp = pow(zz/(4./3.*PI*rhosatp),1./3.);
    drhs = rhs - rhsp;

    aslab = pow((cnlsurf0/clsurf0),1./2.);

    ais = pow((cnlsurf/clsurf),1./2.); // in fm
    agaus2 = ais*ais + pow((PI/(1.-hnm12.ksym/18./hnm12.jsym)),1./2.)*satdata.rhosat0/n0_
        *3.*hnm12.jsym*(ii_ - ii_*ii_)/clsurf*aslab*drhs; 
    agaus = sqrt(agaus2);
    adif = agaus ;
    aratio = adif/r0;

    drskin = drhs*(1. + PI2/3.*adif*adif/rhs/rhsp); // n skin thickness
    ratioskin = drskin/adif ;

    esurfissurf = (clsurf + cnlsurf/adif/adif)*aratio*pow(aa_,2./3.) ;
    esurfiscurv = (clcurv + cnlcurv/adif/adif)*pow(aratio,2.)*pow(aa_,1./3.) ;
    esurfisind = (clind + cnlind/adif/adif)*pow(aratio,3.) ;
    esurfis = esurfissurf + esurfiscurv + esurfisind ; 

    esurfivcoeff = 3.*pow(PI/(1.-hnm12.ksym/18./hnm12.jsym),1./2.)*satdata.rhosat0/n0_*aslab/r0*hnm12.jsym
        *(pow(ratioskin,2.)/4. + (ratioskin - pow(ratioskin,2.)/2.)*ii_
                + ((1. - satdata.jsym0/hnm12.jsym)-ratioskin
                    -(1.+satdata.lsym0*hnm12.lsym/hnm12.jsym/satdata.ksat0)/4.*pow(ratioskin,2.))*pow(ii_,2.));
    esurfivsurf = esurfivcoeff*pow(aa_,2./3.) ; 
    esurfivind = esurfivcoeff*2./(1. - hnm12.ksym/18./hnm12.jsym)*pow(aslab/r0,2.) - 2.*PI2/3.*pow(aratio,2.) ;
    esurfiv = esurfivsurf + esurfivind ;

    esurf = esurfis + esurfiv;

    return esurf;
}

#define N (9)

int be_f (const gsl_vector * x, void *data, gsl_vector * f)
{
    int *zz = ((struct data *)data)->zz;
    int *aa = ((struct data *)data)->aa;
    double *be = ((struct data *)data)->be;
    struct parameters satdata = ((struct data *)data)->satdata;

    struct sf_params prms;
    prms.sigma0 = gsl_vector_get (x, 0);
    prms.b = gsl_vector_get (x, 1);

    size_t i;

    for (i = 0; i < N; i++) // (wc -l spherical_nuclei.data)
    {
        double ii = 1.-2.*zz[i]/aa[i];
        double n0 = satdata.rhosat0*(1.-3.*satdata.lsym0*ii*ii/(satdata.ksat0 + satdata.ksym0*ii*ii));
        // bulk
        struct hnm meta;
        meta = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, n0, ii);
        double ebulk = meta.enpernuc;
        // surface
        double r0 = pow(4.*PI*n0/3.,-1./3.);
        double ypnuc = (1.-ii)/2.;
        double sigma = prms.sigma0*(pow(2.,P_SURF_TENSION+1.) + prms.b)/(pow(ypnuc,-P_SURF_TENSION) 
                + prms.b + pow(1.-ypnuc,-P_SURF_TENSION));
        double esurf = 4.*PI*r0*r0*sigma*pow(aa[i],-1./3.);
        // coulomb
        double Ecoul = calc_coulomb_en(satdata, aa[i], ii, n0, 0.);
        double ecoul = Ecoul/aa[i];
        // total
        double Be = ebulk + esurf + ecoul;
        gsl_vector_set (f, i, Be - be[i]);
    }

    return GSL_SUCCESS;
}

struct sf_params fit_sf_params(struct parameters satdata)
{
    struct sf_params prms;

    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver *s;
    int status, info;
    size_t i;
    const size_t n = N;
    const size_t nb_of_params = 2;

    gsl_matrix *J = gsl_matrix_alloc(n, nb_of_params);
    gsl_matrix *covar = gsl_matrix_alloc (nb_of_params, nb_of_params);
    int zz[n], aa[n];
    double be[n], weights[n];
    struct data d = { zz, aa, be, satdata };
    gsl_multifit_function_fdf f;
    double x_init[2] = { 1.08, 23.0 };
    gsl_vector_view x = gsl_vector_view_array (x_init, nb_of_params);
    gsl_vector_view w = gsl_vector_view_array(weights, n);
    gsl_vector *res_f;
    double chi, chi0;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 0.0;

    f.f = &be_f;
    f.df = NULL; // finite differences
    f.n = n;
    f.p = nb_of_params;
    f.params = &d;

    /* This is the data to be fitted */
    FILE *data = NULL;
    /* data = fopen("../../input/masses/AME2012.data", "r"); */
    data = fopen("../../input/masses/spherical_nuclei.data", "r");
    if (data != NULL)
    {
        for (i = 0; i < n; i++)
        {
            fscanf(data, "%d %d %lf", &zz[i], &aa[i], &be[i]);
            weights[i] = 1.;
        }
        fclose(data);
    } else {
        fprintf(stderr,"ERROR: file issue\n");
        exit(0);
    }

    s = gsl_multifit_fdfsolver_alloc (T, n, nb_of_params);

    /* initialize solver with starting point and weights */
    /* gsl_multifit_fdfsolver_wset (s, &f, &x.vector, &w.vector); */
    gsl_multifit_fdfsolver_wset (s, &f, &x.vector, &w.vector);

    /* compute initial residual norm */
    res_f = gsl_multifit_fdfsolver_residual(s);
    chi0 = gsl_blas_dnrm2(res_f);

    /* solve the system with a maximum of 20 iterations */
    status = gsl_multifit_fdfsolver_driver(s, 20, xtol, gtol, ftol, &info);

    gsl_multifit_fdfsolver_jac(s, J);
    gsl_multifit_covar (J, 0.0, covar);

    /* compute final residual norm */
    chi = gsl_blas_dnrm2(res_f);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

    fprintf(stderr, "summary from method '%s'\n",
            gsl_multifit_fdfsolver_name(s));
    fprintf(stderr, "number of iterations: %zu\n",
            gsl_multifit_fdfsolver_niter(s));
    fprintf(stderr, "function evaluations: %zu\n", f.nevalf);
    fprintf(stderr, "Jacobian evaluations: %zu\n", f.nevaldf);
    fprintf(stderr, "reason for stopping: %s\n",
            (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "initial |f(x)| = %g\n", chi0);
    fprintf(stderr, "final   |f(x)| = %g\n", chi);

    { 
        double dof = n - nb_of_params;
        double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 

        fprintf(stderr, "chisq/dof = %g\n",  pow(chi, 2.0) / dof);

        fprintf (stderr, "sigma0 = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
        fprintf (stderr, "b      = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
    }

    fprintf (stderr, "status = %s\n\n", gsl_strerror (status));
    fprintf (stderr, "==============================================\n\n");

    prms.sigma0 = gsl_vector_get(s->x, 0);
    prms.b = gsl_vector_get(s->x, 1);

    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);
    gsl_matrix_free (J);

    return prms;
}

double calc_ls_surface_en(struct sf_params sparams, double aa_, double ii_, double n0_)
{
    double surf_energy;
    double r0;
    double sigma;
    double ypnuc;

    r0 = pow(4.*PI*n0_/3.,-1./3.);
    ypnuc = (1. - ii_)/2.;
    sigma = sparams.sigma0*(pow(2.,P_SURF_TENSION+1.) + sparams.b)/(pow(ypnuc,-P_SURF_TENSION) 
            + sparams.b + pow(1.-ypnuc,-P_SURF_TENSION));
    surf_energy = 4.*PI*r0*r0*sigma*pow(aa_,2./3.);

    return surf_energy;
}
