#include <math.h>
#include <gsl/gsl_integration.h>

#include "nuclear_matter.h"
#include "nuclear_surface_en.h"

double calc_sly4_ldm_surface_en(double aa_)
{
    float as;

    as = 18.24; // in MeV

    return as*pow(aa_,2./3.);
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

double calc_etf_surface_en(struct parameters satdata, double aa_, double ii_, double n0_)
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

double calc_ls_surface_en(struct parameters satdata, double aa_, double ii_)
{
    double surf_energy;
    double r0;
    double sig;
    double sigs, ss;
    double q;
    double ypnuc;

    r0 = pow(4.*PI*satdata.rhosat0/3.,-1./3.);
    ypnuc = (1. - ii_)/2.;

    sigs = 1.15;                             // !!!!!!!!!!!!!!!!!!!!!
    ss = 45.8;                               // !!!! SkI' values !!!!
    q = 384.*PI*r0*r0*sigs/ss - 16.;         // !!!!!!!!!!!!!!!!!!!!!
    sig = sigs*(16. + q)/(pow(ypnuc,-3.) 
            + q + pow(1.-ypnuc,-3.));

    surf_energy = 4.*PI*r0*r0*sig*pow(aa_,2./3.);

    return surf_energy;
}
