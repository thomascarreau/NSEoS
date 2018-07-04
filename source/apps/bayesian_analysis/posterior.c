#include <math.h>

#include "../../nseos/nuclear_matter.h"
#include "../../nseos/nuclear_surface_en.h"
#include "../../nseos/eos.h"
#include "../../nseos/modeling.h"

void get_low_density_posterior(FILE *prior);

int main(int argc, char* argv[])
{
    if(argc != 2)
    {
        fprintf(stderr, "ERROR: syntax is './posterior prior.in\n");
        return 1;
    }

    FILE *prior = NULL;
    prior = fopen(argv[1], "r");
    if(prior == NULL)
    {
        fprintf(stderr, "ERROR: file issue\n");
        return 1;
    }

    get_low_density_posterior(prior);

    fclose(prior);

    return 0;
}

void get_low_density_posterior(FILE *prior)
{
    struct parameters satdata;
    float m, dm;

    double nn[10] = {0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20}; 
    double e_sm_min[10] = {-4.4921, -7.6542, -10.1682, -12.1990, -13.8226, -15.0815, -16.0017, -16.6000, -16.8879, -16.8733};
    double e_sm_max[10] = {-3.8713, -6.9232, -9.2607, -11.0292, -12.3161, -13.0155, -13.2269, -13.0000, -12.3471, -11.2774};
    double p_sm_min[10] = {-0.0720, -0.2238, -0.4060, -0.5818, -0.7180, -0.7820, -0.7430, -0.5811, -0.2508, 0.2839};
    double p_sm_max[10] = {-0.0688, -0.2067, -0.3503, -0.4479, -0.4512, -0.3138, 0.0099, 0.5654, 1.3974, 2.5511};
    double e_nm_min[10] = {4.2123, 6.0424, 7.4706, 8.8545, 10.3097, 11.8712, 13.5397, 15.2000, 16.8821, 18.5931};
    double e_nm_max[10] = {4.3001, 6.2685, 7.9385, 9.6124, 11.3935, 13.3080, 15.3503, 17.5000, 19.7299, 22.0096};
    double p_nm_min[10] = {0.0463, 0.1225, 0.2476, 0.4510, 0.7375, 1.1104, 1.5754, 2.1269, 2.7534, 3.4383};
    double p_nm_max[10] = {0.0489, 0.1395, 0.2959, 0.5501, 0.9235, 1.4262, 2.0581, 2.8092, 3.6612, 4.5881};
    struct hnm test_hnm_sm;
    struct hnm test_hnm_nm;
    int ld_test;

    double p[3] = {2.5, 3.0, 3.5};
    struct sf_params sparams;
    struct transition_qtt tqtt;
    FILE *posterior = NULL;

    posterior = fopen("posterior_ld.out", "w+"); 

    while(read_table_of_sets(prior, &satdata, &m, &dm) == 0)
    {
        ld_test = 0;
        for(int i = 3; i < 10; i++)
        {
            test_hnm_sm = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, nn[i], 0.);
            test_hnm_nm = calc_meta_model_nuclear_matter(satdata, TAYLOR_EXP_ORDER, nn[i], 1.);

            if(test_hnm_sm.enpernuc < e_sm_min[i] || test_hnm_sm.enpernuc > e_sm_max[i]
                    || test_hnm_sm.p < p_sm_min[i] || test_hnm_sm.p > p_sm_max[i]
                    || test_hnm_nm.enpernuc < e_nm_min[i] || test_hnm_nm.enpernuc > e_nm_max[i]
                    || test_hnm_nm.p < p_nm_min[i] || test_hnm_nm.p > p_nm_max[i])
            {
                ld_test = 1;
                break;
            }
        }

        if(ld_test == 0)
        {
            for(int j = 0; j < 3; j++)
            {
                sparams = fit_sf_params(satdata, p[j]);
                tqtt.nt = 0.0005;
                eval_transition_qtt(satdata, sparams.p, &tqtt);
                if (tqtt.nt > 0.0005)
                    fprintf(posterior, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                            tqtt.nt, tqtt.pt, satdata.rhosat0, satdata.lasat0, satdata.ksat0,
                            satdata.qsat0, satdata.zsat0, satdata.jsym0, satdata.lsym0, 
                            satdata.ksym0, satdata.qsym0, satdata.zsym0, m, dm, satdata.b, 
                            sparams.p, sparams.sigma0, sparams.b, sparams.chi2);
            }
        }
    }

    fclose(posterior);
}
