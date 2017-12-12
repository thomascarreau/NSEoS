#include <stdio.h>

#include "../nseos/nuclear_matter.h"

int main(void)
{
    struct parameters sly4;
    int inn;
    double nn;
    double ii;
    double elfc2_enpernuc;
    double elfc3_enpernuc;
    double elfc4_enpernuc;
    double sly4_enpernuc;

    sly4 = assign_param(sly4);
    ii = 1.;

    for(inn = 0; inn < 61; inn++)
    {
        nn = inn/100.;
        elfc2_enpernuc = calc_meta_model_nuclear_matter_enpernuc(sly4,2,nn,ii);
        elfc3_enpernuc = calc_meta_model_nuclear_matter_enpernuc(sly4,3,nn,ii);
        elfc4_enpernuc = calc_meta_model_nuclear_matter_enpernuc(sly4,4,nn,ii);
        sly4_enpernuc = calc_sly4_nuclear_matter_enpernuc(nn,ii);
        printf("%g %g %g %g %g\n", nn, sly4_enpernuc, elfc2_enpernuc,
                elfc3_enpernuc, elfc4_enpernuc);
    }

    return 0;
}
