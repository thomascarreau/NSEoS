#include "f_icrust4d.h"

int main(void)
{
    int irhob;
    double rhob;
    double guess[4] = {100.,0.25,0.15,1.e-4};

    for(irhob = 3; irhob < 1001; irhob ++)
    {
        rhob = irhob/10000.;
        calc_icrust4d(rhob, guess);
        if (guess[0] != guess[0])
            break;
    }

    return 0;
}
