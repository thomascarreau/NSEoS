#include <math.h>

#include "../../nseos/crust.h"
#include "../../nseos/modeling.h"
#include "../../nseos/phyconst.h"

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "ERROR: syntax is './gs_ocrust mass_table.data outfile'\n");
    return 1;
  }
  struct gs_ocrust gs;
  double           nb = 1.0e-9;

  gs.mun = 0.0;

  FILE *outfile = fopen(argv[2], "w+");
  fprintf(outfile, "nb,aa,zz,pres,mun,mup\n");

  while (gs.mun < RMN && nb < 3.0e-4) {
    gs = calc_ocrust_composition_with_mass_table(nb, 0.0, "sol", argv[1]);
    fprintf(outfile, "%g,%d,%d,%g,%g,%g\n", nb, gs.aa, gs.zz, gs.pres, gs.mun,
        gs.mup);
    nb += nb / 50.0;
  }

  fclose(outfile);

  return 0;
}
