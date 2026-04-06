#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "CalculiX.h"

/*--------------------------------------------------------------
  write_cg_sens_csv

  Purpose:
    Write stress sensitivity arrays to a single CSV file with columns:
      dCGx, dCGy, dCGz

  Parameters:
    path  - output CSV file path
    ne    - number of elements (length of each array)
    dPnorm_drho  - [ne] sensitivities of aggregated (P-norm) stress w.r.t. rho

  Returns:
    0 on success, nonzero on error (also prints a message to stderr).

  CSV format:
    dPnorm_dRho
    val_0
    val_1
    ...
----------------------------------------------------------------*/
int write_Stress_sens(const char *path,
                      size_t ne,
                      const double *dPnorm_drho)
{

    /* Open file for writing (overwrite if exists) */
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "Stress_sens.csv: cannot open '%s': %s\n",
                path, strerror(errno));
        return 2;
    }

    /* Write header */
    if (fprintf(fp, "PNORM GRADIENT\n") < 0) {
        fprintf(stderr, "Stress_sens.csv: write header failed\n");
        fclose(fp);
        return 3;
    }

    /* Write rows (use high precision; CSV uses '.' as decimal separator) */
    for (size_t i = 0; i < ne; ++i) 
    {
        if (fprintf(fp, "%.15g\n",
            dPnorm_drho[i]) < 0) 
            {
                fprintf(stderr, "write_cg_sens_csv: write row %zu failed\n", i);
                fclose(fp);
                return 4;
            }
    }

    /* Flush and close */
    if (fclose(fp) != 0) 
    {
        fprintf(stderr, "write_Stress_sens_csv: close failed for '%s': %s\n",
                path, strerror(errno));
        return 5;
    }

    return 0;
}
