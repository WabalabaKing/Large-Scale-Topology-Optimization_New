#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/**
 * Writes volume sensitivities to a CSV file.
 *
 * @param ne                Number of elements.
 * @param eleVol            Array of original (geometric) element volumes.
 * @param rhoPhys           Array of element densities updated by the optimizer.
 */
void write_objectives(int ne,
                                const double *eleVol,
                                const double *rhoPhys,
                                const double * compliance_sum)
{
    const char *filename = "objectives.csv";

    /* total material volume with rho = 1 */
    double initialVol_sum = 0;

    /* total material volume with optimized rho */
    double  designVol_sum = 0;

    /* Check for old sensitivity file and delete it */
    if (access(filename, F_OK) == 0)
    {
        if (remove(filename) != 0)
        {
            perror("Error deleting existing objectives file");
            exit(EXIT_FAILURE);
        }
    }

    /* Open new sensitivity file */
    FILE *obj_file = fopen(filename, "w");
    if (!obj_file) {
        perror("Error opening objectives.csv");
        exit(EXIT_FAILURE);
    }

    /* Write file header */
    fprintf(obj_file, "COMPLIANCE, INITIAL VOLUME, DESIGN VOLUME\n");

    /* Loop over all elements and compute the initial and current volume*/
    for (int i = 0; i < ne; i++)
    {
        initialVol_sum+= eleVol[i];
        designVol_sum+= (eleVol[i]*rhoPhys[i]);
        //fprintf(sens_file, "%.15f,%.15f,%.15f\n", eleVol[i], eleVol[i] * rhoPhys[i], eleVolFiltered[i]);
    }

    /* Write structure compliance and volume to file */
    fprintf(obj_file, "%.15f, %.15f, %.15f \n", *compliance_sum, initialVol_sum, designVol_sum);
    
    fclose(obj_file);
}
