
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* Function to open a file and handle errors 
FILE *open_file(const char *filename, const char *mode) 
{
/    FILE *file = fopen(filename, mode);
    if (!file) {
        perror(filename);
        exit(EXIT_FAILURE);
    }
    return file;
}
*/
/**
 * Writes compliance sensitivities to a csv file and returns the 
 * total compliance of the structure
 * @param ne                Number of elements
 * @param eleVol            Array of original (geometric) element volume
 * @param rhoPhys           Array of element densities updated by the optimizer
 * @param eleVolFiltered    Array of element volume sensitivities filtered
 */

 void write_volume_sensitivities(int ne,
                                     const double *eleVol,
                                     const double *rhoPhys,
                                     const double *eleVolFiltered)
{
    const char *filename = "volume_sens.csv";

    /* Check for old sensitivity file and delete it */
    if (access(filename, F_OK) == 0)
    {
        if (remove(filename) != 0)
        {
            perror("Error deleting existing compliance sensitibity file! \n");
            exit(EXIT_FAILURE);
        }
    }

    /* Old file purged, write new sensitivity file */
    FILE *sens_file = open_file(filename, "w");

    /* Write file header */
    fprintf(sens_file, "Element volume, Current element volume, volume GRADIENT\n");
    

    /* Loop over all elements and write their sensitivities to file */
    for (int i = 0; i < ne; i++)
    {
        fprintf(sens_file, "%.15f,%.15f\n", eleVol[i], eleVol[i]*rhoPhys[i], eleVolFiltered[i]);
    }

    fclose(sens_file);

}