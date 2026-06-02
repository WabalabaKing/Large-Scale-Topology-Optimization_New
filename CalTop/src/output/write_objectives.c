#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>

/**
 * Writes volume sensitivities to a CSV file.
 *
 * @param ne                Number of elements.
 * @param eleVol            Array of original (geometric) element volumes.
 * @param rhoPhys           Array of filtered element densities.
 * @param mass              Structrual mass.
 * @param cgx               C.G x coordinate
 * @param cgy               C.G y coordinate
 * @param cgz               C.G z coordinate
 * @param Pnorm             Pnorm stress
 */
void write_objectives(int ne,
                                const double *eleVol,
                                const double *rhoPhys,
                                const double * compliance_sum,
                                const double *Mass, 
                                const double *cgx, 
                                const double *cgy, 
                                const double *cgz,
                                const int *passiveIDs,   // ***skin element IDs (from passiveElements)
                                int numPassive,       
                                const double *pnorm)
{
    const char *filename = "objectives.csv";

    /* total material volume with rho = 1 */
    double initialVol_sum = 0;

    /* initial design volume of active domain */
    double initialDesignVol_sum = 0.0;

    /* total material volume with optimized rho */
    double  designVol_sum = 0;

    /* design domain discreteness */
    double discreteness_sum = 0.0;

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

    /* NOTE:
        ORIGINAL VOLUME := Volume of the structure with rho = 1 
        DESIGN VOLUME   := Volume of the structure with filtered rho
    */


    /* Pre-allocate skin element indices here*/
    int *is_skin = calloc(ne, sizeof(int));
    if (!is_skin) 
    {
        perror("calloc failed for is_skin");
        exit(EXIT_FAILURE);
    }
    
    /* passiveIDs are 1-based */
    for (int k =0; k < numPassive; k++)
    {
        int id = passiveIDs[k] -1;
        if (id >= 0 && id < ne) 
        {
            is_skin[id] = 1;
        }
    }

    int nthreads = omp_get_max_threads();

    printf("Using %d CPU(s) for volume and discreteness evaluation\n",
       nthreads);
    fflush(stdout);



    /* Write file header */
    fprintf(obj_file, "COMPLIANCE, ORIGINAL VOLUME, DESIGN VOLUME, VOLUME_FRACTION, DISCRETENESS, MASS, CGx, CGy, CGz, PNORM\n");

    #pragma omp parallel for reduction(+:initialVol_sum,initialDesignVol_sum,designVol_sum,discreteness_sum)
    for (int i = 0; i < ne; i++)
    {
        initialVol_sum += eleVol[i];
        discreteness_sum += rhoPhys[i] * (1.0 - rhoPhys[i]);

        if (!is_skin[i])
        {   initialDesignVol_sum+= eleVol[i];
            designVol_sum += eleVol[i] * rhoPhys[i];
        }
    }

    /* Compute volume fraction for active domain only */
    double volume_fraction = 0.0;
    
    if (initialDesignVol_sum > 0.0)
    {
        volume_fraction = designVol_sum / initialDesignVol_sum;
    }

    /* Compute discreteness */
    double discreteness = (4.0 /ne) * discreteness_sum;

    /* Write structure compliance and volume to file */
    fprintf(obj_file, "%.15f, %.15f, %.15f, %.15f, %.15f, %.15f, %.15f, %.15f, %.15f, %.15f \n", *compliance_sum, initialVol_sum, designVol_sum, volume_fraction, discreteness, *Mass, *cgx, *cgy, *cgz, *pnorm);
    
    fclose(obj_file);

    free(is_skin);
}
