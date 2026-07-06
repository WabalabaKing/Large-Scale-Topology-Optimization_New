
#include "su2_preprocess.h"

int main(void)
{
    char su2file[512];

    printf("\n");
    printf("=====================================\n");
    printf("     SU2 -> CalculiX Preprocessor    \n");
    printf("=====================================\n\n");


    /* Remove old .nam files (if any) */
    remove_existing_nam_files();

    /* Find su2 mesh in CWD */
    find_su2_file(su2file);

    printf("Found SU2 mesh:\n");
    printf("%s\n\n", su2file);


    /* Show all available markers */
    print_su2_markers(su2file);

    /* Get mesh.nam */
    convert_volume_mesh(su2file);

    /* Extract .nam for fixed nodes*/
    extract_marker(su2file,
                   "fixed",
                   "Nfix1.nam");
    
    /* Extract .nam for traction nodes */
    extract_marker(su2file,
                   "surface",
                   "NSurface.nam");

    /* Extract skin elements */
    extract_skin_elements(su2file, "surface", "skinElementList.nam");

    extract_skin_elements(su2file, "tank", "tankElementList.nam");

    printf("\nSU2 preprocessing complete.\n\n");

    return 0;
}