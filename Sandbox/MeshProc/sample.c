/**********************************************************************
 *
 *  SU2 to CalculiX Mesh Preprocessor
 *
 *  What it does:
 *      1. Looks for the first .su2 file in the current directory
 *      2. Reads SU2 nodes and tetrahedral elements
 *      3. Writes CalculiX-compatible mesh.nam
 *      4. Extracts MARKER_TAG= fixed   -> Nfix1.nam
 *      5. Extracts MARKER_TAG= surface -> NSurface.nam
 *
 *  Compile:
 *      gcc su2_to_calculix.c -o su2_to_calculix -lm
 *
 *  Run:
 *      ./su2_to_calculix
 *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>

#define MAXLINE 512


void find_su2_file(char *filename)
{
    DIR *dir;
    struct dirent *entry;

    dir = opendir(".");

    if (dir == NULL) {
        printf("ERROR: Cannot open current directory\n");
        exit(EXIT_FAILURE);
    }

    while ((entry = readdir(dir)) != NULL) {
        if (strstr(entry->d_name, ".su2") != NULL) {
            strcpy(filename, entry->d_name);
            closedir(dir);
            return;
        }
    }

    closedir(dir);

    printf("ERROR: No .su2 mesh file found in current directory\n");
    exit(EXIT_FAILURE);
}


double tet_volume(double a[3],
                  double b[3],
                  double c[3],
                  double d[3])
{
    double bx = b[0] - a[0];
    double by = b[1] - a[1];
    double bz = b[2] - a[2];

    double cx = c[0] - a[0];
    double cy = c[1] - a[1];
    double cz = c[2] - a[2];

    double dx = d[0] - a[0];
    double dy = d[1] - a[1];
    double dz = d[2] - a[2];

    double cross_x = by * cz - bz * cy;
    double cross_y = bz * cx - bx * cz;
    double cross_z = bx * cy - by * cx;

    return (cross_x * dx +
            cross_y * dy +
            cross_z * dz) / 6.0;
}


void convert_volume_mesh(char *su2file)
{
    FILE *fp;
    FILE *out;

    char line[MAXLINE];

    int nnode = 0;
    int nelem = 0;

    printf("Reading SU2 nodes ... ");

    fp = fopen(su2file, "r");

    if (fp == NULL) {
        printf("ERROR: Cannot open SU2 file %s\n", su2file);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, MAXLINE, fp)) {
        if (strstr(line, "NPOIN=") != NULL) {
            sscanf(line, "NPOIN= %d", &nnode);
            break;
        }
    }

    if (nnode <= 0) {
        printf("ERROR: Could not read NPOIN from SU2 file\n");
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    printf("%d nodes\n", nnode);

    double (*xyz)[3];

    xyz = malloc((size_t)nnode * sizeof(*xyz));

    if (xyz == NULL) {
        printf("ERROR: Memory allocation failed for node coordinates\n");
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < nnode; i++) {
        if (fgets(line, MAXLINE, fp) == NULL) {
            printf("ERROR: Unexpected end of file while reading nodes\n");
            free(xyz);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        sscanf(line,
               "%lf %lf %lf",
               &xyz[i][0],
               &xyz[i][1],
               &xyz[i][2]);
    }

    fclose(fp);

    out = fopen("mesh.nam", "w");

    if (out == NULL) {
        printf("ERROR: Cannot open mesh.nam for writing\n");
        free(xyz);
        exit(EXIT_FAILURE);
    }

    fprintf(out, "** Generated from SU2 mesh\n");
    fprintf(out, "*NODE, NSET=NALL\n");

    for (int i = 0; i < nnode; i++) {
        fprintf(out,
                "%d, %.10lf, %.10lf, %.10lf\n",
                i + 1,
                xyz[i][0],
                xyz[i][1],
                xyz[i][2]);
    }

    fp = fopen(su2file, "r");

    if (fp == NULL) {
        printf("ERROR: Cannot reopen SU2 file %s\n", su2file);
        free(xyz);
        fclose(out);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, MAXLINE, fp)) {
        if (strstr(line, "NELEM=") != NULL) {
            sscanf(line, "NELEM= %d", &nelem);
            break;
        }
    }

    if (nelem <= 0) {
        printf("ERROR: Could not read NELEM from SU2 file\n");
        free(xyz);
        fclose(fp);
        fclose(out);
        exit(EXIT_FAILURE);
    }

    printf("Reading SU2 elements ... %d elements\n", nelem);

    fprintf(out, "\n*ELEMENT, TYPE=C3D4, ELSET=EALL\n");

    for (int e = 0; e < nelem; e++) {
        int type;
        int n[4];
        int eid;

        if (fgets(line, MAXLINE, fp) == NULL) {
            printf("ERROR: Unexpected end of file while reading elements\n");
            free(xyz);
            fclose(fp);
            fclose(out);
            exit(EXIT_FAILURE);
        }

        if (sscanf(line,
                   "%d %d %d %d %d %d",
                   &type,
                   &n[0],
                   &n[1],
                   &n[2],
                   &n[3],
                   &eid) != 6) {
            printf("ERROR: Could not parse element line %d\n", e + 1);
            free(xyz);
            fclose(fp);
            fclose(out);
            exit(EXIT_FAILURE);
        }

        if (type != 10) {
            printf("ERROR: Non-tetrahedral element detected. Element type = %d\n", type);
            free(xyz);
            fclose(fp);
            fclose(out);
            exit(EXIT_FAILURE);
        }

        double vol = tet_volume(xyz[n[0]],
                                xyz[n[1]],
                                xyz[n[2]],
                                xyz[n[3]]);

        if (fabs(vol) < 1.0e-30) {
            printf("ERROR: Zero-volume tetrahedron detected at element %d\n", e + 1);
            free(xyz);
            fclose(fp);
            fclose(out);
            exit(EXIT_FAILURE);
        }

        if (vol < 0.0) {
            int tmp = n[1];
            n[1] = n[2];
            n[2] = tmp;
        }

        fprintf(out,
                "%d, %d, %d, %d, %d\n",
                e + 1,
                n[0] + 1,
                n[1] + 1,
                n[2] + 1,
                n[3] + 1);
    }

    fclose(fp);
    fclose(out);

    free(xyz);

    printf("mesh.nam written\n");
}


void extract_marker(char *su2file,
                    char *marker,
                    char *outfile)
{
    FILE *fp;
    FILE *out;

    char line[MAXLINE];
    char key[128];

    int active = 0;
    int found_marker = 0;
    int node_count = 0;

    sprintf(key, "MARKER_TAG= %s", marker);

    fp = fopen(su2file, "r");

    if (fp == NULL) {
        printf("ERROR: Cannot open SU2 file %s\n", su2file);
        exit(EXIT_FAILURE);
    }

    out = fopen(outfile, "w");

    if (out == NULL) {
        printf("ERROR: Cannot open %s for writing\n", outfile);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    if (strcmp(marker, "fixed") == 0) {
        fprintf(out, "*NSET, NSET=Nfix1\n");
    } else {
        fprintf(out, "*NSET, NSET=Nsurface\n");
    }

    while (fgets(line, MAXLINE, fp)) {
        if (strstr(line, key) != NULL) {
            active = 1;
            found_marker = 1;
            continue;
        }

        if (active) {
            if (strstr(line, "MARKER_TAG=") != NULL) {
                break;
            }

            if (strstr(line, "MARKER_ELEMS=") != NULL) {
                continue;
            }

            int type;
            int a, b, c;

            if (sscanf(line,
                       "%d %d %d %d",
                       &type,
                       &a,
                       &b,
                       &c) == 4) {
                fprintf(out, "%d\n", a + 1);
                fprintf(out, "%d\n", b + 1);
                fprintf(out, "%d\n", c + 1);

                node_count += 3;
            }
        }
    }

    fclose(fp);
    fclose(out);

    if (!found_marker) {
        printf("WARNING: MARKER_TAG= %s not found. %s may be empty.\n",
               marker,
               outfile);
    } else {
        printf("%s written using MARKER_TAG= %s, raw node entries = %d\n",
               outfile,
               marker,
               node_count);
    }
}

/**********************************************************************
 *
 * Print all SU2 boundary markers
 *
 * Searches:
 *
 *      MARKER_TAG=
 *
 *********************************************************************/

void print_su2_markers(char *su2file)
{


    FILE *fp;


    char line[MAXLINE];

    int count = 0;



    fp = fopen(su2file,
               "r");



    if(fp == NULL){


        printf("ERROR: Cannot open %s\n",
               su2file);


        exit(EXIT_FAILURE);
    }




    printf("Detected SU2 boundary markers:\n");
    printf("--------------------------------\n");



    while(fgets(line,
                MAXLINE,
                fp)){



        if(strstr(line,
                  "MARKER_TAG=") != NULL){



            char tag[256];



            sscanf(line,
                   "MARKER_TAG= %s",
                   tag);



            printf("  %d) %s\n",
                   count+1,
                   tag);



            count++;

        }

    }



    fclose(fp);



    if(count == 0){

        printf("No MARKER_TAG entries detected\n");

    }



    printf("--------------------------------\n\n");

}



void remove_existing_nam_files(void)
{

    DIR *dir;
    struct dirent *entry;


    int count = 0;


    dir = opendir(".");


    if(dir == NULL){

        printf("ERROR: Cannot open current directory\n");

        exit(EXIT_FAILURE);
    }



    printf("Checking for old .nam files ...\n");



    while((entry = readdir(dir)) != NULL){


        /*
           Check filename extension
        */

        if(strstr(entry->d_name,
                  ".nam") != NULL){



            printf("Removing old file: %s\n",
                   entry->d_name);



            if(remove(entry->d_name) != 0){


                printf("WARNING: Could not remove %s\n",
                       entry->d_name);

            }


            count++;

        }

    }


    closedir(dir);



    if(count == 0){

        printf("No old .nam files found\n");

    }


    printf("\n");

}


void extract_skin_elements(char *su2file,
                           char *markername,
                           char *outfile)
{
    FILE *fp;
    FILE *out;

    char line[MAXLINE];
    char key[128];

    int nnode = 0;
    int nelem = 0;

    sprintf(key, "MARKER_TAG= %s", markername);

    /*
       First pass: read NPOIN and NELEM.
    */

    fp = fopen(su2file, "r");

    if (fp == NULL) {
        printf("ERROR: Cannot open %s\n", su2file);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, MAXLINE, fp)) {
        if (strstr(line, "NPOIN=") != NULL) {
            sscanf(line, "NPOIN= %d", &nnode);
        }

        if (strstr(line, "NELEM=") != NULL) {
            sscanf(line, "NELEM= %d", &nelem);
        }
    }

    fclose(fp);

    if (nnode <= 0 || nelem <= 0) {
        printf("ERROR: Could not read NPOIN or NELEM\n");
        exit(EXIT_FAILURE);
    }

    int *skin_nodes = calloc((size_t)nnode, sizeof(int));
    int *skin_elems = calloc((size_t)nelem, sizeof(int));

    if (skin_nodes == NULL || skin_elems == NULL) {
        printf("ERROR: Memory allocation failed in extract_skin_elements\n");
        exit(EXIT_FAILURE);
    }

    /*
       Second pass: collect nodes on the chosen marker.

       SU2 surface triangle format:
           5 n1 n2 n3

       The node IDs are kept 0-based here.
    */

    fp = fopen(su2file, "r");

    int active = 0;
    int found_marker = 0;

    while (fgets(line, MAXLINE, fp)) {
        if (strstr(line, key) != NULL) {
            active = 1;
            found_marker = 1;
            continue;
        }

        if (active) {
            if (strstr(line, "MARKER_TAG=") != NULL) {
                break;
            }

            if (strstr(line, "MARKER_ELEMS=") != NULL) {
                continue;
            }

            int type;
            int a, b, c;

            if (sscanf(line, "%d %d %d %d",
                       &type, &a, &b, &c) == 4) {
                skin_nodes[a] = 1;
                skin_nodes[b] = 1;
                skin_nodes[c] = 1;
            }
        }
    }

    fclose(fp);

    if (!found_marker) {
        printf("WARNING: MARKER_TAG= %s not found. No skinElementList.nam written.\n",
               markername);

        free(skin_nodes);
        free(skin_elems);

        return;
    }

    /*
       Third pass: scan volume tetrahedra.

       SU2 tetrahedron format:
           10 n1 n2 n3 n4 elementID

       If a tetrahedron contains at least one skin node,
       it is marked as a skin/passive element.
    */

    fp = fopen(su2file, "r");

    int element_section = 0;

    while (fgets(line, MAXLINE, fp)) {
        if (strstr(line, "NELEM=") != NULL) {
            element_section = 1;
            continue;
        }

        if (!element_section) {
            continue;
        }

        int type;
        int n1, n2, n3, n4;
        int eid;

        if (sscanf(line, "%d %d %d %d %d %d",
                   &type, &n1, &n2, &n3, &n4, &eid) == 6) {
            if (type != 10) {
                continue;
            }

            if (skin_nodes[n1] ||
                skin_nodes[n2] ||
                skin_nodes[n3] ||
                skin_nodes[n4]) {
                skin_elems[eid] = 1;
            }
        }
    }

    fclose(fp);

    /*
       Write selected elements as 1-based IDs.
    */

    out = fopen(outfile, "w");

    if (out == NULL) {
        printf("ERROR: Cannot open %s for writing\n", outfile);
        exit(EXIT_FAILURE);
    }

    int count = 0;

    for (int e = 0; e < nelem; e++) {
        if (skin_elems[e]) {
            fprintf(out, "%d\n", e + 1);
            count++;
        }
    }

    fclose(out);

    printf("%s written. Number of skin/passive elements = %d\n",
           outfile,
           count);

    free(skin_nodes);
    free(skin_elems);
}




int main(void)
{
    char su2file[512];

    printf("\n");
    printf("=====================================\n");
    printf("     SU2 -> CalculiX Preprocessor    \n");
    printf("=====================================\n\n");


    /* Remove old .nam files (if any) */
    remove_existing_nam_files();

    find_su2_file(su2file);

    printf("Found SU2 mesh:\n");
    printf("%s\n\n", su2file);


    /* Show all available markers */
    print_su2_markers(su2file);


    convert_volume_mesh(su2file);

    extract_marker(su2file,
                   "fixed",
                   "Nfix1.nam");

    extract_marker(su2file,
                   "surface",
                   "NSurface.nam");

    extract_skin_elements(su2file, "surface", "skinElementList.nam");




    printf("\nSU2 preprocessing complete.\n\n");

    return 0;
}