#ifndef SU2_PREPROCESS_H
#define SU2_PREPROCESS_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>

#define MAXLINE 1024

void find_su2_file(char *filename);


void remove_existing_nam_files(void);


void print_su2_markers(char *su2file);


double tet_volume(double a[3],
                  double b[3],
                  double c[3],
                  double d[3]);


void convert_volume_mesh(char *su2file);

void extract_marker(char *su2file,
                    char *marker,
                    char *outfile);

void extract_skin_elements(char *su2file,
                           char *markername,
                           char *outfile);

void write_deformed_su2_from_v(int nk,
                               double *v);


                               

#endif