// A sample c-script to convert native su2 mesh to CalTop format


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

    if (dir == NULL)
    {
        printf("ERROR: Cannot open current directory\n");
        exit(EXIT_FAILURE);
    }

    


}


int main (void)
{
    /* find su2 file in current directory*/


}

