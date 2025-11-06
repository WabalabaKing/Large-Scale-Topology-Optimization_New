#include <stdio.h>
#include <stdlib.h>

/**
 * Reads skinElementList.nam to and returns the element indices.
 *
 * @param filename The name of the file (e.g., "skinElementList.nam").
 * @param numPassive Pointer to an int where the number of elements will be stored.
 * @return Pointer to an array of integers with passive element IDs.
 */
int* passiveElements(const char *filename, int *numPassive) 
{
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Could not open %s\n", filename);
        *numPassive = 0;
        return NULL;
    }

    /* Set some initial capacity */
    int capacity = 1000;
    int *elements = (int *)malloc(capacity * sizeof(int));
    if (!elements) 
    {
        printf("Memory allocation failed\n");
        fclose(file);
        *numPassive = 0;
        return NULL;
    }

    int id;
    *numPassive = 0;
    /* Start reading the file */
    while (fscanf(file, "%d", &id) == 1) 
    {
        if (*numPassive >= capacity) 
        {
            capacity *= 2;
            elements = (int *)realloc(elements, capacity * sizeof(int));
            if (!elements) 
            {
                printf("Memory reallocation failed\n");
                fclose(file);
                *numPassive = 0;
                return NULL;
            }
        }
        elements[(*numPassive)++] = id;
    }

    fclose(file);
    return elements;
}
