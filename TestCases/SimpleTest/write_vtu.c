#include <stdio.h>
#include <stdlib.h>

void writeVTU(const char *filename, int nk, int ne, double *co, int *kon, int *ipkon) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("Error opening file!\n");
        return;
    }

    // Write XML and VTK headers
    fprintf(file, "<?xml version=\"1.0\"?>\n");
    fprintf(file, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(file, "  <UnstructuredGrid>\n");
    fprintf(file, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nk, ne);

    // Write Points
    fprintf(file, "      <Points>\n");
    fprintf(file, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int i = 0; i < nk; i++) {
        fprintf(file, "        %f %f %f\n", co[3 * i], co[3 * i + 1], co[3 * i + 2]);
    }
    fprintf(file, "        </DataArray>\n");
    fprintf(file, "      </Points>\n");

    // Write Cells
    fprintf(file, "      <Cells>\n");

    // Connectivity
    fprintf(file, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    for (int e = 0; e < ne; e++) {
        int startIdx = ipkon[e];  // Start index in kon array
        for (int j = 0; j < 4; j++) { // Tetrahedron has 4 nodes
            fprintf(file, "        %d ", kon[startIdx + j] - 1); // Convert 1-based to 0-based index
        }
        fprintf(file, "\n");
    }
    fprintf(file, "        </DataArray>\n");

    // Offsets
    fprintf(file, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    for (int e = 0; e < ne; e++) {
        fprintf(file, "        %d\n", (e + 1) * 4); // Offsets increase by 4 for tetrahedra
    }
    fprintf(file, "        </DataArray>\n");

    // Types
    fprintf(file, "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n");
    for (int e = 0; e < ne; e++) {
        fprintf(file, "        10\n"); // VTK_TETRA = 10
    }
    fprintf(file, "        </DataArray>\n");

    fprintf(file, "      </Cells>\n");

    // Close XML structure
    fprintf(file, "    </Piece>\n");
    fprintf(file, "  </UnstructuredGrid>\n");
    fprintf(file, "</VTKFile>\n");

    fclose(file);
    printf("VTU file written successfully: %s\n", filename);
}

int main() {
    // Example: Two tetrahedra sharing three points
    int nk = 5;  // Number of nodes
    int ne = 2;  // Number of tetrahedral elements

    // Node coordinates
    double co[] = {
        0.0, 0.0, 0.0,  // Node 1
        1.0, 0.0, 0.0,  // Node 2
        0.5, 1.0, 0.0,  // Node 3
        0.5, 0.5, 1.0,  // Node 4
        0.5, 0.5, -1.0  // Node 5
    };

    // Tetrahedral connectivity (1-based indexing from CalculiX)
    int kon[] = {
        0, 1, 2, 3,  // First tetrahedron
        0, 1, 2, 4   // Second tetrahedron (shares three nodes)
    };

    // Start index of each element in kon
    int ipkon[] = {0, 4}; 

    // Write the VTU file
    writeVTU("output.vtu", nk, ne, co, kon, ipkon);

    return 0;
}
