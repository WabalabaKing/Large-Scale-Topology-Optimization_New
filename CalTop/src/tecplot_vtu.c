#include <stdio.h>
#include <stdlib.h>

/**
 * @file tecplot_vtu.c
 * @brief Exports an elastic field solution to a VTU (VTK Unstructured Grid) file.
 *
 * This function writes nodal coordinates, connectivity information, and 
 * displacement field data to a VTU file (`elastic_Field.vtu`) for visualization 
 * in ParaView or other VTU-compatible software.
 *
 * @param nk     Number of nodes in the mesh.
 * @param ne     Number of elements in the mesh.
 * @param co     Array of nodal coordinates (size: 3 * nk).
 * @param kon    Element connectivity array, storing node indices per element.
 * @param ipkon  Index array mapping elements to the connectivity array.
 * @param v      Nodal displacement array (size: 3 * nk).
 *
 * @details
 * The function writes the data in VTK Unstructured Grid (.vtu) format, including:
 * - **Nodal coordinates**: Each node's (x, y, z) position.
 * - **Element connectivity**: Defines elements by listing their node indices.
 * - **Offsets**: Specifies the end position of each element's connectivity in the list.
 * - **Cell types**: Assumes quadrilateral elements (VTK type 10).
 * - **Nodal displacements**: Exports displacement vectors at each node.
 *
 * @note 
 * - The function assumes 4-node quadrilateral elements.
 * - The generated file uses ASCII format for readability.
 * - The function does not return any value but writes the data to `elastic_Field.vtu`.
 * - If the file cannot be opened, an error is printed, and the function returns.
 *
 * @usage
 * ```c
 * int nk = 100; // Example number of nodes
 * int ne = 50;  // Example number of elements
 * double co[300];  // Node coordinates (3 per node)
 * int kon[200];  // Connectivity array (4 nodes per element)
 * int ipkon[50]; // Index mapping elements to connectivity array
 * double v[300];  // Nodal displacements (3 per node)
 * 
 * tecplot_vtu(nk, ne, co, kon, ipkon, v);
 * ```
 */

void tecplot_vtu(int nk, int ne, double *co, int *kon, int *ipkon, double *lakon, int mi0, double *v, double *stx, double *rhoPhy) 
    
    {
    FILE *fp = fopen("elastic_Field.vtu", "w");
    if (fp == NULL) {
        perror("Error opening file");
        return;
    }

    // Write VTU XML header
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "  <UnstructuredGrid>\n");
    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nk, ne);

    // Write nodal coordinates
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int node = 0; node < nk; node++) {
        fprintf(fp, "        %.5f %.5f %.5f\n", co[3 * node], co[3 * node + 1], co[3 * node + 2]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");

    // Write connectivity
    fprintf(fp, "      <Cells>\n");
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    for (int ielem = 0; ielem < ne; ielem++) {
        char *lk= &lakon[8*ielem];
        int nope = 4;
        if (lk[3]=='1' && lk[4]=='0') nope = 10;
        for (int j = 0; j < nope; j++) {  // Node # based on nope (C3D4or C3D10)
            fprintf(fp, " %d", kon[ipkon[ielem] + j]-1);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // Write offsets
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    int offset = 0;
    for (int ielem = 0; ielem < ne; ielem++) {
        char *lk = &lakon[8*ielem];
        offset += (lk[3]=='1' && lk[4]=='0') ? 10 : 4;
        fprintf(fp, " %d\n", offset);
    }
    fprintf(fp, "        </DataArray>\n");

    // Write cell types
    fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (int ielem = 0; ielem < ne; ielem++) {
        char *lk = &lakon[8*ielem];
        int vtktype = (lk[3]=='1' && lk[4]=='0') ? 24 : 10;
        fprintf(fp, " %d\n", vtktype);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");

    // Write nodal displacements
    fprintf(fp, "      <PointData Scalars=\"Displacement\">\n");
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int node = 0; node < nk; node++) {
        fprintf(fp, "        %.8f %.8f %.8f\n", v[4*node+1], v[4*node+2], v[4* node+3]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </PointData>\n");

    // Write Cell Stress data
    fprintf(fp, "      <CellData Scalars=\"density\">\n");
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Stress\" NumberOfComponents=\"6\" format=\"ascii\">\n");
    for (int cell = 0; cell < ne; cell++) {
        char *lk = &lakon[8*cell];
        int ngp = (lk[3]=='1' && lk[4]=='0') ? 4 : 1;
        double s[6] = {0,0,0,0,0,0};
        for (int gp = 0; gp < ngp; gp++) {
            for (int c = 0; c < 6; c++)
                s[c] += stx[6*mi0*cell + 6*gp + c];
        }
        for (int c = 0; c < 6; c++) s[c] /= ngp;
        fprintf(fp, "      %.8f %.8f %.8f %.8f %.8f %.8f\n",
                s[0], s[1], s[2], s[3], s[4], s[5]);
    }
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for (int cell = 0; cell < ne; cell++) {
        fprintf(fp, "      %.8f\n",  rhoPhy[cell]);
    }
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "      </CellData>\n");


    // Close XML tags
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");

    fclose(fp);
}
