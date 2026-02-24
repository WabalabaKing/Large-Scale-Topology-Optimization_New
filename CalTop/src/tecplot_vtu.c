#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

void tecplot_vtu(int nk, int ne, double *co, int *kon, int *ipkon, double *v, double *stx, double *rhoPhy) 
    
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
        for (int j = 0; j < 4; j++) {  // Always tetrahedral elements
            fprintf(fp, " %d", kon[ipkon[ielem] + j]-1);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // Write offsets
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    for (int ielem = 0; ielem < ne; ielem++) {
        fprintf(fp, " %d\n", (ielem + 1) * 4);
    }
    fprintf(fp, "        </DataArray>\n");

    // Write cell types
    fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (int ielem = 0; ielem < ne; ielem++) {
        fprintf(fp, " 10\n");  //  VTU type 10 for 4-node tetrahedra
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

    // Write Cell Stress and Density data
    // Set default scalar to VonMises (nice for coloring in ParaView)
    fprintf(fp, "      <CellData Scalars=\"VonMises\">\n");

    // --- Von Mises scalar per cell ---
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"VonMises\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for (int cell = 0; cell < ne; cell++) {
        double sxx = stx[6*cell    ];
        double syy = stx[6*cell + 1];
        double szz = stx[6*cell + 2];
        double sxy = stx[6*cell + 3];
        double syz = stx[6*cell + 4];
        double szx = stx[6*cell + 5];

        double vm = sqrt(
            0.5 * ( (sxx - syy)*(sxx - syy)
                  + (syy - szz)*(syy - szz)
                  + (szz - sxx)*(szz - sxx) )
          + 3.0 * ( sxy*sxy + syz*syz + szx*szx )
        );

        fprintf(fp, "      %.8f\n", vm);
    }
    fprintf(fp, "        </DataArray>\n");

    // --- Full 6-component stress tensor per cell ---
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Stress\" NumberOfComponents=\"6\" format=\"ascii\">\n");
    for (int cell = 0; cell < ne; cell++) 
    {
        fprintf(fp, "      %.8f %.8f %.8f %.8f %.8f %.8f\n",
                stx[6*cell    ], stx[6*cell + 1], stx[6*cell + 2],
                stx[6*cell + 3], stx[6*cell + 4], stx[6*cell + 5]);
    }
    fprintf(fp, "        </DataArray>\n");

    // --- Density per cell ---
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    
    for (int cell = 0; cell < ne; cell++) 
    {
        fprintf(fp, "      %.8f\n", rhoPhy[cell]);
    }

    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "      </CellData>\n");


    // Close XML tags
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");

    fclose(fp);
}

void tecplot_vtu_passive(int nk, int ne,
                         double *co, int *kon, int *ipkon,
                         double *v, double *stx, double *rhoPhy,
                         int *passiveIDs, int numPassive)
{
    FILE *fp = fopen("elastic_Field_passive.vtu", "w");
    if (fp == NULL) {
        perror("Error opening file");
        return;
    }

    // --- VTU header ---
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "  <UnstructuredGrid>\n");
    // Only passive elements are written as cells; keep all nodes
    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nk, numPassive);

    // --- Points (all nodes) ---
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int node = 0; node < nk; node++) {
        fprintf(fp, "        %.5f %.5f %.5f\n",
                co[3 * node    ],
                co[3 * node + 1],
                co[3 * node + 2]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");

    // --- Cells: only passive elements ---
    fprintf(fp, "      <Cells>\n");

    // Connectivity
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    for (int k = 0; k < numPassive; k++) {
        int e = passiveIDs[k] - 1;  // convert 1-based element ID to 0-based index
        for (int j = 0; j < 4; j++) {  // Q4 elements
            fprintf(fp, " %d", kon[ipkon[e] + j] - 1); // nodes are 1-based → 0-based
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // Offsets (still 4 nodes per cell)
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    for (int k = 0; k < numPassive; k++) {
        fprintf(fp, " %d\n", (k + 1) * 4);
    }
    fprintf(fp, "        </DataArray>\n");

    // Types (VTK_TETRA = 10 for 4-node tets)
    fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (int k = 0; k < numPassive; k++) {
        fprintf(fp, " 10\n");
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");

    // --- PointData: Displacement at all nodes ---
    fprintf(fp, "      <PointData Scalars=\"Displacement\">\n");
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int node = 0; node < nk; node++) {
        // CalculiX vold-style: [0]=dummy, [1]=ux, [2]=uy, [3]=uz
        fprintf(fp, "        %.8f %.8f %.8f\n",
                v[4 * node + 1],
                v[4 * node + 2],
                v[4 * node + 3]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </PointData>\n");

    // --- CellData: VonMises, Stress, Density only for passive elements ---
    fprintf(fp, "      <CellData Scalars=\"VonMises\">\n");

    // Von Mises per passive cell
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"VonMises\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for (int k = 0; k < numPassive; k++) {
        int e = passiveIDs[k] - 1;  // 0-based element index

        double sxx = stx[6*e    ];
        double syy = stx[6*e + 1];
        double szz = stx[6*e + 2];
        double sxy = stx[6*e + 3];
        double syz = stx[6*e + 4];
        double szx = stx[6*e + 5];

        double vm = sqrt(
            0.5 * ( (sxx - syy)*(sxx - syy)
                  + (syy - szz)*(syy - szz)
                  + (szz - sxx)*(szz - sxx) )
          + 3.0 * ( sxy*sxy + syz*syz + szx*szx )
        );

        fprintf(fp, "      %.8f\n", vm);
    }
    fprintf(fp, "        </DataArray>\n");

    // Full 6-component stress for passive cells
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Stress\" NumberOfComponents=\"6\" format=\"ascii\">\n");
    for (int k = 0; k < numPassive; k++) {
        int e = passiveIDs[k] - 1;
        fprintf(fp, "      %.8f %.8f %.8f %.8f %.8f %.8f\n",
                stx[6*e    ], stx[6*e + 1], stx[6*e + 2],
                stx[6*e + 3], stx[6*e + 4], stx[6*e + 5]);
    }
    fprintf(fp, "        </DataArray>\n");

    // Density for passive cells
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for (int k = 0; k < numPassive; k++) {
        int e = passiveIDs[k] - 1;
        fprintf(fp, "      %.8f\n", rhoPhy[e]);
    }
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "      </CellData>\n");

    // --- Close ---
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");

    fclose(fp);
}


void tecplot_vtu_active(int nk, int ne,
                        double *co, int *kon, int *ipkon,
                        double *v, double *stx, double *rhoPhy,
                        int *passiveIDs, int numPassive)
{
    // Build a mask of which elements are passive
    int *isPassive = (int*)calloc(ne, sizeof(int));
    if (!isPassive) {
        perror("Error allocating isPassive");
        return;
    }

    for (int i = 0; i < numPassive; i++) {
        int eID = passiveIDs[i];   // 1-based element ID from skinElementList.nam
        int e   = eID - 1;         // convert to 0-based index
        if (e >= 0 && e < ne) {
            isPassive[e] = 1;
        }
    }

    // Count active elements
    int numActive = 0;
    for (int e = 0; e < ne; e++) {
        if (!isPassive[e]) numActive++;
    }

    if (numActive == 0) {
        fprintf(stderr, "Warning: no active elements found (all are passive?). Skipping VTU write.\n");
        free(isPassive);
        return;
    }

    FILE *fp = fopen("elastic_Field_active.vtu", "w");
    if (fp == NULL) {
        perror("Error opening elastic_Field_active.vtu");
        free(isPassive);
        return;
    }

    // --- VTU header ---
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "  <UnstructuredGrid>\n");
    // Only active elements as cells; keep all nodes
    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nk, numActive);

    // --- Points (all nodes) ---
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int node = 0; node < nk; node++) {
        fprintf(fp, "        %.5f %.5f %.5f\n",
                co[3 * node    ],
                co[3 * node + 1],
                co[3 * node + 2]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");

    // --- Cells: only active elements ---
    fprintf(fp, "      <Cells>\n");

    // Connectivity
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    for (int e = 0; e < ne; e++) {
        if (isPassive[e]) continue;
        for (int j = 0; j < 4; j++) {  // Q4 elements
            fprintf(fp, " %d", kon[ipkon[e] + j] - 1); // nodes 1-based → 0-based
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // Offsets: 4 nodes per active cell
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    int offset = 4;
    for (int i = 0; i < numActive; i++) {
        fprintf(fp, " %d\n", offset);
        offset += 4;
    }
    fprintf(fp, "        </DataArray>\n");

    // Types: VTK_TET = 10
    fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (int i = 0; i < numActive; i++) {
        fprintf(fp, " 10\n");
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");

    // --- PointData: Displacement at all nodes ---
    fprintf(fp, "      <PointData Scalars=\"Displacement\">\n");
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int node = 0; node < nk; node++) {
        // CalculiX vold-style: [0]=dummy, [1]=ux, [2]=uy, [3]=uz
        fprintf(fp, "        %.8f %.8f %.8f\n",
                v[4 * node + 1],
                v[4 * node + 2],
                v[4 * node + 3]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </PointData>\n");

    // --- CellData: VonMises, Stress, Density only for active elements ---
    fprintf(fp, "      <CellData Scalars=\"VonMises\">\n");

    // Von Mises per active cell
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"VonMises\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for (int e = 0; e < ne; e++) {
        if (isPassive[e]) continue;

        double sxx = stx[6*e    ];
        double syy = stx[6*e + 1];
        double szz = stx[6*e + 2];
        double sxy = stx[6*e + 3];
        double syz = stx[6*e + 4];
        double szx = stx[6*e + 5];

        double vm = sqrt(
            0.5 * ( (sxx - syy)*(sxx - syy)
                  + (syy - szz)*(syy - szz)
                  + (szz - sxx)*(szz - sxx) )
          + 3.0 * ( sxy*sxy + syz*syz + szx*szx )
        );

        fprintf(fp, "      %.8f\n", vm);
    }
    fprintf(fp, "        </DataArray>\n");

    // Full 6-component stress for active cells
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Stress\" NumberOfComponents=\"6\" format=\"ascii\">\n");
    for (int e = 0; e < ne; e++) {
        if (isPassive[e]) continue;
        fprintf(fp, "      %.8f %.8f %.8f %.8f %.8f %.8f\n",
                stx[6*e    ], stx[6*e + 1], stx[6*e + 2],
                stx[6*e + 3], stx[6*e + 4], stx[6*e + 5]);
    }
    fprintf(fp, "        </DataArray>\n");

    // Density for active cells
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for (int e = 0; e < ne; e++) {
        if (isPassive[e]) continue;
        fprintf(fp, "      %.8f\n", rhoPhy[e]);
    }
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "      </CellData>\n");

    // --- Close ---
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");

    fclose(fp);
    free(isPassive);
}
