#include <stdio.h>

/**
 * @brief Populate the flattened filter matrix and associated row/column index arrays
 *        from sparse data given in coordinate format (drow, dcol, dval).
 *
 * This function assumes that the 2D matrices FilterMatrixs, rowFilters, and colFilters
 * are stored in a flattened 1D array using column-major storage, i.e., (i, j) maps to
 * i + fnnzassumed * j.
 *
 * @param FilterMatrixs  Flattened 2D array of filter weights [fnnzassumed x ne]
 * @param rowFilters     Flattened 2D array of row indices [fnnzassumed x ne]
 * @param colFilters     Flattened 2D array of column indices [fnnzassumed x ne]
 * @param filternnzElems Array storing the number of nonzeros expected per row (size ne)
 * @param drow           Row indices of nonzero entries in COO format (size filternnz)
 * @param dcol           Column indices of nonzero entries in COO format (size filternnz)
 * @param dval           Values of nonzero entries in COO format (size filternnz)
 * @param ne             Number of elements (rows in the filter matrix)
 * @param ttime          Total time (not used here, kept for compatibility)
 * @param time           Current time (not used here, kept for compatibility)
 * @param ne0            Starting element index (not used)
 * @param filternnz      Total number of nonzeros in the filter
 * @param fnnzassumed    Assumed max number of nonzeros per row (used to size flattened arrays)
 */
void assembleFilter(double *FilterMatrixs, int *rowFilters, int *colFilters,
                int *filternnzElems, int *drow, int *dcol, double *dval,
                int ne, int ne0,
                int *filternnz, int *fnnzassumed) 
{

    int i;              // Loop counter for nonzeros
    int rowval, colval; // Row and column indices for the current nonzero
    double value;       // Value of the current nonzero

    int index = 1;      // Shared counter for writing into flattened arrays
                        // (shared across all rows — may cause overwrite)

    for (i = 0; i < *filternnz; i++) 
    {
        rowval = drow[i];    // Target row (0-based index)
        colval = dcol[i];    // Target column
        value  = dval[i];    // Corresponding value

        // Compute flattened array offset for (index, rowval)
        //int offset = index + (*fnnzassumed) * rowval;

        int offset = (index - 1) + (*fnnzassumed) * (rowval - 1);
        //FilterMatrixs[offset] = value;

        // Store data in the filter structure
        rowFilters[offset]    = rowval;
        colFilters[offset]    = colval;
        FilterMatrixs[offset] = value;

        // Update index (shared for all rows — may lead to memory overwrites)
        if (index < filternnzElems[rowval]) 
        {
            index++;
        }
        else 
        {
            index = 1;
        }
    }

    return;
}
