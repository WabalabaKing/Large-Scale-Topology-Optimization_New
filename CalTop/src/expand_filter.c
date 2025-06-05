void mafillsm_expandfilter(double *FilterMatrixs, int *filternnzElems,
                           int *rowFilters, int *colFilters,
                           int ne, int fnnzassumed)
{
    printf("I am trying to expand the filter\n");

    for (int i = 0; i < ne; i++) {
        int row_nnz = filternnzElems[i];

        if (row_nnz > fnnzassumed) {
            printf("[ERROR] filternnzElems[%d] = %d > fnnzassumed = %d\n", i, row_nnz, fnnzassumed);
            continue;
        }

        for (int j = 0; j < row_nnz; j++) {
            int offset = i * fnnzassumed + j;

            int rowval = rowFilters[offset];
            int colval = colFilters[offset];
            double value = FilterMatrixs[offset];

            if (colval >= ne) {
                printf("[ERROR] colval (%d) >= ne (%d)\n", colval, ne);
                continue;
            }

            if (colval > i) {

                printf("I am in colVal \n");
                int new_row = colval;
                int new_col = i;
                int new_idx = filternnzElems[new_row];

                if (new_idx >= fnnzassumed) {
                    printf("[WARNING] Overflow in symmetric expansion at row %d\n", new_row);
                    continue;
                }

                int new_offset = new_row * fnnzassumed + new_idx;

                rowFilters[new_offset] = new_row;
                colFilters[new_offset] = new_col;
                FilterMatrixs[new_offset] = value;

                filternnzElems[new_row]++;
            }
        }
    }
}
