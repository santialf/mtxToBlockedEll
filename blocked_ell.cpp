#include <stdio.h>            // printf
#include <stdlib.h>           // EXIT_FAILURE
#include <zlib.h>
#include <set>

#include <string.h>
#include <time.h>
#include <iostream>

#include "mmio.c"
#include "smsh.c"

float* createRandomArray(int n) {
    float* array = new float[n];

    for (int i = 0; i < n; i++) { 
        array[i] = 1;
    }

    return array;
}

/* Finds the possible amount of column blocks the matrix can have */
int findMaxNnz(int *rowPtr, int *colIndex, int num_rows, int block_size) {

    int max = 0;
    int num_blocks = num_rows / block_size;

    std::set<int> mySet;

    for(int i=0; i < num_blocks; i++) {

        for (int j = 0; j<block_size; j++) {
            int id = block_size*i+j;
            
            for(int k=rowPtr[id]; k<rowPtr[id+1]; k++)
                mySet.insert(colIndex[k]/block_size);
            
            if (mySet.size() > max)
                max = mySet.size();
        }
        mySet.clear();
	}

    return max*block_size;
}

/* Creates the array of block indexes for the blocked ell format */
int *createBlockIndex(int *rowPtr, int *colIndex, int num_rows, int block_size, int ell_cols) {

    long int mb = num_rows/block_size, nb = ell_cols/block_size;
    if (num_rows % block_size != 0)
        mb++;

    int* hA_columns = new int[(long int)nb*mb]();
    int ctr = 0;

    memset(hA_columns, -1, (long int)nb * mb * sizeof(int));
    std::set<int> mySet;

    /* Goes through the blocks of the matrix of block_size */
    for(int i=0; i<mb; i++) {

        /* Iterates through the rows of each block */
        for (int j = 0; j < block_size; j++) {
            int id = block_size*i + j;
            int index = 0;
            if (id >= num_rows)
                break;

            /* Iterates over the nnzs of each row */
            for(int k=rowPtr[id]; k<rowPtr[id+1]; k++) {    
                index = (colIndex[k]/block_size);
                mySet.insert(index);
            }
        }
        for (int elem : mySet)
            hA_columns[ctr++] = elem;
        
        ctr = i*nb+nb;
        mySet.clear();
    }
    return hA_columns; 
}

/* Creates the array of values for the blocked ell format */
float *createValueIndex(int *rowPtr, int *colIndex, float *values, int *hA_columns, int num_rows, int block_size, int ell_cols) {

    /* Allocate enough memory for the array */
    float* hA_values = new float[(long int)num_rows * ell_cols]();
    long int mb = num_rows/block_size, nb = ell_cols/block_size;
    if (num_rows % block_size != 0)
        mb++;

    /* Set all values to 0 */
    memset(hA_values, 0, (long int) num_rows * ell_cols * sizeof(int));

    /* Iterate the blocks in the y axis */
    for (int i=0; i<mb;i++){

        /* Iterate the lines of each block */
        for (int l = 0; l<block_size; l++) {
            int ctr = 0;

            /* Iterate the blocks in the block_id array (x axis) */
            for (int j = 0; j < nb; j++) {
                int id = nb*i + j;
                if (hA_columns[id] == -1)
                    break;

                /* Iterate each line of the matrix */
                for(int k=rowPtr[i*block_size+l]; k<rowPtr[i*block_size+l+1]; k++) {  

                    /* If the element is not in the same block, skip*/
                    if (colIndex[k]/block_size > hA_columns[id])
                        break;
                    else if (colIndex[k]/block_size == hA_columns[id]) 
                        hA_values[(long int)i*ell_cols*block_size+l*ell_cols+j*block_size+(colIndex[k]-(hA_columns[id]*block_size))] = values[k];
                }
            }
        }
    }
    
    return hA_values;
}

int main(int argc, char *argv[]) {

    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int A_num_rows, A_num_cols, nz, A_nnz;
    int i = 0, *I_complete, *J_complete;
    float *V_complete;

    /*******************************************************************/
    if ((f = fopen(argv[1], "r")) == NULL)
    {
        printf("Could not locate the matrix file. Please make sure the pathname is valid.\n");
        exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }
    matcode[4] = '\0';
    
    if ((ret_code = mm_read_mtx_crd_size(f, &A_num_rows, &A_num_cols, &nz)) != 0)
    {
        printf("Could not read matrix dimensions.\n");
        exit(1);
    }
    
    if ((strcmp(matcode, "MCRG") == 0) || (strcmp(matcode, "MCIG") == 0) || (strcmp(matcode, "MCPG") == 0) || (strcmp(matcode, "MCCG") == 0))
    {

        I_complete = (int *)calloc(nz, sizeof(int));
        J_complete = (int *)calloc(nz, sizeof(int));
        V_complete = (float *)calloc(nz, sizeof(float));

        for (i = 0; i < nz; i++)
        {
            if (matcode[2] == 'P') {
                fscanf(f, "%d %d", &I_complete[i], &J_complete[i]);
                V_complete[i] = 1;
            }  
            else {
                fscanf(f, "%d %d %f", &I_complete[i], &J_complete[i], &V_complete[i]);
            } 
            fscanf(f, "%*[^\n]\n");
            /* adjust from 1-based to 0-based */
            I_complete[i]--;
            J_complete[i]--;
        }
    }

    /* If the matrix is symmetric, we need to construct the other half */

    else if ((strcmp(matcode, "MCRS") == 0) || (strcmp(matcode, "MCIS") == 0) || (strcmp(matcode, "MCPS") == 0) || (strcmp(matcode, "MCCS") == 0) || (strcmp(matcode, "MCCH") == 0) || (strcmp(matcode, "MCRK") == 0) || (matcode[0] == 'M' && matcode[1] == 'C' && matcode[2] == 'P' && matcode[3] == 'S'))
    {

        I_complete = (int *)calloc(2 * nz, sizeof(int));
        J_complete = (int *)calloc(2 * nz, sizeof(int));
        V_complete = (float *)calloc(2 * nz, sizeof(float));

        int i_index = 0;

        for (i = 0; i < nz; i++)
        {
            if (matcode[2] == 'P') {
                fscanf(f, "%d %d", &I_complete[i], &J_complete[i]);
                V_complete[i] = 1;
            }
            else {
                fscanf(f, "%d %d %f", &I_complete[i], &J_complete[i], &V_complete[i]);
            }
                
            fscanf(f, "%*[^\n]\n");

            if (I_complete[i] == J_complete[i])
            {
                /* adjust from 1-based to 0-based */
                I_complete[i]--;
                J_complete[i]--;
            }
            else
            {
                /* adjust from 1-based to 0-based */
                I_complete[i]--;
                J_complete[i]--;
                J_complete[nz + i_index] = I_complete[i];
                I_complete[nz + i_index] = J_complete[i];
                V_complete[nz + i_index] = V_complete[i];
                i_index++;
            }
        }
        nz += i_index;
    }
    else
    {
        printf("This matrix type is not supported: %s \n", matcode);
        exit(1);
    }

    /* sort COO array by the rows */
    if (!isSorted(J_complete, I_complete, nz)) {
        quicksort(J_complete, I_complete, V_complete, nz);
    }
    
    /* Convert from COO to CSR */
    int* rowPtr = new int[A_num_rows + 1]();
    int* colIndex = new int[nz]();
    float* values = new float[nz]();

    for (i = 0; i < nz; i++) {
        colIndex[i] = J_complete[i];
        values[i] = V_complete[i];
        rowPtr[I_complete[i] + 1]++;
    }
    for (i = 0; i < A_num_rows; i++) {
        rowPtr[i + 1] += rowPtr[i];
    }
    A_nnz = nz;

    free(I_complete);
    free(J_complete);
    free(V_complete);
    fclose(f);
    /*******************************************************************/

    int A_ell_blocksize = 16;

    int * rowPtr_pad;
    int remainder = A_num_rows % A_ell_blocksize;
    if (remainder != 0) {
        A_num_rows = A_num_rows + (A_ell_blocksize - remainder);
        A_num_cols = A_num_cols + (A_ell_blocksize - remainder);
        rowPtr_pad = new int[A_num_rows + 1];
        for (int i=0; i<A_num_rows - (A_ell_blocksize - remainder); i++)
            rowPtr_pad[i] = rowPtr[i];
        for (int j=A_num_rows - (A_ell_blocksize - remainder); j<A_num_rows + 1; j++)
            rowPtr_pad[j] = nz;
        delete[] rowPtr;
    } else {
        rowPtr_pad = rowPtr;
    }  

    int A_ell_cols = findMaxNnz(rowPtr_pad, colIndex, A_num_rows, A_ell_blocksize);
    int A_num_blocks = A_ell_cols * A_num_rows / (A_ell_blocksize * A_ell_blocksize);
    int *hA_columns = createBlockIndex(rowPtr_pad, colIndex, A_num_rows, A_ell_blocksize, A_ell_cols);
    float *hA_values = createValueIndex(rowPtr_pad, colIndex, values, hA_columns, A_num_rows, A_ell_blocksize, A_ell_cols);

    float mem_ids = (float) A_num_blocks * sizeof(int) / 1000000000;
    float mem_values = (float) A_ell_cols * A_num_rows * sizeof(float) / 1000000000;
    std::cout << "Memory used: " << std::endl;
    std::cout << "Block indexes: " << mem_ids << " Gbytes" << std::endl;
    std::cout << "Block values: " << mem_values << " Gbytes" << std::endl;
    std::cout << "Total: " << mem_ids + mem_values << " Gbytes" << std::endl;

    delete[] rowPtr_pad;
    delete[] colIndex;
    delete[] values;
    delete[] hA_columns;
    delete[] hA_values;
    
    return 0;
}