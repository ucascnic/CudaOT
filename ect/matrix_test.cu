#include "cublas_v2.h"
#include <cusparse.h>
#include"matrixbase.h"
#include<stdio.h>
#include <assert.h>
#include <cuda_runtime.h>
#include<cuda_runtime_api.h>
#include<Common.h>
#include"Common/handler_cuda_error.h"
#include"memory_manage.h"
static cusparseHandle_t cusparse_handle = 0;
static void init_cusparse() {
  if (cusparse_handle == 0) {
    cusparseStatus_t status = cusparseCreate(&cusparse_handle);
    if (status != CUSPARSE_STATUS_SUCCESS) {
      printf("CUSPARSE Library initialization failed");
    }
  }
}

#include<stdlib.h>
void show_kenel(COO_SpareseMatrix_20 *kernel){
    printf("\t\tkernel entries: %d\n",kernel->nonZeros);
    for (int i = 0 ; i < kernel->nonZeros; ++i)
        printf("%d %d %d \n",kernel->row_ptrl[i],kernel->col_ptrl[i],kernel->val[i]);
    printf("\n");
}


__global__ void change_val(COO_SpareseMatrix_20 *kernel_dev,IndPtr_20 *indptr){
    for (int i = 0 ; i  < kernel_dev->nonZeros; ++i)
        ;
}
int main(){

    cusparseStatus_t Status_t;
    init_cusparse();
    COO_SpareseMatrix_20 *kernel_dev = (COO_SpareseMatrix_20 *)malloc(sizeof(COO_SpareseMatrix_20));
    int n = 9;
    kernel_dev->nonZeros = n;
    for (int i = 0 ; i  < n; ++i){
        kernel_dev->val[i] = 1;
    }
    kernel_dev->col_ptrl[0] = 0;
    kernel_dev->col_ptrl[1] = 1;
    kernel_dev->col_ptrl[2] = 1;
    kernel_dev->col_ptrl[3] = 2;
    kernel_dev->col_ptrl[4] = 0;
    kernel_dev->col_ptrl[5] = 3;
    kernel_dev->col_ptrl[6] = 4;
    kernel_dev->col_ptrl[7] = 2;
    kernel_dev->col_ptrl[8] = 4;


    kernel_dev->row_ptrl[0] = 0;
    kernel_dev->row_ptrl[1] = 1;
    kernel_dev->row_ptrl[2] = 0;
    kernel_dev->row_ptrl[3] = 1;
    kernel_dev->row_ptrl[4] = 2;
    kernel_dev->row_ptrl[5] = 2;
    kernel_dev->row_ptrl[6] = 2;
    kernel_dev->row_ptrl[7] = 3;
    kernel_dev->row_ptrl[8] = 3;

    show_kenel(kernel_dev);

    COO_SpareseMatrix_20 *kernel;
    cudaMalloc((void **)&kernel,sizeof(COO_SpareseMatrix_20));
    CHECK(cudaMemcpy(kernel, kernel_dev, sizeof(COO_SpareseMatrix_20),
                     cudaMemcpyHostToDevice));
    // sort row
    size_t workspace_size = 128;
    int num_rows = 4;
    int num_cols = 5;
    Status_t = cusparseXcoosort_bufferSizeExt(
        cusparse_handle,
        num_rows, num_cols,
        n,
        kernel->row_ptrl,
        kernel->col_ptrl,
        &workspace_size);
    assert( CUSPARSE_STATUS_SUCCESS == Status_t);
    void *buffer_for_coo_sort;
    CHECK(cudaMalloc(&buffer_for_coo_sort, sizeof(char) * workspace_size));
    IndPtr_20 *indptrl = (IndPtr_20 *)malloc(sizeof(IndPtr_20));
    for (int i = 0 ;  i <MAXELEMENT_20;++i){
        indptrl->ind[i] = i;
    }
    IndPtr_20 *indptrl_cu = malloc_struct(indptrl);

    Status_t = cusparseXcoosortByRow(
        cusparse_handle,
        num_rows, num_cols,
        n,
        kernel->row_ptrl,
        kernel->col_ptrl,
        indptrl_cu->ind,
        buffer_for_coo_sort);
//    Status_t = cusparseDgthr(
//                cusparse_handle,
//                n,
//                d_cooVals,
//                d_cooVals_sorted,
//                indptrl_cu->ind,
////                CUSPARSE_INDEX_BASE_ZERO);
//    assert( CUSPARSE_STATUS_SUCCESS == status)

    int number_of_rows = 4;
    cusparseXcoo2csr(cusparse_handle,
        kernel->row_ptrl, n, number_of_rows,
        kernel->col_ptrl, CUSPARSE_INDEX_BASE_ZERO);

    CHECK(cudaMemcpy(kernel_dev,kernel, sizeof(COO_SpareseMatrix_20),
                     cudaMemcpyDeviceToHost));


    show_kenel(kernel_dev);

//    show_kenel<<<1,1>>>(kernel);
    return 0;

}
