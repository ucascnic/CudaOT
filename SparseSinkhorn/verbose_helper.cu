#include"verbose_helper.h"
#include<stdio.h>
#include<cuda_runtime.h>
__global__ void show_parameters(CuSinkhornSolverParameters *s){
      printf("%2.f\n",s->maxError);
      printf("%2.f\n",s->absorption_scalingLowerBound);
      printf("%d\n",s->maxIterations);

  }



  __global__ void show_cuCost(CuCost *c){
      printf("%.6f\t%.6f\t%.6f\t",c->weight,c->WFlenscale,c->WFprefactor); // global rescaling parameter for Euclidean distance
      printf("%d\t%d\t\n",c->posDim,c->layerBottom); // global rescaling parameter for Euclidean distance

  }
  __global__ void show_data(double *data,int n, int start ){
      printf("\n");
      int end = start + n;
      for (int i = start ; i< end ;++i)
          printf("%.2f\t",data[i]);
  }
  __global__ void show_kenel(COO_SpareseMatrix_20 *kernel){
      printf("\t\tkernel entries: %d\n",kernel->nonZeros);
      for (int i = 0 ; i < 50; ++i)
          printf("row = %d col = %d %.8f\t",kernel->row_ptrl[i],kernel->col_ptrl[i],kernel->val[i]);
      printf("\n");
  }
  __global__ void show_kenel_sum(COO_SpareseMatrix_20 *kernel){
      double res = 0.0;
      for (int i = 0 ; i < kernel->nonZeros; ++i)
      {
          res += kernel->val[i];
      }
      printf("!!!!!res = %.10f\n",res);
  }
  __global__ void show_crs_kenel(COO_SpareseMatrix_20 *kernel,int xres){
      printf("\t\tkernel entries: %d\n",kernel->nonZeros);
      for (int i = 0 ; i < xres+1 ; ++i){
            printf("crs row = %d\n",kernel->crs_row_ptrl[i]);
      }
      for (int i = 0 ; i < kernel->nonZeros; ++i)
          printf("col = %d %.8f\n",kernel->col_ptrl[i],kernel->val[i]);
      printf("\n");
  }


__global__ void show_double_data(double *data,int *ind,int layer ,int num,double *target){
    int start = ind[layer];
    int end = start + num;
    printf("\n");
    for (int i = start; i < end; ++i){
        printf("%.8f\t",data[i]);
    }
    printf("\n");

}
__global__ void show_dd_data(double **test,int nLayers,int num){
    for (int i = 0 ; i < num ;++i)
        printf("%.8f\t",test[nLayers][i]);
}



__global__ void mu_dd_res(double *u,double *mu,double *res){
    int k = blockIdx.x;
    u[k] = mu[k] / res[k];

}

