#ifndef __SKONDEV__H
#define __SKONDEV__H
#include<THierarchicalPartition.h>

#include"data_struct.h"
#include<stdio.h>
#include<cuda_runtime_api.h>


inline int changeLayer(int newlayer,int *layer,
                 int * kernelValid,
                 THierarchicalPartition *HPX,
                 THierarchicalPartition *HPY,
                 int *xres,int *yres,double *u,double *v,
                       double *origonals_one_cuda,
                       double *muX_cuda,double *muY_cuda,
                       XData *cu_muXH,XData *cu_muYH,
                       XData *muXH,XData *muYH
                 ){
    *layer = newlayer;
    *kernelValid = 0;
    printf("layer=%d\n",*layer);

    *xres=HPX->layers[*layer]->nCells;
    *yres=HPY->layers[*layer]->nCells;
    CHECK(cudaMemcpy(muX_cuda, &cu_muXH->data[muXH->ind[*layer]],
          sizeof(double)*(*xres), cudaMemcpyDeviceToDevice));
    CHECK(cudaMemcpy(muY_cuda, &cu_muYH->data[muYH->ind[*layer]],
          sizeof(double)*(*yres), cudaMemcpyDeviceToDevice));
    CHECK(cudaMemcpy(u, origonals_one_cuda, sizeof(double)* (*xres), cudaMemcpyDeviceToDevice));
    CHECK(cudaMemcpy(v, origonals_one_cuda, sizeof(double)* (*yres), cudaMemcpyDeviceToDevice));


    return 0;
}
#endif
