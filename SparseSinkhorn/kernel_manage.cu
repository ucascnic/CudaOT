#include"kernel_manage.h"
#include <stdio.h>
#include<cuda_runtime_api.h>
#include <cusparse.h>
#include"memory_manage.h"
#include <cuda.h>
#include<cublas_v2.h>
#include<assert.h>




 __device__ double getCostEff(XData *xPos,XData *yPos,
        XData *xRadii,
        XData *yRadii,
        int layerBottom,
        int layerX, int x,
        int layerY, int y){
    double result;

    result=EUCL_lincombSqr(&xPos->data[xPos->ind[layerX]+(x*DIM)],
            &yPos->data[yPos->ind[layerY]+(y *DIM)]);
    result = sqrt(result);
    if((layerX<layerBottom) || (layerY<layerBottom)) {
        // if not at finest layer, need to compute lower bound
        if(layerX<layerBottom) {
            result-=xRadii->data[xRadii->ind[layerX] + x ];
        }
        if(layerY<layerBottom) {
            result-=yRadii->data[yRadii->ind[layerY] + y];
        }
        if(result<0) {
            result=0;
        }
    }
    // if( p_index == 2)
    //    return result;

     return pow(result,((double)p_index));

}

__global__ void kernal_mat_u_v(double *u,double *v,double *mu){

    int k = blockIdx.x;
    double temp = u[k]*v[k];
    u[k] = abs(temp - mu[k]);

}


__global__ void checkAbsorb(int *msg_cu,double *u,double *v,double bound,int xres,int yres){
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if(j<xres)
        if( u[j] > bound){
            *msg_cu = MSG_ABSORB_REITERATE;
            return ;
        }

    if(j<yres)
        if( v[j] > bound){
            *msg_cu = MSG_ABSORB_REITERATE;
            return ;
        }


}




__global__ void  SinkhornAbsorbScaling1(XData *alphaH,
                                       double *u, int layer, double eps,int res)
{

//    int xres=HPX->layer[layer].nCells;
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x >= res)
        return;
    {
        alphaH->data[ alphaH->ind[layer] + x ] += eps*log(u[x]);
        u[x]=1.;
    }


}
__global__ void  SinkhornAbsorbScaling2(XData *alphaH,
                                       double *u, int layer, double eps,int res)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x >= res)return;
//     for(int x=0;x<xres;x++)
     {
        alphaH->data[ alphaH->ind[layer] + x ] += eps * log(u[x]);
        u[x]=1.;

    }


}




