#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include "kernel_manage.h"
#include <stdio.h>
#include <oneapi/mkl.hpp>
#include <dpct/blas_utils.hpp>

#include "memory_manage.h"
#include <assert.h>

 double getCostEff(XData *xPos,XData *yPos,
        XData *xRadii,
        XData *yRadii,
        int layerBottom,
        int layerX, int x,
        int layerY, int y){
    double result;

    result=EUCL_lincombSqr(&xPos->data[xPos->ind[layerX]+(x*DIM)],
            &yPos->data[yPos->ind[layerY]+(y *DIM)]);
    result = sycl::sqrt(result);
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

void kernal_mat_u_v(double *u,double *v,double *mu, sycl::nd_item<3> item_ct1){

    int k = item_ct1.get_group(2);
    double temp = u[k]*v[k];
    u[k] = sycl::fabs(temp - mu[k]);
}


void checkAbsorb(int *msg_cu,double *u,double *v,double bound,int xres,int yres,
                 sycl::nd_item<3> item_ct1){
    int j = item_ct1.get_local_id(2) +
            item_ct1.get_group(2) * item_ct1.get_local_range(2);
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




void  SinkhornAbsorbScaling1(XData *alphaH,
                                       double *u, int layer, double eps,int res,
                                       sycl::nd_item<3> item_ct1)
{

//    int xres=HPX->layer[layer].nCells;
    int x = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
            item_ct1.get_local_id(2);
    if (x >= res)
        return;
    {
        alphaH->data[alphaH->ind[layer] + x] += eps * sycl::log(u[x]);
        u[x]=1.;
    }


}
void  SinkhornAbsorbScaling2(XData *alphaH,
                                       double *u, int layer, double eps,int res,
                                       sycl::nd_item<3> item_ct1)
{
    int x = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
            item_ct1.get_local_id(2);
    if (x >= res)return;
//     for(int x=0;x<xres;x++)
     {
        alphaH->data[alphaH->ind[layer] + x] += eps * sycl::log(u[x]);
        u[x]=1.;

    }


}




