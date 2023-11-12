#ifndef KERNEL_MANAGE_H
#define KERNEL_MANAGE_H
#include"basic_settings.h"
#include"data_struct.h"
#include"family.h"
#include <cusparse.h>
#include<cuda_runtime_api.h>
#include"Common/handler_cuda_error.h"

#define WARP_ 32

inline __device__ double __shfl_down_(double var, unsigned int srcLane, int width=WARP_) {
  int2 a = *reinterpret_cast<int2*>(&var);
  a.x = __shfl_down(a.x, srcLane, width);
  a.y = __shfl_down(a.y, srcLane, width);
  return *reinterpret_cast<double*>(&a);
}





template <int MAXELEMENT>
double getError_alpha_beta(COO_Sparese_Matrix<MAXELEMENT>  *kernel,
                           double *res,
                           int xres,int yres,double *muX);
template <int MAXELEMENT>
__global__ void update_alpha(COO_Sparese_Matrix<MAXELEMENT> *kernel,
                             XData *alphaH,
                                  double eps,int num,double *mu,int layer);
template <int MAXELEMENT>
__global__ void update_beta(COO_Sparese_Matrix<MAXELEMENT> *kernel,
                             XData *alphaH,
                                  double eps,int num,double *mu,int layer);
template <int MAXELEMENT>
__global__ void transport_kernel(COO_Sparese_Matrix<MAXELEMENT> * kernel,
                                 COO_Sparese_Matrix<MAXELEMENT> *kernelT);

template <int MAXELEMENT>
__global__ void check_cells_(int layer,int layerBottom,
                                double eps,
                                double slack,
                                double *muX,double *muY,
                                COO_Sparese_Matrix<MAXELEMENT> *kernel,
                                XData *alpha,XData *beta,
                                XData * xPos, XData *yPos,
                                XData *xRadii, XData *yRadii,
                             int xres,int yres);

template <int MAXELEMENT>
void sortcoo_and_to_csr(int xres,int yres,int N, COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                        Ind_Ptr_ *indptrl_cu_mid,Ind_Ptr_ *indptrl_cu,
                        COO_Buff * buff);


template<int ROW_, unsigned int COL_, int MAXLAYER_,int MAXELEMENT>
void refinekernel_(COO_Sparese_Matrix<MAXELEMENT> *kernel,COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                  double *muX_cuda,
                  double *muY_cuda,
                  int layerBottom,
                  double _eps,int xres,int yres,
                  int layer,
                   Familys<ROW_*COL_,MAXLAYER_> * ,
                   Familys<ROW_*COL_,MAXLAYER_> * ,
                  XData *cu_alphaH,XData *cu_betaH,
                  XData * cu_xposH, XData *cu_yposH,
                  XData *cu_xradii, XData *cu_yradii,
                 Ind_Ptr_ *indptrl_cu_mid,Ind_Ptr_ *indptrl_cu,
                 int *kernelValid,COO_Buff * buff);
template<int MAXELEMENT>
int generateKernel_(COO_Sparese_Matrix<MAXELEMENT> *kernel, COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                   CuSinkhornSolverParameters *cfg, double eps,
                   double *muX_cuda,
                   double *muY_cuda, int layer, XData * alphaH,
                   XData *betaH, XData * xpos, XData *ypos,
                   XData *rax, XData *ray, int, int, Ind_Ptr_ *indptrl_cu_mid,
                   Ind_Ptr_ *indptrl_cu, COO_Buff *buff);


__device__ double getCostEff(XData *xPos,XData *yPos,
        XData *xRadii,
        XData *yRadii,
        int layerBottom,
        int layerX, int x,
        int layerY, int y);

template <int MAXELEMENT>
__global__ void scoreTransportCost(XData *xPos,XData *yPos,
                                   double *res,double *u,double *v,
                                   int N,COO_Sparese_Matrix<MAXELEMENT> *kernel,int layer);


template <int MAXELEMENT>
__global__ void refine_kernel_2(COO_Sparese_Matrix<MAXELEMENT> *newkernel,
                             COO_Sparese_Matrix<MAXELEMENT> *oldKernel,
                             double *muX,
                             double *muY,
                             int layerBottom,
                             double _eps,
                             int layer,
                             XData *alpha,XData *beta,
                             XData * xPos, XData *yPos,
                             XData *xRadii, XData *yRadii,
                                double slack
                             );


template <int MAXELEMENT>
int refinekernel_2(COO_Sparese_Matrix<MAXELEMENT> *kernel,COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                  double *muX_cuda,
                  double *muY_cuda,
                  int layerBottom,
                  double _eps,int xres,int yres,
                  int layer,
                  XData *cu_alphaH,XData *cu_betaH,
                  XData * cu_xposH, XData *cu_yposH,
                  XData *cu_xradii, XData *cu_yradii,
                 Ind_Ptr_ *indptrl_cu_mid,Ind_Ptr_ *indptrl_cu,
                   int *kernelValid,COO_Buff * buff,
                    CuSinkhornSolverParameters *cfg);

inline __device__ double EUCL_lincombSqr(double * a,
                                   double * b) {
    // returns |sa*a+sb*b|^2, dimension given by n
    int i;
    double result,buffer;
    result=0;

    for(i=0;i<DIM;++i) {
//        buffer=sa*a[i]+sb*b[i];
        buffer=a[i]-b[i];
        result+=buffer*buffer;
    }
    return result;
}

#endif // KERNEL_MANAGE_H
