#ifndef SOLVER_H
#define SOLVER_H
#include<cuda_runtime_api.h>
#include<basic_settings.h>
#include<data_struct.h>
template<  int MAXELEMENT>
__global__ void kernal_mat_u(COO_Sparese_Matrix<MAXELEMENT> *kernel,
                             double *u,double *res,
                             int num,double *mu,double *v
                             );

template<  int MAXELEMENT>
__global__ void transport_kernel_2(COO_Sparese_Matrix<MAXELEMENT> * kernel,
                                 COO_Sparese_Matrix<MAXELEMENT> *kernelT,int N);
template<int MAXELEMENT>
int generateKernel_(COO_Sparese_Matrix<MAXELEMENT> *kernel, COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                   CuSinkhornSolverParameters *cfg,double eps,
                   double *muX_cuda,
                   double *muY_cuda,int layer,
                    XData * alphaH,
                   XData *betaH,XData * xpos, XData *ypos,
                   XData *rax, XData *ray,int xres,int yres,
                   Ind_Ptr_ *indptrl_cu_mid,Ind_Ptr_ *indptrl_cu,
                   COO_Buff * buff);
template<  int MAXELEMENT>
__global__ void v_mat_kernal_mat_u_m_mu(COO_Sparese_Matrix<MAXELEMENT> *kernel,
                             double *u,double * res, int num,double *mu,double *v);
__global__ void kernal_mat_u_v(double *u,double *v,double *mu);
template<  int MAXELEMENT>
void iterate(int *msg, int n,int xres, double *muX,
                       int yres,double * muY,
                        COO_Sparese_Matrix<MAXELEMENT> *kernel,
                        COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                        double *u,double *v,double *res) ;

__global__ void checkAbsorb(int *msg_cu,double *u,double *v,double bound,int xres,int yres);

template<  int MAXELEMENT>
double getError(COO_Sparese_Matrix<MAXELEMENT>  *kernel,
                    double *muX,
                    double *u,double *v,double *res,int xres,int yres,
                cublasHandle_t handle);

template<  int MAXELEMENT>
__global__ void get_alpha_error(COO_Sparese_Matrix<MAXELEMENT> *kernel,
                                double *res,
                               int yres,double *muX);


__global__ void  SinkhornAbsorbScaling1(XData *alphaH,
                                       double *u, int layer, double eps,int res);
__global__ void  SinkhornAbsorbScaling2(XData *alphaH,
                                       double *u, int layer, double eps,int res);



template< int ROW_,  int COL_, int MAXLAYER_ ,int MAXELEMENT>
int  solveSingle_(int *, int *kernelValid,
                 COO_Sparese_Matrix<MAXELEMENT> *kernel, COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                  CuSinkhornSolverParameters *cfg,
                 double eps,
                 double *muX_cuda,
                 double *muY_cuda,
                 Familys<ROW_*COL_,MAXLAYER_> * ,
                 Familys<ROW_*COL_,MAXLAYER_> * ,
                 int layerBottom,
                 XData * alphaH, XData *betaH,
                 XData *xpos, XData *ypos,
                 XData *rax, XData *ray,
                 int xres, int yres,
                 double *u, double *v, double *res, int layer,
                 Ind_Ptr_ *indptrl_cu_mid,
                 Ind_Ptr_ *indptrl_cu,
                  COO_Buff * buff, cublasHandle_t handle);


template<  int MAXELEMENT>
double calculate_result(int N,XData *cu_xposH, XData *cu_yposH,
                        double *res,double *u,double *v,
                        COO_Sparese_Matrix<MAXELEMENT> *kernel,int layer);


template<  int MAXELEMENT>
void iterate_2(int xres,int yres, double *muX,double * muY,
                        COO_Sparese_Matrix<MAXELEMENT> *kernel,
                        COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                        XData *alphaH,XData *betaH,double eps,
                        int layer) ;


template <int MAXELEMENT>
int refinekernel_3(COO_Sparese_Matrix<MAXELEMENT>* kernel, COO_Sparese_Matrix<MAXELEMENT>* kernelT,
    double* muX_cuda,
    double* muY_cuda,
    int layerBottom,
    double _eps, int xres, int yres,
    int layer,
    XData* cu_alphaH, XData* cu_betaH,
    XData* cu_xposH, XData* cu_yposH,
    XData* cu_xradii, XData* cu_yradii,
    Ind_Ptr_* indptrl_cu_mid, Ind_Ptr_* indptrl_cu,
    int* kernelValid, COO_Buff* buff,
    CuSinkhornSolverParameters* cfg);

template<int MAXELEMENT>
__global__ void generate_kernel_32_wrap(int layer, int layerBottom,
    double eps,
    double slack,
    double* muX, double* muY,
    COO_Sparese_Matrix<MAXELEMENT>* kernel,
    XData* alpha, XData* beta,
    XData* xPos, XData* yPos,
    XData* xRadii, XData* yRadii,
    int xres, int yres);

#endif // SOLVER_H
