#include<time.h>
#include<../Common/Common.h>
#include<stdio.h>
#include<iostream>
#include"common.h"
#include"cuda_runtime.h"
#include "cublas_v2.h"
static void HandleError( cudaError_t err, const char *file, int line )
{    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n",
                cudaGetErrorString( err ),file, line );
        exit( EXIT_FAILURE );    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#define N 64
#define XROW (N*N)
#define YCOL (N*N)
#define PROBLEMSIZE XROW*YCOL

#define ITERSIZE 20
#define THREAD 128


typedef double datatype;
typedef double costtype;

typedef float mindatatype;





#define ITERSIZE 20
__device__ __constant__  double iter_eps[ITERSIZE];
__device__ __constant__  double iter_eps_[ITERSIZE];

__global__  void iter11(int *iter){
    iter[0]++;
}

__global__  void kernal_row_min(datatype *g,datatype *f,costtype *c,
                                datatype * logx,int *eps){
    unsigned int tid = threadIdx.x;
    unsigned int bid = blockIdx.x;
    unsigned int col = tid + blockIdx.y * THREAD * 8;

    unsigned int j = bid *YCOL + col;
    datatype a,a1,a2,a3,a4,a5,a6,a7;


    a  = ((datatype)c[j           ]) - g[col];
    a1 = ((datatype)c[j +   THREAD]) - g[col+ THREAD];
    a2 = ((datatype)c[j + 2*THREAD]) - g[col+ 2*THREAD];
    a3 = ((datatype)c[j + 3*THREAD]) - g[col+ 3*THREAD];
    a4 = ((datatype)c[j + 4*THREAD]) - g[col+ 4*THREAD];
    a5 = ((datatype)c[j + 5*THREAD]) - g[col+ 5*THREAD];
    a6 = ((datatype)c[j + 6*THREAD]) - g[col+ 6*THREAD];
    a7 = ((datatype)c[j + 7*THREAD]) - g[col+ 7*THREAD];

    a6 = min(a6,a7);
    a4 = min(a4,a5);
    a2 = min(a2,a3);
    a  = min(a ,a1);

    a4 = min(a4,a6);
    a  = min(a,a2);


    __shared__ mindatatype sdata[THREAD];

    sdata[tid] = (mindatatype) min(a,a4);
    __syncthreads();

    if (tid < 64)
        sdata[tid]  = min ( sdata[tid] ,sdata[tid + 64]);
    __syncthreads();

    if (tid < 32) {
        volatile mindatatype *vmem = sdata;
        vmem[tid] = min(vmem[tid],vmem[tid+32]);
        vmem[tid] = min(vmem[tid],vmem[tid+16]);
        vmem[tid] = min(vmem[tid],vmem[tid+8]);
        vmem[tid] = min(vmem[tid],vmem[tid+4]);
        vmem[tid] = min(vmem[tid],vmem[tid+2]);
        vmem[tid] = min(vmem[tid],vmem[tid+1]);
    }
    if (tid == 0)
    {
        f[bid] = sdata[0];
    }
    __syncthreads();
    mindatatype b0,b1,b2,b3,b4,b5,b6,b7;
    b0 =  __expf((float)(__dmul_rn((((datatype)c[j              ]) - g[col            ] - f[bid]) ,iter_eps_[ *eps])));
    b1 =  __expf((float)(__dmul_rn((((datatype)c[j  +     THREAD]) - g[col +    THREAD] - f[bid]),iter_eps_[ *eps])));
    b2 =  __expf((float)(__dmul_rn((((datatype)c[j  +  2* THREAD]) - g[col + 2* THREAD] - f[bid]),iter_eps_[ *eps])));
    b3 =  __expf((float)(__dmul_rn((((datatype)c[j  +  3* THREAD]) - g[col + 3* THREAD] - f[bid]),iter_eps_[ *eps])));
    b4 =  __expf((float)(__dmul_rn((((datatype)c[j  +  4* THREAD]) - g[col + 4* THREAD] - f[bid]),iter_eps_[ *eps])));
    b5 =  __expf((float)(__dmul_rn((((datatype)c[j  +  5* THREAD]) - g[col + 5* THREAD] - f[bid]),iter_eps_[ *eps])));
    b6 =  __expf((float)(__dmul_rn((((datatype)c[j  +  6* THREAD]) - g[col + 6* THREAD] - f[bid]),iter_eps_[ *eps])));
    b7 =  __expf((float)(__dmul_rn((((datatype)c[j  +  7* THREAD]) - g[col + 7* THREAD] - f[bid]),iter_eps_[ *eps])));

    mindatatype mySum = b0+b1+b2+b3+b4+b5+b6+b7;

    sdata[tid] = mySum;
    __syncthreads();

    if (tid < 64)
        sdata[tid] += sdata[tid + 64];
    __syncthreads();

    if (tid < 32) {
    volatile mindatatype *vsmem = sdata;
    vsmem[tid] += vsmem[tid + 32];
    vsmem[tid] += vsmem[tid + 16];
    vsmem[tid] += vsmem[tid + 8];
    vsmem[tid] += vsmem[tid + 4];
    vsmem[tid] += vsmem[tid + 2];
    vsmem[tid] += vsmem[tid + 1];
    }

    if (tid == 0)
    {
        f[bid] +=  (iter_eps[ *eps]) *  (log(sdata[0])- logx[bid] ) ;
    }

}
__global__  void kernal_row_min(datatype *g,costtype *c,
                          mindatatype * mindata){
    unsigned int tid = threadIdx.x;
    unsigned int bid = blockIdx.x;
    unsigned int col = tid + blockIdx.y * THREAD * 8;

    unsigned int j = bid *YCOL + col;
    datatype a,a1,a2,a3,a4,a5,a6,a7;


    a  = ((datatype)c[j           ]) - g[col];
    a1 = ((datatype)c[j +   THREAD]) - g[col+ THREAD];
    a2 = ((datatype)c[j + 2*THREAD]) - g[col+ 2*THREAD];
    a3 = ((datatype)c[j + 3*THREAD]) - g[col+ 3*THREAD];
    a4 = ((datatype)c[j + 4*THREAD]) - g[col+ 4*THREAD];
    a5 = ((datatype)c[j + 5*THREAD]) - g[col+ 5*THREAD];
    a6 = ((datatype)c[j + 6*THREAD]) - g[col+ 6*THREAD];
    a7 = ((datatype)c[j + 7*THREAD]) - g[col+ 7*THREAD];

    a6 = min(a6,a7);
    a4 = min(a4,a5);
    a2 = min(a2,a3);
    a  = min(a ,a1);

    a4 = min(a4,a6);
    a  = min(a,a2);


    __shared__ mindatatype sdata[THREAD];

    sdata[tid] = (mindatatype) min(a,a4);
    __syncthreads();

    if (tid < 64)
        sdata[tid]  = min ( sdata[tid] ,sdata[tid + 64]);
    __syncthreads();

    if (tid < 32) {
    volatile mindatatype *vmem = sdata;
        vmem[tid] = min(vmem[tid],vmem[tid+32]);
        vmem[tid] = min(vmem[tid],vmem[tid+16]);
        vmem[tid] = min(vmem[tid],vmem[tid+8]);
        vmem[tid] = min(vmem[tid],vmem[tid+4]);
        vmem[tid] = min(vmem[tid],vmem[tid+2]);
        vmem[tid] = min(vmem[tid],vmem[tid+1]);
    }
    if (tid == 0)
        mindata[j] = sdata[0];

}


__global__  void kernal_global_row_min(datatype *f,
                                       mindatatype *mindata){
    unsigned int tid = threadIdx.x;
    unsigned int bid = blockIdx.x;
    unsigned int row = bid *THREAD + tid;
    unsigned int j = YCOL * row ;

    mindatatype a1 = min(mindata[j],mindata[j+8*THREAD]);
    mindatatype a2 = min(mindata[j+16*THREAD],mindata[j+24*THREAD]);
    f[row] = min(a1,a2);

}
#include <cooperative_groups.h>

namespace cg = cooperative_groups;

__global__ void kernal_cal_newkernal(datatype *f,datatype *g,
                                     costtype *c,
                                     mindatatype *mindata,int *eps){
    unsigned int tid = threadIdx.x;
    unsigned int bid = blockIdx.x;
    unsigned int col = tid + blockIdx.y * THREAD * 8;
    cg::thread_block cta = cg::this_thread_block();
    unsigned int j = bid *YCOL + col;
    mindatatype a0,a1,a2,a3,a4,a5,a6,a7;
    a0 =  __expf((float)(__dmul_rn((((datatype)c[j              ]) - g[col            ] - f[bid]) ,iter_eps_[ *eps])));
    a1 =  __expf((float)(__dmul_rn((((datatype)c[j  +     THREAD]) - g[col +    THREAD] - f[bid]),iter_eps_[ *eps])));
    a2 =  __expf((float)(__dmul_rn((((datatype)c[j  +  2* THREAD]) - g[col + 2* THREAD] - f[bid]),iter_eps_[ *eps])));
    a3 =  __expf((float)(__dmul_rn((((datatype)c[j  +  3* THREAD]) - g[col + 3* THREAD] - f[bid]),iter_eps_[ *eps])));
    a4 =  __expf((float)(__dmul_rn((((datatype)c[j  +  4* THREAD]) - g[col + 4* THREAD] - f[bid]),iter_eps_[ *eps])));
    a5 =  __expf((float)(__dmul_rn((((datatype)c[j  +  5* THREAD]) - g[col + 5* THREAD] - f[bid]),iter_eps_[ *eps])));
    a6 =  __expf((float)(__dmul_rn((((datatype)c[j  +  6* THREAD]) - g[col + 6* THREAD] - f[bid]),iter_eps_[ *eps])));
    a7 =  __expf((float)(__dmul_rn((((datatype)c[j  +  7* THREAD]) - g[col + 7* THREAD] - f[bid]),iter_eps_[ *eps])));

    __shared__ mindatatype sdata[THREAD];
    mindatatype mySum = a0+a1+a2+a3+a4+a5+a6+a7;

    sdata[tid] = mySum;
    __syncthreads();

    if (tid < 64)
        sdata[tid] += sdata[tid + 64];
    __syncthreads();

    if (tid < 32) {
    volatile mindatatype *vsmem = sdata;
    vsmem[tid] += vsmem[tid + 32];
    vsmem[tid] += vsmem[tid + 16];
    vsmem[tid] += vsmem[tid + 8];
    vsmem[tid] += vsmem[tid + 4];
    vsmem[tid] += vsmem[tid + 2];
    vsmem[tid] += vsmem[tid + 1];


    }
    if (tid == 0)
        mindata[j] = sdata[0];


}
__global__  void kernal_global_sum(datatype *f,datatype *logx,
                                   mindatatype *mindata,int *eps){
    unsigned int tid = threadIdx.x;
    unsigned int bid = blockIdx.x;
    unsigned int row = bid *THREAD + tid;
    unsigned int j = YCOL * row ;
    mindata[j] += mindata[j+8*THREAD] + mindata[j+16*THREAD]  + mindata[j+24*THREAD];
    f[row] +=  (iter_eps[ *eps]) *  (log(mindata[j])- logx[row] ) ;


}


__global__ void get_wdistance(datatype *f,datatype *g,mindatatype *kernal,
                              costtype *cost, int *eps){

    unsigned int tid = threadIdx.x;
    unsigned int bid = blockIdx.x;
    unsigned int col = tid + blockIdx.y * THREAD;

    unsigned int j = bid *YCOL + col;
    double a = f[bid] + g[col];
    kernal[j] = cost[j] * __expf((cost[j] - a) * iter_eps_[ *eps]);

}

__global__ void reduction1(mindatatype *kernal){
    unsigned int tid = threadIdx.x;
    unsigned int bid = blockIdx.x;
    unsigned int col = tid + blockIdx.y * THREAD * 8;

    unsigned int j = bid *YCOL + col;
    kernal[j] +=  kernal[j + THREAD] +
            kernal[j + 2*THREAD] + kernal[j + 3*THREAD] +
            kernal[j + 4*THREAD] + kernal[j + 5*THREAD] +
            kernal[j + 6*THREAD] + kernal[j + 7*THREAD];

   __syncthreads();
   if(tid < 64){
       kernal[j] += kernal[j+64];

   }
   __syncthreads();
   if(tid<32)
       {
           volatile mindatatype *vmem = kernal;
           vmem[j]+=vmem[j+32];
           vmem[j]+=vmem[j+16];
           vmem[j]+=vmem[j+8];
           vmem[j]+=vmem[j+4];
           vmem[j]+=vmem[j+2];
           vmem[j]+=vmem[j+1];
       }

}

__global__ void reduction2(mindatatype *kernal ){

    int bid = blockIdx.x;
    int j = YCOL * bid;
    kernal[j] += kernal[j+8*THREAD] + kernal[j+16*THREAD]  + kernal[j+24*THREAD];
}
__global__ void reduction3(mindatatype *kernal){
    int tid = threadIdx.x;
    int row = tid + blockIdx.x * THREAD;
    int j = row *YCOL ;

    // sum every col!!
    for (int i = blockDim.x/2; i > 0 ; i = i/2){
        if(tid < i){
            kernal[j] += kernal[j+i*YCOL];
        }
        __syncthreads();
    }


}


__global__  void reduction4(mindatatype *kernal){

    int tid = threadIdx.x;
    int j = tid * THREAD * YCOL  ;
    for (int i = blockDim.x/2; i > 0 ; i /= 2){
        if(tid < i){
            kernal[j] += kernal[j+i*THREAD* YCOL];
        }
        __syncthreads();
    }

    if (tid==0)
        printf("%.4f\n",kernal[0]);
}
#include<string>
#include<stdio.h>
int test_64x64(std::vector<double> &muXdat,std::vector<double>&muYdat){



    constexpr int nlist[20] = {1,0,1,0,0,0,1,1,1,1,1,50,50,50,50,30,20,10,10,10};


    // dimensions of grid
    int dim = 2;
    std::vector<int> muXdim = { N ,N };
    std::vector<int> muYdim = { N ,N };

    TDoubleMatrix *posX=GridToolsGetGridMatrix(dim, muXdim.data());
    TDoubleMatrix *posY=GridToolsGetGridMatrix(dim, muYdim.data());
    int xres=posX->dimensions[0];
    int yres=posY->dimensions[0];
    TCostFunctionProvider_Dynamic costFunctionProvider(
        &xres, &yres,
        &(posX->data), &(posY->data),
        1, dim);

    int fmemorysize = sizeof(datatype) * XROW;
    int gmemorysize = sizeof(datatype) * YCOL;
    int costmatrixsize = PROBLEMSIZE * (sizeof(double));
    double *c = (double *)malloc(costmatrixsize);
    costFunctionProvider.getCDense(c);
//    for (int i = 0 ; i < PROBLEMSIZE; ++i){
//        c[i] = sqrt(c[i]);

//    }


    datatype *logx = (datatype *)malloc(fmemorysize);
    datatype *logy = (datatype *)malloc(gmemorysize);
    for (int i = 0 ;i < xres ; ++i)
        logx[i] = (datatype) std::log(muXdat[i]);
    for (int i = 0 ;i < yres ; ++i)
        logy[i] = (datatype) std::log(muYdat[i]);

    datatype *f = (datatype *) malloc(fmemorysize);
    datatype *g = (datatype *) malloc(gmemorysize);


    costtype *cost_matrix_cuda;
    datatype *f_cuda;
    datatype *g_cuda;
    datatype *logx_cuda;
    datatype *logy_cuda;
    datatype *eps_cuda;
    mindatatype *mindata_matrix_cuda;

    int *iter_cuda;
    cudaMalloc((void**) &iter_cuda,sizeof(int));
//    cudaMalloc((void**) &kernal_cuda,kernelsize);
    cudaMalloc((void**) &cost_matrix_cuda,costmatrixsize);

    int minsize = PROBLEMSIZE * (sizeof(mindatatype));
    HANDLE_ERROR(cudaMalloc((void**) &mindata_matrix_cuda,minsize));

    cudaMalloc((void**) &f_cuda,fmemorysize);
    cudaMalloc((void**) &g_cuda,gmemorysize);
    cudaMalloc((void**) &logx_cuda,fmemorysize);
    cudaMalloc((void**) &logy_cuda,gmemorysize);
    cudaMalloc((void**) &eps_cuda,sizeof(datatype));

    cudaMemset(g_cuda, 0, gmemorysize);
    cudaMemset(f_cuda,0,fmemorysize);
//    cudaMemcpy(g_cuda,g,gmemorysize,cudaMemcpyHostToDevice);
//    cudaMemcpy(f_cuda,f,fmemorysize,cudaMemcpyHostToDevice);


    double eps[ITERSIZE];
    double eps_[ITERSIZE];
//    float *eps = (float *)malloc(sizeof(float)*ITERSIZE);
    for (int i= 0 ; i < ITERSIZE;++i){
        eps[i] = - ((1e6 - 1e-5) * exp(-i) + 1e-5);
        eps_[i] = -((float)1.0)/((1e6 - 1e-5) * exp(-i) + 1e-5);
    }


    HANDLE_ERROR(cudaMemcpy(cost_matrix_cuda,c,costmatrixsize,cudaMemcpyHostToDevice));

    HANDLE_ERROR(cudaMemcpy(logx_cuda,logx,fmemorysize,cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(logy_cuda,logy,gmemorysize,cudaMemcpyHostToDevice));

    HANDLE_ERROR(cudaMemcpyToSymbol(iter_eps, eps, sizeof(double)*ITERSIZE,0, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpyToSymbol(iter_eps_, eps_, sizeof(double)*ITERSIZE,0, cudaMemcpyHostToDevice));
    cublasHandle_t handle;
    cublasCreate(&handle);
    float *result = (float *)malloc(sizeof(float));

    unsigned int n = XROW*YCOL;
    int nstreams = 16;
    cudaStream_t *streams = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));

    for (int i = 0; i < nstreams; i++)
    {
        (cudaStreamCreate(&(streams[i])));
    }
    dim3 blockxx(XROW,YCOL/THREAD/8);
    double iStart, iElaps;


    iStart = seconds();
    cudaMemset(iter_cuda,-1,sizeof(int));
    for (int iter = 0 ; iter < ITERSIZE  ; ++iter)
    {


      iter11<<<1,1,0,streams[0]>>>(iter_cuda);
    for(int k=0;k< nlist[iter];k++)
    {
      kernal_row_min<<<blockxx,THREAD,0,streams[0]>>>(g_cuda,cost_matrix_cuda,mindata_matrix_cuda);

      kernal_global_row_min<<<XROW/THREAD,THREAD,0,streams[0]>>>(f_cuda,mindata_matrix_cuda);
      kernal_cal_newkernal<<<blockxx,THREAD,THREAD,streams[0]>>>(f_cuda,g_cuda,
                                           cost_matrix_cuda,mindata_matrix_cuda,iter_cuda);
      kernal_global_sum<<<XROW/THREAD,THREAD,0,streams[0]>>>(f_cuda,logx_cuda,mindata_matrix_cuda,iter_cuda);
      kernal_row_min<<<blockxx,THREAD,0,streams[0]>>>(f_cuda,cost_matrix_cuda,mindata_matrix_cuda);
      kernal_global_row_min<<<XROW/THREAD,THREAD,0,streams[0]>>>(g_cuda,mindata_matrix_cuda);
      kernal_cal_newkernal<<<blockxx,THREAD,THREAD,streams[0]>>>(g_cuda,f_cuda,
                                           cost_matrix_cuda,mindata_matrix_cuda,iter_cuda);
      kernal_global_sum<<<XROW/THREAD,THREAD,0,streams[0]>>>(g_cuda,logy_cuda,mindata_matrix_cuda,iter_cuda);

    }

    }
    iElaps = seconds() - iStart;
    printf("time=%f\n",iElaps);


    get_wdistance<<<dim3(XROW,YCOL/THREAD),THREAD,0,streams[0]>>>(f_cuda,g_cuda,mindata_matrix_cuda,cost_matrix_cuda,iter_cuda);
    cublasSasum(handle,n , mindata_matrix_cuda, 1, result);
    iElaps = seconds() - iStart;
    printf("time=%f\n",iElaps);


    std::cout << *result << std::endl;

    for (int i = 0; i < nstreams; i++)
    {
        (cudaStreamDestroy(streams[i]));
    }

    free(logx);
    free(logy);
    free(c);
    free(posX);
    free(posY);
    free(f);
    free(g);
//    cudaFree(kernal_cuda);
    cudaFree(cost_matrix_cuda);
    cudaFree (f_cuda);
    cudaFree (g_cuda);
    cudaFree (logx_cuda);
    cudaFree (logy_cuda);
    cudaFree (eps_cuda);
    cudaFree(mindata_matrix_cuda);
    cudaFree(iter_cuda);
}

int test_32x32(std::vector<double> &muXdat,std::vector<double>&muYdat){



    constexpr int nlist[20] = {1,0,1,0,0,0,1,1,1,1,1,50,50,50,50,30,20,10,10,10};


    // dimensions of grid
    int dim = 2;
    std::vector<int> muXdim = { N ,N };
    std::vector<int> muYdim = { N ,N };

    TDoubleMatrix *posX=GridToolsGetGridMatrix(dim, muXdim.data());
    TDoubleMatrix *posY=GridToolsGetGridMatrix(dim, muYdim.data());
    int xres=posX->dimensions[0];
    int yres=posY->dimensions[0];
    TCostFunctionProvider_Dynamic costFunctionProvider(
        &xres, &yres,
        &(posX->data), &(posY->data),
        1, dim);


    int fmemorysize = sizeof(datatype) * XROW;
    int gmemorysize = sizeof(datatype) * YCOL;
    int costmatrixsize = PROBLEMSIZE * (sizeof(costtype));
    costtype *c = (costtype *)malloc(costmatrixsize);
    costFunctionProvider.getCDense(c);

    datatype *logx = (datatype *)malloc(fmemorysize);
    datatype *logy = (datatype *)malloc(gmemorysize);
    for (int i = 0 ;i < xres ; ++i)
        logx[i] = (datatype) std::log(muXdat[i]);
    for (int i = 0 ;i < yres ; ++i)
        logy[i] = (datatype) std::log(muYdat[i]);

    datatype *f = (datatype *) malloc(fmemorysize);
    datatype *g = (datatype *) malloc(gmemorysize);


    costtype *cost_matrix_cuda;
    datatype *f_cuda;
    datatype *g_cuda;
    datatype *logx_cuda;
    datatype *logy_cuda;
    datatype *eps_cuda;
    mindatatype *mindata_matrix_cuda;

    int *iter_cuda;
    cudaMalloc((void**) &iter_cuda,sizeof(int));
//    cudaMalloc((void**) &kernal_cuda,kernelsize);
    cudaMalloc((void**) &cost_matrix_cuda,costmatrixsize);

    int minsize = PROBLEMSIZE * (sizeof(mindatatype));
    HANDLE_ERROR(cudaMalloc((void**) &mindata_matrix_cuda,minsize));

    cudaMalloc((void**) &f_cuda,fmemorysize);
    cudaMalloc((void**) &g_cuda,gmemorysize);
    cudaMalloc((void**) &logx_cuda,fmemorysize);
    cudaMalloc((void**) &logy_cuda,gmemorysize);
    cudaMalloc((void**) &eps_cuda,sizeof(datatype));

    cudaMemset(g_cuda, 0, gmemorysize);
    cudaMemset(f_cuda,0,fmemorysize);
//    cudaMemcpy(g_cuda,g,gmemorysize,cudaMemcpyHostToDevice);
//    cudaMemcpy(f_cuda,f,fmemorysize,cudaMemcpyHostToDevice);

    cublasHandle_t handle;
    cublasCreate(&handle);
    float *result = (float *)malloc(sizeof(float));
    unsigned int n = XROW*YCOL;

    double eps[ITERSIZE];
    double eps_[ITERSIZE];
//    float *eps = (float *)malloc(sizeof(float)*ITERSIZE);
    for (int i= 0 ; i < ITERSIZE;++i){
        eps[i] = - ((1e6 - 1e-5) * exp(-i) + 1e-5);
        eps_[i] = -((float)1.0)/((1e6 - 1e-5) * exp(-i) + 1e-5);
    }


    HANDLE_ERROR(cudaMemcpy(cost_matrix_cuda,c,costmatrixsize,cudaMemcpyHostToDevice));

    HANDLE_ERROR(cudaMemcpy(logx_cuda,logx,fmemorysize,cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(logy_cuda,logy,gmemorysize,cudaMemcpyHostToDevice));

    HANDLE_ERROR(cudaMemcpyToSymbol(iter_eps, eps, sizeof(double)*ITERSIZE,0, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpyToSymbol(iter_eps_, eps_, sizeof(double)*ITERSIZE,0, cudaMemcpyHostToDevice));

    int nstreams = 16;
    cudaStream_t *streams = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));

    for (int i = 0; i < nstreams; i++)
    {
        (cudaStreamCreate(&(streams[i])));
    }
    dim3 blockxx(XROW,YCOL/THREAD/8);
    double iStart, iElaps;


    iStart = seconds();
    cudaMemset(iter_cuda,-1,sizeof(int));
    for (int iter = 0 ; iter < ITERSIZE  ; ++iter)
    {


      iter11<<<1,1,0,streams[0]>>>(iter_cuda);
    for(int k=0;k< nlist[iter];k++)
    {
      kernal_row_min<<<blockxx,THREAD,THREAD,streams[0]>>>(g_cuda,f_cuda,cost_matrix_cuda,
                                                      logx_cuda,iter_cuda);


      kernal_row_min<<<blockxx,THREAD,THREAD,streams[0]>>>(f_cuda,g_cuda,cost_matrix_cuda,
                                                      logy_cuda,iter_cuda);

    }

    }
    iElaps = seconds() - iStart;
    printf("time=%f\n",iElaps);
    iStart = seconds();



    get_wdistance<<<dim3(XROW,YCOL/THREAD),THREAD,0,streams[0]>>>(f_cuda,g_cuda,mindata_matrix_cuda,cost_matrix_cuda,iter_cuda);
    cublasSasum(handle,n , mindata_matrix_cuda, 1, result);

    std::cout << *result << std::endl;

    for (int i = 0; i < nstreams; i++)
    {
        (cudaStreamDestroy(streams[i]));
    }

    free(logx);
    free(logy);
    free(c);

    free(posX);
    free(posY);
    free(f);
    free(g);
//    cudaFree(kernal_cuda);
    cudaFree(cost_matrix_cuda);
    cudaFree (f_cuda);
    cudaFree (g_cuda);
    cudaFree (logx_cuda);
    cudaFree (logy_cuda);
    cudaFree (eps_cuda);
    cudaFree(mindata_matrix_cuda);
    cudaFree(iter_cuda);
}
