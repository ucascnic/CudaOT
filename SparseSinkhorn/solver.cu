#ifndef MY_SK_H
#define MY_SK_H
#include<vector>
#include <cusparse_v2.h>
#include<cuda_runtime_api.h>
#include<cublas_v2.h>
#include <math.h>
#include <algorithm>
#include<stdio.h>
#include<iostream>
#include<string>
#include<fstream>   
#include<sstream>
#include<vector>
 
#include"basic_settings.h"
#include<THierarchicalPartition.h>
#include"data_struct.h"
#include"family.h"
#include<cuda_runtime_api.h>
#include<cusparse.h>
#include"eps_handler.h"
#include"memory_manage.h"
#include"layer_manage.h"
#include"kernel_manage.h"
#include"solver.h"
static cusparseHandle_t cusparse_handle = 0;
cusparseStatus_t Status_t = cusparseCreate(&cusparse_handle);
static void init_cusparse() {
  if (cusparse_handle == 0) {
    if (Status_t != CUSPARSE_STATUS_SUCCESS) {
      printf("CUSPARSE Library initialization failed");
    }
  }
}


template<int MAXELEMENT>
__global__ void record_cost(XData *xPos,XData *yPos,
                                   double *res,double *u,double *v,
                                   int N,COO_Sparese_Matrix<MAXELEMENT> *kernel,int layer) {

    unsigned int tid = threadIdx.x;
    unsigned int n = tid + blockIdx.x * blockDim.x;
    if (n>=N) return;
    int x = kernel->row_ptrl[n];
    int y = kernel->col_ptrl[n];

    double result=EUCL_lincombSqr(&xPos->data[xPos->ind[layer]+(x*DIM)],
                                  &yPos->data[yPos->ind[layer]+(y*DIM)]);
    result = pow(result,(double)p_index/2.0);
    res[n] = kernel->val[n] * result * u[x]* v[y];


}

template<int MAXELEMENT>
double record_value(int N,XData *cu_xposH, XData *cu_yposH,
                        double *res,double *u,double *v,
                        COO_Sparese_Matrix<MAXELEMENT> *kernel,int layer,cublasHandle_t handle){

    double result = 0.0;
    int blocks = N/128+1;
 
    record_cost<MAXELEMENT><<<blocks,128>>>(cu_xposH,cu_yposH,res,u,v,
                      N,kernel,layer);
    cublasDasum(handle,N, res, 1, &result);
    return result;
}
/////////////////////////////////////////////////////////////////////////////////
template< int TOTALNUM_, int MAXLAYER_>
__global__ void signal_refine_doublex_(  Familys<TOTALNUM_,MAXLAYER_> *familyX,
                                         XData *alphaH,int lTop) {
    int newlayer = lTop + 1;
    double *signal = &alphaH->data[ alphaH->ind[lTop]];
    double *signalFine = &alphaH->data[ alphaH->ind[newlayer]];

    int x,xFine;
    int nchildren_ind = get_n_children_ind(lTop);
    for(x=0;x < familyX->nCells[lTop]; x++) {
        int total = familyX->nChildren[nchildren_ind + x];
        for(int k=0;k < total; k++) {
            xFine=familyX->children[familyX->nChildren_ind[nchildren_ind + x] +  k];
            signalFine[xFine]=signal[x];

        }
    }
}

template< int TOTALNUM_, int MAXLAYER_>
__global__ void signal_refine_doubley_( Familys<TOTALNUM_,MAXLAYER_> *familyX,XData *alphaH,int lTop) {
    int newlayer = lTop + 1;
    double *signal = &alphaH->data[ alphaH->ind[lTop]];
    double *signalFine = &alphaH->data[ alphaH->ind[newlayer]];

    int x,xFine;
    int nchildren_ind = get_n_children_ind(lTop);
    for(x=0;x < familyX->nCells[lTop];x++) {
        int total = familyX->nChildren[nchildren_ind + x];
        for(int k=0; k<total; k++) {
            xFine=familyX->children[familyX->nChildren_ind[nchildren_ind + x] +  k];
            signalFine[xFine]=signal[x];

        }
    }
}
template<int TOTALNUM_, int MAXLAYER_>
__global__ void signal_propagate_doublex_( Familys<TOTALNUM_,MAXLAYER_> *familyX,
                                           XData *alphaH,int lTop,
                                         int lBottom, int mode){
    double newValue,value;
    int tid = threadIdx.x;

    for(int i=lBottom-1;i>=lTop;i--) {
          int nchildren_ind = get_n_children_ind(i);
          for(int j=tid ;j < familyX->nCells[i]; j += blockDim.x) {

            value=alphaH->data[alphaH->ind[i-lTop+1] +
                    familyX->children[familyX->nChildren_ind[nchildren_ind + j] +  0]];

            for(int k=1; k< familyX->nChildren[nchildren_ind + j];k++) {
                 newValue= alphaH->data[ alphaH->ind[i-lTop+1] +
                     familyX->children[familyX->nChildren_ind[nchildren_ind + j] +  k]];

                if( ((mode==MODE_MAX_) && (newValue>value)) ||
                        ((mode==MODE_MIN_) && (newValue<value))) {
                    value=newValue;

                 }
            }
            alphaH->data[  alphaH->ind[i - lTop] +j  ]=value;
        }
    }
}
template<int TOTALNUM_, int MAXLAYER_>
__global__ void signal_propagate_doubley_( Familys<TOTALNUM_,MAXLAYER_> *familyX,XData *alphaH,int lTop,
                                        int lBottom, int mode){
    double newValue,value;
    int tid = threadIdx.x;
    for(int i=lBottom-1;i>=lTop;i--) {
          int nchildren_ind = get_n_children_ind(i);
          for(int j=tid ;j < familyX->nCells[i]; j += blockDim.x)  {

            value=alphaH->data[alphaH->ind[i-lTop+1] +
                    familyX->children[familyX->nChildren_ind[nchildren_ind + j] +  0]];

            for(int k=1;k< familyX->nChildren[nchildren_ind + j];k++) {
                 newValue= alphaH->data[ alphaH->ind[i-lTop+1] +
                     familyX->children[familyX->nChildren_ind[nchildren_ind + j] +  k]];

                if( ((mode==MODE_MAX_) && (newValue>value)) ||
                        ((mode==MODE_MIN_) && (newValue<value))) {
                    value=newValue;

                 }
            }
            alphaH->data[  alphaH->ind[i-lTop] +j  ]=value;
        }
    }
}
template< int TOTALNUM_, int MAXLAYER_>
void signal_refine_all(Familys<TOTALNUM_,MAXLAYER_> *cu_familyX,
                       Familys<TOTALNUM_,MAXLAYER_> *cu_familyY,
                       XData *cu_alphaH,
                       XData *cu_betaH, int lTop){

    signal_refine_doublex_<TOTALNUM_,MAXLAYER_><<<1,1>>>(cu_familyX,cu_alphaH,lTop);
    signal_refine_doubley_<TOTALNUM_,MAXLAYER_><<<1,1>>>(cu_familyY,cu_betaH,lTop);
}

#ifdef  _RECORD_ 
std::vector<double> optimal_values;
std::vector<int> kernel_entries;
std::vector<int> itertime;
#include <ctime>

#endif


template< int ROW_,  int COL_, int MAXLAYER_, int MAXELEMENT>
int  solveSingle_(int *msg_cu,int *kernelValid,
                 COO_Sparese_Matrix<MAXELEMENT> *kernel,
                 COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                 CuSinkhornSolverParameters *cfg,
                 double eps,
                 double *muX_cuda,
                 double *muY_cuda,
                 Familys<ROW_*COL_,MAXLAYER_> *familyX,
                 Familys<ROW_*COL_,MAXLAYER_> *familyY,
                 int layerBottom,
                 XData * alphaH, XData *betaH,
                 XData *xpos, XData *ypos,
                 XData *rax,XData *ray,
                 int xres,int yres,
                 double *u,double *v,
                 double *res,int layer,
                 Ind_Ptr_ *indptrl_cu_mid,
                 Ind_Ptr_ *indptrl_cu,
                 COO_Buff * buff,
                 cublasHandle_t handle)
{
    int nIterations=0;
    int nAbsorptionLoops=0;
    int msg = 0;
    int N ;

    double error = 0.0;
    cudaMemset(msg_cu,0,sizeof(int));
  
    if(!*kernelValid) {
         if(layer < layerBottom  ){
         // refine the whole region 
           msg = refinekernel_2<MAXELEMENT>(kernel,kernelT,
                                  muX_cuda,
                                  muY_cuda,
                                  layerBottom,
                                  eps,xres,yres,
                                  layer,
                                  alphaH,betaH,xpos,ypos,rax,ray,
                                 indptrl_cu_mid,indptrl_cu,
                                   kernelValid,buff,cfg);
                                   }
        else
        {
             
            msg = refinekernel_3<MAXELEMENT>(kernel,kernelT,
                              muX_cuda,
                              muY_cuda,
                              layerBottom,
                              eps,xres,yres,
                              layer,
                              alphaH,betaH,xpos,ypos,rax,ray,
                             indptrl_cu_mid,indptrl_cu,
                               kernelValid,buff,cfg);

        }

         *kernelValid=1;
        if(msg!=0) return msg;
    }

    cudaMemcpy(&N,&kernel->nonZeros,sizeof(int),cudaMemcpyDeviceToHost);
    printf("\t\tkernel entries: %ld\n",N);

    while(true) {


        iterate(msg_cu,cfg->innerIterations,
                         xres,muX_cuda,yres,muY_cuda,kernel,
                         kernelT,u,v,res
                         );
        // record the optimal value we get
        #ifdef  _RECORD_ 
             cudaMemcpy(&N,&kernel->nonZeros,sizeof(int),cudaMemcpyDeviceToHost);
             kernel_entries.push_back(N);
             double tmp = record_value(N,xpos,ypos,res,u,v,kernel,layer,handle);
             optimal_values.push_back(tmp); 
        #endif


        cudaMemcpy(&msg,msg_cu,sizeof(int),cudaMemcpyDeviceToHost);

        if(msg== MSG_NANSCALING ) return msg;
        checkAbsorb<<<xres/128+1,128>>>(msg_cu,u,v,cfg->absorption_scalingBound,xres,yres);
        cudaMemcpy(&msg,msg_cu,sizeof(int),cudaMemcpyDeviceToHost);

        if(msg == MSG_ABSORB_REITERATE) {

            cudaMemset(msg_cu,0,sizeof(int));
            printf("\t\tabsorbing\n");
            nAbsorptionLoops++;
            if(nAbsorptionLoops>cfg->maxAbsorptionLoops) {
                return MSG_ABSORB_TOOMANYABSORPTIONS;
            }
            SinkhornAbsorbScaling1<<<xres/128+1,128>>>(alphaH, u, layer, eps,xres);
            signal_propagate_doublex_<ROW_*COL_,MAXLAYER_><<<1,128>>>(familyX, alphaH, 0, layer, THierarchicalPartition::MODE_MAX);


            SinkhornAbsorbScaling2<<<yres/128+1,128>>>(betaH, v, layer, eps,yres);
            signal_propagate_doubley_<ROW_*COL_,MAXLAYER_><<<1,128>>>(familyY, betaH, 0, layer, THierarchicalPartition::MODE_MAX);
             
            msg = refinekernel_3<MAXELEMENT>(kernel,kernelT,
                              muX_cuda,
                              muY_cuda,
                              layerBottom,
                              eps,xres,yres,
                              layer,
                              alphaH,betaH,xpos,ypos,rax,ray,
                             indptrl_cu_mid,indptrl_cu,
                               kernelValid,buff,cfg);
 


            cudaMemcpy(&N,&kernel->nonZeros,sizeof(int),cudaMemcpyDeviceToHost);
            printf("\t\tkernel entries: %ld\n",N);


            if(msg!=0) return msg;
            *kernelValid=1;
            continue;
        } else {
            if(msg!=0) return msg;
        }
        nAbsorptionLoops=0;

        error = getError(kernel,
                 muX_cuda,
                 u,v,res,xres,yres,handle);
        if(isnan(error)){
            return -2 ;
        }
        if(error<=cfg->maxError) {
            cudaMemset(msg_cu,0,sizeof(int));
            checkAbsorb<<<xres/128+1,128>>>(msg_cu,u,v,cfg->absorption_scalingLowerBound,xres,yres);
            cudaMemcpy(&msg,msg_cu,sizeof(int),cudaMemcpyDeviceToHost);

            if(msg==MSG_ABSORB_REITERATE) {
                printf("\t\tsafety absorption.\n");

                SinkhornAbsorbScaling1<<<xres/128+1,128>>>(alphaH, u, layer, eps,xres);
                signal_propagate_doublex_<ROW_*COL_,MAXLAYER_><<<1,128>>>(familyY, alphaH, 0, layer, THierarchicalPartition::MODE_MAX);
                SinkhornAbsorbScaling2<<<yres/128+1,128>>>(betaH, v, layer, eps,yres);
                signal_propagate_doubley_<ROW_*COL_,MAXLAYER_><<<1,128>>>(familyY, betaH, 0, layer, THierarchicalPartition::MODE_MAX);

                if(layer < layerBottom)
                    msg=generateKernel_<MAXELEMENT>(kernel,kernelT,cfg,eps,
                                       muX_cuda,muY_cuda,layer,
                                       alphaH,betaH,xpos,ypos,rax,ray,xres,yres,
                                       indptrl_cu_mid,indptrl_cu,buff);
                else{
                    msg = refinekernel_2<MAXELEMENT>(kernel,kernelT,
                                      muX_cuda,
                                      muY_cuda,
                                      layerBottom,
                                      eps,xres,yres,
                                      layer,
                                      alphaH,betaH,xpos,ypos,rax,ray,
                                     indptrl_cu_mid,indptrl_cu,
                                       kernelValid,buff,cfg);

                }

                cudaMemcpy(&N,&kernel->nonZeros,sizeof(int),cudaMemcpyDeviceToHost);
                printf("\t\tkernel entries: %ld\n",N);
                if(msg!=0) return msg;
                *kernelValid=1;

                continue;

            } else {
                // otherwise, check for other errors
                if(msg!=0) return msg;
                printf("\t\tfinal absorption\n");
                SinkhornAbsorbScaling1<<<xres/128+1,128>>>(alphaH, u, layer, eps,xres);
                signal_propagate_doublex_<ROW_*COL_,MAXLAYER_><<<1,128>>>(familyY, alphaH, 0, layer, THierarchicalPartition::MODE_MAX);
                SinkhornAbsorbScaling2<<<yres/128+1,128>>>(betaH, v, layer, eps,yres);
                signal_propagate_doubley_<ROW_*COL_,MAXLAYER_><<<1,128>>>(familyY, betaH, 0, layer, THierarchicalPartition::MODE_MAX);
                *kernelValid=0;
                return 0;
            }
        }
        nIterations+=cfg->innerIterations;
        if(nIterations>cfg->maxIterations) {
            return MSG_EXCEEDMAXITERATIONS;
        }
    }
}





template<unsigned int TOTALNUM_, int MAXLAYER_,int MAXELEMENT>
__global__ void refine_kernel_on_gpu(COO_Sparese_Matrix<MAXELEMENT> *newkernel,
                             COO_Sparese_Matrix<MAXELEMENT> *oldKernel,
                             double *muX,
                             double *muY,
                             int layerBottom,
                             double _eps,
                             int layer,
                             Familys<TOTALNUM_,MAXLAYER_> *familyX,
                             Familys<TOTALNUM_,MAXLAYER_> *familyY,
                             XData *alpha,XData *beta,
                             XData * xPos, XData *yPos,
                             XData *xRadii, XData *yRadii
                             ) {

    int kk =  threadIdx.x + blockIdx.x * blockDim.x;
    if (kk >= oldKernel->nonZeros)return;

        int xOld = oldKernel->row_ptrl[kk];
        int yOld= oldKernel->col_ptrl[kk];
         int nchildren_ind = get_n_children_ind(layer - 1);
        int nChildrenX=familyX->nChildren[nchildren_ind + xOld];
        int nChildrenY=familyY->nChildren[nchildren_ind + yOld];
        int *childrenX=&familyX->children[familyX->nChildren_ind[nchildren_ind + xOld]];
        int *childrenY=&familyY->children[familyY->nChildren_ind[nchildren_ind + yOld]];
        for(int iX=0;iX<nChildrenX;iX++) {
                for(int iY=0;iY<nChildrenY;iY++) {
                    int x = childrenX[iX];
                    int y = childrenY[iY];
                    double value = getCostEff(xPos,yPos,xRadii,yRadii,layer,layer,x,layer,y);
                    value -= (alpha->data[alpha->ind[layer] + x] + beta->data[beta->ind[layer] +y]);

                        double _value= exp( - value / _eps );
                        _value *= muX[x]*muY[y];
                        int rr = atomicAdd(&newkernel->nonZeros, 1);
                        newkernel->row_ptrl[rr] = x;
                        newkernel->col_ptrl[rr] = y;
                        newkernel->val[rr] = _value;
                }
            }


}


template<int ROW_, unsigned int COL_, int MAXLAYER_ ,int MAXELEMENT>
void refinekernel_(COO_Sparese_Matrix<MAXELEMENT> *kernel,COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                  double *muX_cuda,
                  double *muY_cuda,
                  int layerBottom,
                  double _eps,int xres,int yres,
                  int layer,
                  Familys<ROW_*COL_,MAXLAYER_> *familyX,
                   Familys<ROW_*COL_,MAXLAYER_> *familyY,
                  XData *cu_alphaH,XData *cu_betaH,
                  XData * cu_xposH, XData *cu_yposH,
                  XData *cu_xradii, XData *cu_yradii,
                 Ind_Ptr_ *indptrl_cu_mid,Ind_Ptr_ *indptrl_cu,
                   int *kernelValid, COO_Buff * buff){

    int N = 0;
    cudaMemcpy(&N,&kernel->nonZeros,sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(kernelT,kernel,sizeof(COO_Sparese_Matrix<MAXELEMENT>),cudaMemcpyDeviceToDevice);
    cudaMemset(&kernel->nonZeros,0,sizeof(unsigned int));
    refine_kernel_on_gpu<ROW_*COL_,MAXLAYER_,MAXELEMENT><<<N/128+1,128>>>(kernel,kernelT,
                  muX_cuda,
                  muY_cuda,
                  layerBottom,
                  _eps,layer,
                  familyX,familyY,
                  cu_alphaH,cu_betaH,
                  cu_xposH,cu_yposH,
                  cu_xradii,cu_yradii);
//    show_kenel_sum<<<1,1>>>(kernel);
    cudaMemcpy(&N,&kernel->nonZeros,sizeof(int),cudaMemcpyDeviceToHost);
    int blocks = N/128+1;
    transport_kernel<<<blocks,128>>>(kernel,kernelT);
    cudaMemcpy(&kernelT->nonZeros,&kernel->nonZeros,sizeof(int),cudaMemcpyDeviceToDevice);
    printf("\t\tkernel entries: (refinement) %d\n",N);
    sortcoo_and_to_csr<MAXELEMENT>(xres,yres,N, kernel,indptrl_cu_mid,indptrl_cu,buff);
    sortcoo_and_to_csr<MAXELEMENT>(yres,xres,N, kernelT,indptrl_cu_mid,indptrl_cu,buff);

    *kernelValid=1;



}


#include<assert.h>




template<int MAXELEMENT>
void sortcoo_and_to_csr(int xres,int yres,int N, COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                        Ind_Ptr_ *indptrl_cu_mid,Ind_Ptr_ *indptrl_cu,
                         COO_Buff * buff){





    size_t workspace_size = 0;
    Status_t = cusparseXcoosort_bufferSizeExt(
        cusparse_handle,
        xres, yres,
        N,
        kernelT->row_ptrl,
        kernelT->col_ptrl,
        &workspace_size);
    assert( CUSPARSE_STATUS_SUCCESS == Status_t);
    double * buffer_for_coo_sort;

//    printf("i need for %d  size\n",workspace_size);
    CHECK(cudaMalloc(&buffer_for_coo_sort, sizeof(char) * workspace_size));

    cudaMemcpy(indptrl_cu,indptrl_cu_mid,sizeof(Ind_Ptr_) * MAXELEMENT,cudaMemcpyDeviceToDevice);


    Status_t = cusparseXcoosortByRow(
        cusparse_handle,
        xres, yres,
        N,
        kernelT->row_ptrl,
        kernelT->col_ptrl,
        indptrl_cu,
        buffer_for_coo_sort);



   Status_t = cusparseDgthr(cusparse_handle,
                          N,
                          kernelT->val,
                          buffer_for_coo_sort,
                          indptrl_cu,
                          CUSPARSE_INDEX_BASE_ZERO);


   cudaMemcpy(kernelT->val,buffer_for_coo_sort,
              N*sizeof(double),cudaMemcpyDeviceToDevice);
   assert( CUSPARSE_STATUS_SUCCESS == Status_t);

    Status_t = cusparseXcoo2csr(cusparse_handle,
       kernelT->row_ptrl, N, yres,
       kernelT->crs_row_ptrl, CUSPARSE_INDEX_BASE_ZERO);
    assert( CUSPARSE_STATUS_SUCCESS == Status_t);

    cudaFree(buffer_for_coo_sort);




}

template<int MAXELEMENT>
int generateKernel_(COO_Sparese_Matrix<MAXELEMENT> *kernel, COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                   CuSinkhornSolverParameters *cfg,double eps,
                   double *muX_cuda,
                   double *muY_cuda,int layer,
                    XData * alphaH,
                   XData *betaH,XData * xpos, XData *ypos,
                   XData *rax, XData *ray,int xres,int yres,
                   Ind_Ptr_ *indptrl_cu_mid,Ind_Ptr_ *indptrl_cu,
                    COO_Buff * buff) {



    CHECK(cudaMemset(&kernel->nonZeros,0,sizeof(unsigned int)));
    dim3 block(xres,yres/128/16 + 1);

    printf("check cell\n");

    double eps1 = (double )1.0/eps;
    check_cells_<<<block,128>>>(layer,layer,
                              eps1,eps*log(cfg->truncation_thresh),
                               muX_cuda,muY_cuda,kernel,
                              alphaH, betaH, xpos, ypos,rax,ray,xres,yres);


    int N = 0;
    cudaMemcpy(&N,&kernel->nonZeros,sizeof(int),cudaMemcpyDeviceToHost);
    int blocks = N/128+1;
    transport_kernel<<<blocks,128>>>(kernel,kernelT);

    sortcoo_and_to_csr<MAXELEMENT>(xres,yres,N, kernel,indptrl_cu_mid,indptrl_cu,buff);
    sortcoo_and_to_csr<MAXELEMENT>(yres,xres,N, kernelT,indptrl_cu_mid,indptrl_cu,buff);


    return 0;
}


template <int MAXELEMENT>
void sortcoo_and_to_csr2(int xres,int yres,int N, COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                        Ind_Ptr_ *indptrl_cu_mid, Ind_Ptr_ *indptrl_cu,
                         COO_Buff * buff){




    size_t workspace_size = 0;
    Status_t = cusparseXcoosort_bufferSizeExt(
        cusparse_handle,
        xres, yres,
        N,
        kernelT->row_ptrl,
        kernelT->col_ptrl,
        &workspace_size);
    assert( CUSPARSE_STATUS_SUCCESS == Status_t);
    double * buffer_for_coo_sort;


    CHECK(cudaMalloc(&buffer_for_coo_sort, sizeof(char) * workspace_size));


    cudaMemcpy(indptrl_cu,indptrl_cu_mid,sizeof(Ind_Ptr_)*MAXELEMENT,cudaMemcpyDeviceToDevice);
    Status_t = cusparseXcoosortByRow(
        cusparse_handle,
        xres, yres,
        N,
        kernelT->row_ptrl,
        kernelT->col_ptrl,
        indptrl_cu,
        buffer_for_coo_sort);
   Status_t = cusparseDgthr(cusparse_handle,
                          N,
                          kernelT->val,
                          buffer_for_coo_sort,
                          indptrl_cu,
                          CUSPARSE_INDEX_BASE_ZERO);
   cudaMemcpy(kernelT->val,buffer_for_coo_sort,
              N*sizeof(double),cudaMemcpyDeviceToDevice);
   assert( CUSPARSE_STATUS_SUCCESS == Status_t);
    cudaFree(buffer_for_coo_sort);


}
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>

template <int MAXELEMENT>
int refinekernel_3(COO_Sparese_Matrix<MAXELEMENT> *kernel, COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                  double *muX_cuda,
                  double *muY_cuda,
                  int layerBottom,
                  double _eps, int xres, int yres,
                  int layer,
                  XData *cu_alphaH, XData *cu_betaH,
                  XData * cu_xposH, XData *cu_yposH,
                  XData *cu_xradii, XData *cu_yradii,
                    Ind_Ptr_ *indptrl_cu_mid, Ind_Ptr_ *indptrl_cu,
                   int *kernelValid,  COO_Buff * buff,
                    CuSinkhornSolverParameters *cfg){

 
    _eps = (double )1.0/_eps;

 
 
    generate_kernel_32_wrap<MAXELEMENT><<<xres,WARP_>>>(layer,layerBottom,_eps,
                                    0.0,
                                    muX_cuda,muY_cuda,
                                     kernel,
                                    cu_alphaH,cu_betaH,
                                   cu_xposH, cu_yposH,
                                    cu_xradii, cu_yradii,
                                     xres, yres); 
    generate_kernel_32_wrap<MAXELEMENT><<<yres,WARP_>>>(layer,layerBottom,_eps,
                                    0.0,
                                    muY_cuda,muX_cuda,
                                     kernelT,
                                    cu_betaH,cu_alphaH,
                                    cu_yposH,cu_xposH, 
                                    cu_yradii,cu_xradii, 
                                     yres,xres); 


    *kernelValid=1;
    return 0;



}



template <int MAXELEMENT>
int refinekernel_2(COO_Sparese_Matrix<MAXELEMENT> *kernel, COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                  double *muX_cuda,
                  double *muY_cuda,
                  int layerBottom,
                  double _eps, int xres, int yres,
                  int layer,
                  XData *cu_alphaH, XData *cu_betaH,
                  XData * cu_xposH, XData *cu_yposH,
                  XData *cu_xradii, XData *cu_yradii,
                 Ind_Ptr_ *indptrl_cu_mid, Ind_Ptr_ *indptrl_cu,
                   int *kernelValid,  COO_Buff * buff,
                    CuSinkhornSolverParameters *cfg){

    double slack = _eps * log(cfg->truncation_thresh);
    _eps = (double )1.0/_eps;

    int N = 0;
    cudaMemcpy(&N,&kernel->nonZeros,sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(kernelT,kernel,sizeof(COO_Sparese_Matrix<MAXELEMENT>),cudaMemcpyDeviceToDevice);
    cudaMemset(&kernel->nonZeros,0,sizeof(unsigned int));

    int blocks = N/128+1;
    int origionN = N;
    refine_kernel_2<<<N/128/8+1,128>>>(kernel,kernelT,
                  muX_cuda,
                  muY_cuda,
                  layerBottom,
                  _eps,layer,
                  cu_alphaH,cu_betaH,
                  cu_xposH,cu_yposH,
                  cu_xradii,cu_yradii,slack);
    int init = 0;

    // binary operation used to reduce values
    thrust::plus<int> binary_op;

    // compute sum on the device
    thrust::device_ptr<int> dev_ptr(&kernelT->row_ptrl[0]);
    thrust::device_ptr<int> last  = dev_ptr + origionN;
    N = thrust::reduce(dev_ptr, last, init, binary_op);

    cudaMemcpy(&kernel->nonZeros,&N,sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(&kernelT->nonZeros,&kernel->nonZeros,sizeof(int),cudaMemcpyDeviceToDevice);
    transport_kernel_2<<<blocks,128>>>(kernel,kernelT,origionN);
    sortcoo_and_to_csr2<MAXELEMENT>(xres,yres,origionN, kernel,
                            indptrl_cu_mid,indptrl_cu,
                            buff);
    sortcoo_and_to_csr2<MAXELEMENT>(yres,xres,origionN, kernelT,
                            indptrl_cu_mid,indptrl_cu,
                            buff);

    Status_t = cusparseXcoo2csr(cusparse_handle,
       kernel->row_ptrl, N, yres,
       kernel->crs_row_ptrl, CUSPARSE_INDEX_BASE_ZERO);
    assert( CUSPARSE_STATUS_SUCCESS == Status_t);
    Status_t = cusparseXcoo2csr(cusparse_handle,
       kernelT->row_ptrl, N, xres,
       kernelT->crs_row_ptrl, CUSPARSE_INDEX_BASE_ZERO);
    assert( CUSPARSE_STATUS_SUCCESS == Status_t);



    *kernelValid=1;
    return 0;



}


template<int MAXELEMENT>
__global__ void scoreTransportCost(XData *xPos,XData *yPos,
                                   double *res,double *u,double *v,
                                   int N,COO_Sparese_Matrix<MAXELEMENT> *kernel,int layer) {

    unsigned int tid = threadIdx.x;
    unsigned int n = tid + blockIdx.x * blockDim.x;
    if (n>=N) return;
    int x = kernel->row_ptrl[n];
    int y = kernel->col_ptrl[n];

    double result=EUCL_lincombSqr(&xPos->data[xPos->ind[layer]+(x*DIM)],
                                  &yPos->data[yPos->ind[layer]+(y*DIM)]);
    result = pow(result,(double)p_index/2.0);
    kernel->val[n] *= result*u[x]*v[y];


}
template<int MAXELEMENT>
__global__ void check_cells_(int layer,int layerBottom,
                                double eps,
                                double slack,
                                double *muX,double *muY,
                                COO_Sparese_Matrix<MAXELEMENT> *kernel,
                                XData *alpha,XData *beta,
                                XData * xPos, XData *yPos,
                                XData *xRadii, XData *yRadii,
                             int xres,int yres){
    double value = 0.0;
    unsigned int tid = threadIdx.x ;
    unsigned int x = blockIdx.x ;
    unsigned int y = tid   + blockIdx.y * blockDim.x * 16 ;


    if( x >= xres)
        return;


    int totalnozeros = 0;


    double vals[16];
    int row[16];
    int col[16];
    int indx = alpha->ind[layer];
    int indy = beta->ind[layer];
    for (int n = 0 ; n < 16; n++)
    {


        if( y >= yres)
            break;
        value =  (alpha->data[indx + x] + beta->data[ indy + y]) - getCostEff(xPos,yPos,xRadii,yRadii,layerBottom,layer,x,layer,y);
        if( value > slack) {
                double _value= exp( value * eps )*muX[x]*muY[y];

//                unsigned int n = atomicAdd(&kernel->nonZeros, 1);
//                kernel->row_ptrl[n] = x;
//                kernel->col_ptrl[n] = y;
//                kernel->val[n] = _value;
//                unsigned int n = atomicAdd(&kernel->nonZeros, 1);
                vals[totalnozeros] = _value;
                row[totalnozeros] = x;
                col[totalnozeros] = y;
                totalnozeros++;
        }
        y += 128;
    }

    if (totalnozeros){
    unsigned int rr = atomicAdd(&kernel->nonZeros, totalnozeros);

            for (int i = 0 ; i < totalnozeros ;++i){
            kernel->row_ptrl[rr+i] = row[i];
            kernel->col_ptrl[rr+i] = col[i];
            kernel->val[rr+i] = vals[i];

    }
    }


}


#include<float.h>
template<int MAXELEMENT>
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
                             ) {

    int kk =  threadIdx.x + blockIdx.x * blockDim.x * 8;
    int ind_alpha = alpha->ind[layer];
    int ind_beta = beta->ind[layer];
    for (int i = 0 ; i < 8 ;++i){
    if (kk >= oldKernel->nonZeros)return;

        int x = oldKernel->row_ptrl[kk];
        int y= oldKernel->col_ptrl[kk];


        double value = (alpha->data[ind_alpha + x] + beta->data[ind_beta +y]) -
                getCostEff(xPos,yPos,xRadii,yRadii,layer,layer,x,layer,y);

        if( value > slack){
            double _value= exp(  value * _eps )* muX[x]*muY[y];
//            int rr = atomicAdd(&newkernel->nonZeros, 1);
            oldKernel->row_ptrl[kk] = 1;

            newkernel->row_ptrl[kk] = x;
            newkernel->col_ptrl[kk] = y;
            newkernel->val[kk] = _value;

        }else{
            oldKernel->row_ptrl[kk] = 0;
            newkernel->row_ptrl[kk] = DBL_MAX;
            newkernel->col_ptrl[kk] = DBL_MAX;

        }

        kk += blockDim.x;
    }

}

template<unsigned int TOTALNUM_, int MAXLAYER_ ,int MAXELEMENT>
__global__ void refine_kernel_(COO_Sparese_Matrix<MAXELEMENT> *newkernel,
                             COO_Sparese_Matrix<MAXELEMENT> *oldKernel,
                             double *muX,
                             double *muY,
                             int layerBottom,
                             double _eps,
                             int layer,
                             Familys<TOTALNUM_,MAXLAYER_> *familyX,
                             Familys<TOTALNUM_,MAXLAYER_> *familyY,
                             XData *alpha,XData *beta,
                             XData *xPos, XData *yPos,
                             XData *xRadii, XData *yRadii
                             ) {

    int kk =  threadIdx.x + blockIdx.x * blockDim.x;
    if (kk >= oldKernel->nonZeros)return;

        int xOld = oldKernel->row_ptrl[kk];
        int yOld= oldKernel->col_ptrl[kk];
         int nchildren_ind = get_n_children_ind(layer - 1);
        int nChildrenX=familyX->nChildren[nchildren_ind + xOld];
        int nChildrenY=familyY->nChildren[nchildren_ind + yOld];
        int *childrenX=&familyX->children[familyX->nChildren_ind[nchildren_ind + xOld]];
        int *childrenY=&familyY->children[familyY->nChildren_ind[nchildren_ind + yOld]];
        for(int iX=0;iX<nChildrenX;iX++) {
                for(int iY=0;iY<nChildrenY;iY++) {
                    int x = childrenX[iX];
                    int y = childrenY[iY];
                    double value = getCostEff(xPos,yPos,xRadii,yRadii,layer,layer,x,layer,y);
                    value -= (alpha->data[alpha->ind[layer] + x] + beta->data[beta->ind[layer] +y]);

                        double _value= exp( - value / _eps );
                        _value *= muX[x]*muY[y];
                        int rr = atomicAdd(&newkernel->nonZeros, 1);
                        newkernel->row_ptrl[rr] = x;
                        newkernel->col_ptrl[rr] = y;
                        newkernel->val[rr] = _value;
                }
            }



}

template<int MAXELEMENT>
__global__ void transport_kernel(COO_Sparese_Matrix<MAXELEMENT> * kernel,
                                 COO_Sparese_Matrix<MAXELEMENT> *kernelT){

    int i =  threadIdx.x + blockIdx.x * blockDim.x;


    if( i < kernel->nonZeros){
        kernelT->col_ptrl[i] = kernel->row_ptrl[i];
        kernelT->row_ptrl[i] = kernel->col_ptrl[i];
        kernelT->val[i] = kernel->val[i];
    }

}





template<int MAXELEMENT>
__global__ void show_kenel_sum1(COO_Sparese_Matrix<MAXELEMENT> *kernel){
    double res = 0.0;
    for (int i = 0 ; i < kernel->nonZeros; ++i)
    {
        res += kernel->val[i];
    }
    printf("!!!!!res = %.10f\n",res);
}

template<int MAXELEMENT>
__global__ void transport_kernel_2(COO_Sparese_Matrix<MAXELEMENT> * kernel,
                                 COO_Sparese_Matrix<MAXELEMENT> *kernelT,int N){

    int i =  threadIdx.x + blockIdx.x * blockDim.x;


    if( i < N){
        kernelT->col_ptrl[i] = kernel->row_ptrl[i];
        kernelT->row_ptrl[i] = kernel->col_ptrl[i];
        kernelT->val[i] = kernel->val[i];
    }

}



template<int MAXELEMENT>
__global__ void generate_kernel_32_wrap(int layer,int layerBottom,
                                double eps,double slack,
                                double *muX,double *muY,
                                COO_Sparese_Matrix<MAXELEMENT> *kernel,
                                XData *alpha,XData *beta,
                                XData * xPos, XData *yPos,
                                XData *xRadii, XData *yRadii,
                                int xres,int yres
                             ){
    int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    int vector_id = thread_id / WARP_;
    int lane_id = thread_id % WARP_;
    int row = vector_id;
    int indx = alpha->ind[layer];
    int indy = beta->ind[layer];
    if(row < xres){
        int begin_index = kernel->crs_row_ptrl[row];
        int end_index = kernel->crs_row_ptrl[row+1];
        double value;
        for(int i = begin_index + lane_id; i < end_index; i+=WARP_){
                int y  = kernel->col_ptrl[i];
                int x = row;
                value =  (alpha->data[indx + x] + beta->data[ indy + y]) - getCostEff(xPos,yPos,xRadii,yRadii,layerBottom,layer,x,layer,y);
                kernel->val[i] = exp( value * eps )*muX[x]*muY[y];
        }
     }
}



template<int MAXELEMENT>
__global__ void kernal_mat_u_32_wrap(COO_Sparese_Matrix<MAXELEMENT> *kernel,
                             double *u,double *res,
                             int num,double *mu,double *v
                             ){
    int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    int vector_id = thread_id / WARP_;
    int lane_id = thread_id % WARP_;
    int row = vector_id;
    if(row < num){
        int begin_index = kernel->crs_row_ptrl[row];
        int end_index = kernel->crs_row_ptrl[row+1];
        double thread_sum = 0.0;
        for(int i = begin_index + lane_id; i < end_index; i+=WARP_)
            thread_sum += kernel->val[i] * u[kernel->col_ptrl[i]];
        int temp = WARP_/2;
        while(temp >= 1){
            thread_sum += __shfl_down_(thread_sum, temp);
            temp >>= 1;
        }
        if ( lane_id == 0) {
              v[row] = mu[row] / thread_sum;
    }
    }


}
template<int MAXELEMENT>
__global__ void kernal_mat_u(COO_Sparese_Matrix<MAXELEMENT> *kernel,
                             double *u,double *res,
                             int num,double *mu,double *v
                             ){

    int k = blockIdx.x + blockIdx.y * gridDim.x;
    if (k >= num)
        return;
    unsigned int tid = threadIdx.x;
    int start = kernel->crs_row_ptrl[k];
    int end  = kernel->crs_row_ptrl[k+1];
    int delta = end - start;
    double *val = &kernel->val[start];
    int *col_ptrl = &kernel->col_ptrl[start];

    __shared__ double buff[128];
    buff[tid] = 0;
    {
        for (int i = tid ; i < delta ; i += blockDim.x){
             int j = col_ptrl[i];
             buff[tid] += val[i] * u[j];
        }
    }
    __syncthreads();
    for (int stride = blockDim.x/2; stride >0; stride >>=1)
    {
        if (tid <stride)
        {
            buff[tid] += buff[tid + stride];
        }
        __syncthreads();
    }
    if(tid == 0){
         v[k] = mu[k] / buff[0];
    }

}

template<int MAXELEMENT>
__global__ void v_mat_kernal_mat_u_m_mu(COO_Sparese_Matrix<MAXELEMENT> *kernel,
                             double *u,double * res, int num,double *mu,double *v){

    unsigned int tid = threadIdx.x;
    int k = blockIdx.x * blockDim.x +  tid;
    if (k >= num)
        return;
    int start = kernel->crs_row_ptrl[k];
    int end  = kernel->crs_row_ptrl[k+1];
    if (start  == end){
        printf("dangerous!");
    }

    double temp = 0.0;

    for (int i = start ; i < end; ++i){
         int j = kernel->col_ptrl[i];
         temp += (kernel->val[i]) * u[j];

    }
//    printf("v[%d]= %.4f",k,v[k]);
    res[k] =  abs(v[k]*temp - mu[k]);

}





template<int MAXELEMENT>
double calculate_result(int N,XData *cu_xposH, XData *cu_yposH,
                        double *res,double *u,double *v,
                        COO_Sparese_Matrix<MAXELEMENT> *kernel,int layer){
    cublasHandle_t handle;
    cublasCreate(&handle);
    double result = 0.0;
    int blocks = N/128+1;
 


    scoreTransportCost<MAXELEMENT><<<blocks,128>>>(cu_xposH,cu_yposH,res,u,v,
                      N,kernel,layer);
    cublasDasum(handle,N, &kernel->val[0], 1, &result);
    return result;
}

template<int MAXELEMENT>
__global__ void update_beta(COO_Sparese_Matrix<MAXELEMENT> *kernel,
                             XData *alphaH,
                                  double eps,int num,double *mu,int layer){
    unsigned int tid = threadIdx.x;
    int k = blockIdx.x * blockDim.x +  tid;
    if (k >= num)
        return;
    int start = kernel->crs_row_ptrl[k];
    int end  = kernel->crs_row_ptrl[k+1];
    int delta = end - start;
    double temp = 0.0;
    double *alpha =  &alphaH->data[alphaH->ind[layer]];
    double mindata = -1e10;
    for (int i = start ; i < end; ++i){
        mindata = max(mindata,kernel->val[i]);
    }


    for (int i = 0 ; i < delta; ++i){
         temp += (kernel->val[start + i]/mindata) ;
    }
    alpha[k] += eps*log(mu[k]) - eps * log(temp) - eps *log(mindata) ;

}

template<int MAXELEMENT>
__global__ void update_alpha(COO_Sparese_Matrix<MAXELEMENT> *kernel,
                             XData *alphaH,
                                  double eps,int num,double *mu,int layer){
    unsigned int tid = threadIdx.x;
    int k = blockIdx.x * blockDim.x +  tid;
    if (k >= num)
        return;
    int start = kernel->crs_row_ptrl[k];
    int end  = kernel->crs_row_ptrl[k+1];
    int delta = end - start;
    double temp = 0.0;
    double *alpha =  &alphaH->data[alphaH->ind[layer]];
    double mindata = -1e10;
    for (int i = start ; i < end; ++i){
        mindata = max(mindata,kernel->val[i]);
    }


    for (int i = 0 ; i < delta; ++i){
         temp += (kernel->val[start + i]/mindata) ;
    }
    alpha[k] += eps*log(mu[k]) - eps * log(temp) - eps *log(mindata) ;

}

template<int MAXELEMENT>
void iterate_2(int xres,int yres, double *muX,double * muY,
                        COO_Sparese_Matrix<MAXELEMENT> *kernel,
                        COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                        XData *alphaH,XData *betaH,double eps,
                        int layer) {

        update_alpha<MAXELEMENT><<<xres/128+1,128>>>(kernel,alphaH,eps,xres,muX,layer);
        update_beta<MAXELEMENT><<<yres/128+1,128>>>(kernelT,betaH,eps,yres,muY,layer);


}


template<int MAXELEMENT>
void iterate(int *msg, int n,int xres, double *muX,
                       int yres,double * muY,
                        COO_Sparese_Matrix<MAXELEMENT> *kernel,
                        COO_Sparese_Matrix<MAXELEMENT> *kernelT,
                        double *u,double *v,double *res) {


    dim3 block1(xres/2+1,2);
    dim3 block2(yres/2+1,2);

    int i;
    for( i=0; i < n; ++i) {

        kernal_mat_u_32_wrap<MAXELEMENT><<<yres,WARP_>>>(kernel,v,res,xres,muX,u);
        kernal_mat_u_32_wrap<MAXELEMENT><<<xres,WARP_>>>(kernelT,u,res,yres,muY,v);

    }

#ifdef _RECORD_
    itertime.push_back(i);
#endif
}



template<int MAXELEMENT>
double getError(COO_Sparese_Matrix<MAXELEMENT>  *kernel,
                    double *muX,
                    double *u,double *v,double *res,int xres,int yres,
                cublasHandle_t handle) {

    v_mat_kernal_mat_u_m_mu<MAXELEMENT><<<xres/128+1,128>>>(kernel,v,res,yres,muX,u);
    double result = 0.0;
    cublasDasum(handle,xres, res, 1, &result);

    printf("\terror = %.8f\n",result);
    return result;
}

template<int MAXELEMENT>
__global__ void get_alpha_error(COO_Sparese_Matrix<MAXELEMENT> *kernel,
                                double *res,
                               int yres,double *muX){

    unsigned int tid = threadIdx.x;
    int k = blockIdx.x * blockDim.x +  tid;
    if (k >= yres)
        return;
    int start = kernel->crs_row_ptrl[k];
    int end  = kernel->crs_row_ptrl[k+1];
    double temp = 0.0;
    for (int i = start ; i < end; ++i){
         temp += (kernel->val[i]) ;
    }
    res[k] =  abs(temp - muX[k]);

}


template<int MAXELEMENT>
double getError_alpha_beta(COO_Sparese_Matrix<MAXELEMENT>  *kernel,
                           double *res,
                           int xres,int yres,double *muX) {
    get_alpha_error<MAXELEMENT><<<xres/128+1,128>>>(kernel,res,yres,muX);
    cublasHandle_t handle;
    cublasCreate(&handle);
    double result = 0.0;
    cublasDasum(handle,xres, res, 1, &result);
    printf("error = %.8f\n",result);
    return result;
}

 
template<int image_size>
double cuda_image_ot_solver(std::vector<double> &muXdat,
                             std::vector<double> &muYdat){


    #define MAXELEMENT  ( 1 << (ALLOCALTE_MEMORY(image_size)) )


    int dev = 0;
    cudaDeviceProp deviceProp;
    CHECK(cudaGetDeviceProperties(&deviceProp, dev));
    printf("device %d: %s \n", dev, deviceProp.name);
    CHECK(cudaSetDevice(dev));

    std::vector<int> muXdim = {image_size, image_size};
    std::vector<int> muYdim = {image_size, image_size};

    TDoubleMatrix muX,muY;

    muX.data=muXdat.data();
    muX.dimensions=muXdim.data();
    muX.depth=2;
    muY.data = muYdat.data();
    muY.dimensions = muYdim.data();
    muY.depth=2;

    int layerFinest=(GETLOG2(image_size)); // what is the finest layer we want to solve?
    int layerCoarsest = 3; // coarsest layer to solve on. sometimes skip a few layers at the top
    int msg; // store return codes from functions
 
    TMultiScaleSetupSingleGrid MultiScaleSetupX(&muX,layerFinest);
    msg=MultiScaleSetupX.Setup();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupX.SetupRadii();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupX.SetupDuals();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    

    TMultiScaleSetupSingleGrid MultiScaleSetupY(&muY,layerFinest);
    msg=MultiScaleSetupY.Setup();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupY.SetupRadii();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    msg=MultiScaleSetupY.SetupDuals();
    if(msg!=0) { printf("error: %d\n",msg); return msg; }
    int nLayers = MultiScaleSetupX.nLayers;


    Familys<OT_PROBLEM_SIZE(image_size),GETLOG2(image_size)> * familyX =
            createFamily<OT_PROBLEM_SIZE(image_size),GETLOG2(image_size)>(&MultiScaleSetupX);
    Familys<OT_PROBLEM_SIZE(image_size),GETLOG2(image_size)> * familyY =
            createFamily<OT_PROBLEM_SIZE(image_size),GETLOG2(image_size)>(&MultiScaleSetupY);
    familyX->check_for_layers(&MultiScaleSetupX);

    
//  epsScaling ///////////////////////////////////////////////////////////////////////////////////////
    double epsStart=1E8;
    double epsTarget=1E-1;
    int epsSteps=20;
    TEpsScalingHandler epsScalingHandler(epsStart,epsTarget,epsSteps); // basic eps scaling
    epsScalingHandler.getEpsScalesFromBox(muXdim[0],2,MultiScaleSetupX.nLayers); // eps scales for each layer
    epsScalingHandler.getEpsScalingSplit(layerCoarsest,1); // create sub eps lists for each layer
    printf("total layers = %d\n",MultiScaleSetupX.nLayers);
    CuEspsScalingHandler_ h_eps;
    h_eps.init_eps(epsScalingHandler);


    CuSinkhornSolverParameters cfg_={
        1E-4, // maxError
        100000, // maxIterations
        1000, // innerIterations
        20, // maxAbsorptionLoops
        1E3, // absorption_scalingBound
        1E3, // absorption_scalingLowerBound
        1E-1, // truncation_thresh
        1 // refineKernel when refining layer (as opposed to attempting to estimate directly)
        };



     

     int *numsX = (int *)malloc(sizeof(int)*nLayers);
     int *numsY = (int *)malloc(sizeof(int)*nLayers);
     for (int i = 0 ; i < nLayers ; ++i){
         numsX[i] = MultiScaleSetupX.HP->layers[i]->nCells;
         numsY[i] = MultiScaleSetupY.HP->layers[i]->nCells;
     }

     int *numsradiiX = (int *)malloc(sizeof(int)*nLayers);
     int *numsradiiY = (int *)malloc(sizeof(int)*nLayers);
     for (int i = 0 ; i < nLayers -2 ; ++i){
         numsradiiX[i] = MultiScaleSetupX.HB->layers[i].nodes.size();
         numsradiiY[i] = MultiScaleSetupY.HB->layers[i].nodes.size();
     }

     int *numposX = (int *)malloc(sizeof(int)*nLayers);
     int *numposY = (int *)malloc(sizeof(int)*nLayers);
     for (int i = 0 ; i < nLayers ; ++i){
         numposX[i] = MultiScaleSetupX.HB->layers[i].nodes.size()*2;
         numposY[i] = MultiScaleSetupY.HB->layers[i].nodes.size()*2;
     }

     XData *xposH = create_dev_struct<XData>(numposX,nLayers,MultiScaleSetupX.posH);
     XData *cu_xposH =  malloc_struct<XData>(xposH);
     XData *yposH = create_dev_struct<XData>(numposY,nLayers,MultiScaleSetupY.posH);
     XData *cu_yposH = malloc_struct<XData>(yposH);

     XData *xradii = create_dev_struct<XData>(numsradiiX,nLayers-2,MultiScaleSetupX.radii);
     XData *cu_xradii =  malloc_struct<XData>(xradii);

     XData *yradii = create_dev_struct<XData>(numsradiiY,nLayers-2,MultiScaleSetupY.radii);
     XData *cu_yradii = malloc_struct<XData>(yradii);


     XData *muXH = create_dev_struct<XData>(numsX,nLayers,MultiScaleSetupX.muH);
     XData *cu_muXH =  malloc_struct<XData>(muXH);
     XData *muYH = create_dev_struct<XData>(numsY,nLayers,MultiScaleSetupY.muH);
     XData *cu_muYH = malloc_struct<XData>(muYH);

     XData *alphaH = create_dev_struct<XData>(numsX,nLayers,MultiScaleSetupX.alphaH);
     XData *cu_alphaH =  malloc_struct<XData>(alphaH);
     XData *betaH = create_dev_struct<XData>(numsY,nLayers,MultiScaleSetupY.alphaH);
     XData *cu_betaH =  malloc_struct<XData>(betaH);




     Ind_Ptr_ *indptrl = (Ind_Ptr_ *)malloc(sizeof(Ind_Ptr_) * MAXELEMENT);
     for (int i = 0 ;  i < MAXELEMENT; ++i){
         indptrl[i] = i;
     }
     Ind_Ptr_ *indptrl_cu_mid = malloc_struct(indptrl,MAXELEMENT);
     Ind_Ptr_ *indptrl_cu;
     CHECK(cudaMalloc((void **)&indptrl_cu, sizeof(Ind_Ptr_) * MAXELEMENT ));

    Familys<OT_PROBLEM_SIZE(image_size),GETLOG2(image_size)> *cu_familyX,*cu_familyY;
    size_t cufamily_size = sizeof(Familys<OT_PROBLEM_SIZE(image_size),GETLOG2(image_size)>);
    CHECK(cudaMalloc((void **)&cu_familyX, cufamily_size));
    CHECK(cudaMalloc((void **)&cu_familyY, cufamily_size));
    CHECK(cudaMemcpy(cu_familyX, familyX, cufamily_size, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(cu_familyY, familyY, cufamily_size, cudaMemcpyHostToDevice));

    COO_Buff * buff;
    //CHECK(cudaMalloc((void **)&buff, sizeof( COO_Buff) * 100));
    double *muX_cuda ;
    double *muY_cuda ;

    double *u;
    double *v;
    double *res;
    double * origonals_one_cuda;
    CHECK(cudaMalloc((void **)&origonals_one_cuda,sizeof(double)*OT_PROBLEM_SIZE(image_size)));

    CHECK(cudaMalloc((void **)&muX_cuda,sizeof(double)*OT_PROBLEM_SIZE(image_size)));
    CHECK(cudaMalloc((void **)&muY_cuda,sizeof(double)*OT_PROBLEM_SIZE(image_size)));
    CHECK(cudaMalloc((void **)&u,sizeof(double)*OT_PROBLEM_SIZE(image_size)));
    CHECK(cudaMalloc((void **)&v,sizeof(double)*OT_PROBLEM_SIZE(image_size)));
    #ifdef _RECORD_
    
    CHECK(cudaMalloc((void **)&res,sizeof(double)* MAXELEMENT));
    #else
    CHECK(cudaMalloc((void **)&res,sizeof(double)*OT_PROBLEM_SIZE(image_size)));
    #endif
    

    int layer = layerCoarsest;
    int kernelValid = 0;
    int xres = 0;
    int yres = 0;
    int layerBottom = layerFinest;
    double * origionals = (double *)malloc(sizeof(double) * OT_PROBLEM_SIZE(image_size));
    for (int i = 0 ; i < OT_PROBLEM_SIZE(image_size); ++i)
        origionals[i] = 1.0;
    CHECK(cudaMemcpy(origonals_one_cuda, origionals, sizeof(double)*OT_PROBLEM_SIZE(image_size),
                     cudaMemcpyHostToDevice));
    int *msg_cu;
    cudaMalloc((void **)&msg_cu,sizeof(int));
    msg = changeLayer(layerCoarsest,&layer,&kernelValid,
                MultiScaleSetupX.HP,MultiScaleSetupY.HP,
                &xres,&yres,u,v,origonals_one_cuda,
                           muX_cuda,muY_cuda,
                           cu_muXH,cu_muYH,muXH,muYH);
    COO_Sparese_Matrix<MAXELEMENT> *kernel,*kernelT;
    CHECK(cudaMalloc((void **)&kernel,sizeof(COO_Sparese_Matrix<MAXELEMENT>)));
    CHECK(cudaMalloc((void **)&kernelT,sizeof(COO_Sparese_Matrix<MAXELEMENT>)));
    double eps = 1e10;
    if(msg!=0) return msg;

    cublasHandle_t handle;
    cublasCreate(&handle);


    #ifdef _RECORD_
    clock_t begin = clock();
    #endif
    //generate the whole region
    msg=generateKernel_<MAXELEMENT>(kernel,kernelT,&cfg_,eps,
                    muX_cuda,muY_cuda,layer,
                    cu_alphaH,cu_betaH,cu_xposH,cu_yposH,cu_xradii,cu_yradii,xres,yres,
                    indptrl_cu_mid,indptrl_cu,buff);
    kernelValid = 1;




    while(true) {
     
        if( layer >= 8){
            // for large problems
            cfg_.maxError = 1E-4;
        }
        for(int nEps=0; nEps<h_eps.nEpsLists[layer]; nEps++) {
            eps = h_eps.epsLists[layer][nEps];
            printf("\teps=%e\n",eps);
            kernelValid = 0;
            if(cfg_.refineKernel && (nEps==0) && (layer>layerCoarsest)) {
                refinekernel_<image_size,image_size,GETLOG2(image_size),MAXELEMENT>(kernel,kernelT,muX_cuda,
                                  muY_cuda,
                                  layerBottom,
                                  eps,xres,yres,
                                  layer,
                                  cu_familyX,cu_familyY,
                                  cu_alphaH,cu_betaH,
                                  cu_xposH, cu_yposH,
                                  cu_xradii, cu_yradii,
                                 indptrl_cu_mid,indptrl_cu,
                                 &kernelValid,buff);


            }
            int msg=solveSingle_<image_size,image_size,GETLOG2(image_size),MAXELEMENT>(msg_cu,&kernelValid,
                                kernel,kernelT,
                                &cfg_,eps,
                                muX_cuda,muY_cuda,
                                cu_familyX,cu_familyY,
                                layerBottom,
                                cu_alphaH,cu_betaH,
                                cu_xposH,cu_yposH,
                                cu_xradii,cu_yradii,
                                xres,yres,u,v,res,layer,
                                indptrl_cu_mid,
                                indptrl_cu,buff,
                                 handle
                                );

            if(msg!=0) return msg;
        }



        if(msg!=0) return msg;
        if(layer < layerFinest) {

 
            int newLayer = layer + 1;
            if(newLayer>0) {
                signal_refine_all<OT_PROBLEM_SIZE(image_size),GETLOG2(image_size)>(cu_familyX,cu_familyY,cu_alphaH,cu_betaH,layer);
            }
            kernelValid=0;
            if(msg!=0) return msg;

            msg = changeLayer(newLayer,&layer,&kernelValid,
                        MultiScaleSetupX.HP,MultiScaleSetupY.HP,
                        &xres,&yres,u,v,origonals_one_cuda,
                                   muX_cuda,muY_cuda,
                                   cu_muXH,cu_muYH,muXH,muYH
                             );

            if(msg!=0) return msg;

        } else {
            break;
        }

    }




    int N = 0;
    msg = refinekernel_2<MAXELEMENT>(kernel,kernelT,
                      muX_cuda,
                      muY_cuda,
                      layerBottom,
                      eps,xres,yres,
                      layer,
                      cu_alphaH,cu_betaH,cu_xposH,cu_yposH,cu_xradii,cu_yradii,
                     indptrl_cu_mid,indptrl_cu,
                      &kernelValid,buff,&cfg_);

    cudaMemcpy(&N,&kernel->nonZeros,sizeof(int),cudaMemcpyDeviceToHost);

    eprintf("\t\tkernel entries: %ld\n",N);
    printf("return code: %d\n",msg);
        // save kernel to the file

    COO_Sparese_Matrix<MAXELEMENT> *kernel_host = (COO_Sparese_Matrix<MAXELEMENT> *)  malloc(  sizeof(COO_Sparese_Matrix<MAXELEMENT>) );
    cudaMemcpy(kernel_host, kernel,sizeof(COO_Sparese_Matrix<MAXELEMENT>),cudaMemcpyDeviceToHost);
    std::ofstream kernel_outFile;
 
    std::string filename_kernel =   "output.csv";
    kernel_outFile.open(filename_kernel, std::ios::out);
    for (int i = 0 ; i < kernel_host->nonZeros; ++i){
        kernel_outFile << kernel_host->row_ptrl[i] << "," << kernel_host->col_ptrl[i] << "," << kernel_host->val[i] <<   std::endl;
    } 
    kernel_outFile.close();      
    free(kernel_host);
    double result = calculate_result(N,cu_xposH,cu_yposH,res,u,v,kernel,layer);
    #ifdef  _RECORD_ 
            cudaMemcpy(&N,&kernel->nonZeros,sizeof(int),cudaMemcpyDeviceToHost);
            kernel_entries.push_back(N);
            auto maxPosition = std::max_element(kernel_entries.begin(), kernel_entries.end());
            optimal_values.push_back(result);     
             std::ofstream outFile;
 
            std::string s =  "./output/" + std::to_string(image_size) +  "_" + std::to_string(int(p_index*10)) + ".csv";
            outFile.open(s, std::ios::out);
            for (int i = 0 ; i < optimal_values.size();++i){
                outFile << optimal_values[i] <<  std::endl;
            } 
            outFile.close();   
    #endif

    #ifdef _RECORD_
    // record the number of
            std::ofstream outFile_iter;
 
            std::string s_iter =  "./output/" + std::to_string(image_size) +  "_iteration_size_" + std::to_string(int(p_index*10)) + ".csv";
            outFile_iter.open(s_iter, std::ios::out);
            for (int i = 0; i < itertime.size(); ++i){
                outFile_iter << itertime[i] <<  std::endl;
            }
            outFile_iter.close(); 
    #endif


    #ifdef _RECORD_
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        double t =  elapsed_secs;

        std::ofstream outFile_time;
 
        std::string s_time =  "./output/" + std::to_string(image_size) +  "_iteration_time_" + std::to_string(int(p_index*10)) + ".csv";
        outFile_time.open(s_time, std::ios::app);

        outFile_time << t <<  std::endl;

        outFile_time.close(); 

    #endif


     
 

    cublasDestroy(handle);
    cudaFree(msg_cu);
    cudaFree(muX_cuda);
    cudaFree(muY_cuda);
    cudaFree(u);
    cudaFree(v);
    cudaFree(kernel);
    cudaFree(kernelT);
    //cudaFree(buff);
    cudaFree(indptrl_cu);
    cudaFree(indptrl_cu_mid);
    cudaFree(cu_xposH);
    cudaFree(cu_yposH);
    cudaFree(cu_betaH);
    cudaFree(cu_alphaH);
    cudaFree(cu_xradii);
    cudaFree(cu_yradii);
    cudaFree(res);
    cudaFree(cu_muXH);
    cudaFree(cu_muYH);
    cudaFree(origonals_one_cuda);



    cudaFree(cu_familyX);
    cudaFree(cu_familyY);



    return result;


}


#endif // MY_SK_H

