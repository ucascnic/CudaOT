#ifndef __MEMEMORY_M__H__
#define __MEMEMORY_M__H__
#include"cuda_runtime_api.h"
#include"Common/handler_cuda_error.h"



template <typename T>
inline T* malloc_struct(T *data){
    T * cudata;
    size_t nBytesmuXH_cuda = sizeof(T);
    CHECK(cudaMalloc((void **)&cudata, nBytesmuXH_cuda));
    CHECK(cudaMemcpy(cudata, data, nBytesmuXH_cuda, cudaMemcpyHostToDevice));

    return cudata;
}

template <typename T>
inline T* malloc_struct(T *data,size_t size){
    T * cudata;
    size_t nBytesmuXH_cuda = sizeof(T)*size;
    CHECK(cudaMalloc((void **)&cudata, nBytesmuXH_cuda));
    CHECK(cudaMemcpy(cudata, data, nBytesmuXH_cuda, cudaMemcpyHostToDevice));

    return cudata;
}

template <typename T>
inline T * create_dev_struct(int *nums,int total,double **data){
    T *mu = (T*) malloc(sizeof(T));
    mu->ind[0] = 0;
    int mu_ind = 0;
    for (int i = 0 ; i < total ; ++i){
        int num = nums[i];
        memcpy(&mu->data[mu->ind[i]],&data[i][0],num*sizeof(double));
        mu_ind += num;
        mu->ind[i+1] = mu_ind;
    }
    return mu;
}

#endif
