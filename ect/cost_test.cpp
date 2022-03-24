#include"Common/handler_cuda_error.h"
#include<cuda_runtime_api.h>
#include <stdio.h>
#include"vshelp.h"
#include<cublas_v2.h>
#include<Common.h>

#include <iostream>
#include"my_sk.h"
#include<string>
#include"memory_manage.h"
#include"family.h"



int main(int argc, char **argv){

    int dev = 0;
    cudaDeviceProp deviceProp;
    CHECK(cudaGetDeviceProperties(&deviceProp, dev));
     printf("device %d: %s \n", dev, deviceProp.name);
    CHECK(cudaSetDevice(dev));
    std::string filex,filey;
    filex = "../data/ot_benchmark/CauchyDensity/data"+ std::to_string(ROW)+"_1001.dat";
    filey = "../data/ot_benchmark/CauchyDensity/data"+ std::to_string(COL)+"_1002.dat";
    std::vector<double> muXdat=readFile<double>(filex.c_str());
    std::vector<double> muYdat=readFile<double>(filey.c_str());
    std::vector<int> muXdim = {ROW, COL};
    std::vector<int> muYdim = {ROW, COL};
    TDoubleMatrix muX,muY;
    muX.data=muXdat.data();
    muX.dimensions=muXdim.data();
    muX.depth=2;
    muY.data = muYdat.data();
    muY.dimensions = muYdim.data();
    muY.depth=2;

    int msg; // store return codes from functions
    int layerFinest=MAXLAYER; // what is the finest layer we want to solve?

    Family    familyx(muXdim,muXdat,layerFinest);
    familyx.init();

    for (int i = 0 ; i < 10; ++i){
    int n = get_children_ind(i);

        printf("n = %d\n",n);

    }
    for (int i = 0 ; i < 10; ++i){
        int k = get_n_children_ind(i);
         printf("k = %d\n",k);

    }


    int n = TOTALNUM + ((1 << (DIM * ( MAXLAYER  ))) - (1 << DIM ))/((1 << DIM) - 1);

    printf("n = %d\n",n);
    n = 1 + ((1 << (DIM * (MAXLAYER ))) - (1 << DIM ))/((1 << DIM) - 1) + TOTALNUM ;
    printf("n = %d\n",n);
    return 0;
}
