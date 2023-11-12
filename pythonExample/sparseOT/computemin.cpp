


#include <torch/extension.h>
#include "computemin.h"
#include<cuda_runtime_api.h>
#include<cublas_v2.h>
#include<stdio.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <cmath>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/extrema.h>

 
void launch_computemin(int *,int *,int,int,double*,double*);
void computemin(torch::Tensor &ind, 
			torch::Tensor &col,int m,int n,
			torch::Tensor &out,torch::Tensor &data) {


	
	int * ind_ = (int *)ind.data_ptr();
	int * col_ = (int *)col.data_ptr();
	double * out_ = (double *)out.data_ptr();
	double * data_ = (double *)data.data_ptr();
	
	launch_computemin(ind_,col_,m,n,out_,data_);
	return ;
}

void computesum(torch::Tensor &ind, 
			torch::Tensor &col,int m,int n,
			torch::Tensor &out,torch::Tensor &data) {


	
	int * ind_ = (int *)ind.data_ptr();
	int * col_ = (int *)col.data_ptr();
	double * out_ = (double *)out.data_ptr();
	double * data_ = (double *)data.data_ptr();
	
	launch_computesum(ind_,col_,m,n,out_,data_);
	return ;
}


PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
   
    m.def("computemin",
          &computemin,
          "computemin");
		      m.def("computesum",
          &computesum,
          "computesum");
		
}