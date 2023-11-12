

__global__ void computemin_kernel(int *ind,int *col,
						int m,int n,double* out,double* data) {
						
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; \
            i < m; i += gridDim.x * blockDim.x) {

		int local_row = ind[i];
		int local_row_end = ind[i+1];
		double out_temp = data[local_row];
		//double * temp = data + local_row;
		for (int j = local_row + 1; j < local_row_end; j++){
			
			out_temp =  min(data[j],out_temp);				
		}
		out[i] = out_temp;
    }
	//out[0] = 13;
}


__global__ void computesum_kernel(int *ind,int *col,
						int m,int n,double* out,double* data) {
						
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; \
            i < m; i += gridDim.x * blockDim.x) {

		int local_row = ind[i];
		int local_row_end = ind[i+1];
		double out_temp = 0.0;
		//double * temp = data + local_row;
		for (int j = local_row; j < local_row_end; j++){
			
			//if (out_temp < data[j])
				{
					out_temp += data[j];
				}
		}
		out[i] = out_temp;
    }
}

#include<stdio.h>
void launch_computemin(int *ind,int *col,int m,int n,double* out,double* data){


	
    dim3 grid((m + 1023) / 1024);
    dim3 block(1024);
	//printf("hello world");
    computemin_kernel<<<grid, block>>>(ind,col,m,n,out,data);
	
	
	return ;
  
  
}

void launch_computesum(int *ind,int *col,int m,int n,double* out,double* data){


	
    dim3 grid((m + 1023) / 1024);
    dim3 block(1024);
	//printf("hello world");
    computesum_kernel<<<grid, block>>>(ind,col,m,n,out,data);
	
	
	return ;
  
}