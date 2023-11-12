
#include <torch/extension.h>
void launch_computesum(int *ind,int *col,int m,int n,double* out,double* data);

void launch_computemin(int *,int *,int,int,double*,double*);
