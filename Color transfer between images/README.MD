We recommend running our code on a GPU with more than 8GB of memory for satisfactory acceleration
To run the jupyter notebook, the following libraries are required
CMake 3.18 or higher
pandas 
torch
opencv-python  
pykeops   please reffer to http://www.kernel-operations.io/keops 



```python
import os.path
import sys
sys.path.append('/home/chenyidong/keops')
import numpy as np
import torch
from random import choices
import imageio
from pykeops.torch import LazyTensor
from pykeops.torch.cluster import grid_cluster
import pandas as pd
import numpy as np
from pykeops.torch.cluster import cluster_ranges_centroids
from pykeops.torch import generic_argmin
from pykeops.torch import generic_argmin
from pykeops.torch import generic_sum,generic_logsumexp,generic_argmin
from pykeops.torch import generic_argmin
from pykeops.torch import generic_sum,generic_logsumexp,generic_argmin
use_cuda = torch.cuda.is_available()
dtype = torch.cuda.DoubleTensor if use_cuda else torch.FloatTensor
from mytools import load_image,RGB_cloud,get_measure,combine_measure,squared_distances,sort_clusters,progate
from SparseSinkhorn import  sinkhorn_on_log_domain_GPU_keops,sinkhorn_on_log_domain_GPU
```

 



```python
d = 3
n = 500
X_i_cpu = np.random.randn(n,d) 
Y_j_cpu = np.random.randn(n,d)
X_i = torch.from_numpy(X_i_cpu).type(dtype).view(-1, 3)
Y_j = torch.from_numpy(Y_j_cpu).type(dtype).view(-1, 3)
```


```python
# solve by the pot !pip install POTa
import ot
a = np.ones((len(X_i_cpu),1))/X_i_cpu.shape[0]
b = np.ones((len(Y_j_cpu),1))/Y_j_cpu.shape[0]
M = np.sum((X_i_cpu[:, None, :] - Y_j_cpu[None, :, :]) ** 2, 2)
T = ot.emd(np.squeeze(a), np.squeeze(b), M)  # exact linear program
distance = np.sum(np.multiply(T,M))
```


```python
# solve by the sinkhorn on log domain
reg=1e-3
f,g,regk = sinkhorn_on_log_domain_GPU(a, b, M, reg, epsilon0 = 1e6, numItermax=30,stopThr=1e-9)
def get_K(alpha, beta,reg,M):
    """log space computation"""
    device = torch.device('cuda')
    M_ = torch.as_tensor(M, dtype=torch.float64,device=device)
    return torch.exp(-(M_ - alpha - beta) / reg),M_
T_log_domain,M_ = get_K(f, g,regk,M)
distance_log = torch.sum(torch.mul(T_log_domain,M_))
```

    tensor(1.2109e-30, device='cuda:0', dtype=torch.float64)
    tensor(1., device='cuda:0', dtype=torch.float64)
    tensor(6.3862e-31, device='cuda:0', dtype=torch.float64)
    tensor(2., device='cuda:0', dtype=torch.float64)
    tensor(1.1816e-30, device='cuda:0', dtype=torch.float64)
    tensor(3., device='cuda:0', dtype=torch.float64)
    tensor(5.6668e-31, device='cuda:0', dtype=torch.float64)
    tensor(4., device='cuda:0', dtype=torch.float64)
    tensor(1.4085e-30, device='cuda:0', dtype=torch.float64)
    tensor(5., device='cuda:0', dtype=torch.float64)
    tensor(5.6809e-31, device='cuda:0', dtype=torch.float64)
    tensor(6., device='cuda:0', dtype=torch.float64)
    tensor(6.0627e-31, device='cuda:0', dtype=torch.float64)
    tensor(7., device='cuda:0', dtype=torch.float64)
    tensor(1.1605e-30, device='cuda:0', dtype=torch.float64)
    tensor(8., device='cuda:0', dtype=torch.float64)
    tensor(5.5991e-31, device='cuda:0', dtype=torch.float64)
    tensor(9., device='cuda:0', dtype=torch.float64)
    tensor(1.2439e-30, device='cuda:0', dtype=torch.float64)
    tensor(10., device='cuda:0', dtype=torch.float64)
    tensor(6.2141e-31, device='cuda:0', dtype=torch.float64)
    tensor(11., device='cuda:0', dtype=torch.float64)
    tensor(7.5824e-31, device='cuda:0', dtype=torch.float64)
    tensor(12., device='cuda:0', dtype=torch.float64)
    tensor(4.2393e-31, device='cuda:0', dtype=torch.float64)
    tensor(13., device='cuda:0', dtype=torch.float64)
    tensor(4.2515e-31, device='cuda:0', dtype=torch.float64)
    tensor(14., device='cuda:0', dtype=torch.float64)
    tensor(3.1221e-31, device='cuda:0', dtype=torch.float64)
    tensor(15., device='cuda:0', dtype=torch.float64)
    tensor(2.3764e-31, device='cuda:0', dtype=torch.float64)
    tensor(16., device='cuda:0', dtype=torch.float64)
    tensor(8.5670e-31, device='cuda:0', dtype=torch.float64)
    tensor(17., device='cuda:0', dtype=torch.float64)
    tensor(2.4636e-30, device='cuda:0', dtype=torch.float64)
    tensor(18., device='cuda:0', dtype=torch.float64)
    tensor(1.4544e-29, device='cuda:0', dtype=torch.float64)
    tensor(19., device='cuda:0', dtype=torch.float64)
    tensor(5.3121e-29, device='cuda:0', dtype=torch.float64)
    tensor(20., device='cuda:0', dtype=torch.float64)
    tensor(1.7076e-28, device='cuda:0', dtype=torch.float64)
    tensor(21., device='cuda:0', dtype=torch.float64)
    tensor(2.7141e-28, device='cuda:0', dtype=torch.float64)
    tensor(22., device='cuda:0', dtype=torch.float64)
    tensor(4.3559e-28, device='cuda:0', dtype=torch.float64)
    tensor(23., device='cuda:0', dtype=torch.float64)
    tensor(4.8757e-28, device='cuda:0', dtype=torch.float64)
    tensor(24., device='cuda:0', dtype=torch.float64)
    tensor(3.5295e-28, device='cuda:0', dtype=torch.float64)
    tensor(25., device='cuda:0', dtype=torch.float64)
    tensor(4.5008e-28, device='cuda:0', dtype=torch.float64)
    tensor(26., device='cuda:0', dtype=torch.float64)
    tensor(5.0937e-28, device='cuda:0', dtype=torch.float64)
    tensor(27., device='cuda:0', dtype=torch.float64)
    tensor(5.1905e-28, device='cuda:0', dtype=torch.float64)
    tensor(28., device='cuda:0', dtype=torch.float64)
    tensor(4.0012e-28, device='cuda:0', dtype=torch.float64)
    tensor(29., device='cuda:0', dtype=torch.float64)
    tensor(3.6679e-28, device='cuda:0', dtype=torch.float64)
    tensor(30., device='cuda:0', dtype=torch.float64)



```python
# solve by the sinkhorn on log domain with keops
reg=1e-3
measurex = torch.as_tensor(a).type(dtype)
measurey = torch.as_tensor(b).type(dtype) 
f_keops,g_keops,regk = sinkhorn_on_log_domain_GPU_keops(measurex, measurey, X_i, Y_j, reg, epsilon0 = 1e6, numItermax=30,stopThr=1e-9)
def get_K(alpha, beta,reg,M):
    """log space computation"""
    device = torch.device('cuda')
    M_ = torch.as_tensor(M, dtype=torch.float64,device=device)
    return torch.exp(-(M_ - alpha - beta) / reg),M_
T_keops,M_ = get_K(f_keops, g_keops.T, regk,M)
distance_keops = torch.sum(torch.mul(T_keops,M_))
```

    tensor(1.5119e-16, device='cuda:0', dtype=torch.float64)
    tensor(0.0002, device='cuda:0', dtype=torch.float64)
    tensor(7.5715e-05, device='cuda:0', dtype=torch.float64)



```python
print(distance_log)
print(distance)
print(distance_keops)
```

    tensor(0.2506, device='cuda:0', dtype=torch.float64)
    0.2504988107340873
    tensor(0.2507, device='cuda:0', dtype=torch.float64)