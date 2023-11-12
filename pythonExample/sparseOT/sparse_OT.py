import time
import torch
import torch.nn.functional as F

import numpy as np
import os

from log_sinkhorn_sparse import *

os.environ['CUDA_LAUNCH_BLOCKING'] = '3'

torch.set_printoptions(precision=8, threshold=np.inf, sci_mode=False)

torch.cuda.set_device(3)


dist_matrix = torch.load('dist_matrix_15000',map_location=torch.device('cuda:3')) #dense matrix
indices = torch.load('indices_15000',map_location=torch.device('cuda:3'))

print('indices.device:', indices.device)

[N,m] = torch.load('mn_15000')

print('dist_matrix.shape:', dist_matrix.shape)
# print('dist_matrix[:3,:]: ', dist_matrix[:3,:])  


print('N: ', N)  
print('m: ', m)  

sp_matrix  = torch.sparse_coo_tensor(indices, dist_matrix.contiguous().view(-1), size=[N, m])
abs_matrix = sp_matrix.to_dense().contiguous() 

# row_1 = torch.nonzero(abs_matrix[0,:],as_tuple=False) 
# print('row_1: ', row_1)

# value_row_1 = torch.index_select(abs_matrix[0,:], dim=0, index=row_1.squeeze())
# print('value_row_1: ', value_row_1)


sp_matrix.shape

# get_ipython().run_line_magic('time', '')
reg=1e-6
# a = torch.ones((N,1),dtype=torch.float64)/N
# b = torch.ones((m,1),dtype=torch.float64)/m
a = torch.ones((N,1),dtype=torch.float64)
b = torch.ones((m,1),dtype=torch.float64)*N/float(m)

print(torch.cuda.memory_allocated() / (1024*1024))

start = time.time()
sparse_coo_transport = log_sinkhorn_sparse(sp_matrix,indices,a, b)
end = time.time()

print('whole time: ', end-start)

print(torch.cuda.memory_allocated() / (1024*1024))
# print(sparse_coo_transport)
dense_transport = sparse_coo_transport.to_dense()
# print(dense_transport[:2,:])

print('dist_matrix[0,:6]:', dist_matrix[0,:6])


row_1 = torch.nonzero(dense_transport[0,:],as_tuple=False) 
print('row_1: ', row_1)
value_row_1 = torch.index_select(dense_transport[0,:], dim=0, index=row_1.squeeze())
print('value_row_1: ', value_row_1)

print('indices.shape:', indices.shape)
print('indices[1,:6]: ', indices[1,:6]) 
print('-'*20)
 
# dense_transport = F.normalize(dense_transport, p=2, dim=-1)
# row_1 = torch.nonzero(dense_transport[0,:],as_tuple=False) 
# print('row_1: ', row_1)
# value_row_1 = torch.index_select(dense_transport[0,:], dim=0, index=row_1.squeeze())
# print('value_row_1: ', value_row_1)

# indices = indices.T
# print('indices.shape:', indices.shape)
# print('indices[:6,:]: ', indices[:6,:]) 




    
dense_N = torch.sum(dense_transport,dim=0) 
print('dense_N.shape: ', dense_N.shape)
print(dense_N[:10])


dense_m = torch.sum(dense_transport,dim=1) 
print('dense_m.shape: ', dense_m.shape)
print(dense_m[:10])



torch.sum( (sparse_coo_transport * sp_matrix).to_dense() )