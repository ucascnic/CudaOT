# -*- coding: utf-8 -*-
"""

@author: Administrator
"""

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
from pykeops.torch.cluster import cluster_ranges_centroids
from pykeops.torch import generic_argmin
from pykeops.torch import generic_sum,generic_logsumexp,generic_argmin
use_cuda = torch.cuda.is_available()

from collections import Counter
dtype = torch.cuda.DoubleTensor if use_cuda else torch.FloatTensor
#use_cuda = 1
def load_image(fname):
    img = imageio.imread(fname)  # RGB
    return img / 255.0  # Normalized to [0,1]


def RGB_cloud(fname, sampling, dtype=torch.FloatTensor):
    A = load_image(fname)
    #A = A[::sampling, ::sampling, :]
    
    data = torch.from_numpy(A).type(dtype).view(-1, 3)
    
    
    X_i = (
    torch.stack((data[:,0], data[:,1], data[:,2]))
    .t()
    .contiguous()
    )
    return X_i

 


from collections import Counter
def get_measure(lab,N,allN):
    lab_i_cpu = list(np.array(lab.cpu())) 
    xdict = Counter(lab_i_cpu)
    a = np.ones(shape=(N,),dtype=np.float64)/N
    alla = 0
    for key in xdict.keys():
        alla += xdict[key]

    for key in xdict.keys():
        a[key] = xdict[key]

    a = a/sum(a)
    return a

def squared_distances(x, y):
    if x.dim() == 2:
        D_xx = (x * x).sum(-1).unsqueeze(1)  # (N,1)
        D_xy = torch.matmul(x, y.permute(1, 0))  # (N,D) @ (D,M) = (N,M)
        D_yy = (y * y).sum(-1).unsqueeze(0)  # (1,M)
    elif x.dim() == 3:  # Batch computation
        D_xx = (x * x).sum(-1).unsqueeze(2)  # (B,N,1)
        D_xy = torch.matmul(x, y.permute(0, 2, 1))  # (B,N,D) @ (B,D,M) = (B,N,M)
        D_yy = (y * y).sum(-1).unsqueeze(1)  # (B,1,M)
    else:
        print("x.shape : ", x.shape)
        raise ValueError("Incorrect number of dimensions")

    return D_xx - 2 * D_xy + D_yy
def sort_clusters(x, lab):
    lab, perm = torch.sort(lab.view(-1))
    if type(x) is tuple:
        x_sorted = tuple(a[perm] for a in x)
    elif type(x) is list:
        x_sorted = list(a[perm] for a in x)
    else:
        x_sorted = x[perm]

    return x_sorted, lab, perm


def combine_measure(x_labels,measure):
    lab_i_cpu = list(np.array(x_labels.cpu())) 
    data = pd.DataFrame()
    data['index'] = lab_i_cpu
    data['value'] = measure
    xx = data.groupby('index').sum()
    xx = torch.as_tensor(xx.value).type(dtype).unsqueeze(1)
    return xx


# broadcast the solution
def progate(f,x_label):
    nn = x_label.shape[0]
    x_labelscpu = x_label.cpu()
    temp = pd.DataFrame(np.array(f.cpu())) 
    fnew = torch.as_tensor(np.array(temp.iloc[x_labelscpu])).type(dtype)
    return fnew

def generate_the_optimal_transport_map(mu, nv, x,y,finit,ginit,   epsilon0 = 1e6,ranges = None):
    
    x_i, y_j = LazyTensor(x[:, None, :]), LazyTensor(y[:, None, :])
    
    logss = generic_argmin(
       '(  (SqNorm2(x-y) - f - g))  ', # Formula
        
        'a = Vi(1)',           
        'x = Vi(3)',           

        'y = Vj(3)',  'f = Vi(1)', 
        
        'g = Vj(1)',
        
        dtype  = 'float64'  )          

    K =   logss(x,y,finit,ginit,ranges = ranges)

 
    return  K

def generate_the_optimal_transport_map_smooth(mu, nv, x,y,finit,ginit,   epsilon0 = 1e6,ranges = None):
    
    x_i, y_j = LazyTensor(x[:, None, :]), LazyTensor(y[:, None, :])
    
    logss2 = generic_logsumexp(
       '(   ( (  f + g - SqNorm2(x-y)) * eps) + z - mu ) ', # Formula
        
        'a = Vi(1)',           
        'x = Vi(3)',           

        'y = Vj(3)',  'f = Vi(1)', 
        
        'g = Vj(1)','eps = Pm(1)','z = Vj(1)', 'mu = Vi(1)', 
        
        dtype  = 'float64'  )          

    K1 =   logss2(x,y,finit,ginit,epsilon0,torch.log(y[:,0]).unsqueeze(1),torch.log(mu),ranges = ranges)
    K2 =   logss2(x,y,finit,ginit,epsilon0,torch.log(y[:,1]).unsqueeze(1),torch.log(mu),ranges = ranges)
    K3 =   logss2(x,y,finit,ginit,epsilon0,torch.log(y[:,2]).unsqueeze(1),torch.log(mu),ranges = ranges)

 
    return  torch.exp(K1),torch.exp(K2),torch.exp(K3)


def show_sorted_image(orginal,xdata,sizex,sizey):
    data = pd.DataFrame()
    data['index'] = np.array(orginal.cpu())
    temp =  np.array(xdata.cpu()) 
    data['r1'] = temp[:,0]
    data['r2'] = temp[:,1]
    data['r3'] = temp[:,2]
    d = torch.tensor(np.array(data.sort_values('index')))
    need = (
        torch.stack((d[:,1], d[:,2], d[:,3]))
    .t()
    .contiguous()
    )
    return np.array(need.reshape((sizex,sizey,3)).cpu())
 