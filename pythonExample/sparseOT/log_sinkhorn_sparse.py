import os.path
import sys
import time

import numpy as np
import torch
from random import choices
import pandas as pd
use_cuda = torch.cuda.is_available()
dtype = torch.cuda.DoubleTensor if use_cuda else torch.FloatTensor
from torch.utils.cpp_extension import load
cuda_module = load(name="computemin",
                   sources=["computemin.cpp","computemin.cu"],
                   verbose=True)

def log_sinkhorn_sparse(sp_matrix, indices,a, b):
    def get_K( alpha, beta,reg,csr_value, col_, row_, crow, col, n,m):
        """log space computation"""
        f_ = torch.squeeze(alpha)
        g_ = torch.squeeze(beta)        
        data = torch.clone(csr_value)
        data -= g_[col_]
        data -= f_[row_]

        data /= -regk
        data = torch.exp(data)
        data = data.contiguous()
        sumf = torch.clone(f_)
        cuda_module.computesum(crow,col,n,m,sumf,data) 
        return torch.squeeze(sumf)      
    reg = 1e-6
    # epsilon0 = 1e6
    epsilon0 = 1e-4
    # numItermax=10
    numItermax=1
    stopThr=1e-9
    
    # print('indices.dtype: ', indices.dtype)
    
    xind = indices[0]
    yind = indices[1]
    # print('xind.dtype: ', xind.dtype)
    # print('yind.dtype: ', yind.dtype)
    
    # print('xind.shape:', xind.shape)

        
    indicesT = torch.clone(indices)
    temp = torch.clone(indicesT[0])
    indicesT[0] = torch.clone(indicesT[1])
    indicesT[1] = torch.clone(temp)

    # print('a.shape:', a.shape)
    # print('b.shape:', b.shape)
    # print('len(a):', len(a))
    # print('len(b):', len(b))

    n = len(a)
    m = len(b)
    
    device = torch.device('cuda:3')
    reg = torch.as_tensor(reg, dtype=torch.float64,device=device)
    a = torch.as_tensor(a, dtype=torch.float64,device=device)
    b = torch.as_tensor(b, dtype=torch.float64,device=device)
    a = torch.reshape(a,(n,1))
    b = torch.reshape(b,(1,m))
    
    
    regk = torch.as_tensor(epsilon0, dtype=torch.float64,device=device)
 
    loop = 1
    cpt = torch.as_tensor(1, dtype=torch.float64,device=device)
    cptt = 0
    err = 1
    
    f = torch.zeros_like(a , dtype=torch.float64,device=device)
    g = torch.zeros_like(b , dtype=torch.float64,device=device)
    

    def get_reg(n):  # exponential decreasing
        return (epsilon0 - reg) * torch.exp(-n) + reg
    # numIterinner = [200,500,500,500,800,1000]
    numIterinner = [20,500,500,500,800,1000]
    nn = len(numIterinner) - 1
    
    
    data = sp_matrix.coalesce().values()
    # print('sp_matrix.device: ', sp_matrix.device)
    # print('data.device: ', data.device)
    
    sp_matrix_data = torch.as_tensor(data, dtype=torch.float64,device=device)
    csr_data = sp_matrix.to_sparse_csr()
    csr_v = csr_data.values().contiguous()
    csr_v =  torch.as_tensor(csr_v, dtype=torch.float64,device=device)
    crow = torch.as_tensor( csr_data.crow_indices(), dtype=torch.int32,device=device)
    col = torch.as_tensor( csr_data.col_indices(), dtype=torch.int32,device=device)
    col_ = torch.as_tensor(col, dtype=xind.dtype,device=device)
    # print('crow.shape:', crow.shape)
    # print('crow[:10]:', crow[:10])
    # print('col.shape:', col.shape)
    # print('col_[:12]:', col_[:12])
    # print('csr_v.shape:', csr_v.shape)
    # print('g.shape:', g.shape)
    # input()
 
    sp_matrixT  = torch.sparse_coo_tensor(indicesT, sp_matrix_data, size=[m, n])

    csc_dataT = sp_matrixT.to_sparse_csr()
    csc_v = csc_dataT.values().contiguous()
    csc_v =  torch.as_tensor(csc_v, dtype=torch.float64,device=device)
    crowT = torch.as_tensor( csc_dataT.crow_indices(), dtype=torch.int32,device=device)
    colT = torch.as_tensor( csc_dataT.col_indices(), dtype=torch.int32,device=device)   
    colT_ = torch.as_tensor(colT, dtype=xind.dtype,device=device)
    # print('crowT.shape:', crowT.shape)
    # print('crowT[:10]:', crowT[:10])
    # print('colT_.shape:', colT_.shape)
    # print('colT_[:12]:', colT_[:12])
    # print('csc_v.shape:', csc_v.shape)
    # input()
    

    row_ = torch.zeros((len(col_),),dtype=torch.int64,device='cpu')
    assert (n+1) == len(crow)
    count = 0
    crowcpu = crow.cpu()
    for i in range(n):
        start = crowcpu[i]
        end = crowcpu[i+1]
        row_[start:end] = count
        count += 1
        
    rowT_ = torch.zeros((len(colT_),),dtype=torch.int64,device='cpu')
    assert (m+1) == len(crowT)
    count = 0
    crowTcpu = crowT.cpu()
    for i in range(m):
        start = crowTcpu[i]
        end = crowTcpu[i+1]

        rowT_[start:end] = count
        count += 1  
        



    def get_coo( alpha, beta,reg):
        """log space computation"""
        f_ = torch.squeeze(alpha)
        g_ = torch.squeeze(beta)        
        data = torch.clone(csr_v)
        data -= g_[col_]
        data -= f_[row_]
        data /= -regk
        data = torch.exp(data)
        out = torch.sparse_csr_tensor(crow,col,data,(n,m),device=device)
        out = out.to_sparse_coo()
        return out
    
    start0 = time.time()
    
    while loop:

        regk = get_reg(cpt)
        ind = [nn,cptt][cptt<=nn]
        # print('regk:', regk)
        print('numIterinner[ind]:', numIterinner[ind])
        start = time.time()
        
        for i in range(numIterinner[ind]):
            # print('-'*10, i, '-'*10)

            torch.cuda.empty_cache()
            g_ = torch.squeeze(g)
            data = torch.clone(csr_v)
            data -= g_[col_]
            f_ = torch.squeeze(f)
            min_f = torch.clone(f_)
            # print('min_f.shape:', min_f.shape)        
            cuda_module.computemin(crow,col,n,m,min_f,data) 
            # print('min_f[:10]:', min_f[:10])
            
            data -= min_f[row_]
            data /= -regk      
            data = torch.exp(data)
            # print('data[:10]:', data[:10])
            
            sumf = torch.clone(f_)                    
            cuda_module.computesum(crow,col,n,m,sumf,data) 
            # print('sumf[:10]:', sumf[:10])
            
            f = regk * torch.log(a) -\
            regk * torch.log(sumf.reshape(f.shape)) + min_f.reshape(f.shape) 
            # print('f[:10]:', f[:10])
 
            data = torch.clone(csc_v)
            f_ = torch.squeeze(f)
            data -= f_[colT_]  
            g_ = torch.squeeze(g)
            min_g = torch.clone(g_)
            cuda_module.computemin(crowT,colT,m,n,min_g,data) 
            # print('min_g[:10]:', min_g[:10])

            data -= min_g[rowT_]
            data /= -regk
            data = torch.exp(data)
            # print('data[:10]:', data[:10])
            
            sumg = torch.clone(g_)
            cuda_module.computesum(crowT,colT,m,n,sumg,data) 
            # print('sumg[:10]:', sumg[:10])
            
            g = regk * torch.log(b) -\
            regk * torch.log(sumg.reshape(g.shape)) +  min_g.reshape(g.shape)
            # print('g.shape:', g.shape)         
            # print('g[0,:10]:', g[0,:10])
            
            # input()
        end = time.time()
        print('time: ', int(cpt.item()), end-start)  
            
        if cpt % 1 == 0:
            print('---')
            res = get_K( f, g,regk ,csr_v, col_, row_, crow, col, n,m)
            res = torch.squeeze(res)
            res -= torch.squeeze(a)
            # print('res:', torch.norm(res ,2) **2)
            
            res2 = get_K( g, f,regk , csc_v, colT_, rowT_, crowT, colT, m,n)
            res2 = torch.squeeze(res2)
            res2 -= torch.squeeze(b)
            # print('res2:', torch.norm(res2 ,2) **2)
            
            err =  torch.norm(res ,2) **2 + torch.norm(res2 ,2) **2
            print('err: ', err)
            pass
        
        if err <= stopThr and cpt >= 50:
            loop = False
        if cpt >= numItermax:
            loop = False 
        #print(cpt)
        cpt += 1
        cptt += 1
        
        # input()
    
    end0 = time.time()
    print('time of while loop: ', end0-start0)  
    
    out = get_coo(f,g,regk)
    out = torch.detach(out)
   
    torch.cuda.empty_cache()
    return out