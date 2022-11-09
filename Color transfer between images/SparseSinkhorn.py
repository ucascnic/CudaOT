# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 11:08:23 2021

@author: Administrator
"""
import os.path
import sys
sys.path.append('/home/chenyidong/keops')
import numpy as np
import torch
from random import choices

from pykeops.torch import LazyTensor
from pykeops.torch.cluster import grid_cluster
import pandas as pd
from pykeops.torch.cluster import cluster_ranges_centroids
from pykeops.torch import generic_argmin
from pykeops.torch import generic_sum,generic_logsumexp,generic_argmin
use_cuda = torch.cuda.is_available()
dtype = torch.cuda.DoubleTensor if use_cuda else torch.FloatTensor
def sinkhorn_on_log_domain_GPU_keops(mu, nv, x,y,  reg, epsilon0 = 1e6, numItermax=10,stopThr=1e-9,debug=0,ranges = None,
                                    warmf=None,warmg=None):
    
    

    device = torch.device('cuda')
    reg = torch.as_tensor(reg, device=device).type(dtype)


    regk = torch.as_tensor(epsilon0,  device=device).type(dtype)
 
    loop = 1
    cpt = torch.as_tensor([1],  device=device).type(dtype)
    
    one = torch.as_tensor([1.0],  device=device).type(dtype)
    cptt = 0
    err = 1
    if warmf is not None:
        f = torch.as_tensor(warmf  ).type(dtype)
        g = torch.as_tensor(warmg ).type(dtype)
    else:
        
        f = torch.zeros_like(mu ,device=device).type(dtype)
        g = torch.zeros_like(nv ,device=device).type(dtype)
    
 
    #if debug:
    #    print(f)
    #    return f,f,f
    def get_reg(n):  # exponential decreasing
        return  (epsilon0 - reg) * torch.exp(-n) + reg 
    
    
    nearest_neighbor = generic_argmin(
    'SqNorm2(x-y)-g',   # Formula
    'a = Vi(1)',      # Output: 1 scalar per line
    'x = Vi(3)',    
    'y = Vj(3)',
    'g = Vj(1)',dtype  = 'float64'  )   
    
    nearest_neighbor2 = generic_argmin(
    'SqNorm2(x-y)-f',   # Formula
    'a = Vj(1)',      # Output: 1 scalar per line
    'x = Vi(3)',    
    'y = Vj(3)',
    'f = Vi(1)',dtype  = 'float64'  )    
    
    log_likelihood = generic_logsumexp(
       '(-eps * (SqNorm2(x-y)-g))  ', # Formula
        
        'a = Vi(1)',              # Output: 1 scalar per line

        'x = Vi(3)',              # 1st input: dim-1 vector per line

        'y = Vj(3)', 
        
        'g = Vj(1)',
        
        'eps = Pm(1)',dtype  = 'float64'  )               

    log_likelihood2 = generic_logsumexp(
       '(-eps * (SqNorm2(x-y)-g-minf))  ', # Formula
        
        'a = Vi(1)',              # Output: 1 scalar per line

        'x = Vi(3)',              # 1st input: dim-1 vector per line

        'y = Vj(3)', 
        
        'g = Vj(1)',
        
        'minf = Vi(1)',
        
        'eps = Pm(1)',dtype  = 'float64'  )  
        
    log_likelihood_trans = generic_logsumexp(
       '(-eps * (SqNorm2(x-y)-f))  ', # Formula
        
        'a = Vj(1)',              # Output: 1 scalar per line

        'x = Vi(3)',              # 1st input: dim-1 vector per line

        'y = Vj(3)', 
        
        'f = Vi(1)',
        
        'eps = Pm(1)',dtype  = 'float64'  )   
    
    log_likelihood_trans2 = generic_logsumexp(
       '(-eps * (SqNorm2(x-y)-f-ming))  ', # Formula
        
        'a = Vj(1)',              # Output: 1 scalar per line

        'x = Vi(3)',              # 1st input: dim-1 vector per line

        'y = Vj(3)', 
        
        'f = Vi(1)',         
        'ming = Vj(1)',
        
        'eps = Pm(1)',dtype  = 'float64'  )    
    sums = generic_logsumexp(
       ' ( (-eps * (SqNorm2(x-y)-f-g)))  ', # Formula
        
        'a = Vi(1)',              # Output: 1 scalar per line

        'x = Vi(3)',              # 1st input: dim-1 vector per line

        'y = Vj(3)',
        
        'f = Vi(1)', 
        
        'g = Vj(1)',
        
        'eps = Pm(1)',dtype  = 'float64'  ) 
    sumsy = generic_logsumexp(
       ' ( (-eps * (SqNorm2(x-y)-f-g)))  ', # Formula
        
        'a = Vj(1)',              # Output: 1 scalar per line

        'x = Vi(3)',              # 1st input: dim-1 vector per line

        'y = Vj(3)',
        
        'f = Vi(1)', 
        
        'g = Vj(1)',
        
        'eps = Pm(1)',dtype  = 'float64' ) 
    numIterinner = [5,10,10,50,50,100]
    nn = len(numIterinner) - 1
    
 
    while loop:
 
        
        regk = torch.max(regk*0.1,reg)
        regk = get_reg(cpt)

        regk1 = one/regk
        ind = [nn,cptt][cptt<=nn]
 
        for i in range(numIterinner[ind]): 
             
            a_ = nearest_neighbor(x,y,g, ranges = ranges )
            temp = x - y[ a_.view(-1).long()]
            #print(x.shape)
            #print(y.shape)
            #print(f.shape)
            #print(temp.shape)
            min_f = torch.multiply(temp,temp).sum(dim=1, keepdim = True)  - g[ a_.view(-1).long()]
            a = log_likelihood2(x,y,g,min_f,regk1, ranges = ranges )
 
            f = regk * torch.log(mu) - regk * a  + min_f
            #f = f*0
            b_ = nearest_neighbor2(x,y,f, ranges = ranges ) 
            
            temp_ = x[ b_.view(-1).long()] -  y
            min_g = (temp_ ** 2).sum(1, keepdim = True) - f[ b_.view(-1).long()]
            #min_g = torch.multiply(temp,temp).sum(dim=0, keepdim = True) 
            #min_g = min_g.T
            #return  f,min_g,regk,temp_
            b = log_likelihood_trans2(x,y,f,min_g, regk1, ranges = ranges )   
 
            g = regk * torch.log(nv) - regk * b + min_g
 
        #if debug:
            #print(regk * torch.log(nv)- regk * b)
            #break
        if cpt % 10 == 0:
            
            res = torch.exp(sums(x,y,f,g,regk1, ranges = ranges )) 
            res2 = torch.exp(sumsy(x,y,f,g,regk1, ranges = ranges )) 
            print(torch.norm(res - mu)+torch.norm(res2 - nv))
            if debug:
                #print(res)
                print(sum(res))
                print(sum(res2))
                #print(torch.norm(res - mu))
                #print(torch.norm(res2 - nv))
                break 
            
    
            pass
        if err <= stopThr and cpt >= 50:

            loop = False

        if cpt >= numItermax:
            loop = False 
        #print(cpt)
        cpt += 1
        cptt += 1

  
    return  f,g,regk


def sinkhorn_on_log_domain_GPU(a, b, M, reg, epsilon0 = 1e6, numItermax=10,stopThr=1e-9):
    
    
    n = len(a)
    m = len(b)
    device = torch.device('cuda')
    reg = torch.as_tensor(reg, dtype=torch.float64,device=device)
    a = torch.as_tensor(a, dtype=torch.float64,device=device)
    b = torch.as_tensor(b, dtype=torch.float64,device=device)
    a = torch.reshape(a,(n,1))
    b = torch.reshape(b,(1,m))
    M = torch.as_tensor(M, dtype=torch.float64,device=device)
    regk = torch.as_tensor(epsilon0, dtype=torch.float64,device=device)
 
    loop = 1
    cpt = torch.as_tensor(1, dtype=torch.float64,device=device)
    cptt = 0
    err = 1
    
    f = torch.zeros_like(a , dtype=torch.float64,device=device)
    g = torch.zeros_like(b , dtype=torch.float64,device=device)
    
    def get_reg(n):  # exponential decreasing
        return (epsilon0 - reg) * torch.exp(-n) + reg
    
    def get_K(alpha, beta,reg):
        """log space computation"""
        return torch.exp(-(M - alpha - beta) / reg)
 
    
    numIterinner = [5,10,10,50,50,200]
    nn = len(numIterinner) - 1
    while loop:
 
        regk = torch.max(regk*0.1,reg)
        regk = get_reg(cpt)
        ind = [nn,cptt][cptt<=nn]
 
        for i in range(numIterinner[ind]): 
            min_f,_ = torch.min(M - g, dim = 1,  keepdim = True)
            #print(min_f)
            #break
            f = regk * torch.log(a) -\
            regk * torch.log (torch.sum (torch.exp( -( M - g - min_f ) /regk ) ,dim = 1,keepdim = True ) )+ min_f
            
      
            min_g,_ = torch.min(M - f, dim = 0, keepdim = True)
            g = regk * torch.log(b) -\
            regk * torch.log (torch.sum (torch.exp( -( M - f  - min_g ) /regk) ,dim = 0,keepdim = True ) ) + min_g
            #print(g)
        #print(b.shape)
        #break      
        if cpt % 1 == 0:
 
            res = get_K(f,g,regk)
            m = res.sum(0)
            err2 =  torch.norm(res.sum(0) - a ,2) **2
            print(err2)
            #break
            #return res,res,res
            pass
        if err <= stopThr and cpt >= 50:

            loop = False

        if cpt >= numItermax:
            loop = False 
        print(cpt)
        cpt += 1
        cptt += 1

    #res = get_K(f,g,regk)
    return f,g,regk
def sinkhorn_on_log_domain_GPU_keops_lazyM(M,mu, nv, x,y,finit,ginit,  reg, epsilon0 = 1e6, numItermax=10,stopThr=1e-9):
    
    
    n = len(mu)
    m = len(nv)
    device = torch.device('cuda')
    reg = torch.as_tensor(reg, device=device).type(dtype)
    regk = torch.as_tensor(epsilon0,  device=device).type(dtype)
    loop = 1
    cpt = torch.as_tensor([1],  device=device).type(dtype)
    one = torch.as_tensor([1.0],  device=device).type(dtype)
    cptt = 0
    err = 1
    
    finit = torch.zeros_like(mu ,device=device).type(dtype)
    ginit = torch.zeros_like(nv ,device=device).type(dtype)
    
    f_ = torch.as_tensor(finit ,device=device).type(dtype)
    g_ = torch.as_tensor(ginit ,device=device).type(dtype)
    f =  f_
    g =  g_
     

    x_i, y_j = LazyTensor(x[:, None, :]), LazyTensor(y[:, None, :])

    D_ij = torch.sum((x[:, None, :] - y[None, :, :]) ** 2, 2)
    #D_ij =  torch.as_tensor(M ,device=device).type(dtype)
    sums = generic_logsumexp(
       ' ( (-eps * (SqNorm2(x-y)- f - g )))  ', # Formula
        
        'a = Vi(1)',              # Output: 1 scalar per line

        'x = Vi(3)',              # 1st input: dim-1 vector per line

        'y = Vj(3)',
        
        'f = Vi(1)', 
        
        'g = Vj(1)',
        
        'eps = Pm(1)',dtype  = 'float64'  ) 
 

    
    
    def get_reg(n):  # exponential decreasing
        return  (epsilon0 - reg) * torch.exp(-n) + reg 
 
 
    numIterinner = [5,10,10,50,50,100]
    nn = len(numIterinner) - 1
    
    N = f.shape[0]
    bb =   torch.ones(N, 1).type(dtype)
 
    while loop:
 
        regk = torch.max(regk*0.1,reg)
        regk = get_reg(cpt)
        regk1 = one/regk
        ind = [nn,cptt][cptt<=nn]
        #print(regk)
        print(cpt)
        for i in range(numIterinner[ind]): 
            #print(D_ij.shape)
            #print(g.shape)
            K1 =   ( -regk1 * (  D_ij  -  g.t()    ) ).exp()
            #K1.ranges = ranges_ij
            print((  D_ij  -  g.t()    ).min())
            break
            temp = K1.sum(dim=1)
            a =  torch.log(temp).unsqueeze(1)  
            
             
            f_ =  regk * torch.log(mu) - regk * (a)
 
            f =  f_

            K2 =   ( -regk1 * ( (D_ij) - f   ) ).exp()
            #K2.ranges = ranges_ij
            #b =  K2.sum(dim=0) 
            b =  torch.log(K2.sum(dim=0)).unsqueeze(1)
            g_ =    regk * torch.log(nv) - regk * b
            #print(b.shape)
            #print( regk * torch.log(nv)[0:3])
            #break
            g =  g_            
            #print(g)
            cpt += 1
            print(cpt)
    
        #break
        print(reg)
        if cpt % 1 == 0:
            #print(f_.shape)
            res1 = torch.exp(sums(x,y,f_,g_,regk1))
            res2 = torch.exp(sums(y,x,g_,f_,regk1))
            m = res2-nv 
            print(torch.norm(m))
            #break
            #print(res1)
            #print(g_)             
 
            #break
            
            pass
        if err <= stopThr and cpt >= 50:

            loop = False

        if cpt >= numItermax:
            loop = False 
        
        
        cptt += 1

  
    return  f_,g_,regk
def sinkhorn_on_log_domain_GPU_keops_lazy(ranges_ij,mu, nv, x,y,finit,ginit,  reg, epsilon0 = 1e6, numItermax=10,stopThr=1e-9,debug=0):
    
    
    n = len(mu)
    m = len(nv)
    device = torch.device('cuda')
    reg = torch.as_tensor(reg, device=device).type(dtype)
    regk = torch.as_tensor(epsilon0,  device=device).type(dtype)
    loop = 1
    cpt = torch.as_tensor([1],  device=device).type(dtype)
    one = torch.as_tensor([1.0],  device=device).type(dtype)
    cptt = 0
    err = 1
    
    

    
    
    if len(ranges_ij):
        f_ = torch.as_tensor(finit ,device=device).type(dtype)
        g_ = torch.as_tensor(ginit ,device=device).type(dtype)
    else:
        f_ = torch.zeros_like(mu ,device=device).type(dtype)
        g_ = torch.zeros_like(nv ,device=device).type(dtype)        
    f = LazyTensor(f_[:, None])
    g = LazyTensor(g_[:, None])  
 

    x_i, y_j = LazyTensor(x[:, None, :]), LazyTensor(y[:, None, :])

 
    D_ij = ( ( x_i - y_j.t() ) ** 2 ).sum(dim=2)  
    sums = generic_logsumexp(
       ' ( (-eps * (SqNorm2(x-y)- f - g )))  ', # Formula
        
        'a = Vi(1)',              # Output: 1 scalar per line

        'x = Vi(3)',              # 1st input: dim-1 vector per line

        'y = Vj(3)',
        
        'f = Vi(1)', 
        
        'g = Vj(1)',
        
        'eps = Pm(1)',dtype  = 'float64'  ) 
 

    
    
    def get_reg(n):  # exponential decreasing
        return  (epsilon0 - reg) * torch.exp(-n) + reg 
 
 
    numIterinner = [5,10,10,10,10,150]
    nn = len(numIterinner) - 1
    
    N = f.shape[0]
    bb =   torch.ones(N, 1).type(dtype)
 
    while loop:
 
        regk = torch.max(regk*0.1,reg)
        regk = get_reg(cpt)
        regk1 = one/regk
        ind = [nn,cptt][cptt<=nn]
        print("reg = ")
        print(regk)
        for i in range(numIterinner[ind]): 
            
            #fmin_ =  (D_ij - g.t()).min(dim=1)
            #print(fmin_)
            #fmin = LazyTensor(fmin_[:, None])
            #K1 =   (  regk1 * ( fmin + g.t() - D_ij ) ).exp()
            K1 =    regk1 * (   g.t() - D_ij ) 
            if len(ranges_ij):
                K1.ranges = ranges_ij
            
            a =  K1.logsumexp(dim=1)
            
            f_ =  regk * torch.log(mu) - regk * (a)
            f = LazyTensor(f_[:, None])
            K2 =   (  regk1 * (   f - D_ij ) ).exp()
            if len(ranges_ij):                
                K2.ranges = ranges_ij
            b =  torch.log(K2.sum(dim=0))
            g_ =    regk * torch.log(nv) - regk * b

            g = LazyTensor(g_[:, None])               
        
     
        if debug:
            print(a)
            break
        
        if cpt % 5 == 0:
            res1 = torch.exp(sums(x,y,f_,g_,regk1)) 
            res2 = torch.exp(sums(y,x,g_,f_,regk1)) 
            print("residual = ")
            print(torch.norm(res2 - nv))             
 
            #break
            pass
        if err <= stopThr and cpt >= 50:

            loop = False

        if cpt >= numItermax:
            loop = False 
        print(cpt)
        cpt += 1
        cptt += 1

  
    return  f_,g_,regk
