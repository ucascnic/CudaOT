# CudaOT
Repository for solving discrete optimal transport problems via  Cplex, Sinkhorn, FastEMD and the GPU implementation of the Multi-Scale Sparse Sinkhorn.

This library is designed to solve a very large instance of discrete optimal transport problems.

## Environment Requirements

**Programming Language:** CUDA C/C++ (tested using cuda/11.1)

**eigen3** library headers: These can be installed via the package manager on most distributions, e.g. via package libeigen3-dev on Ubuntu.

## Installation Instructions

For unified memory implementation

(1) In the `CMakeLists.txt`, edit the variable `CUDA_INSTALL_PATH` to match the CUDA installation directory. Edit the variable `CPLEX_LIBRARY` to match the CPLEX installation directory for building TransportNetwork.

(2) Type `cmake .` and  `make`  to compile.

## Examples

(1)Cplex: The Cplex solver computes the unregularized OT problems. 

(2)FastEMD: The details of FastEMD can be found at: https://github.com/tillhainbach/FastEMD

(3)Multi-scale OT: The CPU version of the multi-scale Sinkhorn Algorithm. The Document could be found at: https://github.com/bernhard-schmitzer/MultiScaleOT 

(4)M3S: The implementation of the multi-scale sparse Sinkhorn. More details will be released later. 

(5)LogDomainSK: Implemented the Sinkhorn algorithm on log domain.

(6)FastTransport: Fast Network Simplex for Optimal Transport: https://github.com/nbonneel/network_simplex

(7)3D Wasserstein Distance: code/3D Wasserstein Distance/plot3D_Wasserstein.m

(8)Color transfer between images: code/Color transfer between images/Multi_scale_Color_Transfer.ipynb


## Reference

[1] I. Corporation, Ibm ilog cplex optimizer, (2013), http://www.ilog.com/products/cplex/.

[2] R. Flamary and N. Courty, Pot python optimal transport library, 2017, https://pythonot.github.io/.

[3] B. Schmitzer and C. Schnorr, A hierarchical approach to optimal transport, vol. 7893, 06 2013, pp. 452–571464, https://doi.org/10.1007/978-3-642-38267-3 38.

[4] B. Schmitzer, Stabilized sparse scaling algorithms for entropy regularized transport problems, SIAM  Journal on Scientific Computing, 41 (2016).

[5] Computational optimal transport: Complexity by accelerated gradient descent is better than by sinkhorn’s algorithm, in Proceedings of the 35th International Conference on Machine Learning(ICML) in PMLR, 02 2018, pp. 80:1367–1376.

[6] N. Bonneel, M. van de Panne, S. Paris, and W. Heidrich, Displacement Interpolation Using Lagrangian Mass Transport, ACM Transactions on Graphics (SIGGRAPH ASIA 2011), 30 (2011).

[7] Pele O., Werman M. (2008) A Linear Time Histogram Metric for Improved SIFT Matching. In: Forsyth D., Torr P., Zisserman A. (eds) Computer Vision – ECCV 2008. ECCV 2008. Lecture Notes in Computer Science, vol 5304. Springer, Berlin, Heidelberg

[8] O. Pele and M. Werman, "Fast and robust Earth Mover's Distances," 2009 IEEE 12th International Conference on Computer Vision, Kyoto, 2009, pp. 460-467.

[9] Charlier, B., Feydy, J., Glaunès, J. A., Collin, F.-D. & Durif, G. Kernel Operations on the GPU, with Autodiff, without Memory Overflows. Journal of Machine Learning Research 22, 1–6 (2021).

[10] Jean Feydy, Thibault Séjourné, Franc¸ois-Xavier Vialard, Shun-ichi Amari, Alain Trouve, and Gabriel Peyré. Interpolating between optimal transport and mmd using sinkhorn divergences. In The 22nd International Conference on Artificial Intelligence and Statistics, pages 2681–2690, 2019.