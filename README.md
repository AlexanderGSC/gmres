# GMRES Algorithm from Scratch

This project implements the Generalized Minimal Residual (GMRES) algorithm from its fundamentals, exploring computational efficiency and numerical stability. It focuses on solving large and sparse linear systems, which are fundamental in High-Performance Computing (HPC) simulations.

## ✨ Key Features
* Matrix-Free Implementation: Does not require storing matrix $A$; operates using matrix-vector products (user-defined operator).
* Robust Stability: Implementation of Householder transformations following Walker's approach, ensuring a robust orthonormal basis.
* Implementation of the numerical core in C Fortran 2008 and OpenMP. Subsequent migration to CUDA/C++ (under construction).
* Use of a Chebyshev polynomial preconditioner to reduce the number of iterations until convergence. 

## 🚀 Project status
Currently in development. Once a functional and validated version of the algorithm has been obtained, GPU parallelization via CUDA will be explored.

[x] Modified Gram Schmidt base implementation with reorthogonalization.

[x] Householder implementation (Walker '84).

[x] Chebyshev Matrix-Free Preconditioner.

[ ] CUDA kernel for Block-Householder (WY Representation) and Stencil Vector products.

[ ] Integration in the Rosenbrock-Krylov scheme for reaction-diffusion.

## 📊 Validation
The solver is validated with classic problems known for their ill-conditioning:

* Hilbert matrices: famous for their extreme conditioning.
* Poisson 2D: typical stencil operators in computational fluid dynamics (CFD).

## 📊 Performance Analysis: The Hilbert Stress-Test
The Hilbert matrix is a classic example of ill-conditioning ($rank(H) \approx N$ but eigenvalues decay exponentially). It is used here to test the limits of the Krylov subspace construction.Stability: Householder GMRES maintains an orthonormal basis $V$ such that $\|I - V^T V\| \approx 10^{-29}$, effectively reaching the theoretical limit of numerical stability.MGS Weakness: Even with reorthogonalization, MGS remains 15 orders of magnitude less precise in basis maintenance.Efficiency: The computational overhead of Householder reflections is $<5\%$ compared to MGS, making it the superior choice for high-reliability scientific computing.
All tests were performed with a convergence tolerance of $1.0 \times 10^{-15}$.Matrix Size

## 📁 Structure
* `src/`: Contains the GMRES core modules (Householder and MGS versions).
* `preconds/`: contains the code for implemented preconditioners.
* `problems/`: problems defined for testing the algorithm.
* `tests/`: Unitary tests for GMRES (Poisson 2D, Hilbert matrix).
* `CMakeLists.txt`: Configuration for compiling with CMake.

## 🛠️ Compilation If you have `gfortran` and `cmake`, you can compile it like this: 
```
bash mkdir build && cd build 
cmake .. 
make
```
## 📚 References & Further Reading
This implementation is built upon the theoretical foundations of Numerical Linear Algebra and Iterative Methods. The following works are the primary references for the stability analysis and algorithmic structures used:

### Core Algorithm & Stability
* Walker, H. F. (1984). Implementation of the GMRES method using Householder transformations. SIAM Journal on Scientific and Statistical Computing, 5(4), 942-953.

(The primary basis for the Householder implementation in this project).

* Saad, Y., & Schultz, M. H. (1986). GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems. SIAM Journal on Scientific and Statistical Computing, 7(3), 856-869.

(The original paper introducing the GMRES method).

* Saad, Y. (2003). Iterative Methods for Sparse Linear Systems (2nd Ed.). Society for Industrial and Applied Mathematics.

(The comprehensive guide for Krylov subspace methods and preconditioning).

### Numerical Linear Algebra Foundations
* Trefethen, L. N., & Bau III, D. (1997). Numerical Linear Algebra. Society for Industrial and Applied Mathematics.

(Essential reference for the stability of Householder reflections vs. Gram-Schmidt).

* Golub, G. H., & Van Loan, C. F. (2013). Matrix Computations (4th Ed.). Johns Hopkins University Press.

(The definitive "Bible" for matrix algorithms and Block-Householder representations).
