# GMRES Algorithm from Scratch

This project implements the Generalized Minimal Residual (GMRES) algorithm from its fundamentals, exploring computational efficiency and numerical stability. It focuses on solving large and sparse linear systems, which are fundamental in High-Performance Computing (HPC) simulations.

## ✨ Key Features
* Matrix-Free Implementation: Does not require storing matrix $A$; operates using matrix-vector products (user-defined operator).
* Robust Stability: Implementation of Householder transformations following Walker's approach, ensuring a robust orthonormal basis.
* Implementation of the numerical core in C Fortran 2008 and OpenMP. Subsequent migration to CUDA/C++ (under construction).
* Use of a Chebyshev polynomial preconditioner to reduce the number of iterations until convergence. 

## 🚀 Project status
Currently in development. Once a functional and validated version of the algorithm has been obtained, GPU parallelization via CUDA will be explored.

## 🚀 Roadmap
[x] Modified Gram Schmidt base implementation with reorthogonalization.

[x] Householder implementation (Walker '84).

[x] Chebyshev Matrix-Free Preconditioner.

[ ] CUDA kernel for Block-Householder (WY Representation) and Stencil Vector products.

[ ] Integration in the Rosenbrock-Krylov scheme for reaction-diffusion.

## 📁 Structure
* `src/`: Contains the GMRES core modules (Householder and MGS versions).
* `preconds/`: contains the code for implemented preconditioners.
* `problems/`: problems defined for testing the algorithm.
* `tests/`: Unitary tests for GMRES (Poisson 2D, Hilbert matrix).
* `CMakeLists.txt`: Configuration for compiling with CMake.

## 📊 Benchmarks de Validación
The solver is validated with classic problems known for their ill-conditioning:

* Hilbert matrices: famous for their extreme conditioning.
* Poisson 2D: typical stencil operators in computational fluid dynamics (CFD).

## 🛠️ Compilation If you have `gfortran` and `cmake`, you can compile it like this: 
```
bash mkdir build && cd build 
cmake .. 
make
```
## 📚 References
* Walker, H. F. (1984). Implementation of the GMRES method using Householder transformations.
* Golub, G. H., & Van Loan, C. F. Matrix Computations.
