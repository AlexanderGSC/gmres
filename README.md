# 🚀 Krylov Lab: High Performance Iterative Methods from Scratch

This repository is a laboratory for Krylov Subspace Methods, implementing state-of-the-art iterative solvers from the ground up. The project focuses on the intersection of numerical stability, algorithmic efficiency, and hardware-aware optimization (targeting future FPGA/HLS implementation).

* Project Goal: To bridge the gap between theoretical numerical stability and hardware-efficient implementations of Krylov subspace methods.

## ✨ Key Features
* Diverse Solver Suite: GMRES: Robust nonsymmetric solver with Householder and MGS-R variants.BiCGSTAB: Efficient, short-memory nonsymmetric solver. Conjugate Gradient (CG): Optimal performance for Symmetric Positive Definite (SPD) systems.
* Matrix-Free Architecture: Operates exclusively via user-defined stencil operators, eliminating $\mathcal{O}(N^2)$ storage requirements.
* Numerical Precision: Householder-based GMRES achieving orthogonality limits of $\approx 10^{-30}$.
* Spectral Preconditioning: Adaptive Chebyshev polynomial preconditioner with spectral radius estimation via Lanczos iteration.

## 📊  Performance & Scaling Analysis
The laboratory includes a comprehensive benchmarking suite for Strong and Weak Scaling.
1. The Power of Preconditioning (GMRES 78k vars)
2. GMRES restart tuning
3. Weak and Strong Scaling of GMRES/BICGSTAB/CG.

## 🤖 Smart Restart Tuning (Experimental)
GMRES performance is highly sensitive to the restart parameter ($n$). This project explores an adaptive tuning approach:The "Efficiency Valley": Empirical data shows that for a 90k Poisson problem, the optimal $m$ lies around 95. Smaller values suffer from "stagnation," while larger values incur quadratic orthogonalization costs. 
* Future Work: Integration of a lightweight Neural Network (MLP) to predict the optimal $m$ based on initial residual decay and grid dimensions.

## 🚀 Roadmap
[x] Modified Gram Schmidt base implementation with reorthogonalization.

[x] Householder implementation (Walker '84).

[x] Chebyshev Matrix-Free Preconditioner.

[x] BICGSTAB and Conjugate Gradient Implementations.

[ ] Porting to CUDA/C++ (WIP).

## 📁 Structure
* `src/`: Contains the solvers core modules (Householder and MGS versions).
* `preconds/`: contains the code for implemented preconditioners.
* `problems/`: problems defined for testing the algorithm.
* `tests/`: Unitary tests (Poisson 2D, Hilbert matrix).
* `CMakeLists.txt`: Configuration for compiling with CMake.

## 📊 Validation Benchmarks
The solver is validated with classic problems known for their ill-conditioning:

* Hilbert matrices: famous for their extreme conditioning.
* Poisson 2D: typical stencil operators in computational fluid dynamics (CFD).
* Anisotropic Diffusion Equation (2D) (WIP)

## 🛠️ Compilation If you have `gfortran` and `cmake`, you can compile it like this: 
```
bash mkdir build && cd build 
cmake .. 
make
```
## 📚 References
* Trefethen, L. N., & Bau III, D. (1997). Numerical Linear Algebra. Society for Industrial and Applied Mathematics (SIAM).
* Golub, G. H., & Van Loan, C. F. Matrix Computations.
* Walker, H. F. (1984): Implementation of the GMRES method using Householder transformations.
* Saad, Y. (2003): Iterative Methods for Sparse Linear Systems.
* Van der Vorst, H. A. (1992): Bi-CGSTAB: A Fast and Smoothly Convergent Variant of Bi-CG.

