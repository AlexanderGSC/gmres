# GMRES Algorithm from Scratch

Implementation from scratch of the GMRES algorithm, using MGS with reorthogonalization and Householder reflections.

## 🚀 Project status
Currently in development. The goal is to have a clean and functional implementation from scratch.

## 📁 Structure
* `src/`: Contains the Arnoldi, Givens, and GMRES core modules.
* `CMakeLists.txt`: Configuration for compiling with CMake.

## 🛠️ Compilation If you have `gfortran` and `cmake`, you can compile it like this: 
```
bash mkdir build && cd build 
cmake .. 
make
