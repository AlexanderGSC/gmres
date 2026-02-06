!!>This module implements the generate_matrix and stencil_vec methods 
!!> (defined in interfaces.f90):
!!> - generate_matrix: generates a dense matrix of size N^2 x N^2 for a Poisson problem 
!!>   with the classic five-point stencil.
!!> - stencil_vec: defines the stencil_vec operator for the Poisson problem 
!!>   (used in the matrix-free version)
MODULE poisson
    implicit none
    real(8)::diag = 4.0d0
    real(8)::ndiag=-1.0d0
contains

    subroutine generate_matrix(A,nsize)
        real(8),allocatable,intent(out)::A(:,:)
        integer,intent(in)::nsize
        integer::i,j,row

        allocate(A(nsize*nsize,nsize*nsize))
        A = 0.0d0;
        do i=1,nsize
            do j=1,nsize
                row = i+(j-1)*nsize
                A(row,row) =  diag
                if (i > 1) A(row,row-1)         = ndiag
                if (i < nsize) A(row,row+1)     = ndiag
                if (j > 1) A(row-nsize,row)     = ndiag
                if (j < nsize) A(row+nsize,row) = ndiag
            end do
        end do
    end subroutine generate_matrix


    SUBROUTINE stvec(x,y,n)   !y = A * x
        real(8), intent(in):: x(:)
        real(8), intent(out):: y(:)
        integer, intent(in)::n
        integer i,j,idx
        !$omp do collapse(2) private(idx)
        do j=2,n-1 !col
            do i=2,n-1 !row
                idx = i + (j-1)*n
                y(idx) = 4.0d0*x(idx)-1.0d0*(x(idx-1)+x(idx+1)+x(idx+n)+x(idx-n))
            end do
        end do
        !$omp end do
        !$omp do nowait
        do i=2,n-1 !borders col=1
            y(i)   = 4.0d0*x(i) - 1.0d0*(x(i-1)+x(i+1)+x(i+n)) !col=1
        end do
        !$omp do nowait private(idx)
        do i=2,n-1 !borders col=n
            idx = n*n-n+i
            y(idx) = 4.0d0*x(idx)-1.0d0*(x(idx-1)+x(idx+1)+x(idx-n))
        end do
        !$omp end do
        !$omp do nowait private(idx) !borders row=1
        do i=2,n-1
            idx  = (i-1)*n+1 
            y(idx) = 4.0d0*x(idx)-1.0d0*(x(idx+1)+x(idx+n)+x(idx-n))
        end do
        !$omp end do
        !$omp do nowait private(idx)
        do i=2,n-1
            idx  = (i-1)*n+n !row=n
            y(idx) = 4.0d0*x(idx)-1.0d0*(x(idx-1)+x(idx+n)+x(idx-n))
        end do
        !$omp end do
        !corners
        !$omp single
        y(1)  = 4.0d0*x(1) - 1.0d0*(x(2)+x(1+n))     !row=1, col=1
        y(n)  = 4.0d0*x(n) - 1.0d0*(x(n-1)+x(n+n))   !row=n, col=1
        idx   = n*(n-1)+1                                       
        y(idx)= 4.0d0*x(idx)-1.0d0*(x(idx+1)+x(idx-n)) !row=1, col=n
        idx   = n*(n-1)+n                                    
        y(idx)= 4.0d0*x(idx)-1.0d0*(x(idx-1)+x(idx-n)) !row=n, col=n                             
        !$omp end single
    END SUBROUTINE stvec

    SUBROUTINE stv_poisson(x,y,n)   !y = A * x
        real(8), intent(in):: x(:)
        real(8), intent(out):: y(:)
        integer, intent(in)::n
        integer i,j,idx
        !$omp do collapse(2), private(idx)
        do j=1,n
            do i=1,n
                idx = i + (j-1)*n
                y(idx) = 4.0 * x(idx)
                if (i > 1) y(idx) = y(idx) -1.0 * x(idx-1)
                if (i < n) y(idx) = y(idx) -1.0 * x(idx+1)
                if (j > 1) y(idx) = y(idx) -1.0 * x(idx-n)
                if (j < n) y(idx) = y(idx) -1.0 * x(idx+n)
            end do
        end do
        !$omp end do
    END SUBROUTINE stv_poisson
END MODULE poisson