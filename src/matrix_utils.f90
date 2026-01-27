MODULE MATRIX_UTILS
    implicit none
    private
    public :: generate_poisson
    public :: generate_hilbert
    real(8)::diag = 4.0d0
    real(8)::ndiag=-1.0d0
CONTAINS
    subroutine generate_poisson(A,b,nsize)
        real(8),allocatable,intent(out)::A(:,:)
        real(8),allocatable,intent(out)::b(:)
        integer,intent(in)::nsize
        integer::i,j,row

        allocate(A(nsize*nsize,nsize*nsize), b(nsize*nsize))
        A = 0.0d0; b= 1.0d0
        do i=1,nsize
            do j=1,nsize
                row = i+(j-1)*nsize
                A(row,row) =  diag
                if (i > 1) A(row,row-1)       = ndiag
                if (i < nsize) A(row,row+1)   = ndiag
                if (j > 1) A(row-nsize,row)     = ndiag
                if (j < nsize) A(row+nsize,row) = ndiag
            end do
        end do
        b = matmul(A,b)
    end subroutine generate_poisson
    
    subroutine generate_hilbert(H,b,n)
        real(8),allocatable,intent(out)::H(:,:)
        real(8),allocatable,intent(out)::b(:)
        integer,intent(in)::n
        integer::i,j

        allocate(H(n,n),b(n))
        H = 0.0d0; b = 1.0d0
        do concurrent (i=1:n)
            do concurrent (j=1:n)
                H(i,j) = 1 / real(i+j-1)
            end do
        end do
        b = matmul(H,b)
    end subroutine generate_hilbert
END MODULE MATRIX_UTILS