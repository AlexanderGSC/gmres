MODULE hilbert
    implicit none
    private
    public :: generate_matrix
CONTAINS
    subroutine generate_matrix(H,n)
        real(8),allocatable,intent(out)::H(:,:) !Hilbert's matrix
        integer,intent(in)::n                   !size of the matrix
        integer::i,j

        allocate(H(n,n))
        H = 0.0d0;
        do concurrent (i=1:n)
            do concurrent (j=1:n)
                H(i,j) = 1 / real(i+j-1)
            end do
        end do
    end subroutine generate_matrix
END MODULE hilbert