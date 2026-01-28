MODULE PRECOND
    implicit none
    private
    public :: jacobi
CONTAINS
    SUBROUTINE jacobi(A,b)
        real(8),intent(inout)::A(:,:)
        real(8),intent(inout)::b(:)
        integer::i
        real(8)::d
        do i=1,size(A,1)
            d      = A(i,i)
            b(i)   = b(i) / d
            A(i,:) = A(i,:) / d
            A(i,i) = 1.0d0
        end do
    END SUBROUTINE jacobi
END MODULE PRECOND   