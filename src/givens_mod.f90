!This module performs Givens rotations of a Matrix. The Givens is performed
!on i-th, i+1-th row row and j-th column and returns
!Rotated Matrix
!
MODULE Givens_Mod
    IMPLICIT none
    CONTAINS
    SUBROUTINE givens(A,b,i,j)
        integer, intent(in)::i,j
        real(8), intent(inout)::A(:,:),b(:)
        real(8)::c,s,d !sin, cos, distance
        real(8)::tmp_a, tmp_b
        integer::k
        d = hypot(A(i,j),A(i+1,j))
        c = A(i,j)/d
        s = A(i+1,j)/d
        do concurrent (k=j:size(A,2)) !bycol
            tmp_a   = A(i,k)
            A(i,k)  = c*tmp_a+s*A(i+1,k)
            A(i+1,k)=-s*tmp_a+c*A(i+1,k)
        end do
        A(i+1,j) = 0.0
        tmp_b = b(i)
        b(i)  = c*tmp_b+s*b(i+1)
        b(i+1)=-s*tmp_b+c*b(i+1)
    END SUBROUTINE
END MODULE Givens_Mod
