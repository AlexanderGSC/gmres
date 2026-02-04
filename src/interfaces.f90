MODULE Interfaces
    implicit none
        abstract interface 
        subroutine stencil_vector(x, y, n)
            real(8), intent(in) :: x(:)   
            real(8), intent(out):: y(:)  !y = A x
            integer, intent(in) :: n     !Size of the operator A
        end subroutine stencil_vector
    end interface
    abstract interface
        subroutine precond(A_x, r, z, aux, params, n)
            procedure(stencil_vector) :: A_x    !Operator A x
            real(8), intent(in)  :: r(:)        !Residual
            real(8), intent(out) :: z(:)        !Preconditioned residual
            real(8), intent(out) :: aux(:)      !buffer for performance
            real(8), intent(in)  :: params(:)   !parameters of the algorithm
            integer, intent(in)  :: n           !Size of operator A
        end subroutine
    end interface
END MODULE