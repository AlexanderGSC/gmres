!!> This module contains all the interfaces used in this repository. 
!!> Specifically, it contains the following interfaces:
!!>   - stencil_vector: this interface defines the operation y = Ax. 
!!>     Each problem must offer a different implementation of this interface in 
!!>     order to work in GMRES Matrix free.
!!>   - precond: this interface defines the operation M_inv, which is the 
!!>     preconditioner on the left side of the system M_inv A x = M_inv b. 
!!>     Each preconditioner must implement this interface in order to work in 
!!>     preconditioned GMRES. 
MODULE Interfaces
    implicit none  !interface of the stencil vector prod 
        abstract interface 
        subroutine stencil_vector(x, y, n)
            real(8), intent(in) :: x(:)   
            real(8), intent(out):: y(:)  !y = A x
            integer, intent(in) :: n     !Size of the operator A
        end subroutine stencil_vector
    end interface
    abstract interface
        subroutine precond(A_x, r, z, aux, params, n) !interface of preconditioner
            procedure(stencil_vector) :: A_x    !Operator A x
            real(8), intent(in)  :: r(:)        !Residual
            real(8), intent(out) :: z(:)        !Preconditioned residual
            real(8), intent(out) :: aux(:)      !buffer for performance
            real(8), intent(in)  :: params(:)   !parameters of the algorithm
            integer, intent(in)  :: n           !Size of operator A
        end subroutine
    end interface
END MODULE