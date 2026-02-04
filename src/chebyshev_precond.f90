MODULE CHEBYSHEV_PRECOND
    use interfaces
    implicit none
    private
    public :: cbpr2

CONTAINS
    SUBROUTINE cbpr2(A_x, r, z, aux, params, n)
        procedure(stencil_vector) :: A_x
        real(8), intent(in)    :: r(:) !residual
        real(8), intent(out)   :: z(:) !precond residual
        real(8), intent(out)   :: aux(:) !buffer as param for better perfeormance
        real(8), intent(in)    :: params(:) !eigenvalues
        integer, intent(in)    :: n
        !other variables
        real(8) :: c, d, alpha, beta, eigen_min, eigen_max
        !-----------------------------------------------------------
        eigen_min = params(1); eigen_max = params(2)
        c     = (eigen_max - eigen_min) / 2.0d0
        d     = (eigen_max + eigen_min) / 2.0d0
        alpha = 1.0d0 / d
        beta  = (c * alpha / 2.0d0)**2
        alpha = 1.0d0 / (d - beta)
        z = r / d
        call A_x(z,aux,n)
        z = z + alpha*(r - aux)
    END SUBROUTINE cbpr2
END MODULE Chebyshev_PRECOND   