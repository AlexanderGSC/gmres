MODULE conjugate_gradient
    use interfaces
    use omp_lib
    implicit none
    private
    public :: cg
    public :: pcg
    public :: cg_omp
    public :: pcg_omp
CONTAINS
    subroutine cg(Ax_op, b, x, tol, iter, res)  !Solver Ax=b
        procedure(stencil_vector) :: Ax_op          !Stencil vector operator
        real(8), intent(in) :: b(:)                 !rhs of the system
        real(8), allocatable, intent(out):: x(:)    !Solution
        real(8), intent(in) :: tol                  !tolerance
        integer, intent(inout) :: iter              !Maximun number of iterations
        real(8), intent(out):: res                  !Residual
        !Other variables
        real(8), allocatable :: ax(:)       ! ax stores the Ax product 
        real(8), allocatable :: p(:)        ! gradient direction p
        real(8), allocatable :: r(:)
        real(8) :: alpha, beta, rr
        integer :: i,n,grid_size
        n = size(b)
        grid_size = int(sqrt(real(n)))
        allocate(x(n),ax(n),p(n),r(n))
        x = 0.0d0; r = b; p = r;    !x0 = 0; r0 = b; p0 = r0
        do i=1,iter
            call Ax_op(p,ax,grid_size)
            rr    = dot_product(r,r)        !rr = (r.T,r)      !
            alpha = rr / dot_product(ax,p)  !step length
            x     = x + alpha * p           !approximate solution
            r     = r - alpha * ax          !residual vector
            res   = norm2(r)                !residual norm as the error
            beta  = dot_product(r,r) / rr   !improvement this step
            p     = r + beta * p            !search direction
            if (res < tol) then
                iter = i
                exit
            endif
        end do
    end subroutine cg

    subroutine pcg(Ax_op, b, x, tol, iter, res, M_inv, params)  !Solver Ax=b
        procedure(stencil_vector) :: Ax_op          !Stencil vector operator
        real(8), intent(in) :: b(:)                 !rhs of the system
        real(8), allocatable, intent(out):: x(:)    !Solution
        real(8), intent(in) :: tol                  !tolerance
        integer, intent(inout) :: iter              !Maximun number of iterations
        real(8), intent(out):: res                  !Residual
        procedure(precond) :: M_inv                 !Preconditioner
        real(8), intent(in) :: params(:)            !Precond params
        !Other variables
        real(8), allocatable :: ax(:)       ! ax stores the Ax product 
        real(8), allocatable :: p(:)        ! gradient direction p
        real(8), allocatable :: r(:)        ! residual vector
        real(8), allocatable :: z(:), aux(:)! precond vectors
        real(8) :: alpha, beta, rr
        integer :: i,n,grid_size
        n = size(b)
        grid_size = int(sqrt(real(n)))
        allocate(x(n),ax(n),p(n),r(n),z(n),aux(n))
        x = 0.0d0; r = b;    !x0 = 0; r0 = b; p0 = r0
        call M_inv(Ax_op,r,z,aux,params,grid_size)   !M^-1 r = z
        p = z;
        do i=1,iter
            call Ax_op(p,ax,grid_size)
            rr    = dot_product(r,z)        !rr = (r.T,z)      !
            alpha = rr / dot_product(ax,p)  !step length
            x     = x + alpha * p           !approximate solution
            r     = r - alpha * ax          !residual vector
            res   = norm2(r)                !residual norm as the error
            call M_inv(Ax_op,r,z,aux,params,grid_size) !precond M^-1 r = z
            beta  = dot_product(r,z) / rr   !improvement this step
            p     = z + beta * p            !search direction
            if (res < tol) then
                iter = i
                exit
            endif
        end do
    end subroutine pcg

    subroutine cg_omp(Ax_op, b, x, tol, iter, res)  !Solver Ax=b
        procedure(stencil_vector) :: Ax_op          !Stencil vector operator
        real(8), intent(in) :: b(:)                 !rhs of the system
        real(8), allocatable, intent(out):: x(:)    !Solution
        real(8), intent(in) :: tol                  !tolerance
        integer, intent(inout) :: iter              !Maximun number of iterations
        real(8), intent(out):: res                  !Residual
        !Other variables
        real(8), allocatable :: ax(:)       ! ax stores the Ax product 
        real(8), allocatable :: p(:)        ! gradient direction p
        real(8), allocatable :: r(:)        ! residual vector
        logical :: converged
        real(8) :: alpha, beta, rr
        integer :: i,j,n,grid_size
        n = size(b)
        grid_size = int(sqrt(real(n)))
        allocate(x(n),ax(n),p(n),r(n))
        converged = .false.
        !$omp parallel
        !$omp do
            do j=1,n !x0 = 0; r0 = b; p0 = r0
                x(j) = 0.0d0
                r(j) = b(j)
                p(j) = b(j)
            end do
        !$omp end do   
        do i=1,iter
            if (converged) cycle
            call Ax_op(p,ax,grid_size)
            !$omp single
                rr    = 0.0d0
                alpha = 0.0d0
                res   = 0.0d0
                beta  = 0.0d0
            !$omp end single
            !$omp do reduction(+:rr,alpha) 
                do j=1,n
                    rr    = rr + r(j) * r(j)
                    alpha = alpha + ax(j) * p(j)
                end do 
            !$omp end do
            !$omp single
                alpha = rr / alpha
            !$omp end single
            !$omp do reduction(+:res,beta)
                do j=1,n
                    x(j) = x(j) + alpha * p(j)
                    r(j) = r(j) - alpha * ax(j)
                    res  = res  + r(j) * r(j)
                    beta = beta + r(j) * r(j)
                end do
            !$omp end do
            !$omp single
                res  = sqrt(res)
                beta = beta / rr
            !$omp end single
            !$omp do
                do j=1,n
                    p(j) = r(j) + beta * p(j)
                end do
            !$omp end do
            !$omp single
                if (res < tol) then
                    converged = .true.
                    iter = i
                endif
            !$omp end single
        end do
        !$omp end parallel 
    end subroutine cg_omp

    subroutine pcg_omp(Ax_op, b, x, tol, iter, res, M_inv, params)  
        procedure(stencil_vector) :: Ax_op          !Stencil vector operator
        real(8), intent(in) :: b(:)                 !rhs of the system
        real(8), allocatable, intent(out):: x(:)    !Solution
        real(8), intent(in) :: tol                  !tolerance
        integer, intent(inout) :: iter              !Maximun number of iterations
        real(8), intent(out):: res                  !Residual
        procedure(precond) :: M_inv                 !Preconditioner
        real(8), intent(in) :: params(:)            !Preconditioner params
        !Other variables
        real(8), allocatable :: ax(:)       ! ax stores the Ax product 
        real(8), allocatable :: p(:)        ! gradient direction p
        real(8), allocatable :: r(:)        ! residual vector
        real(8), allocatable :: z(:),aux(:) ! Auxiliary vectors for precond
        logical :: converged
        real(8) :: alpha, beta, rr
        integer :: i,j,n,grid_size
        n = size(b)
        grid_size = int(sqrt(real(n)))
        allocate(x(n),ax(n),p(n),r(n),z(n),aux(n))
        converged = .false.
        !$omp parallel
            !$omp do
            do j=1,n !x0 = 0; r0 = b; p0 = r0
                x(j) = 0.0d0
                r(j) = b(j)
            end do
            !$omp end do
            call M_inv(Ax_op,r,z,aux,params,grid_size)   !M^-1 r = z
            !$omp do 
            do j=1,n
                p(j) = z(j)
            end do
            !$omp end do
            do i=1,iter
                if (converged) cycle
                call Ax_op(p,ax,grid_size)    
                !$omp single
                    rr    = 0.0d0
                    alpha = 0.0d0
                    res   = 0.0d0
                    beta  = 0.0d0
                !$omp end single
                !$omp do reduction(+:rr,alpha) 
                    do j=1,n
                        rr    = rr + r(j) * z(j)
                        alpha = alpha + ax(j) * p(j)
                    end do 
                !$omp end do
                !$omp single
                    alpha = rr / alpha
                !$omp end single
                !$omp do reduction(+:res)
                    do j=1,n
                        x(j) = x(j) + alpha * p(j)
                        r(j) = r(j) - alpha * ax(j)
                        res  = res  + r(j) * r(j)
                    end do
                !$omp end do
                call M_inv(Ax_op,r,z,aux,params,grid_size) !precond M^-1 r = z
                !$omp do reduction(+:beta)
                    do j=1,n
                        beta = beta + r(j) * z(j)
                    end do
                !$omp end do
                !$omp single
                    res  = sqrt(res)
                    beta = beta / rr
                    if (res < tol) then
                        converged = .true.
                        iter = i
                    endif
                !$omp end single
                !$omp do
                    do j=1,n
                        p(j) = z(j) + beta * p(j)
                    end do
                !$omp end do
            end do
        !$omp end parallel 
    end subroutine pcg_omp
END MODULE conjugate_gradient
