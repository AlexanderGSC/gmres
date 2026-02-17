!The biconjugated gradient stabilized algorith implemented from scratch
module bicgstab_mod
    use interfaces
    use omp_lib
    implicit none
    private
    public :: bicgstab
    public :: pbicgstab
    public :: pbicgstab_omp

CONTAINS 
    subroutine bicgstab(ax_op,b,x,tol,iter,res) 
        procedure(stencil_vector) :: ax_op          !A matrix defined as operator
        real(8), intent(in) :: b(:)                 !rhs
        real(8), allocatable, intent(out) :: x(:)   !solution vector
        real(8), intent(in) :: tol                  !tolerance
        integer, intent(inout) :: iter              !number of iterations done
        real(8), intent(out):: res                  !residual reached
        !other vars
        real(8) :: rr0  !dot product  (r,r0)
        real(8) :: alpha, beta, omega
        real(8), allocatable :: ap(:), r(:), r0(:),s(:),as(:),p(:)
        integer :: i,n,grid_size
        n = size(b)
        grid_size = int(sqrt(real(n)))
        allocate(x(n),r(n),r0(n),ap(n),s(n),as(n),p(n))
        !compute r = b - Ax, r0 arbitrary r0 = r
        !we start with x = 0, so r = b
        x = 0.0d0; r=b; r0 = r; p = r0
        do i=1,iter
            call ax_op(p,ap,grid_size)   !ap = A*p 
            rr0   = dot_product(r,r0)    !rr0 = (r,r0)
            alpha = rr0 / dot_product(ap,r0) 
            s = r - alpha*ap
            call ax_op(s,as,grid_size)
            omega = dot_product(as,s) / dot_product(as,as)
            x = x + alpha*p + omega*s
            r = s - omega*as
            res = norm2(r)
            if (res < tol) then
                iter = i
                exit
            end if
            beta = (dot_product(r,r0) / rr0) * (alpha / omega)
            p = r + beta*(p - omega * ap)
        end do
    end subroutine bicgstab

    subroutine pbicgstab(ax_op,b,x,tol,iter,res,m_inv,params) 
        procedure(stencil_vector) :: ax_op          !A matrix defined as operator
        real(8), intent(in) :: b(:)                 !rhs
        real(8), allocatable, intent(out) :: x(:)   !solution vector
        real(8), intent(in) :: tol                  !tolerance
        integer, intent(inout) :: iter              !number of iterations done
        real(8), intent(out):: res                  !residual reached
        procedure(precond) :: m_inv                 !preconditioner
        real(8), intent(in):: params(:)             !Precond params
        !other vars
        real(8) :: rr0  !dot product  (r,r0)
        real(8) :: alpha, beta, omega
        real(8), allocatable :: ap(:), r(:), r0(:),s(:),as(:),p(:)
        real(8), allocatable :: z1(:),z2(:),aux(:)
        integer :: i,n,grid_size
        n = size(b)
        grid_size = int(sqrt(real(n)))
        allocate(x(n),r(n),r0(n),ap(n),s(n),as(n),p(n),z1(n),z2(n),aux(n))
        !compute r = b - Ax, r0 arbitrary r0 = r
        !we start with x = 0, so r = b
        x = 0.0d0; r=b; r0 = r; p = r0
        do i=1,iter
            call m_inv(ax_op, p, z1, aux, params, grid_size)  !z = M^-1 r0
            call ax_op(z1,ap,grid_size)   !ap = A*p 
            rr0   = dot_product(r,r0)     !rr0 = (r,z=M^-1r0)
            alpha = rr0 / dot_product(ap,r0) 
            s = r - alpha*ap
            call m_inv(ax_op, s, z2, aux, params, grid_size)  !z = M^-1 r0
            call ax_op(z2,as,grid_size)
            omega = dot_product(as,s) / dot_product(as,as)
            x = x + alpha*z1 + omega*z2
            r = s - omega*as
            res = norm2(r)
            if (res < tol) then
                iter = i
                exit
            end if
            beta = (dot_product(r,r0) / rr0) * (alpha / omega)
            p = r + beta*(p - omega * ap)
        end do
    end subroutine pbicgstab

    subroutine pbicgstab_omp(ax_op,b,x,tol,max_iter,res,m_inv,params) 
        procedure(stencil_vector) :: ax_op          !A matrix defined as operator
        real(8), intent(in) :: b(:)                 !rhs
        real(8), allocatable, intent(out) :: x(:)   !solution vector
        real(8), intent(in) :: tol                  !tolerance
        integer, intent(inout) :: max_iter          !number of iterations done
        real(8), intent(out):: res                  !residual reached
        procedure(precond) :: m_inv                 !preconditioner
        real(8), intent(in):: params(:)             !Precond params
        !other vars
        logical :: converged
        real(8) :: rr0, ap_r0, as_s, as_as,r_r0_new !dot products accumulators
        real(8) :: alpha, beta, omega
        real(8), allocatable :: ap(:), r(:), r0(:),s(:),as(:),p(:)
        real(8), allocatable :: z1(:),z2(:),aux(:)
        integer :: i,j,n,grid_size,iters
        converged = .false.
        n = size(b)
        grid_size = int(sqrt(real(n)))
        allocate(x(n),r(n),r0(n),ap(n),s(n),as(n),p(n),z1(n),z2(n),aux(n))
        !compute r = b - Ax, r0 arbitrary r0 = r
        !we start with x = 0, so r = b
        !$omp parallel
        !$omp do
            do j=1,n
                x(j)=0.0d0; r(j)=b(j); r0(j)=r(j);p(j)=r0(j)
            end do
        !$omp end do
        do i=1,max_iter
            if (converged) cycle
            call m_inv(ax_op, p, z1, aux, params, grid_size)  !z = M^-1 r0
            call ax_op(z1,ap,grid_size)   !ap = A*p 
            !$omp do reduction(+:rr0,ap_r0)
            do j=1,n
                rr0   = rr0 + r(j)*r0(j)
                ap_r0 = ap_r0 + ap(j)*r0(j)
            end do
            !$omp end do
            !$omp single
                alpha = rr0 / ap_r0
            !$omp end single
            !$omp do
                do j=1,n
                    s(j) = r(j) - alpha*ap(j)
                end do
            !$omp end do
            call m_inv(ax_op, s, z2, aux, params, grid_size)  !z = M^-1 r0
            call ax_op(z2,as,grid_size)
            !$omp do reduction(+:as_s,as_as)
                do j=1,n
                    as_s = as_s + as(j)*s(j)
                    as_as= as_as+ as(j)*as(j)
                end do
            !$omp end do
            !$omp single 
                omega = as_s / as_as
            !$omp end single
            !$omp do
                do j=1,n
                    x(j) = x(j) + alpha*z1(j)+omega*z2(j)
                    r(j) = s(j) - omega*as(j)
                end do
            !$omp end do
            !$omp single
                res = norm2(r)
                if (res < tol) then
                    iters = i
                    converged = .true.
                end if
            !$omp end single
            !$omp do reduction(+:r_r0_new)
            do j=1,n
                r_r0_new = r_r0_new + r(j)*r0(j)
            end do
            !$omp end do
            !$omp single
                beta = (r_r0_new / rr0) * (alpha / omega)
                r_r0_new = 0.0d0
                as_s = 0.0d0
                as_as= 0.0d0
                rr0  = 0.0d0
                ap_r0= 0.0d0
            !$omp end single
            !$omp do
                do j=1,n
                    p(j) = r(j) + beta*(p(j) - omega*ap(j))
                end do
            !$omp end do 
        end do
        !$omp end parallel
        max_iter = iters !set max_iter as output parameter with the number of iterations
    end subroutine pbicgstab_omp
    
end module bicgstab_mod
