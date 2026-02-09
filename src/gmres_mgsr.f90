MODULE GMRES_MGSR_MOD
    use interfaces
    use omp_lib
    implicit none
    private 
    integer :: max_restarts = 1000
    public :: gmres_mgsr_dense
    public :: gmres_mgsr_mf
    public :: gmres_mgsr_omp
    contains
    subroutine gmres_mgsr_dense(A, b, x, m, tol,final_err,v_err,n_out,restart_out)
        real(8), intent(in) :: A(:,:)   !Original matrix
        real(8), intent(in) :: b(:)     !Initial vector
        real(8), allocatable, intent(out):: x(:) !Sol. vector
        integer, intent(in) :: m        !Max iterations
        real(8), intent(in) :: tol      !Max tolerance
        real(8), allocatable, intent(out):: final_err(:), v_err(:) !size of the error of the subspace
        integer, intent(out):: n_out    !Iterations done
        integer, intent(out):: restart_out !Num 
        real(8), allocatable:: V(:,:), H(:,:)!V and Hessemberg
        ! Other variables
        integer :: i, j, k, n, st
        real(8), allocatable:: w(:), g(:), y(:)
        real(8) :: tmp, ds, h_val, h_tmp, beta, beta0
        real(8), allocatable:: cs(:), sn(:)
        n = size(A,1)
        allocate(V(n,m+1),y(m),H(m+1,m), x(n), final_err(m), v_err(m+1), w(n), g(m+1))
        allocate(cs(m),sn(m))
        V = 0.0d0;H=0.0d0;final_err=0.0d0;v_err=0.0d0;g=0.0d0;x=0.0d0
        beta0 = norm2(b)
        do st=1,max_restarts
            g = 0.0d0; H = 0.0d0; V = 0.0d0
            w = b - matmul(A,x)
            beta = norm2(w)
            V(:,1) = w / beta
            g(1)   = beta
            do j=1,m
                n_out = j
                w = matmul(A,V(:,j))
                ! ----------- Modified Gram Schmidt (MGSR) -------------
                ! the sequential process is required in order to converge
                ! better orthogonalization.
                do k=1,2 ! Twice is enough for reorthogonalization
                    do i=1,j
                        h_tmp = dot_product(w, V(:,i))
                        H(i,j) = H(i,j) + h_tmp
                        w = w - h_tmp*V(:,i)
                    end do
                end do
                h_val = norm2(w)
                H(j+1,j) = h_val
                !------ GIVENS ------------------------
                do i=1,j-1
                    tmp     = H(i,j)
                    H(i,j)  = cs(i)*tmp + sn(i)*H(i+1,j)
                    H(i+1,j)=-sn(i)*tmp + cs(i)*H(i+1,j)
                end do
                ds    = hypot(H(j+1,j),H(j,j))
                cs(j) = H(j,j) / ds
                sn(j) = H(j+1,j) / ds
                H(j,j)= cs(j)*H(j,j) + sn(j)*H(j+1,j)
                H(j+1,j) = 0.0d0
                tmp=g(j)
                !Apply the last rotation to g
                g(j)  = cs(j)*tmp + sn(j)*g(j+1)
                g(j+1)=-sn(j)*tmp + cs(j)*g(j+1)
                !----------END GIVENS-------------------------
                final_err(j) = abs(g(j+1)) / beta0
                if (h_val < tol .or. final_err(j) < tol) then
                    n_out = j
                    exit
                end if
                V(:,j+1) = w / h_val
            end do !Arnoldi j end loop 
            y = 0.0d0
            y(n_out) = g(n_out) / H(n_out,n_out)
            do i=n_out-1,1,-1
                y(i) = (g(i) - dot_product(H(i,i+1:n_out),y(i+1:n_out))) / H(i,i)
            end do
        
            x = x + matmul(V(:,1:n_out),y(1:n_out))
            !print *, "STAGE=",st,"IT=",(st-1)*m+n_out, "ERROR=", final_err(n_out)
            if (h_val < tol .or. final_err(n_out) < tol) then
                restart_out = st
                exit
            end if
        end do !restart k loop
        do j=1,n_out !Sizing the orthogonality 
            do i=1,j
                v_err(j+1) = v_err(j+1) + 2.0*(dot_product(V(:,i),V(:,j+1))**2)
            end do
            v_err(j+1) = v_err(j+1) + (dot_product(V(:,j+1),V(:,j+1))-1.0d0)**2
            v_err(j+1) = sqrt(v_err(j)**2 + v_err(j+1))
        end do
    end subroutine gmres_mgsr_dense


    subroutine gmres_mgsr_mf(Ax_vec, b, x, m, tol,final_err,v_err,n_out,restart_out,M_inv,params)
        procedure(stencil_vector) :: Ax_vec  !Operator A b = x
        real(8), intent(in) :: b(:)          !Initial vector
        real(8), allocatable, intent(out):: x(:) !Sol. vector
        integer, intent(in) :: m             !Max iterations
        real(8), intent(in) :: tol           !Max tolerance
        real(8), allocatable, intent(out):: final_err(:) !residual of the solution 
        real(8), allocatable, intent(out):: v_err(:) !orthogonality error
        integer, intent(out):: n_out         !Iterations done
        integer, intent(out):: restart_out   !Num of restarts done
        procedure(precond) :: M_inv          !Preconditioner
        real(8), intent(in) :: params(:)     !Preconditioner params 
        real(8), allocatable:: V(:,:), H(:,:)!V and Hessemberg matrixes
        ! Other variables
        integer :: i, j, k, n, st, nsize
        real(8), allocatable:: w(:), g(:), y(:), z(:), aux(:)
        real(8) :: tmp, ds, h_val, h_tmp, beta, beta0
        real(8), allocatable:: cs(:), sn(:) !givens rotations
        n = size(b,1)               !n defines the size of the problem
        nsize = int(sqrt(real(n)))  !nsize defines the size of the grid
        !--------------------------------------------------------------------------
        !Allocating and initialization
        allocate(V(n,m+1),y(m),H(m+1,m), x(n), z(n), aux(n), final_err(m), v_err(m+1), w(n), g(m+1))
        allocate(cs(m),sn(m))
        V = 0.0d0;H=0.0d0;final_err=0.0d0;v_err=0.0d0;g=0.0d0;x=0.0d0
        !--------------------------------------------------------------------------
        ! beta0 = ||b - A x0||, usually we take x0 = 0
        beta0 = norm2(b)
        !the outer loop is the number of stages we compute
        do st=1,max_restarts
            g = 0.0d0; H = 0.0d0; V = 0.0d0
            call Ax_vec(x,w,nsize) !w = Ax
            z = b - w              !z = b - Ax
            call M_inv(Ax_vec,z,w,aux,params,nsize)
            beta = norm2(w)     !beta = ||w||
            V(:,1) = w / beta   !stores the first vector of the orthogonal basis
            g(1)   = beta       !g is rhs of Hessemberg's system
            !the inner loop: Arnoli's Iteration
            do j=1,m
                n_out = j
                call Ax_vec(V(:,j), z, nsize)
                call M_inv(Ax_vec,z,w,aux,params,nsize)
                ! ----------- Modified Gram Schmidt (MGSR) -------------
                ! the sequential process is required in order to converge
                ! better orthogonalization.
                do k=1,2 ! Twice is enough for reorthogonalization
                    do i=1,j
                        h_tmp = dot_product(w, V(:,i))
                        H(i,j) = H(i,j) + h_tmp
                        w = w - h_tmp*V(:,i)
                    end do
                end do
                h_val = norm2(w)
                H(j+1,j) = h_val
                !------ GIVENS ROTATIONS ------------------------
                do i=1,j-1          !performs all previoous rotations 
                    tmp     = H(i,j)
                    H(i,j)  = cs(i)*tmp + sn(i)*H(i+1,j)
                    H(i+1,j)=-sn(i)*tmp + cs(i)*H(i+1,j)
                end do
                ds    = hypot(H(j+1,j),H(j,j))
                cs(j) = H(j,j) / ds !calculate the new rotation
                sn(j) = H(j+1,j) / ds
                H(j,j)= cs(j)*H(j,j) + sn(j)*H(j+1,j)
                H(j+1,j) = 0.0d0
                tmp=g(j)
                !Apply the last rotation to g
                !In order to solver the system H y = g, we need to perform all
                !rotations to rhs g
                g(j)  = cs(j)*tmp + sn(j)*g(j+1)
                g(j+1)=-sn(j)*tmp + cs(j)*g(j+1)
                !----------END GIVENS-------------------------
                !residual error is computed as g(j+1)
                final_err(j) = abs(g(j+1)) / beta0
                if (h_val < tol .or. final_err(j) < tol) then
                    n_out = j
                    exit
                end if
                V(:,j+1) = w / h_val
            end do !Arnoldi j end loop
            !Solving the triangular system H y = g, using backward sustitution
            y = 0.0d0
            y(n_out) = g(n_out) / H(n_out,n_out)
            do i=n_out-1,1,-1
                y(i) = (g(i) - dot_product(H(i,i+1:n_out),y(i+1:n_out))) / H(i,i)
            end do
            !update the solution vector 
            x = x + matmul(V(:,1:n_out),y(1:n_out))
            !print *, "STAGE=",st,"IT=",(st-1)*m+n_out, "ERROR=", final_err(n_out)
            if (h_val < tol .or. final_err(n_out) < tol) then
                restart_out = st
                exit
            end if
        end do !restart k loop
        do j=1,n_out !Sizing the orthogonality 
            do i=1,j
                v_err(j+1) = v_err(j+1) + 2.0*(dot_product(V(:,i),V(:,j+1))**2)
            end do
            v_err(j+1) = v_err(j+1) + (dot_product(V(:,j+1),V(:,j+1))-1.0d0)**2
            v_err(j+1) = sqrt(v_err(j)**2 + v_err(j+1))
        end do
    end subroutine gmres_mgsr_mf


    subroutine gmres_mgsr(A, b, x, m, tol, final_err, v_err, n_out)
        real(8), intent(in) :: A(:,:)   !Original matrix
        real(8), intent(in) :: b(:)     !Initial vector
        real(8), allocatable, intent(out):: x(:) !Sol. vector
        integer, intent(in) :: m        !Max iterations
        real(8), intent(in) :: tol      !Max tolerance
        real(8), allocatable, intent(out):: final_err(:), v_err(:) !size of the error of the subspace
        integer, intent(out):: n_out    !Iterations done
        real(8), allocatable:: V(:,:), H(:,:)!V and Hessemberg
        ! Other variables
        integer :: i, j, k, n
        real(8), allocatable:: w(:), g(:), y(:)
        real(8) :: tmp, ds, h_val, h_tmp, beta
        real(8), allocatable:: cs(:), sn(:)
        n = size(A,1)
        allocate(V(n,m+1), H(m+1,m), final_err(m), v_err(m+1), w(n), g(m+1))
        allocate(cs(m),sn(m))
        V = 0.0d0;H=0.0d0;final_err=0.0d0;v_err=0.0d0;g=0.0d0
        beta = norm2(b)
        g(1) = beta
        V(:,1) = b / beta

        do j=1,m
            n_out = j
            w = matmul(A,V(:,j))
            ! ----------- Modified Gram Schmidt (MGSR) -------------
            ! the sequential process is required in order to converge
            ! better orthogonalization.
            do k=1,2 ! Twice is enough for reorthogonalization
                do i=1,j
                    h_tmp = dot_product(w, V(:,i))
                    H(i,j) = H(i,j) + h_tmp
                    w = w - h_tmp*V(:,i)
                end do
            end do
            h_val = norm2(w)
            H(j+1,j) = h_val
            !------ GIVENS ------------------------
            do i=1,j-1
                tmp     = H(i,j)
                H(i,j)  = cs(i)*tmp + sn(i)*H(i+1,j)
                H(i+1,j)=-sn(i)*tmp + cs(i)*H(i+1,j)
            end do
            ds    = hypot(H(j+1,j),H(j,j))
            cs(j) = H(j,j) / ds
            sn(j) = H(j+1,j) / ds
            H(j,j)= cs(j)*H(j,j) + sn(j)*H(j+1,j)
            H(j+1,j) = 0.0d0
            tmp=g(j)
            !Apply the last rotation to g
            g(j)  = cs(j)*tmp + sn(j)*g(j+1)
            g(j+1)=-sn(j)*tmp + cs(j)*g(j+1)
            !----------END GIVENS-------------------------
            final_err(j) = abs(g(j+1)) / beta
            if (h_val < tol .or. final_err(j) < tol) then
                n_out = j
                exit
            else
                V(:,j+1) = w / h_val
                do i=1,j
                    v_err(j+1) = v_err(j+1) + 2.0*(dot_product(V(:,i),V(:,j+1))**2)
                end do
                v_err(j+1) = v_err(j+1) + (dot_product(V(:,j+1),V(:,j+1))-1.0d0)**2
                v_err(j+1) = sqrt(v_err(j)**2 + v_err(j+1))
            end if
        end do
        allocate(y(n_out))
        y(n_out) = g(n_out) / H(n_out,n_out)
        do i=n_out-1,1,-1
            y(i) = (g(i) - dot_product(H(i,i+1:n_out),y(i+1:n_out))) / H(i,i)
        end do
        allocate(x(n))
        x = matmul(V(:,1:n_out),y)
    end subroutine gmres_mgsr

    subroutine gmres_mgsr_omp(Ax_vec, b, x, m, tol,final_err,v_err,n_out,restart_out,M_inv,params)
        procedure(stencil_vector) :: Ax_vec !Operator A x 
        real(8), intent(in) :: b(:)         !Initial vector
        real(8), allocatable, intent(out):: x(:) !Sol. vector
        integer, intent(in) :: m            !Max iterations
        real(8), intent(in) :: tol          !Max tolerance
        real(8), allocatable, intent(out):: final_err(:) !residual of the solution 
        real(8), allocatable, intent(out):: v_err(:) !orthogonality error
        integer, intent(out):: n_out        !Iterations done
        integer, intent(out):: restart_out  !Num of restarts done
        procedure(precond) :: M_inv          !Preconditioner
        real(8), intent(in) :: params(:)     !Preconditioner params 
        real(8), allocatable:: V(:,:), H(:,:)!V and Hessemberg matrixes
        ! Other variables
        logical :: converged
        integer :: i, j, k, n, st, idx, nsize
        real(8), allocatable:: w(:), g(:), y(:), z(:), aux(:)
        real(8) :: tmp, ds, h_val, h_tmp, beta, beta0
        real(8), allocatable:: cs(:), sn(:) !givens rotations
        !--------------------------------------------------------------------------
        n = size(b,1)               !n defines the size of the problem
        nsize = int(sqrt(real(n)))  !nsize defines the size of the grid
        converged = .false.
        !--------------------------------------------------------------------------
        !Allocating and initialization
        allocate(V(n,m+1),y(m),H(m+1,m), x(n), z(n),aux(n), final_err(m), v_err(m+1), w(n), g(m+1))
        allocate(cs(m),sn(m))
        V = 0.0d0;H=0.0d0;final_err=0.0d0;v_err=0.0d0;g=0.0d0;x=0.0d0
        !--------------------------------------------------------------------------
        ! beta0 = ||b - A x0||, usually we take x0 = 0
        beta0 = norm2(b)
        !the outer loop is the number of stages we compute
        do st=1,max_restarts
            !$omp parallel
                !$omp workshare 
                g = 0.0d0; H = 0.0d0; V = 0.0d0
                !$omp end workshare
                call Ax_vec(x,w,nsize)  !w = Ax
                !$omp do
                    do idx=1,n            !z = b - w  
                        z(idx) = b(idx) - w(idx)
                    end do
                !$omp end do 
                call M_inv(Ax_vec,z,w,aux,params,nsize) !precond M^-1 z = w
                !$omp single
                beta = norm2(w)     !beta = ||w||
                g(1)   = beta       !g is rhs of Hessemberg's system
                !$omp end single
                !$omp do
                    do idx=1,n
                        V(idx,1) = w(idx) / beta
                    end do
                !$omp end do
                !V(:,1) = w / beta   !stores the first vector of the orthogonal basis
            !$omp end parallel 
            !the inner loop: Arnoli's Iteration
            !$omp parallel
            do j=1,m
                if (converged) cycle
                call Ax_vec(V(:,j), z, nsize) !z = A * V(:,j)
                call M_inv(Ax_vec,z,w,aux,params,nsize) !precond w = M^-1 z
                ! ----------- Modified Gram Schmidt (MGSR) -------------
                ! the sequential process is required in order to converge
                ! better orthogonalization.
                do k=1,2 ! Twice is enough for reorthogonalization
                    do i=1,j
                        !$omp single
                        h_tmp = 0.0d0
                        !$omp end single
                        !$omp do reduction(+:h_tmp)
                        do idx=1,n !h_tmp = dot_product(w, V(:,i))
                            h_tmp = h_tmp + w(idx)*V(idx,i)
                        end do 
                        !$omp end do
                        !$omp master
                        H(i,j) = H(i,j) + h_tmp
                        !$omp end master
                        !$omp do
                        do idx=1,n !w = w - h_tmp*V(:,i)
                            w(idx) = w(idx) - h_tmp*V(idx,i)
                        end do
                        !$omp end do
                    end do
                end do
                !$omp single
                h_val = norm2(w)
                H(j+1,j) = h_val
                !------ GIVENS ROTATIONS ------------------------
                do i=1,j-1          !performs all previous rotations 
                    tmp     = H(i,j)
                    H(i,j)  = cs(i)*tmp + sn(i)*H(i+1,j)
                    H(i+1,j)=-sn(i)*tmp + cs(i)*H(i+1,j)
                end do
                ds    = hypot(H(j+1,j),H(j,j))
                cs(j) = H(j,j) / ds !calculate the new rotation
                sn(j) = H(j+1,j) / ds
                H(j,j)= cs(j)*H(j,j) + sn(j)*H(j+1,j)
                H(j+1,j) = 0.0d0             
                !Apply the last rotation to g
                !In order to solver the system H y = g, we need to perform all
                !rotations to rhs g
                tmp=g(j)
                g(j)  = cs(j)*tmp + sn(j)*g(j+1)
                g(j+1)=-sn(j)*tmp + cs(j)*g(j+1)
                !----------END GIVENS-------------------------
                !residual error is computed as g(j+1)
                final_err(j) = abs(g(j+1)) / beta0
                V(:,j+1) = w / h_val
                if (final_err(j) < tol) then
                    restart_out = st
                    converged = .true.
                end if
                n_out = j
                !$omp end single
            end do !Arnoldi j end loop
            !$omp end parallel
            !Solving the triangular system H y = g, using backward sustitution
            y = 0.0d0
            y(n_out) = g(n_out) / H(n_out,n_out)
            do i=n_out-1,1,-1
                y(i) = (g(i) - dot_product(H(i,i+1:n_out),y(i+1:n_out))) / H(i,i)
            end do
            !update the solution vector
            !$omp parallel
            !$omp do 
                do idx=1,n
                    x(idx) = x(idx)+dot_product(V(idx,1:n_out),y(1:n_out))
                end do
            !$omp end do 
            !$omp end parallel 
            !x = x + matmul(V(:,1:n_out),y(1:n_out))
            !print *, "STAGE=",st,"IT=",(st-1)*m+n_out, "ERROR=", final_err(n_out)
            if (h_val < tol .or. final_err(n_out) < tol) then
                restart_out = st
                exit
            end if
        end do !restart k loop
        do j=1,n_out !Sizing the orthogonality 
            do i=1,j
                v_err(j+1) = v_err(j+1) + 2.0*(dot_product(V(:,i),V(:,j+1))**2)
            end do
            v_err(j+1) = v_err(j+1) + (dot_product(V(:,j+1),V(:,j+1))-1.0d0)**2
            v_err(j+1) = sqrt(v_err(j)**2 + v_err(j+1))
        end do
    end subroutine gmres_mgsr_omp
END MODULE GMRES_MGSR_MOD
