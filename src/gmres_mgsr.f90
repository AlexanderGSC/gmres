MODULE GMRES_MGSR_MOD
    implicit none
    private 
    integer :: max_restarts = 100
    public :: gmres_mgsr_restarted
    public :: gmres_mgsr
    contains
        subroutine gmres_mgsr_restarted(A, b, x, m, tol,final_err,v_err,n_out,restart_out)
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
    end subroutine gmres_mgsr_restarted


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
END MODULE GMRES_MGSR_MOD
