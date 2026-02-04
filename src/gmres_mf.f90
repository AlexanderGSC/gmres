MODULE GMRES_MF
    use interfaces
    use omp_lib
!GMRES matrix free implementation
!This subroutine takes a stencil vector function to calculate the 
!Krylov basis vectors
    integer :: stages = 1000
    public :: gmres_hh_mf
    public :: stencil_vector

CONTAINS

    subroutine gmres_mgsr_mf_mpo(Ax_vec, b, x, m, tol,final_err,v_err,n_out,restart_out,M_inv,params)
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
        integer :: i, j, k, n, st, idx, nsize
        real(8), allocatable:: w(:), g(:), y(:), z(:), aux(:)
        real(8) :: tmp, ds, h_val, h_tmp, beta, beta0
        real(8), allocatable:: cs(:), sn(:) !givens rotations
        !--------------------------------------------------------------------------
        n = size(b,1)               !n defines the size of the problem
        nsize = int(sqrt(real(n)))  !nsize defines the size of the grid
        !--------------------------------------------------------------------------
        !Allocating and initialization
        allocate(V(n,m+1),y(m),H(m+1,m), x(n), z(n),aux(n), final_err(m), v_err(m+1), w(n), g(m+1))
        allocate(cs(m),sn(m))
        V = 0.0d0;H=0.0d0;final_err=0.0d0;v_err=0.0d0;g=0.0d0;x=0.0d0
        !--------------------------------------------------------------------------
        ! beta0 = ||b - A x0||, usually we take x0 = 0
        beta0 = norm2(b)
        !the outer loop is the number of stages we compute
        do st=1,stages
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
    end subroutine gmres_mgsr_mf_mpo

    subroutine gmres_hh_mf(Ax_vec,b,x,m,tol,final_err,v_err,n_out,stages_out)
        procedure(stencil_vector) :: Ax_vec !Operator A b = x
        real(8), intent(in) :: b(:)     !Initial vector
        real(8), allocatable, intent(out):: x(:) !Sol. vector
        integer, intent(in) :: m   !Max number of stages
        real(8), intent(in) :: tol      !Max tolerance
        real(8), allocatable, intent(out):: final_err(:), v_err(:) !size of the error of the subspace
        integer, intent(out):: n_out    !Iterations done
        integer, intent(out):: stages_out !Number of stages done
        real(8), allocatable:: P(:,:), H(:,:)!Reflections and Hessemberg
        ! Other variables
        integer :: i, j, k, n ! indexes
        integer :: omp_idx
        real(8) :: omp_dotsum
        real(8), allocatable:: w(:), g(:), y(:), v_j(:)!aux vectors
        real(8) :: tmp, ds, h_val, beta, beta0
        real(8), allocatable:: cs(:), sn(:) !sin/cos for Givens rotations
        
        !Allocating all memory
        n = size(b,1)
        nsize = int(sqrt(real(n)))
        allocate(P(n,m+1), H(m+1,m), final_err(m), v_err(m+1), v_j(n), w(n), g(m+1))
        allocate(cs(m),sn(m))
        allocate(x(n),y(m))
        x=0.0d0;g=0.0d0;
        final_err=0.0d0;v_err=0.0d0;
        beta0 = norm2(b)
        do k=1,stages
            !$omp parallel default(shared), private(omp_idx)
            !$omp workshare 
            g=0.0d0;P = 0.0d0;H=0.0d0
            !$omp end workshare
            call Ax_vec(x, w, nsize)
            !$omp do     
            do omp_idx=1,n
                w(omp_idx) = b(omp_idx) - w(omp_idx)
            end do
            !$omp end do
            !$omp end parallel 
            beta = norm2(w)
            g(1) = -sign(beta,w(1))
            w(1) = sign(beta,w(1)) + w(1)
            P(:,1) = w / norm2(w) !new reflector
            !$omp parallel
            do j=1,m   
                !v_j=0.0d0;v_j(j)=1.0d0 !canonical vector v = ej (canonical))
                !$omp do     
                do omp_idx=1,n
                    v_j(omp_idx) = 0.0d0
                end do
                !$omp end do
                !$omp single
                n_out=j
                v_j(j) = 1.0d0
                !$omp end single
                !do i=j,1,-1
                !    v_j = v_j - 2.0d0*P(:,i)*dot_product(v_j,P(:,i))
                !end do
                do i=j,1,-1    
                    !$omp single
                    omp_dotsum = 0.0d0
                    !$omp end single 
                    !$omp do reduction(+:omp_dotsum)  
                    do omp_idx=1,n
                        omp_dotsum = omp_dotsum + v_j(omp_idx) * P(omp_idx,i)
                    end do
                    !$omp end do
                    !$omp do
                    do omp_idx=1,n
                        v_j(omp_idx) = v_j(omp_idx) - 2.0d0 * P(omp_idx,i) * omp_dotsum
                    end do
                    !$omp end do 
                end do
                
                call Ax_vec(v_j,w,nsize)
                
                !do i=1,j
                !    w = w - 2.0d0*P(:,i)*dot_product(w,P(:,i)) !apply reflections
                !end do
                do i=1,j
                    !$omp single
                    omp_dotsum = 0.0d0
                    !$omp end single
                    !$omp do reduction(+:omp_dotsum)  
                    do omp_idx=1,n
                        omp_dotsum = omp_dotsum + w(omp_idx) * P(omp_idx,i)
                    end do
                    !$omp end do
                    !$omp do
                    do omp_idx=1,n
                        w(omp_idx) = w(omp_idx) - 2.0d0*P(omp_idx, i) * omp_dotsum
                    end do
                    !$omp end do
                end do
                !$omp single
                H(1:j,j) = w(1:j)
                if (j < n) then
                    tmp    = norm2(w(j+1:n))
                    if (w(j+1) > 0.0d0) then
                        H(j+1,j) =-tmp
                    else
                        H(j+1,j) = tmp
                    end if
                    h_val    = abs(H(j+1,j))
                    w(1:j)   = 0.0d0
                    w(j+1)   = w(j+1) - H(j+1,j)
                    w        = w / norm2(w)
                    P(:,j+1) = w   !Stores the new reflector
                else
                    H(j+1,j) = 0.0d0
                end if
                !------ GIVENS ------------------------
                do i=1,j-1
                    tmp     = H(i,j)
                    H(i,j)  = cs(i)*tmp + sn(i)*H(i+1,j)
                    H(i+1,j)=-sn(i)*tmp + cs(i)*H(i+1,j)
                end do
                ds    = hypot(H(j+1,j),H(j,j))
                cs(j) = H(j,j) / ds
                sn(j) = H(j+1,j) / ds
                H(j,j)= cs(j)*H(j,j)+sn(j)*H(j+1,j)
                H(j+1,j) = 0.0d0
                !Now H is triangular
                tmp=g(j)
                !Apply the last rotation to g
                g(j)  = cs(j)*tmp + sn(j)*g(j+1)
                g(j+1)=-sn(j)*tmp + cs(j)*g(j+1)
                !------ END GIVENS -------------------------------------
                final_err(j) = abs(g(j+1)) / beta0
                !if (h_val < tol .or. final_err(j) < tol) then
                !    n_out = j
                !    stages_out = k
                !    exit
                !end if
                !$omp end single
            end do !j loop
            !$omp end parallel 
            !solving system g = Hy, where Hessenberg is triangular

            y = 0.0d0;
            y(n_out) = g(n_out) / H(n_out,n_out)  
            do i=n_out-1,1,-1
                y(i) = (g(i) - dot_product(H(i,i+1:n_out),y(i+1:n_out))) / H(i,i)
            end do
            !using w vector for calculating increment of x
            w = 0.0d0
            w(1:n_out) = y(1:n_out)
            !do i=n_out,1,-1
            !    w = w - 2.0d0*P(:,i)*dot_product(P(:,i),w)
            !end do
            do i=n_out,1,-1
                omp_dotsum = 0.0d0
                !$omp parallel do reduction(+:omp_dotsum)
                do omp_idx=1,n
                    omp_dotsum = omp_dotsum + w(omp_idx) * P(omp_idx,i)
                end do
                !$omp end parallel do
                !$omp parallel do
                do omp_idx=1,n
                    w(omp_idx) = w(omp_idx) - 2.0d0*P(omp_idx, i) * omp_dotsum
                end do
                !$omp end parallel do
            end do
            !$omp parallel do
            do omp_idx=1,n
                x(omp_idx) = x(omp_idx) + w(omp_idx)
            end do 
            !$omp end parallel do
            !x = x + w   !x updateS
            !print *, "STAGE=",k,"IT=",(k-1)*m+n_out, "ERROR=", final_err(n_out)
            stages_out = k
            if (final_err(n_out) < tol) exit
        end do  !restart loop
        call calculate_verr(P,w,y,v_err,n_out)
    end subroutine gmres_hh_mf
    
    subroutine calculate_verr(P, x, y, v_err, n_iter)
        real(8), intent(in):: P(:,:)    !reflectors
        real(8), intent(inout)::x(:)    ! solution vector
        real(8), intent(in)::y(:)       !vector y from Hessembreg
        real(8), intent(inout)::v_err(:)!ortogonality
        integer, intent(in)::n_iter     !number of iterations done
        real(8), allocatable::V(:,:)    !Ortogonal base
        real(8) :: dot
        integer::i,j,k,n
        n = size(P,1)
        allocate(V(n,n_iter))
        !$omp parallel
        !$omp workshare 
        V = 0.0d0
        x = 0.0d0
        !$omp end workshare
        !$omp do
        do i=1,n_iter
            V(i,i) = 1.0d0
        end do
        !$omp end do
        do i=1,n_iter

            do j=i,1,-1
                !$omp single
                dot = 0.0d0
                !$omp end single
                !$omp do reduction(+:dot)
                do k=1,n
                    dot = dot + V(k,i) * P(k,j)
                end do
                !$omp end do
                !$omp do
                do k=1,n
                    V(k,i) = V(k,i) - 2.0d0*P(k,j)*dot
                end do
                !$omp end do

                !V(:,i) = V(:,i) - 2.0d0*P(:,j) * dot_product(V(:,i),P(:,j))
            end do
        end do
        
        !$omp do private(j)
        do i=2,n_iter
            do j=1,i-1
                v_err(i) = v_err(i) + 2.0d0*(dot_product(V(:,i),V(:,j))**2)
            end do
        end do
        !$omp end do 
        !$omp do private(j)
        do i=1,n
            do j=1,n_iter
                x(j) = x(j) + V(i,j) * y(j)
            end do
        end do
        !$omp end do
        !x = matmul(V,y(1:n_iter))
        !$omp end parallel 
    end subroutine calculate_verr


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
        do st=1,stages
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
END MODULE GMRES_MF