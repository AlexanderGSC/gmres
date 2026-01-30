MODULE GMRES_MF
    use omp_lib
!GMRES matrix free implementation
!This subroutine takes a stencil vector function to calculate the 
!Krylov basis vectors
    integer :: stages = 100
    public :: gmres_hh_mf
    public :: stencil_vector
    abstract interface 
        subroutine stencil_vector(x, y, n)
            real(8), intent(in) :: x(:)   
            real(8), intent(out):: y(:)  !y = A x
            integer, intent(in) :: n
        end subroutine stencil_vector
    end interface
CONTAINS

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
            g=0.0d0;P = 0.0d0;H=0.0d0
            call Ax_vec(x, w, nsize)
            !w    = b - w
            !$omp parallel do     
            do omp_idx=1,n
                w(omp_idx) = b(omp_idx) - w(omp_idx)
            end do
            !$omp end parallel do
            beta = norm2(w)
            g(1) = -sign(beta,w(1))
            w(1) = sign(beta,w(1)) + w(1)
            P(:,1) = w / norm2(w) !new reflector 
            do j=1,m
                n_out=j
                !v_j=0.0d0;v_j(j)=1.0d0 !canonical vector v = ej (canonical))
                !$omp parallel do     
                do omp_idx=1,n
                    v_j(omp_idx) = 0.0d0
                end do
                !$omp end parallel do
                v_j(j) = 1.0d0
                !do i=j,1,-1
                !    v_j = v_j - 2.0d0*P(:,i)*dot_product(v_j,P(:,i))
                !end do
                do i=j,1,-1     
                    omp_dotsum = 0.0d0    
                    !$omp parallel do reduction(+:omp_dotsum)
                    do omp_idx=1,n
                        omp_dotsum = omp_dotsum + v_j(omp_idx) * P(omp_idx,i)
                    end do
                    !$omp end parallel do
                    !$omp parallel do
                    do omp_idx=1,n
                        v_j(omp_idx) = v_j(omp_idx) - 2.0d0 * P(omp_idx,i) * omp_dotsum
                    end do
                    !$omp end parallel do 
                end do

                call Ax_vec(v_j,w,nsize)

                !do i=1,j
                !    w = w - 2.0d0*P(:,i)*dot_product(w,P(:,i)) !apply reflections
                !end do
                do i=1,j
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
                if (h_val < tol .or. final_err(j) < tol) then
                    n_out = j
                    stages_out = k
                    exit
                end if
            end do !j loop
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
            if (h_val < tol .or. final_err(n_out) < tol) then
                stages_out = k
                exit
            end if
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
        integer::i,j
        allocate(V(size(P,1),n_iter))
        V = 0.0d0;
        do concurrent (i=1:n_iter)
            V(i,i) = 1.0d0
        end do
        do concurrent (i=1:n_iter)
            do j=i,1,-1
                V(:,i) = V(:,i) - 2.0d0*P(:,j) * dot_product(V(:,i),P(:,j))
            end do
        end do
        
        do concurrent (i=2:n_iter)
            do j=1,i-1
                v_err(i) = v_err(i) + 2.0*(dot_product(V(:,i),V(:,j))**2)
            end do
        end do
        x = matmul(V,y(1:n_iter))
    end subroutine calculate_verr
END MODULE GMRES_MF