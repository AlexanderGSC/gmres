MODULE GMRES_HH_MOD
    implicit none
    private
    public :: gmres_hh
CONTAINS
    subroutine gmres_hh(A, b, x, m, tol, final_err, v_err, n_out)
        real(8), intent(in) :: A(:,:)   !Original matrix
        real(8), intent(in) :: b(:)     !Initial vector
        real(8), allocatable, intent(out):: x(:) !Sol. vector
        integer, intent(in) :: m        !Max iterations
        real(8), intent(in) :: tol      !Max tolerance
        real(8), allocatable, intent(out):: final_err(:), v_err(:) !size of the error of the subspace
        integer, intent(out):: n_out    !Iterations done
        real(8), allocatable:: P(:,:), H(:,:)!Reflections and Hessemberg
        ! Other variables
        integer :: i, j, n ! indexes
        real(8), allocatable:: w(:), g(:), y(:), v_j(:) !aux vectors
        real(8) :: tmp, ds, h_val, beta
        real(8), allocatable:: cs(:), sn(:) !sin/cos for Givens rotations
        real(8), allocatable:: x2(:)
        n = size(A,1)
        allocate(P(n,m+1), H(m+1,m), final_err(m), v_err(m+1), v_j(n), w(n), g(m+1))
        allocate(cs(m),sn(m))
        P = 0.0d0;H=0.0d0;final_err=0.0d0;v_err=0.0d0;g=0.0d0
        beta = norm2(b)
        g(1) = -sign(beta,b(1))
        !print *,"g=",g
        w    = b
        w(1) = sign(beta,b(1)) + b(1)
        P(:,1) = w / norm2(w) !new reflector 
    
        do j=1,m
            n_out=j
            v_j=0.0d0;v_j(j)=1.0d0 !canonical vector v = ej (canonical))
            do i=j,1,-1
                v_j = v_j - 2.0d0*P(:,i)*dot_product(v_j,P(:,i))
            end do
            w = matmul(A,v_j)
            do i=1,j
                w = w - 2.0d0*P(:,i)*dot_product(w,P(:,i)) !apply reflections
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
                !print *,"NEW REFLECTOR: ", w/norm2(w)
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
            final_err(j) = abs(g(j+1)) / beta
            if (h_val < 1d-15 .or. final_err(j) < tol) then
                n_out = j
                exit
            end if
        end do

        allocate(y(n_out))
        y(n_out) = g(n_out) / H(n_out,n_out)
        
        do i=n_out-1,1,-1
            y(i) = (g(i) - dot_product(H(i,i+1:n_out),y(i+1:n_out))) / H(i,i)
        end do
        
        allocate(x(n),x2(n))
        x=0.0d0
        x(1:n_out) = y
        
        do i=n_out,1,-1
            x = x - 2.0d0*P(:,i)*dot_product(P(:,i),x)
        end do
        
        call calculate_verr(P,x2,y,v_err,n_out)
        x = x2
    end subroutine gmres_hh

    subroutine calculate_verr(P,x,y, v_err, n_iter)
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
        x = matmul(V,y)
    end subroutine calculate_verr
end module GMRES_HH_MOD