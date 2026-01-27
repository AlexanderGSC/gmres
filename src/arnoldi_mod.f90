MODULE Arnoldi_Mod
    use Givens_Mod
    implicit none
    contains
    subroutine arnoldi_mgs(A, v1, m, tol, V, H, final_err, n_out)
        real(8), intent(in) :: A(:,:)   !Original matrix
        real(8), intent(in) :: v1(:)    !Initial vector
        integer, intent(in) :: m        !Max iterations
        real(8), intent(in) :: tol      !Max tolerance
        real(8), allocatable, intent(out) :: V(:,:), H(:,:) !V and Hessemberg
        real(8), intent(out):: final_err!size of the error of the subspace
        integer, intent(out):: n_out    !Iterations done
        ! Other variables
        integer :: i, j, k, n
        real(8), allocatable:: w(:)
        real(8) :: h_val, h_tmp
        real(8), allocatable:: b(:)

        n = size(A,1)
        allocate(V(n,m+1), H(m+1,m),b(m+1))
        V = 0.0d0; H = 0.0d0; b= 1.0d0
        allocate(w(n))
        V(:,1) = v1 / norm2(v1)
        do j=1,m
            n_out = j
            w = matmul(A,V(:,j))
            ! ----------- Modified Gram Schmidt -------------
            ! the sequential process is required in order to converge
            ! better orthogonalization
            do k=1,2 ! Twice is enough
                do i=1,j
                    h_tmp = dot_product(w, V(:,i))
                    H(i,j) = H(i,j) + h_tmp
                    w = w - h_tmp*V(:,i)
                end do
            end do
            h_val = norm2(w)
            final_err = h_val
            H(j+1,j) = h_val
            call givens(H,b,j,j)
            if (h_val < tol) then
                H = H(1:j+1,1:j)
                V = V(:,1:j)
                return
            end if
            V(:,j+1) = w / h_val
        end do
    end subroutine arnoldi_mgs

    subroutine arnoldi_hh(A, b, m, tol, P, H, final_err, n_out)
        real(8), intent(in):: A(:,:)   !Original matrix
        real(8), intent(in):: b(:)   !Initial vector
        integer, intent(in):: m        !Max iterations
        real(8), intent(in):: tol      !Max tolerance
        real(8), allocatable, intent(out)::P(:,:), H(:,:) !Reflections and Hessemberg
        real(8), intent(out):: final_err!size of the error of the subspace
        integer, intent(out):: n_out    !Iterations done
        real(8), allocatable::cs(:),sn(:)
        ! Other variables
        integer :: i, j, n
        real(8), allocatable:: w(:), v_j(:)
        real(8) ::n_tmp, ds, tmp

        n = size(A,1)
        allocate(P(n,m+1), H(m+1,m),w(n),v_j(n))
        allocate(cs(m),sn(m))
        P=0.0d0;H=0.0d0;w=0.0d0;cs=0.0d0;sn=0.0d0
        w = b
        w(1) = sign(norm2(b),b(1)) + b(1)
        print *,"NEW REFLECTOR: ", w/norm2(w)
        P(:,1) = w / norm2(w)
        do j=1,m
            print *,"IT=",j
            n_out=j
            v_j=0.0d0;v_j(j)=1.0d0 !canonical vector v = ej (canonical))
            do i=j,1,-1
                v_j = v_j - 2.0d0*P(:,i)*dot_product(v_j,P(:,i))
            end do
            print *,"vj=",v_j
            w = matmul(A,v_j)
            do i=1,j
                w = w - 2.0d0*P(:,i)*dot_product(w,P(:,i)) !apply reflections
            end do
            print *,"w=",w
            H(1:j,j) = w(1:j)
            if (j < n) then
                n_tmp    = sign(norm2(w(j+1:n)),w(j+1))
                H(j+1,j) = n_tmp
                final_err = abs(H(j+1,j))
                w(1:j)   = 0.0d0
                w(j+1)   = n_tmp + w(j+1)
                w        = w / norm2(w)
                print *,"NEW REFLECTOR: ", w/norm2(w)
                P(:,j+1) = w   !Stores the new reflector
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
            else
                H(j+1,j) = 0.0d0
                final_err= 0.0d0
            end if
            ! ----------- Householder Reflections -------------
            ! this variation of the algorithm uses the reflections stored
            ! on P for making zeros on the Hessenberg matrix
            ! step j -> makes zeros from j+1 to m
            if (final_err < tol) then
                H = H(1:j+1,1:j)
                P = P(:,1:j+1)
                exit
            end if
        end do
        H = H(1:n_out+1,1:n_out)
        P = P(:,1:n_out)
    end subroutine arnoldi_hh
END MODULE Arnoldi_Mod

