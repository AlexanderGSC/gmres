MODULE Householder
    IMPLICIT none
    CONTAINS
    SUBROUTINE householder_sb(A,H)
        real(8),intent(inout)::A(:,:)
        real(8),allocatable,intent(out)::H(:,:)
        real(8),allocatable::v(:),x(:)
        integer::i,j,n
        integer::alpha,beta
        n= size(A,2)
        allocate(H(n,n),v(n),x(n))
        v=0.0d0;x=0.0d0;H = 0.0d0
        do j=1,n
            v(j:n)    = A(j:n,j)
            v(j)      = sign(norm2(v(j:n)),v(j)) + v(j)
            v(j:n)    = v(j:n) / norm2(v(j:n))
            H(j:n,j)  = v(j:n) !stores the reflection vector
            do i=j,n
                A(j:n,i) = A(j:n,i) - 2.0d0 * v(j:n) * dot_product(v(j:n),A(j:n,i))
            end do
        end do
    END SUBROUTINE
    SUBROUTINE HH_buildQ(H,Q)
        real(8),intent(in)::H(:,:)
        real(8),allocatable, intent(out)::Q(:,:)
        integer::nx,ny,i,j
        nx = size(H,1); ny = size(H,2)
        allocate(Q(nx,ny))
        Q = 0.0d0
        print *,"H="
        do i=1,nx
            print "(5F12.8)", H(i,1:5)
        end do
        do i=1,min(nx,ny)
            Q(i,i) = 1.0d0
        end do
        do j=1,nx !col n (n-1)...1
            do concurrent (i=1:ny)!row 1,2,...n
                Q(i,j:nx) = Q(i,j:nx) - 2.0d0*H(j:nx,j)*(dot_product(H(j:nx,j),Q(i,j:nx)))
            end do
        end do
    END SUBROUTINE
END MODULE

