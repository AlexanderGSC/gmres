PROGRAM main
    use GMRES_Mod
    use Householder
    use Arnoldi_Mod
    IMPLICIt none
    integer::nsize,max_iter
    real(8)::diag, ndiag

    nsize    = 140
    max_iter = 2000
    diag     = 4.0
    ndiag    =-1.0
    print *,"---------------------------------------------------------------------------------------------"
    !call test_Poisson_2DHH(nsize,max_iter,diag,ndiag)
    print *,"---------------------------------------------------------------------------------------------"
    print *,"---------------------------------------------------------------------------------------------"
    call test_Poisson_2DMGSR(nsize,max_iter,diag,ndiag)
    print *,"---------------------------------------------------------------------------------------------"
        call test_Poisson_2DHH(nsize,max_iter,diag,ndiag)
    print *,"---------------------------------------------------------------------------------------------"
    CONTAINS
    SUBROUTINE test_Poisson_2DHH(nsize, max_iter, diag, ndiag)
        IMPLICIT none
        integer, intent(in)::nsize
        integer, intent(in)::max_iter
        real(8), intent(in)::diag, ndiag  !diag->main diag value / ndiag->neighbours
        real(8), allocatable::A(:,:),b(:),x(:),errn(:),verr(:)
        real(8)::tol
        integer::row,i,j,N, n_iter
        tol = 1.d-15
        N = nsize * nsize
        print *, "GMRES POISSON HH 2D TEST"
        print *,"PARAMS ","N=",N," MAX ITER=",max_iter," TOL=",tol," DIAG=",diag," NEIGHBOURS=",ndiag

        allocate(A(N,N),b(N),x(N))
        allocate(errn(max_iter))
        do i=1,nsize
            do j=1,nsize
                row = i+(j-1)*nsize
                A(row,row) =  diag
                if (row > 1) A(row,row-1)   = ndiag
                if (row < N) A(row,row+1)   = ndiag
                if (j > 1) A(row-nsize,row) = ndiag
                if (j < N) A(row+nsize,row) = ndiag
            end do
        end do
        x = 1.0d0
        b = matmul(A,x)


        print *, "Solving Poisson 2D (N =", N, ")..."
        call gmres_hh(A,b,x,max_iter,tol,errn,verr,n_iter)
        print *, "Num iterations util convergence:", n_iter, "(MAX=",max_iter,")"
        print *, "Final residual:", errn(n_iter)
        !print *, "Final ||I - V.t * V||:", verr(n_iter)
        print *, "Max error:", maxval(abs(x - 1.0d0))
        print *, "Showing first 10 elements of the solution (must be 1.0)..."
        print "(10F14.10)", x(1:10)
        print *, "L2 Sol. Norm: ", norm2(x-1.0d0)
        deallocate(A)
        deallocate(b,x,errn,verr)
    END SUBROUTINE
    SUBROUTINE test_Poisson_2DMGSR(nsize, max_iter, diag, ndiag)
        IMPLICIT none
        integer, intent(in)::nsize
        integer, intent(in)::max_iter
        real(8), intent(in)::diag, ndiag  !diag->main diag value / ndiag->neighbours
        real(8), allocatable::A(:,:),b(:),x(:),errn(:),verr(:)
        real(8)::tol
        integer::row,i,j,N, n_iter
        tol = 1.d-15
        N = nsize * nsize
        print *, "GMRES POISSON MGSR 2D TEST"
        print *,"PARAMS ","N=",N," MAX ITER=",max_iter," TOL=",tol," DIAG=",diag," NEIGHBOURS=",ndiag

        allocate(A(N,N),b(N),x(N),errn(max_iter))

        do i=1,nsize
            do j=1,nsize
                row = i+(j-1)*nsize
                A(row,row) =  diag
                if (row > 1) A(row,row-1)   = ndiag
                if (row < N) A(row,row+1)   = ndiag
                if (j > 1) A(row-nsize,row) = ndiag
                if (j < N) A(row+nsize,row) = ndiag
            end do
        end do
        x = 1.0d0
        b = matmul(A,x)

        print *, "Solving Poisson 2D (N =", N, ")..."
        call gmres_mgsr(A,b,x,max_iter,tol,errn,verr,n_iter)
        print *, "Num iterations util convergence:", n_iter, "(MAX=",max_iter,")"
        print *, "Final residual:", errn(n_iter)
        print *, "Final ||I - V.t * V||:", verr(n_iter)
        print *, "Max error:", maxval(abs(x - 1.0d0))
        print *, "Showing first 10 elements of the solution (must be 1.0)..."
        print "(10F12.8)", x(1:10)
        print *, "L2 Sol. Norm: ", norm2(x-1.0d0)
    END SUBROUTINE
    SUBROUTINE test_QR_HH(n)
        integer, intent(in)::n
        real(8),allocatable::A(:,:),H(:,:),Q(:,:)
        integer::i,j

        allocate(A(n,n),H(n,n))
        do i=1,n
            do j=i,n
                A(i,j) = 1.0d0
                if (i==j) A(i,j)=3.0d0
                if (i/=j) A(j,i)=-1.0d0
            end do
        end do
        print *,"A="
        do i=1,n
            print "(5F12.8)", A(i,1:5)
        end do
        call householder_sb(A,H)
        print *,"R="
        do i=1,n
            print "(5F12.8)", A(i,1:5)
        end do
        call HH_buildQ(H,Q)
        print *,"Q="
        do i=1,n
            print "(5F12.8)", Q(i,1:5)
        end do
        A = matmul(Q,A)
        print *,"QR="
        do i=1,n
            print "(5F12.8)", A(i,1:5)
        end do
    END SUBROUTINE

    SUBROUTINE test_Arnoldi(M)
        integer, intent(in)::M
        real(8), allocatable::A(:,:)    !input matrix
        real(8), allocatable::b(:)      !w vector
        real(8), allocatable::H(:,:)    !Hessenberg
        real(8), allocatable::P(:,:)    !reflections
        real(8)::tol,f_err
        integer::i,n_iter,max_iter!number of iterations
        allocate(A(M,M))
        A(1,1)=4;A(1,2)=1;A(3,1)=0
        A(2,1)=1;A(2,2)=3;A(2,3)=1
        A(3,1)=0;A(3,2)=1;A(3,3)=2
        max_iter=M;tol=1.e-6
        print *,"MAX ITER=",max_iter
        allocate(H(M+1,M),b(M))
        b(1)=1.0d0;b(2)=0.0d0
        call arnoldi_hh(A,b,max_iter,tol,P,H,f_err,n_iter)
        print *, "A="
        do i=1,M
            print "(3F12.8)", A(i,1:3)
        end do
        print *, "P="
        do i=1,M
            print "(3F12.8)", P(i,1:3)
        end do
        print *, "H="
        do i=1,M+1
            print "(3F12.8)", H(i,1:3)
        end do
    END SUBROUTINE
    SUBROUTINE TEST_GMRES_HH()
        real(8), allocatable::A(:,:),b(:),x(:),final_err(:),v_err(:)
        real(8)::tol
        integer::max_iter,n_iter,n
        max_iter = 3
        n=3
        tol=1e-6
        allocate(A(n,n),b(n))
        A(1,1)=10;A(1,2)=1;A(1,3)=2
        A(2,1)= 1;A(2,2)=10;A(2,3)=-1
        A(3,1)=1;A(3,2)=1;A(3,3)=10
        b(1)=13;b(2)=10;b(3)=12
        call gmres_hh(A,b,x,max_iter,tol,final_err,v_err,n_iter)
        print *, "Num iterations util convergence:", n_iter, "(MAX=",max_iter,")"
        print *, "Final residual:", final_err(n_iter)
        print *, "Showing solution..."
        print "(3F12.8)", x(1:3)
    END SUBROUTINE
end program
