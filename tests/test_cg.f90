PROGRAM test_cg
    use poisson
    use omp_lib
    use conjugate_gradient
    use chebyshev_precond
    use utils
    implicit none
    character(len=30)::desc    !Description of the test
    integer::ntests     !Num tests performed
    integer::nsize      !size of the grid
    integer::iter       !Maximun number of iterations per restart
    real(8), allocatable::b(:),x(:)     !vectors
    real(8)::tol        !tolerance of the test
    real(8)::err        !residual of the test
    real(8)::l2, linf   !l2 and linf norms of the solution
    integer::num_threads!num_iterations, num_stages performed
    real(8)::start_time, end_time !clock vars
    real(8), allocatable :: params(:) !params for Chebyshev preconditioner max/min eigens
    integer:: i !loop
    tol    = 1.d-9
    nsize  = 300 !grid size
    ntests = 15
    iter = 2000
    call omp_set_dynamic(.false.)
    call omp_set_num_threads(6)
    !$omp parallel 
    if (omp_get_thread_num() == 0) num_threads = omp_get_num_threads()
    !$omp end parallel 
    allocate(params(2))
    params(1) = 8.2d0; params(2) = 0.2d0 
    desc = " "
    print *,'Conjugate Gradient Convergence Test'
    write(*,'(A25 I4 A25 I2)') "Number of Tests:",ntests, "Number of Threads: ",num_threads
    write(*, '(150("-"))')
    write(*,'(A14 A14 A14 A14 A14 A14 A14 A14)') '# Test', 'Grid Size', 'Num Iters', 'Tol.', 'Error','L2','LINF','Time'
    !call print_header(header)
    do i=1, ntests
        iter = 10000
        allocate(b(nsize*nsize), x(nsize*nsize))
        x = 1.0d0
        call stvec(x,b,nsize)! b = A*1 all solutions must be 1.0
        start_time = omp_get_wtime()
        call pcg_omp(stvec,b,x,tol,iter,err,cbpr2,params)
        !call cg_omp(stvec,b,x,tol,iter,err)
        end_time = omp_get_wtime()
        l2 = norm2(x - 1.0d0)
        linf = maxval(abs(x - 1.0d0))
        write(*,'(I14 I14 I14 ES14.2 ES14.4 ES14.4 ES14.4 F14.4 F14.8)') i, nsize*nsize,iter,tol, err,l2, &
        & linf, end_time - start_time, x(1)
        deallocate(x,b)
        nsize = nsize+50
    end do
    write(*, '(150("-"))')
END PROGRAM test_cg