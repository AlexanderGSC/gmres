!!> TEST1: testing restart size for fixed size problem
PROGRAM strong_scaling
    use poisson
    use omp_lib
    use gmres_mgsr_mod
    use gmres_hh_mod
    use chebyshev_precond
    use utils
    implicit none
    character(len=110)::header !Header print info
    character(len=30)::desc !Description of the test
    character(len=32)::arg_str !Argument string
    integer::ntests     !Num tests performed
    integer::nsize      !size of the grid
    integer::max_iter   !Maximun number of iterations per restart
    integer::n_args
    real(8), allocatable::b(:),x(:),errn(:),verr(:) !gmres params
    real(8)::tol        !tolerance of the test
    integer::n_iter,n_stages,num_threads !num_iterations, num_stages performed
    real(8)::start_time, end_time !clock vars
    real(8), allocatable :: params(:) !params for Chebyshev preconditioner max/min eigens
    integer:: i !loop
    tol    = 1.d-15
    ntests = 6
    max_iter = 50
    call omp_set_dynamic(.false.)
    n_args = command_argument_count()
    if (n_args < 1) then
        print *, "usage ./strong_scaling <grid size>"
        stop
    end if
    call get_command_argument(1, arg_str)
    read(arg_str, *) nsize
    allocate(params(2))
    params(1) = 8.2d0; params(2) = 0.2d0

    allocate(b(nsize*nsize), x(nsize*nsize))
    x = 1.0d0
    call stvec(x,b,nsize)! b = A*1 all solutions must be 1.0

    print *,'GMRES Strong Scaling Test (MGSR with Chebyshev precond)'
    write(header,'(A25 I4 A25 I2)') "Number of Tests:",ntests
    call print_header(header)
    do i=1, ntests
        call omp_set_num_threads(i)
        !$omp parallel 
        if (omp_get_thread_num() == 0) num_threads = omp_get_num_threads()
        !$omp end parallel 
        write (desc,'(A15 I2)') "threads=",num_threads
        start_time = omp_get_wtime()
        call gmres_mgsr_omp(stvec,b,x,max_iter,tol,errn,verr,n_iter,n_stages,cbpr2,params)
        end_time = omp_get_wtime()
        call print_line(i, nsize*nsize,end_time-start_time, (n_stages-1)*max_iter+n_iter, n_stages, max_iter,tol, errn(n_iter), &
            & verr(n_iter), norm2(x - 1.0d0), maxval(abs(x - 1.0d0)), desc)
    end do
    write(*, '(150("-"))')
END PROGRAM strong_scaling