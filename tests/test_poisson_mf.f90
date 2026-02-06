PROGRAM TEST_POISSON_MF
    use poisson
    use omp_lib
    use gmres_mgsr_mod
    use gmres_hh_mod
    use chebyshev_precond
    implicit none
    integer::nsize,max_iter,n_args
    character(len=32)::arg_str

    n_args = command_argument_count()
    if (n_args < 1) then
        print *, "usage ./test_poisson <grid size> <iterations per restart>"
        stop
    end if
    call get_command_argument(1, arg_str)
    read(arg_str, *) nsize
    call get_command_argument(2, arg_str)
    read(arg_str, *) max_iter
    write(*, '(60("-"))') 
    call test(nsize,max_iter)
    write(*, '(60("-"))')
    call test2(nsize,max_iter)
    write(*, '(60("-"))')
CONTAINS

    SUBROUTINE test(nsize, max_iter) 
        integer, intent(in)::nsize
        integer, intent(inout)::max_iter
        real(8), allocatable::b(:),x(:),errn(:),verr(:)
        real(8)::tol
        integer::n_iter,n_stages,num_threads
        real(8)::start_time, end_time
        real(8), allocatable :: params(:)
        tol = 1.d-15
        if (omp_get_thread_num() == 0) num_threads = omp_get_num_threads()
        allocate(b(nsize*nsize), x(nsize*nsize),params(2))
        params(1) = 8.2d0; params(2) = 0.2d0
        x = 1.0d0
        call stvec(x,b,nsize)! b = A*1 all solutions must be 1.0
        write(*,'(A)') 'GMRES Poisson 2D Test Matrix Free (MGSR Chebyshev)'
        write(*,'(A I8 A18 I5 A8 ES10.2 A10 I2)') "N VARS=", nsize*nsize, " MAX ITERS/STAGE=", max_iter, & 
              & " TOL=", tol, "THREADS=",num_threads
        start_time = omp_get_wtime()
        call gmres_mgsr_mf(stvec,b,x,max_iter,tol,errn,verr,n_iter,n_stages,cbpr2,params)
        end_time = omp_get_wtime()
        write(*,'(A30, I8, A10, I4)') 'Iterations until convergence:', (n_stages-1)*max_iter+n_iter, ' Stages=', n_stages
        write(*,'(A30, ES12.4)') "Final ||I - V.t * V||:", verr(n_iter)
        write(*,'(A30, ES12.4)') 'Final residual:', errn(n_iter) 
        write(*,'(A30, ES12.4)') 'Max error L_max:', maxval(abs(x - 1.0d0))
        write(*,'(A30, ES12.4)') 'L2 norm:', norm2(x - 1.0d0)
        write(*,'(A30, 10F10.4)') 'First 10 solution elements', x(1:10)
        write(*,'(A30, F12.4, A)') 'Elapsed time:', end_time-start_time, ' secs.'
    END SUBROUTINE test

    SUBROUTINE test2(nsize, max_iter) 
        integer, intent(in)::nsize
        integer, intent(inout)::max_iter
        real(8), allocatable::b(:),x(:),errn(:),verr(:)
        real(8)::tol
        integer::n_iter,n_stages,num_threads
        real(8)::start_time, end_time
        real(8), allocatable :: params(:)
        tol = 1.d-15
        !$omp parallel
            if (omp_get_thread_num() == 0) num_threads = omp_get_num_threads()
        !$omp end parallel
        allocate(b(nsize*nsize), x(nsize*nsize), params(2))
        params(1) = 8.2d0; params(2) = 0.2d0
        x = 1.0d0
        call stvec(x,b,nsize)! b = A*1 all solutions must be 1.0
        write(*,'(A)') 'GMRES Poisson 2D Test Matrix Free (MGSR MPOpen Chebyshev version)'
        write(*,'(A I8 A18 I5 A8 ES10.2 A10 I2)') "N VARS=", nsize*nsize, " MAX ITERS/STAGE=", max_iter, & 
              & " TOL=", tol, "THREADS=",num_threads
        start_time = omp_get_wtime()
        call gmres_mgsr_omp(stvec,b,x,max_iter,tol,errn,verr,n_iter,n_stages,cbpr2,params)
        end_time = omp_get_wtime()
        write(*,'(A30, I8, A10, I4)') 'Iterations until convergence:', (n_stages-1)*max_iter+n_iter, ' Stages=', n_stages
        write(*,'(A30, ES12.4)') "Final ||I - V.t * V||:", verr(n_iter)
        write(*,'(A30, ES12.4)') 'Final residual:', errn(n_iter) 
        write(*,'(A30, ES12.4)') 'Max error L_max:', maxval(abs(x - 1.0d0))
        write(*,'(A30, ES12.4)') 'L2 norm:', norm2(x - 1.0d0)
        write(*,'(A30, 10F10.4)') 'First 10 solution elements', x(1:10)
        write(*,'(A30, F12.4, A)') 'Elapsed time:', end_time-start_time, ' secs.'
    END SUBROUTINE test2


END PROGRAM TEST_POISSON_MF