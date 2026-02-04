PROGRAM TEST_POISSON_MF
    use omp_lib
    use gmres_mf
    use matrix_utils
    implicit none
    integer::nsize,max_iter,n_args
    character(len=32)::arg_str

    n_args = command_argument_count()
    if (n_args < 1) then
        print *, "usage ./test_poisson <size> <max_iterations>"
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
        tol = 1.d-15
        !$omp parallel
            if (omp_get_thread_num() == 0) num_threads = omp_get_num_threads()
        !$omp end parallel
        allocate(b(nsize*nsize), x(nsize*nsize))
        x = 1.0d0
        call stv_poisson2(x,b,nsize)! b = A*1 all solutions must be 1.0
        write(*,'(A)') 'GMRES Poisson 2D Test Matrix Free (Householder version)'
        write(*,'(A I8 A18 I5 A8 ES10.2 A10 I2)') "N VARS=", nsize*nsize, " MAX ITERS/STAGE=", max_iter, & 
              & " TOL=", tol, "THREADS=",num_threads
        start_time = omp_get_wtime()
        call gmres_hh_mf(stv_poisson2,b,x,max_iter,tol,errn,verr,n_iter,n_stages)
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
        tol = 1.d-15
        !$omp parallel
            if (omp_get_thread_num() == 0) num_threads = omp_get_num_threads()
        !$omp end parallel
        allocate(b(nsize*nsize), x(nsize*nsize))
        x = 1.0d0
        call stv_poisson2(x,b,nsize)! b = A*1 all solutions must be 1.0
        write(*,'(A)') 'GMRES Poisson 2D Test Matrix Free (MGSR version)'
        write(*,'(A I8 A18 I5 A8 ES10.2 A10 I2)') "N VARS=", nsize*nsize, " MAX ITERS/STAGE=", max_iter, & 
              & " TOL=", tol, "THREADS=",num_threads
        start_time = omp_get_wtime()
        call gmres_mgsr_mf_mpo(stv_poisson2,b,x,max_iter,tol,errn,verr,n_iter,n_stages)
        end_time = omp_get_wtime()
        write(*,'(A30, I8, A10, I4)') 'Iterations until convergence:', (n_stages-1)*max_iter+n_iter, ' Stages=', n_stages
        write(*,'(A30, ES12.4)') "Final ||I - V.t * V||:", verr(n_iter)
        write(*,'(A30, ES12.4)') 'Final residual:', errn(n_iter) 
        write(*,'(A30, ES12.4)') 'Max error L_max:', maxval(abs(x - 1.0d0))
        write(*,'(A30, ES12.4)') 'L2 norm:', norm2(x - 1.0d0)
        write(*,'(A30, 10F10.4)') 'First 10 solution elements', x(1:10)
        write(*,'(A30, F12.4, A)') 'Elapsed time:', end_time-start_time, ' secs.'
    END SUBROUTINE test2

    SUBROUTINE stv_poisson(x,y,n)   !y = A * x
        real(8), intent(in):: x(:)
        real(8), intent(out):: y(:)
        integer, intent(in)::n
        integer i,j,idx
        !$omp do collapse(2), private(idx)
        do j=1,n
            do i=1,n
                idx = i + (j-1)*n
                y(idx) = 4.0 * x(idx)
                if (i > 1) y(idx) = y(idx) -1.0 * x(idx-1)
                if (i < n) y(idx) = y(idx) -1.0 * x(idx+1)
                if (j > 1) y(idx) = y(idx) -1.0 * x(idx-n)
                if (j < n) y(idx) = y(idx) -1.0 * x(idx+n)
            end do
        end do
        !$omp end do
    END SUBROUTINE stv_poisson

    SUBROUTINE stv_poisson2(x,y,n)   !y = A * x
        real(8), intent(in):: x(:)
        real(8), intent(out):: y(:)
        integer, intent(in)::n
        integer i,j,idx
        !$omp do collapse(2) private(idx)
        do j=2,n-1 !col
            do i=2,n-1 !row
                idx = i + (j-1)*n
                y(idx) = 4.0d0*x(idx)-1.0d0*(x(idx-1)+x(idx+1)+x(idx+n)+x(idx-n))
            end do
        end do
        !$omp end do
        !$omp do nowait
        do i=2,n-1 !borders col=1
            y(i)   = 4.0d0*x(i) - 1.0d0*(x(i-1)+x(i+1)+x(i+n)) !col=1
        end do
        !$omp do private(idx)
        do i=2,n-1 !borders col=n
            idx = n*n-n+i
            y(idx) = 4.0d0*x(idx)-1.0d0*(x(idx-1)+x(idx+1)+x(idx-n))
        end do
        !$omp end do
        !$omp do private(idx) !borders row=1
        do i=2,n-1
            idx  = (i-1)*n+1 
            y(idx) = 4.0d0*x(idx)-1.0d0*(x(idx+1)+x(idx+n)+x(idx-n))
        end do
        !$omp end do
        !$omp do private(idx)
        do i=2,n-1
            idx  = (i-1)*n+n !row=n
            y(idx) = 4.0d0*x(idx)-1.0d0*(x(idx-1)+x(idx+n)+x(idx-n))
        end do
        !$omp end do
        !corners
        !$omp single
        y(1)  = 4.0d0*x(1) - 1.0d0*(x(2)+x(1+n))     !row=1, col=1
        y(n)  = 4.0d0*x(n) - 1.0d0*(x(n-1)+x(n+n))   !row=n, col=1
        idx   = n*(n-1)+1                                       
        y(idx)= 4.0d0*x(idx)-1.0d0*(x(idx+1)+x(idx-n)) !row=1, col=n
        idx   = n*(n-1)+n                                    
        y(idx)= 4.0d0*x(idx)-1.0d0*(x(idx-1)+x(idx-n)) !row=n, col=n                             
        !$omp end single
    END SUBROUTINE stv_poisson2

END PROGRAM TEST_POISSON_MF