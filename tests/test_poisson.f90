PROGRAM TEST_POISSON
    use gmres_mgsr_mod
    use gmres_hh_mod
    use poisson
    implicit none
    integer::nsize,max_iter,n_args
    character(len=32)::arg_str

    n_args = command_argument_count()
    if (n_args < 1) then
        print *, "usage ./test_poisson <size> <max iterations per stage>"
        stop
    end if
    call get_command_argument(1, arg_str)
    read(arg_str, *) nsize
    call get_command_argument(2, arg_str)
    read(arg_str, *) max_iter
    write(*, '(60("-"))') 
    call test_Poisson_2DMGSR(nsize,max_iter)
    write(*, '(60("-"))')
    call test_Poisson_2DHH(nsize,max_iter)
    write(*, '(60("-"))')
CONTAINS
    SUBROUTINE test_Poisson_2DHH(nsize, max_iter)
        IMPLICIT none
        integer, intent(in)::nsize
        integer, intent(inout)::max_iter
        real(8), allocatable:: A(:,:), b(:),x(:),errn(:),verr(:)
        real(8)::tol
        integer::n_iter,n_stages
        real(8)::start_time, end_time
        tol = 1.d-15

        write(*,'(A)') 'GMRES Poisson 2D Test (Householder Restarted version)'
        write(*,'(A I5 A9 I5 A5 ES10.2)') "N=", nsize*nsize, " ITER/STAGE=", max_iter, " TOL=", tol
        call generate_matrix(A,nsize)
        allocate(b(nsize*nsize))
        b = 1.0d0
        b = matmul(A,b)
        !call jacobi(A,b)
        !allocate(errn(max_iter))
        call cpu_time(start_time)
        call gmres_hh_dense(A,b,x,max_iter,tol,errn,verr,n_iter,n_stages)
        call cpu_time(end_time)
        write(*,'(A30, I6, A10, I3)') 'Iterations until convergence:', (n_stages-1)*max_iter+n_iter, ' Stages=', n_stages
        write(*,'(A30, ES12.4)') "Final ||I - V.t * V||:", verr(n_iter)
        write(*,'(A30, ES12.4)') 'Final residual:', errn(n_iter) 
        write(*,'(A30, ES12.4)') 'Max error L_max:', maxval(abs(x - 1.0d0))
        write(*,'(A30, ES12.4)') 'L2 norm:', norm2(x - 1.0d0)
        write(*,'(A30, 10F10.4)') 'First 10 solution elements', x(1:10)
        write(*,'(A30, F10.6, A)') 'Elapsed time:', end_time-start_time, ' secs.'
    END SUBROUTINE test_Poisson_2DHH

    SUBROUTINE test_Poisson_2DMGSR(nsize, max_iter)
        integer, intent(in)::nsize
        integer, intent(in)::max_iter
        real(8), allocatable::A(:,:),b(:),x(:),errn(:),verr(:)
        real(8)::tol
        integer::n_iter, n_stages
        real(8)::start_time, end_time
        tol = 1.d-15
        
        write(*,'(A)') 'GMRES Poisson 2D Test (MGSR restarted version)'
        write(*,'(A I5 A I6 A ES10.2)') "N=", nsize*nsize, " ITER/STAGE=", max_iter, " TOL=", tol

        !allocate(errn(max_iter),verr(max_iter))

        call generate_matrix(A,nsize)
        allocate(b(nsize*nsize))
        b = 1.0d0
        b = matmul(A,b)
    
        call cpu_time(start_time)
        call gmres_mgsr_dense(A,b,x,max_iter,tol,errn,verr,n_iter,n_stages)
        call cpu_time(end_time)
        write(*,'(A30, I6, A10, I3)') 'Iterations until convergence:', (n_stages-1)*max_iter+n_iter, ' Stages=', n_stages
        write(*,'(A30, ES12.4)') "Final ||I - V.t * V||:", verr(n_iter)
        write(*,'(A30, ES12.4)') 'Final residual:', errn(n_iter) 
        write(*,'(A30, ES12.4)') 'Max error L_max:', maxval(abs(x - 1.0d0))
        write(*,'(A30, ES12.4)') 'L2 norm:', norm2(x - 1.0d0)
        write(*,'(A30, 10F10.4)') 'First 10 solution elements', x(1:10)
        write(*,'(A30, F10.6, A)') 'Elapsed time:', end_time-start_time, ' secs.'
    END SUBROUTINE
END PROGRAM TEST_POISSON