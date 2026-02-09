MODULE utils
    implicit none
    private
    public :: print_results
    public :: print_table
    public :: print_line
    public :: print_header
CONTAINS
    subroutine print_results(num_vars,iterations,max_iter, stages, verr, errn, tol, lmax, l2, x, cpu_time)
        integer, intent(in) :: num_vars, stages, iterations, max_iter
        real(8), intent(in) :: verr, errn, tol, lmax, l2, cpu_time, x(:)
        write(*,'(A I5 A I6 A ES10.2)') "N=", num_vars, " ITER/STAGE=", max_iter, " TOL=", tol
        write(*,'(A30, I8, A10, I4)') 'Iterations until convergence:', (stages-1)*max_iter+iterations, ' Stages=', stages
        write(*,'(A30, ES12.4)') "Final ||I - V.t * V||:", verr
        write(*,'(A30, ES12.4)') 'Final residual:', errn
        write(*,'(A30, ES12.4)') 'Max error L_max:', lmax
        write(*,'(A30, ES12.4)') 'L2 norm:', l2
        write(*,'(A30, 10F10.4)') 'First 10 solution elements', x(1:10)
        write(*,'(A30, F12.4, A)') 'Elapsed time:', cpu_time, ' secs.'
        write(*, '(60("-"))') 
    END subroutine print_results

    subroutine print_table(header, nvars, cpu_time, iterations, restarts, tol, errn, verr, l2, linf)
        character(len=60), intent(in) :: header
        integer, intent(in) :: nvars(:), iterations(:), restarts(:)
        real(8), intent(in):: cpu_time(:), errn(:), verr(:), l2(:), linf(:), tol(:)
        integer :: i
        print *, header
        write(*,'(A10, A10, A10, A14, A14, A14, A14, A14, A10)') "# Vars", "# Iters","# Rest", "Tol.", & 
        "L2 Norm","L_inf Norm", "Residual", "||I-V.t*V||", "Time"
        do i=1, size(nvars)
            write(*,'(I10, I10, I10, ES14.4, ES14.4, ES14.4, ES14.4, ES14.4, F10.4)') nvars(i), iterations(i), &
            & restarts(i), tol(i), l2(i), linf(i), errn(i), verr(i), cpu_time(i)
        end do
    end subroutine print_table

    subroutine print_header(header)
        character(len=110), intent(in) :: header
        print *, header
        write(*,'(A3, A10, A10, A10, A10, A14, A14, A14, A14, A14, A10, A15)') "#", "Vars", "Iters", &
        & "Restarts", "gmres(n)","Tol.", "L2 Norm","L_inf Norm", "Residual", "||I-V.t*V||", "Time", "Info"
        write(*, '(150("-"))')
    end subroutine print_header

    subroutine print_line(test, nvars, cpu_time, iterations, restarts, max_iters, tol, errn, verr, l2, linf,desc)
        integer, intent(in) :: test, nvars, iterations, restarts, max_iters
        real(8), intent(in):: cpu_time, errn, verr, l2, linf, tol
        character(len=30)::desc
        write(*,'(I3 I10, I10, I10, I10, ES14.2, ES14.4, ES14.4, ES14.4, ES14.4, F10.4 A20)') test, nvars, iterations, &
            & restarts, max_iters, tol, l2, linf, errn, verr, cpu_time, desc
    end subroutine print_line
END MODULE 