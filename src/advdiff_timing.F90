module advdiff_timing

  use advdiff_precision
  use omp_lib

  implicit none

  private

  public :: timer, reset, start, stop, print

  type timer
    real(kind = dp) :: time = 0.0_dp
    logical :: running = .false.
    real(kind = dp) :: start_cpu_time
    
!     real(kind = dp) :: wtime = 0.0_dp
!     real(kind = dp) :: start_wtime
  end type timer

  interface reset
    module procedure reset_timer
  end interface reset

  interface start
    module procedure start_timer
  end interface start

  interface stop
    module procedure stop_timer
  end interface stop

  interface print
    module procedure print_timer, print_timer_prefix
  end interface print

contains

  subroutine reset_timer(tmr)
    type(timer), intent(out) :: tmr

    tmr%time = 0.0_dp
    tmr%running = .false.
    
!     tmr%wtime = 0.0_dp
    
  end subroutine reset_timer

  subroutine start_timer(tmr)
    type(timer), intent(inout) :: tmr
    
    !assert(.not. tmr%running)

    call cpu_time(tmr%start_cpu_time)
    tmr%running = .true.
    
!     tmr%start_wtime = omp_get_wtime()

  end subroutine start_timer

  subroutine stop_timer(tmr)
    type(timer), intent(inout) :: tmr

    real(kind = dp) :: end_cpu_time
    real(kind = dp) :: end_wtime

    !assert(tmr%running)
    
    call cpu_time(end_cpu_time)
    tmr%time = tmr%time + end_cpu_time - tmr%start_cpu_time
    
    tmr%running = .false.
    
!     tmr%wtime = tmr%wtime + omp_get_wtime() - tmr%start_wtime
    
  end subroutine stop_timer

  subroutine print_timer(tmr, name)
    type(timer), intent(in) :: tmr
    character(len = *), intent(in) :: name

    call print(tmr, name, prefix = "")

  end subroutine print_timer

  subroutine print_timer_prefix(tmr, name, prefix)
    type(timer), intent(in) :: tmr
    character(len = *), intent(in) :: name
    character(len = *), intent(in) :: prefix

    character(len = *), parameter :: name_pad = "                    "
    real(kind = dp) :: end_cpu_time, time
    real(kind = dp) :: end_wtime, wtime

    if(tmr%running) then
      call cpu_time(end_cpu_time)
      time = tmr%time + end_cpu_time - tmr%start_cpu_time
      
!       end_wtime = omp_get_wtime()
!       wtime = tmr%wtime + end_wtime - tmr%start_wtime
    else
      time = tmr%time
      
!       wtime = tmr%wtime
    end if

    write(6, "(a,a,a,a,f15.6)") prefix, trim(name), &
      & name_pad(1:max(len(name_pad) - len_trim(name), 0)), &
      & " cpu time (s) = ", time

!     write(6, "(a,a,a,a,f15.6)") prefix, trim(name), &
!       & name_pad(1:max(len(name_pad) - len_trim(name), 0)), &
!       & " wall time (s) = ", wtime
      
  end subroutine print_timer_prefix

end module advdiff_timing