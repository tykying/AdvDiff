#include "advdiff_configuration.h"

module advdiff_timing

  use advdiff_precision
#if OMP0MPI1 == 0
  use omp_lib
#endif

  implicit none

  private

  public :: timer, reset, start, stop, print

  type timer
    real(kind = dp) :: time = 0.0_dp
    logical :: running = .false.
    real(kind = dp) :: start_cpu_time

#if OMP0MPI1 == 0
    real(kind = dp) :: wtime = 0.0_dp
    real(kind = dp) :: start_wtime
#endif
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
    
#if OMP0MPI1 == 0
    tmr%wtime = 0.0_dp
#endif
  end subroutine reset_timer

  subroutine start_timer(tmr)
    type(timer), intent(inout) :: tmr
    
    !assert(.not. tmr%running)
    
    call cpu_time(tmr%start_cpu_time)
    tmr%running = .true.
    
#if OMP0MPI1 == 0
    tmr%start_wtime = omp_get_wtime()
#endif
  end subroutine start_timer

  subroutine stop_timer(tmr)
    type(timer), intent(inout) :: tmr

    real(kind = dp) :: end_cpu_time
#if OMP0MPI1 == 0
    real(kind = dp) :: end_wtime
#endif

    !assert(tmr%running)
    
    call cpu_time(end_cpu_time)
    tmr%time = tmr%time + end_cpu_time - tmr%start_cpu_time
    
    tmr%running = .false.
    
#if OMP0MPI1 == 0
    tmr%wtime = tmr%wtime + omp_get_wtime() - tmr%start_wtime
#endif
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
#if OMP0MPI1 == 0
    real(kind = dp) :: end_wtime, wtime
#endif

    if(tmr%running) then
      call cpu_time(end_cpu_time)
      time = tmr%time + end_cpu_time - tmr%start_cpu_time
      
#if OMP0MPI1 == 0
      end_wtime = omp_get_wtime()
      wtime = tmr%wtime + end_wtime - tmr%start_wtime
#endif
    else
      time = tmr%time
      
#if OMP0MPI1 == 0
      wtime = tmr%wtime
#endif
    end if

    write(6, "(a,a,a,a,f15.6)") prefix, trim(name), &
      & name_pad(1:max(len(name_pad) - len_trim(name), 0)), &
      & " cpu time (s) = ", time

#if OMP0MPI1 == 0
    write(6, "(a,a,a,a,f15.6)") prefix, trim(name), &
      & name_pad(1:max(len(name_pad) - len_trim(name), 0)), &
      & " wall time (s) = ", wtime
#endif
  end subroutine print_timer_prefix

end module advdiff_timing
