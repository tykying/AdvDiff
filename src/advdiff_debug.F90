#include "advdiff_configuration.h"

module advdiff_debug
  use advdiff_precision
  use mpi

  implicit none

  private
  
  public :: abort_handle
  public :: print_array
  
  interface print_array
    module procedure print_array_2D, print_array_1D
  end interface print_array
  
contains
  
  subroutine abort_handle(msg, file, line)
    character(len = *), intent(in) :: msg
    character(len = *), intent(in) :: file
    integer, intent(in) :: line

    integer :: ierr

    write(0, "(a,a,a,i0,a,a)") "ERROR on line ", trim(file), ":", line, ": ", &
      & trim(msg)
    flush(0)
    call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
    
  end subroutine abort_handle

  subroutine print_array_2D(array, array_name)
    real(kind=dp), dimension(:,:), intent(in) :: array
    character(len = *), optional, intent(in) :: array_name
    character(len = 20) :: collen
    integer :: i
    
    write(collen, "(i0)") size(array, 2)
    
    if (present(array_name)) then
      write(6, "(a, i0, a, i0)") trim(array_name)//": " , size(array, 1), " by ", size(array, 2)
    else
      write(6, "(i0, a, i0)") size(array, 1), " by ", size(array, 2)
    end if

    do i = 1, size(array, 1)
      write(6, "("//trim(collen)//sp_chr//")") array(i, :)
    end do
    
  end subroutine print_array_2D
  
  subroutine print_array_1D(array, array_name)
    real(kind=dp), dimension(:), intent(in) :: array
    character(len = *), optional, intent(in) :: array_name
    character(len = 20) :: collen
    
    write(collen, "(i0)") size(array, 1)
    
    if (present(array_name)) then
      write(6, "(a, i0, a, i0)") trim(array_name)//": " , size(array, 1), " by ", 1
    else
      write(6, "(i0, a, i0)") size(array, 1), " by ", 1
    end if

    write(6, "("//trim(collen)//sp_chr//")") array(:)
    
  end subroutine print_array_1D

  
end module advdiff_debug
