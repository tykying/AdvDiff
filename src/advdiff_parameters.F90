module advdiff_parameters

  use advdiff_precision

  implicit none

  private

  public :: L, ncell, dt, mtime, nd_scale, input_unit, output_unit, &
    & non_dimensionalise_parameters, &
    & SCR_DT, DAT_DT

  ! Physical parameters from 
  ! 1) Karabasov, Berloff and Goloviznin, Ocean Modelling, 2009, pp. 155 - 168
  ! 2) Maddison, Marshallb and Shiptonc, Ocean Modelling, 2015 
  real(kind = dp) :: L = 3840.0D5    ! Domain size (physical length, in cm or kmD5)
 
  integer, parameter :: ncell = 32   ! Number of control cell
  real(kind = dp) :: dt = 10800.0_dp      ! Timestep (in second)
  real(kind = dp) :: mtime = 3600.0_dp * 24.0_dp * 365.25_dp * 64.0_dp
  
  ! Scale factor: to convert non-dimensional quantities back into dimensional form
  ! e.g. in qg.F90, dt*nd_scale = physical time step 
  ! Will be altered in the non_dimensionalise_parameters subrountine
  real(kind = dp) :: nd_scale = 1.0_dp

  real(kind=dp) :: SCR_DT = 3600.0_dp * 24.0_dp * 120.00_dp * 1.0_dp    ! Time intv to print to screen
  real(kind=dp) :: DAT_DT = 3600.0_dp * 24.0_dp * 365.25_dp * 1.0_dp  ! Time intv to write binary to hdd

  integer, parameter :: input_unit = 10, output_unit = 11
  
contains

  subroutine non_dimensionalise_parameters(sc)  ! sc = L/ncell = 3840D5/512
    real(kind = dp), intent(in) :: sc

    real(kind = dp) :: lsc

    ! Copy the input. This allows a parameter to be passed in.
    lsc = sc

    L = L / lsc
    dt = dt / lsc
    mtime = mtime / lsc
    
    SCR_DT = SCR_DT / lsc
    DAT_DT = DAT_DT / lsc
    
    nd_scale = nd_scale * lsc
    
  end subroutine non_dimensionalise_parameters

  pure real(kind=dp) function x(i, m)
    integer, intent(in) :: i
    integer, intent(in) :: m

    x = L * real(i - 1, kind = dp) / real(m - 1, kind = dp)
    
  end function x

  pure real(kind=dp) function y(j, n)
    integer, intent(in) :: j
    integer, intent(in) :: n

    y = L * real(j - 1, kind = dp) / real(n - 1, kind = dp)

  end function y

end module advdiff_parameters
