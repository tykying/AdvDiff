module advdiff_precision
  implicit none

  private

  public :: dp, dp_chr, int_chr, max_int_len, field_name_len, max_line_len
  public :: sp_chr

  integer, parameter :: dp = kind(0.0D0)
  character(len = *), parameter :: dp_chr = "es24.16e3"
  character(len = *), parameter :: sp_chr = "es12.3e2"
  character(len = *), parameter :: int_chr = "i11"
  integer, parameter :: max_int_len = &
    & int(log10(real(huge(0), kind = dp)) + 0.5D0) + 1
  integer, parameter :: field_name_len = 256
  integer, parameter :: max_line_len = 1024

  
end module advdiff_precision

!  s =    3.000E+01 days ; dt =    1.200E+01 hours
! 1-th step: Current log-likelihood = -1.9460586552253601E+005
! 2-th step: Current log-likelihood = -1.9438510069490981E+005

!  s =    3.000E+01 days ; dt =    2.000E+00 hours
! 1-th step: Current log-likelihood = -1.9460586191927362E+005
! 2-th step: Current log-likelihood = -1.9438508487030680E+005
! 3-th step: Current log-likelihood = -1.9415750535022139E+005

! dt = 12 hours, 250th optim: -1.8611737302370396E+005, max = 0.0032
