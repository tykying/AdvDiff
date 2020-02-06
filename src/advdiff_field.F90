module advdiff_field
  use advdiff_precision
  use advdiff_debug

  implicit none

  private

  public :: field, field_ptr, allocate, deallocate, zeros, &
    & scale, set, addto, addto_product, &
    & dx_field, dy_field, &
    & psi_to_uv, int_field, int_sq_field, &
    & iszero, &
    & dx_field_HoriEdge, dy_field_VertEdge, &
    & dx_field_HoriEdge_BC, dy_field_VertEdge_BC
    
  public :: indicator_field
  public :: uv_signed, K_flux
  public :: print_info, diff_int
  public :: testcase_fields, fld_x, assign_q_periodic, S_rhs_nonstat
  public :: bilinear_intpl, pwc_intpl, sine_filter, linear_intpl
  public :: neg_int
  public :: dt_CFL
  
  ! Note: field and field_ptr are different
  type field
    integer :: m, n
    real(kind = dp), dimension(:, :), pointer :: data  ! It is an array, not a pointer
    character(len = field_name_len) :: name
    integer :: type_id
  end type field

  type field_ptr
    type(field), pointer :: ptr
  end type field_ptr

  type uv_signed
    real(kind = dp), dimension(:, :), pointer :: up, um, vp, vm
    real(kind = dp), dimension(:, :), pointer :: u_abs, v_abs
    integer, dimension(:, :), pointer :: u_sgn, v_sgn
  end type uv_signed

  type K_flux
    real(kind = dp), dimension(:, :), pointer :: K11e, K21e, K12e, K22e
    real(kind = dp), dimension(:, :), pointer :: K11c, K12c, K22c
  end type K_flux

  interface allocate
    module procedure allocate_field, allocate_uv, allocate_Kflux
  end interface allocate

  interface deallocate
    module procedure deallocate_field, deallocate_uv, deallocate_Kflux
  end interface deallocate

  interface zeros
    module procedure zero_field
  end interface zeros

  interface scale
    module procedure scale_field
  end interface scale

  interface set
    module procedure set_field_scalar, set_field_field
  end interface set

  interface addto
    module procedure addto_field_field_scale
  end interface addto

  interface addto_product
    module procedure addto_product_field_fields_scale
  end interface addto_product

  interface print_info
    module procedure print_info_field
  end interface print_info

contains
  subroutine allocate_field(fld, m, n, name, type_id)
    type(field), intent(out) :: fld
    integer, intent(in) :: m
    integer, intent(in) :: n
    character(len = *), intent(in) :: name
    integer, intent(in) :: type_id
    
    fld%m = m
    fld%n = n
    fld%type_id = type_id
    
    if (fld%type_id .eq. 2) then
      ! Defined at cell centre
      allocate(fld%data(m, n))
    elseif (fld%type_id .eq. 1) then
      ! Defined at corner
      allocate(fld%data(m+1, n+1))
    else
      ! DOF type
      allocate(fld%data(m, n))
    end if

    if(len_trim(name) > field_name_len) then
      write(0, "(a,a,a)") 'WARNING: Field name "', trim(name), '" truncated'
    end if
    fld%name = trim(name)

  end subroutine allocate_field

  subroutine deallocate_field(fld)
    type(field), intent(inout) :: fld

    deallocate(fld%data)

  end subroutine deallocate_field

#include "advdiff_configuration.h"
  subroutine sine_filter(out, in)
#if FFTW3 == 1
    ! in, out: no boundaries => size(in) = [2^m+1-2, 2^n+1-2]
    use, intrinsic :: iso_c_binding 
    include 'fftw3.f03'

    real(kind=dp), dimension(:, :), intent(out) :: out
    real(kind=dp), dimension(:, :), intent(in) :: in
  
    type(C_PTR) :: plan
    real(kind=dp), dimension(:, :), allocatable :: Bpq, Bpq_pad
    real(kind=dp), dimension(:, :), allocatable :: in_copy
    
    integer :: i, j
    
    allocate(Bpq(size(in,1), size(in,2)))
    allocate(in_copy(size(in,1), size(in,2)))
    allocate(Bpq_pad(size(out,1), size(out,2)))

    in_copy = in
    
    ! Forward DCT
    plan = fftw_plan_r2r_2d(size(in,2),size(in,1), in_copy, Bpq,  &
                            FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE)
    call fftw_execute_r2r(plan, in_copy, Bpq)
    call fftw_destroy_plan(plan)
    
!     call print_array(Bpq, "Bpq")
    
    ! Zero padding: only need to pad zeros beyond Bpq
    !! Relationship with matlab dct2
    !! http://www.voidcn.com/article/p-pduypgxe-du.html
    Bpq_pad = Bpq(1:size(Bpq_pad,1), 1:size(Bpq_pad,2))
    
    plan = fftw_plan_r2r_2d(size(out,2),size(out,1), Bpq_pad, &
                          out, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE)
    call fftw_execute_r2r(plan, Bpq_pad, out)
    call fftw_destroy_plan(plan)

    out = out/(4*(size(in,1)+1)*(size(in,2)+1))
    
    deallocate(Bpq)
    deallocate(Bpq_pad)
    deallocate(in_copy)
    call fftw_cleanup()
#else
    real(kind=dp), dimension(:, :), intent(out) :: out
    real(kind=dp), dimension(:, :), intent(in) :: in
    
    call abort_handle("No FFTW3 library", __FILE__, __LINE__)
    out = 0.0_dp * in(1,1)
#endif
  end subroutine sine_filter

  subroutine linear_intpl(out, in)
    real(kind=dp), dimension(:), intent(out) :: out
    real(kind=dp), dimension(:), intent(in) :: in
    
    integer :: i, n_in, n_out
    real(kind=dp) :: x, alpha
    
    n_in = size(in, 1)
    n_out = size(out, 1)
    
    out(1) = in(1)
    out(n_out) = in(n_in)
    
    do i = 2, (n_out-1)
      x = real((i-1)*(n_in-1), kind=dp)/real(n_out-1, kind=dp) + 1.0_dp
      alpha = x - int(x)
      out(i) = alpha * in(int(x)+1) + (1.0_dp - alpha) * in(int(x))
    end do
    
  end subroutine linear_intpl
  
  subroutine bilinear_intpl(out, in, Mx, My)
    real(kind=dp), dimension(:, :), intent(out) :: out
    real(kind=dp), dimension(:, :), intent(in) :: in
    integer, intent(in) :: Mx, My  ! Factor of scaling up
    
    integer :: i, j, i_r, j_r
    
    if ((size(out,1)-1) .ne. (size(in,1)-1) * Mx) then
       write(6, *) "Mx = ", Mx
       write(6, *) "size(out) = ", size(out, 1), size(out, 2), "size(in) = ", size(in, 1), size(in, 2)
      call abort_handle("bilinear_intpl: size mismatch", __FILE__, __LINE__)
    end if
    if ((size(out,2)-1) .ne. (size(in,2)-1) * My) then
      write(6, *) "My = ", My
      write(6, *) "size(out) = ", size(out, 1), size(out, 2), "size(in) = ", size(in, 1), size(in, 2)
      call abort_handle("bilinear_intpl: size mismatch", __FILE__, __LINE__)
    end if
    
    ! Bilinear interpolations: for K
    ! Assumed type_id = 1 (at corners)
    out = 0.0_dp
    
    ! Along x
    do j = 1, size(in, 2)
      j_r = (j-1)*My + 1
      call linear_intpl(out(:, j_r), in(:, j))
    end do
    
    ! Along y
    do i = 1, size(in, 1)
      i_r = (i-1)*Mx + 1
      call linear_intpl(out(i_r, :), in(i, :))
    end do
    
    ! Along y; fill in intermediate space
    do i = 1, size(out, 1)
      ! stride = My
      call linear_intpl(out(i, :), out(i, 1:size(out,2):My))
    end do

  end subroutine bilinear_intpl
  
  subroutine pwc_intpl(rfld, fld)
    ! Assumed type_id = 2 (at centres)
    real(kind=dp), dimension(:,:), intent(inout) :: rfld
    real(kind=dp), dimension(:,:), intent(in) :: fld
    
    integer :: i_reflv, j_reflv, i_mesh, j_mesh
    integer :: i, j, i_lv, j_lv
    
    i_reflv = size(rfld, 1)/size(fld, 1)
    j_reflv = size(rfld, 2)/size(fld, 2)
    
    do j = 1, size(fld, 2)
      do i = 1, size(fld, 1)
        do j_lv = 1, j_reflv
          do i_lv = 1, i_reflv
            i_mesh = i_reflv*(i-1)+i_lv
            j_mesh = j_reflv*(j-1)+j_lv
            rfld(i_mesh, j_mesh) = fld(i, j)
          end do
        end do
      end do
    end do
    
  end subroutine pwc_intpl
  
  subroutine allocate_uv(uv_fld, psi)
    type(uv_signed), intent(out) :: uv_fld
    type(field), intent(in) :: psi
    integer :: i, j

    if (psi%type_id .ne. 1) &
      call abort_handle("psi not defined at corners", __FILE__, __LINE__)

    allocate(uv_fld%up(psi%m+1, psi%n))
    allocate(uv_fld%um(psi%m+1, psi%n))
    allocate(uv_fld%vp(psi%m, psi%n+1))
    allocate(uv_fld%vm(psi%m, psi%n+1))
    allocate(uv_fld%u_sgn(psi%m+1, psi%n))
    allocate(uv_fld%v_sgn(psi%m, psi%n+1))
    allocate(uv_fld%u_abs(psi%m+1, psi%n))
    allocate(uv_fld%v_abs(psi%m, psi%n+1))
    
    ! Temporary values
    call psi_to_uv(uv_fld%up, uv_fld%vp, psi)

    ! maxmin u and v wrt 0
    do j = 1, size(uv_fld%up, 2)
      do i = 1, size(uv_fld%up, 1)
        if (dabs(uv_fld%up(i,j)) .gt. 1D-14) then
          uv_fld%u_sgn(i,j) = int(sign(1.0_dp, uv_fld%up(i,j)))
        else
          uv_fld%u_sgn(i,j) = 0
        end if

        uv_fld%u_abs(i,j) = dabs(uv_fld%up(i,j))
        uv_fld%um(i,j) = min(uv_fld%up(i, j), 0.0_dp)
        uv_fld%up(i,j) = max(uv_fld%up(i, j), 0.0_dp)
      end do
    end do

    do j = 1, size(uv_fld%vp, 2)
      do i = 1, size(uv_fld%vp, 1)
        if (dabs(uv_fld%vp(i,j)) .gt. 1D-14) then
          uv_fld%v_sgn(i,j) = int(sign(1.0_dp, uv_fld%vp(i,j)))
        else
          uv_fld%v_sgn(i,j) = 0
        end if
  
        uv_fld%v_abs(i,j) = dabs(uv_fld%vp(i,j))
        uv_fld%vm(i,j) = min(uv_fld%vp(i, j), 0.0_dp)
        uv_fld%vp(i,j) = max(uv_fld%vp(i, j), 0.0_dp)
      end do
    end do

  end subroutine allocate_uv

  subroutine deallocate_uv(uv_fld)
    type(uv_signed), intent(inout) :: uv_fld

    deallocate(uv_fld%up)
    deallocate(uv_fld%um)
    deallocate(uv_fld%vp)
    deallocate(uv_fld%vm)
    deallocate(uv_fld%u_sgn)
    deallocate(uv_fld%v_sgn)
    deallocate(uv_fld%u_abs)
    deallocate(uv_fld%v_abs)

  end subroutine deallocate_uv

  subroutine allocate_Kflux(K_e, K11, K22, K12)
    type(K_flux), intent(inout) :: K_e
    type(field), intent(in) :: K11, K22, K12

    integer :: m, n, i, j

    real(kind=dp), dimension(:,:), pointer :: K11e, K12e, K21e, K22e

    m = K11%m
    n = K11%n

    if (K11%type_id .ne. 1) then
      call abort_handle("Wrong K11%type_id", __FILE__, __LINE__)
    end if 
    
    allocate(K_e%K11e(m+1, n))  ! Define at vertical edges
    allocate(K_e%K12e(m+1, n))  ! Define at vertical edges
    allocate(K_e%K21e(m, n+1))  ! Define at horizonal edges
    allocate(K_e%K22e(m, n+1))  ! Define at horizonal edges
    
    allocate(K_e%K11c(m+1, n+1))  ! Define at corners
    allocate(K_e%K22c(m+1, n+1))  ! Define at corners
    allocate(K_e%K12c(m+1, n+1))  ! Define at corners
    
    K11e => K_e%K11e
    K12e => K_e%K12e
    K21e => K_e%K21e
    K22e => K_e%K22e

    K_e%K11c = K11%data
    K_e%K22c = K22%data
    K_e%K12c = K12%data
    
    ! K: at corners
    ! For K11 and K12
    do j = 1, n
      do i = 1, m+1
        K11e(i,j) = 0.5_dp*(K11%data(i,j) + K11%data(i,j+1))
        K12e(i,j) = 0.5_dp*(K12%data(i,j) + K12%data(i,j+1))
      end do
    end do
    ! For K21 and K22
    do j = 1, n+1
      do i = 1, m
        K21e(i,j) = 0.5_dp*(K12%data(i,j) + K12%data(i+1,j))
        K22e(i,j) = 0.5_dp*(K22%data(i,j) + K22%data(i+1,j))
      end do
    end do
    
  end subroutine allocate_Kflux

  subroutine deallocate_Kflux(K_e)
    type(K_flux), intent(inout) :: K_e

    deallocate(K_e%K11e)
    deallocate(K_e%K21e)
    deallocate(K_e%K12e)
    deallocate(K_e%K22e)
    
    deallocate(K_e%K11c)
    deallocate(K_e%K22c)
    deallocate(K_e%K12c)
    
  end subroutine deallocate_Kflux

  subroutine zero_field(fld)
    type(field), intent(inout) :: fld

    fld%data = 0.0_dp
  end subroutine zero_field

  subroutine scale_field(fld, scale)
    type(field), intent(inout) :: fld
    real(kind = dp), intent(in) :: scale

    fld%data = fld%data * scale

  end subroutine scale_field

  subroutine set_field_scalar(fld, value)
    type(field), intent(inout) :: fld
    real(kind = dp), intent(in) :: value

    fld%data = value

  end subroutine set_field_scalar

  subroutine set_field_field(fld1, fld2)
    type(field), intent(inout) :: fld1
    type(field), intent(in) :: fld2

    fld1%data = fld2%data

  end subroutine set_field_field

  subroutine addto_field_field_scale(fld1, fld2, scale)
    type(field), intent(inout) :: fld1
    type(field), intent(in) :: fld2
    real(kind = dp), intent(in) :: scale

    integer :: i, j

    real(kind = dp), dimension(:, :), pointer :: lfld1
    real(kind = dp), dimension(:, :), pointer :: lfld2

    lfld1 => fld1%data
    lfld2 => fld2%data

    do j = 1, size(lfld1, 2)
      do i = 1, size(lfld1, 1)
        lfld1(i, j) = lfld1(i, j) + scale * lfld2(i, j)
      end do
    end do

  end subroutine addto_field_field_scale

  subroutine addto_product_field_fields_scale(fld1, fld2, fld3, scale)
    type(field), intent(inout) :: fld1
    type(field), intent(in) :: fld2
    type(field), intent(in) :: fld3
    real(kind = dp), intent(in) :: scale

    integer :: i, j

    real(kind = dp), dimension(:, :), pointer :: lfld1
    real(kind = dp), dimension(:, :), pointer :: lfld2
    real(kind = dp), dimension(:, :), pointer :: lfld3

    lfld1 => fld1%data
    lfld2 => fld2%data
    lfld3 => fld3%data

    do j = 1, size(lfld1, 2)
      do i = 1, size(lfld1, 1)
        lfld1(i, j) = lfld1(i, j) + &
                   & scale * lfld2(i, j) * lfld3(i, j)
      end do
    end do

  end subroutine addto_product_field_fields_scale

  subroutine psi_to_uv(u, v, psi_fld)
    real(kind=dp), dimension(:,:), intent(out) :: u, v
    type(field), intent(in) :: psi_fld
    integer :: i, j

    if (psi_fld%type_id .ne. 1) &
      call abort_handle("psi_to_uv: psi is not at corners", __FILE__, __LINE__)
      
    do j = 1, psi_fld%n
      do i = 1, psi_fld%m+1
        u(i, j) = psi_fld%data(i, j) - psi_fld%data(i, j+1)
      end do
    end do

    do j = 1, psi_fld%n+1
      do i = 1, psi_fld%m
        v(i, j) = psi_fld%data(i+1, j) - psi_fld%data(i, j)
      end do
    end do

  end subroutine psi_to_uv

  subroutine dx_field(d_fld, fld)
    real(kind=dp), dimension(:,:), intent(out) :: d_fld
    type(field), intent(in) :: fld
    integer :: i, j

    real(kind = dp), dimension(:, :), pointer :: lfld

    if (fld%type_id .ne. 2) &
      call abort_handle("dx_field: fld is not at centres", __FILE__, __LINE__)
      
    lfld => fld%data

    ! Interior
    do j = 1, fld%n
      do i = 2, fld%m
        d_fld(i, j) = lfld(i, j) - lfld(i-1, j)
      end do
    end do
    
    ! Boundary: 2nd extrapolation
    ! f(x) = 2f(x+h) - f(x+2h) + O(h^2)
    do j = 1, fld%n
      d_fld(1, j) = 2.0_dp*d_fld(2, j) - d_fld(3, j)
      d_fld(fld%m+1, j) = 2.0_dp*d_fld(fld%m, j) - d_fld(fld%m-1, j)
    end do

  end subroutine dx_field

  subroutine dy_field(d_fld, fld)
    real(kind=dp), dimension(:,:), intent(out) :: d_fld
    type(field), intent(in) :: fld
    integer :: i, j

    real(kind = dp), dimension(:, :), pointer :: lfld

    lfld => fld%data

    ! Interior
    do j = 2, fld%n
      do i = 1, fld%m
        d_fld(i, j) = lfld(i, j) - lfld(i, j-1)
      end do
    end do
    
    ! Boundary: 2nd extrapolation
    ! f(x) = 2f(x+h) - f(x+2h) + O(h^2)
    do i = 1, fld%m
      d_fld(i, 1) = 2.0_dp*d_fld(i, 2) - d_fld(i, 3)
      d_fld(i, fld%n+1) = 2.0_dp*d_fld(i, fld%n) - d_fld(i, fld%n-1)
    end do

  end subroutine dy_field

  subroutine dx_field_HoriEdge(d_fld, fld)
    real(kind=dp), dimension(:,:), intent(out) :: d_fld
    type(field), intent(in) :: fld
    integer :: i, j

    real(kind = dp), dimension(:, :), pointer :: lfld

    lfld => fld%data

    ! Left & Right, Near Boundary: 2nd extrapolation
    ! f'(x) = 1/h*(-3/2*f(x) +2f(x+h) -1/2*f(x+2h)) + O(h^2)
    do j = 2, fld%n
      d_fld(1, j) = -0.75_dp*(lfld(1, j-1)+lfld(1, j)) &
                    + (lfld(2, j-1)+lfld(2, j)) &
                    - 0.25_dp*(lfld(3, j-1)+lfld(3, j))  ! Left
      d_fld(fld%m, j) = 0.75_dp*(lfld(fld%m, j-1)+lfld(fld%m, j)) &
                    - (lfld(fld%m-1, j-1)+lfld(fld%m-1, j)) &
                    + 0.25_dp*(lfld(fld%m-2, j-1)+lfld(fld%m-2, j))  ! Right
    end do

    ! Interior
    do j = 2, fld%n
      do i = 2, (fld%m-1)
        d_fld(i, j) = 0.25_dp*(lfld(i+1, j) + lfld(i+1, j-1) &
                    - lfld(i-1, j) - lfld(i-1, j-1))
      end do
    end do
    
    ! Boundary: Bottom & Top
    ! 2nd extrapolation: f(x) = 2f(x+h) - f(x+2h) + O(h^2)
    do i = 1, fld%m
      d_fld(i, 1) = 2.0_dp*d_fld(i, 2) - d_fld(i, 3)
      d_fld(i, fld%n+1) = 2.0_dp*d_fld(i, fld%n) - d_fld(i, fld%n-1)
    end do

  end subroutine dx_field_HoriEdge

  subroutine dy_field_VertEdge(d_fld, fld)
    real(kind=dp), dimension(:,:), intent(out) :: d_fld
    type(field), intent(in) :: fld
    integer :: i, j

    real(kind = dp), dimension(:, :), pointer :: lfld

    lfld => fld%data

    ! Top & Bottom, Near Boundary: 2nd extrapolation
    ! f'(x) = 1/h*(-3/2*f(x) +2f(x+h) -1/2*f(x+2h)) + O(h^2)
    do i = 2, fld%m
      d_fld(i, 1) = -0.75_dp*(lfld(i-1, 1)+lfld(i, 1)) &
                    + (lfld(i-1, 2)+lfld(i, 2)) &
                    - 0.25_dp*(lfld(i-1, 3)+lfld(i, 3))   ! Bottom
      d_fld(i, fld%n) = 0.75_dp*(lfld(i-1, fld%n)+lfld(i, fld%n)) &
                    - (lfld(i-1, fld%n-1)+lfld(i, fld%n-1)) &
                    + 0.25_dp*(lfld(i-1, fld%n-2)+lfld(i, fld%n-2))  ! Top
    end do

    ! Interior
    do j = 2, (fld%n-1)
      do i = 2, fld%m
        d_fld(i, j) = 0.25_dp*(lfld(i, j+1) + lfld(i-1, j+1) &
                    - lfld(i, j-1) - lfld(i-1, j-1))
      end do
    end do
    
    ! Boundary: Left & Right
    ! 2nd extrapolation: f(x) = 2f(x+h) - f(x+2h) + O(h^2)
    do j = 1, fld%n
      d_fld(1, j) = 2.0_dp*d_fld(2, j) - d_fld(3, j)
      d_fld(fld%m+1, j) = 2.0_dp*d_fld(fld%m, j) - d_fld(fld%m-1, j)
    end do
    
  end subroutine dy_field_VertEdge
  
  subroutine dx_field_HoriEdge_BC(dqx_G, dqx_F, dqy_G, K11, K12)
    real(kind=dp), dimension(:,:), intent(inout) :: dqx_G
    real(kind=dp), dimension(:,:), intent(in) :: dqx_F, dqy_G, K11, K12
    integer :: j, m, n

    m = size(dqy_G, 1)
    n = size(dqx_F, 2)

    ! Left & Right, Near Boundary; Use BC
    do j = 2, n
      dqx_G(1, j) = 0.5_dp*( &
        - K12(1,j)/K11(1,j)*(1.5_dp*dqy_G(1,j) - 0.5_dp*dqy_G(2,j)) &
        + 0.5_dp* dqx_F(2,j) + 0.5_dp* dqx_F(2,j-1) )
      dqx_G(m, j) = 0.5_dp*( &
        - K12(m+1,j)/K11(m+1,j)*(1.5_dp*dqy_G(m,j) - 0.5_dp*dqy_G(m-1,j)) &
        + 0.5_dp* dqx_F(m,j) + 0.5_dp* dqx_F(m,j-1) )
    end do

  end subroutine dx_field_HoriEdge_BC
  
  subroutine dy_field_VertEdge_BC(dqy_F, dqx_F, dqy_G, K12, K22)
    real(kind=dp), dimension(:,:), intent(inout) :: dqy_F
    real(kind=dp), dimension(:,:), intent(in) :: dqx_F, dqy_G, K12, K22
    
    integer :: i, m, n
    
    m = size(dqy_G, 1)
    n = size(dqx_F, 2)
    
    ! Bottom & Top, Near Boundary; Use BC
    do i = 2, m
      dqy_F(i, 1) = 0.5_dp*( &
        - K12(i,1)/K22(i,1)*(1.5_dp*dqx_F(i,1) - 0.5_dp*dqx_F(i,2)) &
        + 0.5_dp* dqy_G(i,2) + 0.5_dp* dqy_G(i-1,2) )
      dqy_F(i, n) = 0.5_dp*( &
        - K12(i,n+1)/K22(i,n+1)*(1.5_dp*dqx_F(i,n) - 0.5_dp*dqx_F(i,n-1)) &
        + 0.5_dp* dqy_G(i,n) + 0.5_dp* dqy_G(i-1,n) )
    end do
    
  end subroutine dy_field_VertEdge_BC
  
  subroutine testcase_fields(psi, K11, K22, K12, q)
    type(field), intent(inout) :: psi, K11, K22, K12, q
    integer :: i, j, m, n
    
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    real(kind = dp) :: x, y
    
    m = q%m
    n = q%n

    ! test case = 2-Period q
    ! psi
    if (psi%type_id .ne. 1) &
      call abort_handle("psi not defined at corners", __FILE__, __LINE__)
    
    do j = 1, size(psi%data, 2)
      do i = 1, size(psi%data, 1)
        ! Analytic formula: Defined on [0,1]^2 -> Need a prefactor to scale to [0,m] x [0,n]
        x = fld_x(i, m, psi%type_id)  ! in [0, 1]
        y = fld_x(j, n, psi%type_id)  ! in [0, 1]
        psi%data(i, j) = m*n/((2.0_dp*Pi)**2) *dsin(2.0_dp*Pi*x)*dsin(Pi*y)
        
        if (dabs(psi%data(i, j)) .le. 1D-15) psi%data(i, j) = 0.0_dp
      end do
    end do

    ! K11, K22, K12
    if (K11%type_id .ne. 1) &
      call abort_handle("K11 not defined at corners", __FILE__, __LINE__)
    
    do j = 1, size(K11%data, 2)
      do i = 1, size(K11%data, 1)
        ! Analytic formula: Defined on [0,1]^2 -> Need a prefactor to scale to [0,m] x [0,n]
        x = fld_x(i, m, K11%type_id)
        y = fld_x(j, n, K11%type_id)
        K11%data(i, j) = m*m/(1000.0_dp)* 4.0_dp*dexp(-0.5_dp*((x-0.5_dp)**2+(y-0.75_dp)**2))
        K22%data(i, j) = n*n/(1000.0_dp)* 2.0_dp*dexp(-0.5_dp*((x-0.5_dp)**2+(y-0.75_dp)**2))
        K12%data(i, j) = m*n/(1000.0_dp)* dexp(-0.5_dp*((x-0.5_dp)**2+(y-0.75_dp)**2))
        
        if (dabs(K11%data(i, j)) .le. 1D-15) K11%data(i, j) = 0.0_dp
        if (dabs(K22%data(i, j)) .le. 1D-15) K22%data(i, j) = 0.0_dp
        if (dabs(K12%data(i, j)) .le. 1D-15) K12%data(i, j) = 0.0_dp
      end do
    end do
    
    ! q0
    call assign_q_periodic(q, 0.0_dp)
  end subroutine testcase_fields
  
  subroutine assign_q_periodic(q, t)
    real(kind=dp), intent(in) :: t
    type(field), intent(inout) :: q
    
    integer :: i, j
    
    real(kind = dp) :: x, y
    
    do j = 1, q%n
      do i = 1, q%m
        x = fld_x(i, q%m, q%type_id)
        y = fld_x(j, q%n, q%type_id)
        q%data(i, j) = q_periodic(x,y,t)
      end do
    end do
    
  end subroutine assign_q_periodic
  
  pure real(kind=dp) function q_periodic(x, y, t)
    ! Define on [0, 1] times [0, 1]
    real(kind = dp), intent(in) :: x, y, t
    real(kind=dp), parameter :: Pi =  4.0_dp* datan (1.0_dp)
    
    q_periodic = 1.0_dp + dcos(Pi*t) * Pi**2 * (4.0_dp*x*(1.0_dp-y))**2*(4.0_dp*y*(1.0_dp-x))**3
    
  end function q_periodic
  
    pure real(kind=dp) function S_rhs_nonstat(x, y, t)
    ! Define on [0, m] times [0, n]
    real(kind = dp), intent(in) :: x, y, t
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    
    real(kind = dp) :: dqt, udqx, vdqy, divKgradq
    
    dqt = -1024.0_dp*Pi**3*x**2*(1.0_dp-y)**2*y**3*(1.0_dp-x)**3*dsin(Pi*t)
    udqx = 1280.0_dp*dcos(Pi*y)*Pi*y**3*dsin(2.0_dp*Pi*x)*(-2.0_dp/5.0_dp+x)*(-1.0_dp+y)**2*dcos(Pi*t)*x*(-1.0_dp+x)**2
    vdqy = -512.0_dp*Pi*dcos(2.0_dp*Pi*x)*dsin(Pi*y)*x**2*y**2*(-1.0_dp+x)**3*dcos(Pi*t)*(5.0_dp*y**2-8.0_dp*y+3)
    divKgradq = 0.01_dp* ( (512.0_dp*(-1.0_dp+x))*exp(-(0.5_dp)*x**2+(0.5_dp)*x-13.0_dp/32.0_dp-(0.5_dp)*y**2+0.75_dp*y)*Pi**2 &
                * cos(Pi*t)*y*((y**3-8.0_dp/5.0_dp*(y**2)+3.0_dp/5.0_dp*y)*x**5 &
                + (-12.0_dp/5.0_dp+6.0_dp*y**4-76.0_dp/5.0_dp*(y**3)+18.0_dp/5.0_dp*(y**2)+36.0_dp/5.0_dp*y)*x**4 &
                + (24.0_dp/5.0_dp+y**5-287.0_dp/20.0_dp*(y**4)+191.0_dp/10.0_dp*(y**3)+53.0_dp/4.0_dp*(y**2)-111.0_dp/5.0_dp*y)*x**3 &
                +(-7.0_dp/5.0_dp*(y**5)-23.0_dp/4.0_dp*(y**4)+57.0_dp/2.0_dp*(y**3)-731.0_dp/20.0_dp*(y**2)+84.0_dp/5.0_dp*y-12.0_dp/5.0_dp)*x**2 &
                +(0.2_dp)*(2.0_dp*(-1.0_dp+y))*(y**3+113.0_dp/4.0_dp*(y**2)-157.0_dp/4.0_dp*y+6.0_dp)*y*x-8.0_dp*y**2*(-1.0_dp+y)**2*(0.2_dp)) )
    
    S_rhs_nonstat = dqt + udqx + vdqy - divKgradq
    
  end function S_rhs_nonstat
  
  pure real(kind=dp) function fld_x(i, m, type_id)
    integer, intent(in) :: i, m, type_id
    
    if (type_id .eq. 1) then
      ! Defined at corners
      ! m = number of cells, not number of corners
      fld_x = (real(i,kind=dp)-1.0_dp)/real(m,kind=dp)
    else
      ! Defined at centres
      fld_x = (real(i,kind=dp)-0.5_dp)/real(m,kind=dp)
    end if
  end function fld_x
  
  pure real(kind=dp) function int_field(fld)
    type(field), intent(in) :: fld
    integer :: i, j

    int_field = 0.0_dp

    do j = 1, fld%n
      do i = 1, fld%m
        int_field = int_field + fld%data(i,j)
      end do
    end do

  end function int_field
  
  pure real(kind=dp) function int_sq_field(fld)
    type(field), intent(in) :: fld
    integer :: i, j


    int_sq_field = 0.0_dp

    do j = 1, fld%n
      do i = 1, fld%m
        int_sq_field = int_sq_field + fld%data(i,j) ** 2
      end do
    end do

  end function int_sq_field

  subroutine indicator_field(fld, i, j)
    type(field), intent(inout) :: fld
    integer, intent(in) :: i, j

    fld%data(i,j) = 1.0_dp

  end subroutine indicator_field

  pure integer function count_neg(fld)
    type(field), intent(in) :: fld
    integer :: i, j

    count_neg = 0

    do j = 1, fld%n
      do i = 1, fld%m
        if (fld%data(i, j) .lt. 0.0_dp) count_neg = count_neg + 1
      end do
    end do

  end function count_neg

  pure real(kind=dp) function max_neg(fld)
    type(field), intent(in) :: fld
    integer :: i, j

    max_neg = 0.0_dp

    do j = 1, fld%n
      do i = 1, fld%m
        if (fld%data(i,j) .lt. 0.0_dp) then
          max_neg = fld%data(i,j)
          !exit
        end if
      end do
    end do

    do j = 1, fld%n
      do i = 1, fld%m
        if ((fld%data(i,j) .lt. 0.0_dp) .and. (fld%data(i,j) >= max_neg)) then
          max_neg = fld%data(i,j)
        end if
      end do
    end do

  end function max_neg

  pure real(kind=dp) function neg_int(fld)
    type(field), intent(in) :: fld
    integer :: i, j

    neg_int = 0.0_dp

    do j = 1, fld%n
      do i = 1, fld%m
        if (fld%data(i,j) .lt. 0.0_dp) neg_int = neg_int + fld%data(i,j)
      end do
    end do

  end function neg_int

  subroutine print_info_field(fld)
    type(field), intent(in) :: fld

    write(6, "(a)") ""
    write(6, "(a,a)") " Display info: ", trim(fld%name)
    write(6, "(a,i0)") "|  type_id: ", fld%type_id
    write(6, "(a,i0,a,i0)") "|  size(fld%data, 1), size(fld%data, 2): ", size(fld%data, 1), ", ", size(fld%data, 2)
    write(6, "(a,i0,a,i0,a,i0)") "|  m, n, dof: ", fld%m, ", ", fld%n, ", ", fld%m*fld%n
    write(6, "(a,i0,a,i0)") "|  Negative dof = ", count_neg(fld), " out of ", fld%m*fld%n
    write(6, "(a,"//dp_chr//",a,"//dp_chr//")") "|  Least -ve dof = ", max_neg(fld), "; -ve_int ", neg_int(fld)
    write(6, "(a,"//dp_chr//",a,"//dp_chr//")") "|  Max = ", maxval(fld%data), "; Min = ", minval(fld%data)
    write(6, "(a,"//dp_chr//")") "|  fld_int = ", int_field(fld)
    !write(6, "(a)") "---------------------------------- "
    write(6, "(a)") ""
  end subroutine print_info_field

  pure real(kind=dp) function diff_int(fld1, fld2)
    type(field), intent(in) :: fld1, fld2

    ! Test conservation of mass
    diff_int = dabs(int_field(fld1)-int_field(fld2))
  end function diff_int

  pure logical function iszero(fld)
    real(kind=dp), dimension(:,:), intent(in) :: fld
    integer :: i, j

    iszero = .true.
    do j = 1, size(fld, 2)
      do i = 1, size(fld, 1)
        if (abs(fld(i,j)) .gt. 1e-13) then
          iszero = .false.
        end  if
      end do
    end do

   end function iszero
   
  real(kind=dp) function dt_CFL(psi, CFL)
    type(field), intent(in) :: psi
    real(kind=dp), intent(in) :: CFL
    
    real(kind=dp) :: u_max, v_max, uv_max
    real(kind=dp), dimension(:,:), allocatable :: u, v
    
    allocate(u(psi%m+1,psi%n))
    allocate(v(psi%m,psi%n+1))
    
    call psi_to_uv(u, v, psi)
    
    u_max = maxval(dabs(u))
    v_max = maxval(dabs(v))
    uv_max = u_max + v_max
    
    deallocate(u)
    deallocate(v)
    
    ! CFL = u*dt/dx
    dt_CFL = CFL/uv_max
    
  end function dt_CFL

end module advdiff_field
