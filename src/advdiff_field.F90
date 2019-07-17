module advdiff_field
  use advdiff_precision
  use advdiff_debug

  implicit none

  private

  public :: field, field_ptr, allocate, deallocate, zeros, &
    & scale, set, addto, addto_product, &
    & dx_field, dy_field, manual_field, &
    & psi_to_uv, int_field, int_sq_field, &
    & eval_field, iszero, &
    & assemble_psi, assemble_K, &
    & psi_fn, K_fn, &
    & assemble_Ke

  public :: indicator_field
  public :: uv_signed, K_flux, T_operator
  public :: print_info, diff_int
  public :: testcase_fields, fld_x, assign_q_periodic
  public :: q_stat, S_rhs_nonstat, S_rhs_stat
  public :: fourier_intpl, cosine_intpl
  
  ! Note: field and field_ptr are different
  type field
    integer :: m, n
    real(kind = dp), dimension(:, :), pointer :: data  ! It is an array, not a pointer
    character(len = field_name_len) :: name
    integer :: glayer
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
  end type K_flux

  type T_operator
    real(kind = dp), dimension(4, 4) :: T
  end type T_operator

  interface allocate
    module procedure allocate_field, allocate_uv, &
                     allocate_Kflux, allocate_KfluxE, &
                     allocate_refined_field
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
  subroutine allocate_field(fld, m, n, name, glayer, type_id)
    type(field), intent(out) :: fld
    integer, intent(in) :: m
    integer, intent(in) :: n
    character(len = *), intent(in) :: name
    integer, optional, intent(in) :: glayer
    integer, optional, intent(in) :: type_id

    integer :: gl = 0
    integer :: tid = 0


    if (present(glayer)) then
      gl = glayer
    end if

    if (present(type_id)) then
      tid = type_id
    end if

    fld%m = m
    fld%n = n
    fld%glayer = gl
    fld%type_id = tid

    if (fld%type_id .eq. 0) then
      ! Defined at cell centre
      allocate(fld%data(fld%m+2*gl, fld%n+2*gl))
    elseif (fld%type_id .eq. 1) then
      ! Defined at corner
      allocate(fld%data(fld%m+1+2*gl, fld%n+1+2*gl))
    elseif (fld%type_id .eq. 2) then
      ! Defined at corner, without right and top boundaries for Fourier transform
      allocate(fld%data(fld%m+2*gl, fld%n+2*gl))
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
  
  ! Fourier interpolations
  subroutine fourier_intpl(out, in, Mx, My)
    use, intrinsic :: iso_c_binding 
    include 'fftw3.f03'

    real(kind=dp), dimension(:, :), intent(out) :: out
    real(kind=dp), dimension(:, :), intent(in) :: in
    integer, intent(in) :: Mx, My  ! Factor of scaling up
  
    type(C_PTR) :: plan
    complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: Zk, Zk_pad
    real(kind=dp), dimension(:, :), allocatable :: in_copy
    
    integer :: i, j
    
    if (size(out, 1) .ne. size(in, 1) * Mx) &
      call abort_handle("fourier_intpl: size mismatch", __FILE__, __LINE__)
    if (size(out, 2) .ne. size(in, 2) * My) &
      call abort_handle("fourier_intpl: size mismatch", __FILE__, __LINE__)
    
    allocate(in_copy(size(in,1), size(in,2)))
    allocate(Zk(size(in,1)/2+1, size(in,2)))

    in_copy = in
    
    ! Forward FFT
    plan = fftw_plan_dft_r2c_2d(size(in_copy,2), size(in_copy,1), in_copy, Zk, FFTW_ESTIMATE)
    call fftw_execute_dft_r2c(plan, in_copy, Zk)
    call fftw_destroy_plan(plan)

    allocate(Zk_pad(size(out,1)/2+1, size(out,2)))
    Zk_pad = 0.0_dp
    
    !! Algorithm: https://web.eecs.umich.edu/~fessler/course/451/l/pdf/c5.pdf
    !! FFT-based approach: zero-pad in the frequency domain
    
    ! Zero padding along row (x-direction)
    ! Note r2c reprentation: last row = N/2+1 term
    do j = 1, size(Zk, 2)
      do i = 1, size(Zk, 1)
        Zk_pad(i, j) = Zk(i,j)
      end do
    end do
    
    ! Halve the N/2+1 term
    i = size(Zk, 1)
    do j = 1, size(Zk, 2)
      Zk_pad(i, j) = 0.5_dp*Zk_pad(i, j)
    end do
    
    ! Zero padding along column (y-direction)
    ! Shift along column
    do j = (size(Zk, 2)/2+2), size(Zk, 2)
      do i = 1, size(Zk_pad, 1)-1
        Zk_pad(i, j+(My-1)*size(Zk, 2)) = Zk_pad(i, j)
        Zk_pad(i, j) = 0.0_dp
      end do
    end do

    ! Halve the N/2+1 term
    j = size(Zk, 2)/2 + 1
    do i = 1, size(Zk, 1)
      Zk_pad(i, j+(My-1)*size(Zk, 2)) = 0.5_dp*Zk_pad(i, j)
      Zk_pad(i, j) = 0.5_dp*Zk_pad(i, j)
    end do
    
   ! Inverse FFT
    plan = fftw_plan_dft_c2r_2d(size(out,2),size(out,1), Zk_pad, out, FFTW_ESTIMATE)
    call fftw_execute_dft_c2r(plan, Zk_pad, out)
    call fftw_destroy_plan(plan)

    out = out/real(size(in,1)*size(in,2),kind=dp)
    
    deallocate(Zk)
    deallocate(Zk_pad)
    deallocate(in_copy)
    
  end subroutine fourier_intpl

  subroutine cosine_intpl(out, in, Mx, My)
    use, intrinsic :: iso_c_binding 
    include 'fftw3.f03'

    real(kind=dp), dimension(:, :), intent(out) :: out
    real(kind=dp), dimension(:, :), intent(in) :: in
    integer, intent(in) :: Mx, My  ! Factor of scaling up
  
    type(C_PTR) :: plan
    real(kind=dp), dimension(:, :), allocatable :: Bpq, Bpq_pad
    real(kind=dp), dimension(:, :), allocatable :: in_copy
    
    integer :: i, j
    
    if (size(out, 1) .ne. size(in, 1) * Mx) &
      call abort_handle("cosine_intpl: size mismatch", __FILE__, __LINE__)
    if (size(out, 2) .ne. size(in, 2) * My) &
      call abort_handle("cosine_intpl: size mismatch", __FILE__, __LINE__)
    
    allocate(Bpq(size(in,1), size(in,2)))
    allocate(in_copy(size(in,1), size(in,2)))
    allocate(Bpq_pad(size(out,1), size(out,2)))

    in_copy = in
    
    ! Forward DCT
    plan = fftw_plan_r2r_2d(size(in,2),size(in,1), in_copy, Bpq,  &
                            FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE)
    call fftw_execute_r2r(plan, in_copy, Bpq)
    call fftw_destroy_plan(plan)

    ! Zero padding: only need to pad zeros beyond Bpq
    !! Relationship with matlab dct2
    !! http://www.voidcn.com/article/p-pduypgxe-du.html
    Bpq_pad = 0.0_dp
    do j = 1, size(Bpq, 2)
      do i = 1, size(Bpq, 1)
        Bpq_pad(i, j) = Bpq(i, j)
      end do
    end do
    
    plan = fftw_plan_r2r_2d(size(out,2),size(out,1), Bpq_pad, &
                          out, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
    call fftw_execute_r2r(plan, Bpq_pad, out)
    call fftw_destroy_plan(plan)

    out = out/(4*size(in,1)*size(in,2))
    
    deallocate(Bpq)
    deallocate(Bpq_pad)
    deallocate(in_copy)
    
  end subroutine cosine_intpl
  
  
  subroutine allocate_refined_field(rfld, fld, i_reflv, j_reflv)
    type(field), intent(out) :: rfld
    type(field), intent(in) :: fld
    integer, intent(in) :: i_reflv, j_reflv

    call allocate(rfld, fld%m*i_reflv, fld%n*j_reflv, "r"//trim(fld%name), fld%glayer, fld%type_id)

    if (fld%type_id .eq. 0) then
      ! Cosine interpolations
      call cosine_intpl(rfld%data, fld%data, i_reflv, j_reflv)
      
    elseif (fld%type_id .eq. 1) then
      ! Fourier interpolations
      call fourier_intpl(rfld%data(1:fld%m*i_reflv, 1:fld%n*j_reflv), &
                  fld%data(1:fld%m, 1:fld%n), i_reflv, j_reflv)
    end if

  end subroutine allocate_refined_field

  subroutine allocate_uv(uv_fld, psi)
    type(uv_signed), intent(out) :: uv_fld
    type(field), intent(in) :: psi

    real(kind = dp), dimension(:, :), pointer :: lup, lum, lvp, lvm
    real(kind = dp), dimension(:, :), pointer :: lu_abs, lv_abs
    integer, dimension(:, :), pointer :: lu_sgn, lv_sgn
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


    ! Assign values to u and v
    lup => uv_fld%up
    lum => uv_fld%um
    lvp => uv_fld%vp
    lvm => uv_fld%vm

    lu_sgn => uv_fld%u_sgn
    lv_sgn => uv_fld%v_sgn
    lu_abs => uv_fld%u_abs
    lv_abs => uv_fld%v_abs

    
    ! Temporary values
    call psi_to_uv(lup, lvp, psi)

    ! maxmin u and v wrt 0
    do j = 1, size(lup, 2)
      do i = 1, size(lup, 1)
        if (dabs(lup(i,j)) > 1e-14) then
          lu_sgn(i,j) = int(sign(1.0_dp, lup(i,j)))
        else
          lu_sgn(i,j) = 0
        end if

        lum(i,j) = min(lup(i, j), 0.0_dp)
        lup(i,j) = max(lup(i, j), 0.0_dp)
        lu_abs(i,j) = lup(i,j) - lum(i,j)
      end do
    end do

    do j = 1, size(lvp, 2)
      do i = 1, size(lvp, 1)
        if (dabs(lvp(i,j)) > 1e-14) then
          lv_sgn(i,j) = int(sign(1.0_dp, lvp(i,j)))
        else
          lv_sgn(i,j) = 0
        end if

        lvm(i,j) = min(lvp(i, j), 0.0_dp)
        lvp(i,j) = max(lvp(i, j), 0.0_dp)
        lv_abs(i,j) = lvp(i,j) - lvm(i,j)
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
    n = K22%n

    allocate(K_e%K11e(m+1, n))  ! Define at vertical edges
    allocate(K_e%K12e(m+1, n))  ! Define at vertical edges
    allocate(K_e%K21e(m, n+1))  ! Define at horizonal edges
    allocate(K_e%K22e(m, n+1))  ! Define at horizonal edges

    K11e => K_e%K11e
    K12e => K_e%K12e
    K21e => K_e%K21e
    K22e => K_e%K22e

    ! K: cell-centred
    if ( (K11%type_id .eq. 0) .and. &
         (K22%type_id .eq. 0) .and. (K12%type_id .eq. 0) ) then
      ! For K11 and K12
      do j = 1, n
        do i = 2, m
          K11e(i,j) = 0.5_dp*(K11%data(i,j) + K11%data(i-1,j))
          K12e(i,j) = 0.5_dp*(K12%data(i,j) + K12%data(i-1,j))
        end do
      end do

      ! Set diffusivity at domain boundaries at the nearest value
      ! (To be replaced by BC later)
      do j = 1, n
        K11e(1,j) = K11%data(1,j)
        K12e(1,j) = K12%data(1,j)
        K11e(m+1,j) = K11%data(m,j)
        K12e(m+1,j) = K12%data(m,j)
      end do

      ! For K21 and K22
      do j = 2, n
        do i = 1, m
          K21e(i,j) = 0.5_dp*(K12%data(i,j) + K12%data(i,j-1))
          K22e(i,j) = 0.5_dp*(K22%data(i,j) + K22%data(i,j-1))
        end do
      end do

      ! Set diffusivity at domain boundaries at the nearest value
      ! (To be replaced by BC later)
      do i = 1, m
        K21e(i,1) = K12%data(i,1)
        K22e(i,1) = K22%data(i,1)
        K21e(i,n+1) = K12%data(i,n)
        K22e(i,n+1) = K22%data(i,n)
      end do

    ! K: at corners
    elseif ( (K11%type_id .eq. 1) .and. &
             (K22%type_id .eq. 1) .and. (K12%type_id .eq. 1) ) then
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
    else
      write(6, *) "Wrong type_id in allocate_Kflux"
      stop 1
    end if
  end subroutine allocate_Kflux

  subroutine allocate_KfluxE(K_e, K11e, K12e, K21e, K22e)
    type(K_flux), intent(inout) :: K_e
    type(field), intent(in) :: K11e, K12e, K21e, K22e

    allocate(K_e%K11e(K11e%m, K11e%n))  ! Define at vertical edges
    allocate(K_e%K12e(K12e%m, K12e%n))  ! Define at vertical edges
    allocate(K_e%K21e(K21e%m, K21e%n))  ! Define at horizonal edges
    allocate(K_e%K22e(K22e%m, K22e%n))  ! Define at horizonal edges

    K_e%K11e = K11e%data
    K_e%K12e = K12e%data
    K_e%K21e = K21e%data
    K_e%K22e = K22e%data

  end subroutine allocate_KfluxE

  subroutine deallocate_Kflux(K_e)
    type(K_flux), intent(inout) :: K_e

    deallocate(K_e%K11e)
    deallocate(K_e%K21e)
    deallocate(K_e%K12e)
    deallocate(K_e%K22e)

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

    real(kind = dp), dimension(:, :), pointer :: lpsi_fld

    lpsi_fld => psi_fld%data

    if (psi_fld%type_id .ne. 1) &
      call abort_handle("psi_to_uv: psi is not at corners", __FILE__, __LINE__)

    do j = 1, psi_fld%n
      do i = 1, psi_fld%m+1
        u(i, j) = lpsi_fld(i, j) - lpsi_fld(i, j+1)
        
        if (dabs(u(i, j)) < 1e-14) then
          u(i, j) = 0.0_dp
        end if
      end do
    end do

    do j = 1, psi_fld%n+1
      do i = 1, psi_fld%m
        v(i, j) = lpsi_fld(i+1, j) - lpsi_fld(i, j)
        
        if (dabs(v(i, j)) < 1e-14) then
          v(i, j) = 0.0_dp
        end if
      end do
    end do

  end subroutine psi_to_uv

  subroutine dx_field(d_fld, fld)
    real(kind=dp), dimension(:,:), intent(out) :: d_fld
    type(field), intent(in) :: fld
    integer :: i, j, gl

    real(kind = dp), dimension(:, :), pointer :: lfld

    ! Exclude ghost cell in dx_field and dy_field
    lfld => fld%data
    gl = fld%glayer

!     ! Dimension check, assumed ghost layer exist
!     if ((size(d_fld,1) .ne. fld%m+1) .or. (size(d_fld,2) .ne. fld%n)) then
!       call abort_handle("dx_field: Invalid size", __FILE__, __LINE__)
!     end if
! 
!     if (gl .le. 0) then
!       call abort_handle("dx_field: Need to impose ghost layer", __FILE__, __LINE__)
!     end if

    do j = 1, fld%n
      do i = 1, (fld%m+1)
        d_fld(i, j) = lfld(i+gl, j+gl) - lfld(i+gl-1, j+gl)
      end do
    end do

  end subroutine dx_field

  subroutine dy_field(d_fld, fld)
    real(kind=dp), dimension(:,:), intent(out) :: d_fld
    type(field), intent(in) :: fld
    integer :: i, j, gl

    real(kind = dp), dimension(:, :), pointer :: lfld

    lfld => fld%data
    gl = fld%glayer

!     ! Dimension check
!     if ((size(d_fld,1) .ne. fld%m) .or. (size(d_fld,2) .ne. fld%n+1)) then
!       call abort_handle("dy_field: Invalid size", __FILE__, __LINE__)
!     end if
! 
!     if (gl .le. 0) then
!       call abort_handle("dy_field: Need to impose ghost layer", __FILE__, __LINE__)
!     end if

    do j = 1, (fld%n+1)
      do i = 1, fld%m
        d_fld(i, j) = lfld(i+gl, j+gl) - lfld(i+gl, j+gl-1)
      end do
    end do

  end subroutine dy_field

  subroutine manual_field(fld, param)
    type(field), intent(out) :: fld
    integer, intent(in) :: param
    integer :: i, j
    real(kind = dp) :: im, jm, x, y

    real(kind = dp), dimension(:, :), pointer :: lfld
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)

    lfld => fld%data

    do j = 1, size(fld%data,2)
      do i = 1, size(fld%data,1)
        if (param .eq. 1) then
          lfld(i, j) = 2.0_dp * real(i)
        elseif (param .eq. 2) then
          lfld(i, j) = real(j)*real(i)
          !lfld(i, j) = (fld%n-real(j))*(fld%m-real(i))
        elseif (param .eq. 3) then
          lfld(i, j) = (1.0_dp/(2.0_dp*Pi*(real(fld%m)/16)**2)) &
                      * dexp( ((real(i)-0.5_dp*real(fld%m))**  2+(real(j)-0.5_dp*real(fld%n))**2) &
                      /(-2.0_dp*((real(fld%m)/16)**2)) )
        elseif (param .eq. 4) then
          if ((i .eq. fld%m/2) .and. (j .eq. fld%n/2)) then
            lfld(i, j) = 1.0_dp
          else
            lfld(i, j) = 0.0_dp
          end if
        elseif (param .eq. 5) then
            lfld(i, j) = 2.0_dp *i + 3.0_dp *j - 1.0* i*j
            !lfld(i, j) = 2.0_dp *i
        elseif (param .eq. 6) then
            ! Ridge body rotation for psi
            im = (fld%m+1+fld%type_id)*0.5_dp
            jm = (fld%n+1+fld%type_id)*0.5_dp
            lfld(i, j) = 0.5_dp*( (i-im)*(i-im) + (j-jm)*(j-jm))
        elseif (param .eq. 7) then
            ! TG
            x = fld_x(i, fld%m, fld%type_id)  ! in [0, 1]
            y = fld_x(j, fld%n, fld%type_id)  ! in [0, 1]
            fld%data(i, j) = fld%m*fld%n/((2.0_dp*Pi)**2) *dsin(2.0_dp*Pi*x)*dsin(Pi*y)
        end if
      end do
    end do

  end subroutine manual_field
  
  subroutine testcase_fields(psi, K11, K22, K12, q, test_id)
    type(field), intent(inout) :: psi, K11, K22, K12, q
    integer, intent(in) :: test_id
    integer :: i, j, m, n
    
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    real(kind = dp) :: x, y
    integer :: gl

    m = q%m
    n = q%n

    ! test_id = 1: 2-Period q
    ! test_id = 2: steady state    
    if ((test_id .eq. 1) .or. (test_id .eq. 2)) then     
      ! psi
      gl = psi%glayer
      
      if (psi%type_id .ne. 1) &
        call abort_handle("psi not defined at corners", __FILE__, __LINE__)
      
      do j = 1, psi%n + 1
        do i = 1, psi%m + 1
          ! Analytic formula: Defined on [0,1]^2 -> Need a prefactor to scale to [0,m] x [0,n]
          x = fld_x(i, m, psi%type_id)  ! in [0, 1]
          y = fld_x(j, n, psi%type_id)  ! in [0, 1]
          psi%data(i+gl, j+gl) = m*n/((2.0_dp*Pi)**2) *dsin(2.0_dp*Pi*x)*dsin(Pi*y)
          
          ! Ad-hoc
          !psi%data(i+gl, j+gl) = 0.01_dp*psi%data(i+gl, j+gl)

        end do
      end do

      ! K11, K22, K12
      gl = K11%glayer
      
      do j = 1, K11%n
        do i = 1, K11%m
          ! Analytic formula: Defined on [0,1]^2 -> Need a prefactor to scale to [0,m] x [0,n]
          x = fld_x(i, m, K11%type_id)
          y = fld_x(j, n, K11%type_id)
          K11%data(i+gl, j+gl) = m*m/((16.0*Pi)**2)* 4.0_dp*dsin(Pi*x)*dsin(Pi*y)
          K22%data(i+gl, j+gl) = n*n/((16.0*Pi)**2)* 2.0_dp*dsin(Pi*x)*dsin(Pi*y)
          K12%data(i+gl, j+gl) = m*n/((16.0*Pi)**2)* dsin(Pi*x)*dsin(Pi*y)
        end do
      end do
    end if
    
    
    if (test_id .eq. 1) then
      ! q0
      call assign_q_periodic(q, 0.0_dp)
      
    elseif (test_id .eq. 2) then
      ! Steady state solution
      gl = q%glayer
      
      do j = 1, q%n
        do i = 1, q%m
          x = fld_x(i, m, q%type_id)
          y = fld_x(j, n, q%type_id)
          q%data(i+gl, j+gl) = q_stat(x,y)
        end do
      end do
    end if

  end subroutine testcase_fields
  
  subroutine assign_q_periodic(q, t)
    real(kind=dp), intent(in) :: t
    type(field), intent(inout) :: q
    
    integer :: i, j
    
    real(kind = dp) :: x, y
    integer :: gl
    
    ! q0
    gl = q%glayer
      
    do j = 1, q%n
      do i = 1, q%m
        x = fld_x(i, q%m, q%type_id)
        y = fld_x(j, q%n, q%type_id)
        q%data(i+gl, j+gl) = q_periodic(x,y,t)
      end do
    end do
    
  end subroutine assign_q_periodic
  
  pure real(kind=dp) function q_periodic(x, y, t)
    ! Define on [0, 1] times [0, 1]
    real(kind = dp), intent(in) :: x, y, t
    real(kind=dp), parameter :: Pi =  4.0_dp* datan (1.0_dp)
    
    q_periodic = 1.0_dp + dcos(Pi*t) * Pi**2 * (4.0_dp*x*(1.0_dp-y))**2*(4.0_dp*y*(1.0_dp-x))**3
    
  end function q_periodic
  
  pure real(kind=dp) function q_stat(x, y)
    ! Define on [0, 1] times [0, 1]
    real(kind = dp), intent(in) :: x, y
    real(kind=dp), parameter :: Pi =  4.0_dp* datan (1.0_dp)
    
    q_stat = (x*(1.0_dp-y))**2 * (y*(1.0_dp-x))**3
    !q_stat = (dsin(2.0_dp*Pi*x)*dsin(Pi*y))**2
    !q_stat = (dsin(2.0_dp*Pi*x)*dsin(Pi*y))**3 - dsin(2.0_dp*Pi*x)*dsin(Pi*y)
  end function q_stat
  
  pure real(kind=dp) function S_rhs_nonstat(x, y, t)
    ! Define on [0, m] times [0, n]
    real(kind = dp), intent(in) :: x, y, t
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    
    real(kind = dp) :: dqt, udqx, vdqy, divKgradq
    
    dqt = -1024.0_dp*Pi**3*x**2*(1.0_dp-y)**2*y**3*(1.0_dp-x)**3*dsin(Pi*t)
    udqx = 1280.0_dp*dcos(Pi*y)*Pi*y**3*dsin(2.0_dp*Pi*x)*(-2.0_dp/5.0_dp+x)*(-1.0_dp+y)**2*dcos(Pi*t)*x*(-1.0_dp+x)**2
    vdqy = -512.0_dp*Pi*dcos(2.0_dp*Pi*x)*dsin(Pi*y)*x**2*y**2*(-1.0_dp+x)**3*dcos(Pi*t)*(5.0_dp*y**2-8.0_dp*y+3)
    divKgradq = dcos(Pi*x)*dsin(Pi*y)*(2048.0_dp*Pi**2*x*(1.0_dp-y)**2*y**3*(1.0_dp-x)**3*dcos(Pi*t) &
                -3072.0_dp*Pi**2*x**2*(1.0_dp-y)**2*y**3*(1.0_dp-x)**2*dcos(Pi*t))/(64.0_dp*Pi) &
                +dsin(Pi*x)*dsin(Pi*y)*(2048.0_dp*Pi**2*(1.0_dp-y)**2*y**3*(1.0_dp-x)**3*dcos(Pi*t) &
                -12288.0_dp*Pi**2*x*(1.0_dp-y)**2*y**3*(1.0_dp-x)**2*dcos(Pi*t) &
                +6144.0_dp*Pi**2*x**2*(1.0_dp-y)**2*y**3*(1.0_dp-x)*dcos(Pi*t))/(64.0_dp*Pi**2) &
                +dcos(Pi*x)*dsin(Pi*y)*(-2048.0_dp*Pi**2*x**2*(1.0_dp-y)*y**3*(1.0_dp-x)**3*dcos(Pi*t) &
                +3072.0_dp*Pi**2*x**2*(1.0_dp-y)**2*y**2*(1.0_dp-x)**3*dcos(Pi*t))/(256.0_dp*Pi) &
                +dsin(Pi*x)*dsin(Pi*y)*(-4096.0_dp*Pi**2*x*(1.0_dp-y)*y**3*(1.0_dp-x)**3*dcos(Pi*t) &
                +6144.0_dp*Pi**2*x**2*(1.0_dp-y)*y**3*(1.0_dp-x)**2*dcos(Pi*t) &
                +6144.0_dp*Pi**2*x*(1.0_dp-y)**2*y**2*(1.0_dp-x)**3*dcos(Pi*t) &
                -9216.0_dp*Pi**2*x**2*(1.0_dp-y)**2*y**2*(1.0_dp-x)**2*dcos(Pi*t))/(128.0_dp*Pi**2) &
                +dsin(Pi*x)*dcos(Pi*y)*(2048.0_dp*Pi**2*x*(1.0_dp-y)**2*y**3*(1.0_dp-x)**3*dcos(Pi*t) &
                -3072.0_dp*Pi**2*x**2*(1.0_dp-y)**2*y**3*(1.0_dp-x)**2*dcos(Pi*t))/(256.0_dp*Pi) &
                +dsin(Pi*x)*dcos(Pi*y)*(-2048.0_dp*Pi**2*x**2*(1.0_dp-y)*y**3*(1.0_dp-x)**3*dcos(Pi*t) &
                +3072.0_dp*Pi**2*x**2*(1.0_dp-y)**2*y**2*(1.0_dp-x)**3*dcos(Pi*t))/(128.0_dp*Pi) &
                +dsin(Pi*x)*dsin(Pi*y)*(2048.0_dp*Pi**2*x**2*y**3*(1.0_dp-x)**3*dcos(Pi*t) &
                -12288.0_dp*Pi**2*x**2*(1.0_dp-y)*y**2*(1.0_dp-x)**3*dcos(Pi*t) &
                +6144.0_dp*Pi**2*x**2*(1.0_dp-y)**2*y*(1.0_dp-x)**3*dcos(Pi*t))/(128.0_dp*Pi**2)
    
    S_rhs_nonstat = dqt + udqx + vdqy - divKgradq
    !S_rhs_nonstat = dqt - divKgradq
    !S_rhs_nonstat = dqt + udqx + vdqy
    
  end function S_rhs_nonstat
  
  pure real(kind=dp) function S_rhs_stat(x, y, qxy, lambda)
    ! Define on [0, 1] times [0, 1]
    real(kind = dp), intent(in) :: x, y, qxy, lambda
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    
    real(kind = dp) :: dqt, udqx, vdqy, divKgradq
    
    ! Explicit source term
    dqt = 0.0_dp
    udqx = (5.0_dp*x-2.0_dp)*dcos(Pi*y)*(y-1.0_dp)**2*y**3*dsin(2.0_dp*Pi*x)*(x-1.0_dp)**2*x/(4.0_dp*Pi)
    vdqy = -(5.0_dp*y-3.0_dp)*dsin(Pi*y)*(y-1.0_dp)*y**2*dcos(2.0_dp*Pi*x)*(x-1.0_dp)**3.0_dp*x**2/(2.0_dp*Pi)
    divKgradq = (5.0_dp*y*(x-1.0_dp)*((((8.0_dp*y**2-48.0_dp/5.0_dp*y+12.0_dp/5.0_dp)*x**4 &
                +(10.0_dp*y**3-32.0_dp*y**2+126.0_dp/5.0_dp*y-24.0_dp/5.0_dp)*x**3+(16.0_dp*y**4-46.0_dp*y**3 &
                +232.0_dp/5.0_dp*(y**2)-18.0_dp*y+12.0_dp/5.0_dp)*x**2+(-96.0_dp/5.0_dp*(y**2)+12.0_dp/5.0_dp*y &
                -64.0_dp/5.0_dp*(y**4)+148.0_dp/5.0_dp*(y**3))*x+8.0_dp*y**2*(y-1.0_dp)**2*(1.0_dp/5.0_dp))*dsin(Pi*y) &
                +(2.0_dp*(y-1.0_dp))*y*dcos(Pi*y)*(x-1.0_dp)*((y-3.0_dp/5.0_dp)*x**2+((1.0_dp/2.0_dp)*y**2 &
                -3.0_dp/2.0_dp*y+3.0_dp/5.0_dp)*x-(1.0_dp/5.0_dp)*y**2+(1.0_dp/5.0_dp)*y)*x*Pi)*dsin(Pi*x) &
                +dcos(Pi*x)*(y-1.0_dp)*y*dsin(Pi*y)*((y-3.0_dp/5.0_dp)*x**2+(4.0_dp*y**2-5.0_dp*y+3.0_dp/5.0_dp)*x &
                -8.0_dp*y**2*(1.0_dp/5.0_dp)+8.0_dp*y*(1.0_dp/5.0_dp))*(x-1.0_dp)*x*Pi))/(-256.0_dp*Pi*Pi)
    
    
!     S_rhs_stat = dqt + udqx + vdqy - divKgradq
    S_rhs_stat = lambda*(q_stat(x,y) - qxy) + dqt + udqx + vdqy - divKgradq
!     S_rhs_stat = lambda*(q_stat(x,y) - qxy) + dqt - divKgradq
!     S_rhs_stat = lambda*(q_stat(x,y) - qxy) + dqt + udqx + vdqy
  end function S_rhs_stat

  pure real(kind=dp) function fld_x(i, m, type_id)
    integer, intent(in) :: i, m, type_id
    
    if (type_id .eq. 1) then
      ! Defined at corners
      fld_x = (real(i,kind=dp)-1.0_dp)/real(m,kind=dp)
    else
      ! Defined at centres
      fld_x = (real(i,kind=dp)-0.5_dp)/real(m,kind=dp)
    end if
  end function fld_x
  
  pure real(kind=dp) function psi_fn(xl, yl, param)
    real(kind=dp), intent(in) :: xl, yl  ! in [0, 1]
    integer, intent(in) :: param

    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)

    psi_fn = 0.0_dp

    ! Taylor green
    if (param .eq. 1) psi_fn = dsin(Pi*xl) * dsin(Pi*2.0_dp*yl)
    ! Tilted Taylor green
    if (param .eq. 2) psi_fn = (xl + 2.0_dp*yl)*dsin(Pi*xl)*dsin(Pi*2.0_dp*yl)
  end function psi_fn

  subroutine assemble_psi(fld, param)
    type(field), intent(out) :: fld
    integer, intent(in) :: param
    integer :: i, j, m, n

    real(kind = dp), dimension(:, :), pointer :: lfld
    real(kind = dp) :: xl, yl

    lfld => fld%data

    ! Assumed m = n
    m = size(fld%data,1)-1  ! =number of cells
    n = size(fld%data,2)-1  ! =number of cells

    ! Functions defined on logical space [0, 1]
    ! xl: [0, 1] -> x: [1, nc+1], at edges
    do j = 1, size(fld%data,2)
      do i = 1, size(fld%data,1)
        xl = real(i-1)/m
        yl = real(j-1)/n
        fld%data(i, j) = psi_fn(xl, yl, param) * m  ! Rescale to computational domain [1, m+1]

        ! Assign 0.0
        if (dabs(fld%data(i, j)) < 1.0e-12) fld%data(i, j) = 0.0_dp
      end do
    end do
  end subroutine assemble_psi

  pure real(kind=dp) function K_fn(xl, yl, param)
    real(kind=dp), intent(in) :: xl, yl  ! in [0, 1]
    integer, intent(in) :: param

    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)

    if (param .eq. 0) K_fn = 0.0_dp
    if (param .eq. 1) K_fn = 1.0_dp
    if (param .eq. 2) K_fn = (dsin(2.0_dp*Pi*xl) * dsin(2.0_dp*Pi*yl))**2
    if (param .eq. 3) K_fn = (dsin(2.0_dp*Pi*xl) * dsin(2.0_dp*Pi*yl))**2 + dsin(Pi*xl) * dcos(Pi*yl)
    if (param .eq. 11) K_fn = 1.0_dp*xl + 2.0_dp*yl
    if (param .eq. 22) K_fn = 3.0_dp*xl - 4.0_dp*yl
    if (param .eq. 12) K_fn = 2.0_dp*xl + 5.0_dp*yl
  end function K_fn

  subroutine assemble_K(fld, param)
    type(field), intent(out) :: fld
    integer, intent(in) :: param
    integer :: i, j, m, n

    real(kind = dp), dimension(:, :), pointer :: lfld
    real(kind = dp) :: xl, yl

    lfld => fld%data

    m = fld%m ! =number of cells
    n = fld%n ! =number of cells

    ! Functions defined on logical space [0, 1]
    ! xl: [0, 1] -> x: [1, nc+1], at centre
    do j = 1, size(fld%data,2)
      do i = 1, size(fld%data,1)
        xl = (real(i) -0.5_dp)/m
        yl = (real(j) -0.5_dp)/n

        fld%data(i, j) = K_fn(xl, yl, param)

        ! Rescale to computational domain [1, m+1], assumed m = n
        !fld%data(i, j) = fld%data(i, j) * (m-1)
      end do
    end do
  end subroutine assemble_K

  subroutine assemble_Ke(lfld, param)
    real(kind=dp), dimension(:,:), intent(inout):: lfld
    integer, intent(in) :: param
    integer :: i, j, m, n

    real(kind = dp) :: xl, yl


    m = size(lfld,1)
    n = size(lfld,2)
    
    ! Functions defined on logical space [0, 1]
    ! xl: [0, 1] -> x: [1, nc+1], at edges
    if (m .eq. (n+1)) then
      do j = 1, n
        do i = 1, m
          xl = (real(i) -1.0_dp)/(m-1)
          yl = (real(j) -0.5_dp)/n

          lfld(i, j) = K_fn(xl, yl, param)

          ! Rescale to computational domain [1, m+1], assumed m = n
          !lfld(i, j) = lfld(i, j) * (m-1)
        end do
      end do
    elseif ((m+1) .eq. n) then
      do j = 1, n
        do i = 1, m
          xl = (real(i) -0.5_dp)/m
          yl = (real(j) -1.0_dp)/(n-1)

          lfld(i, j) = K_fn(xl, yl, param)

          ! Rescale to computational domain [1, m+1], assumed m = n
          !fld%data(i, j) = fld%data(i, j) * (m-1)
        end do
      end do
    end if
    
  end subroutine assemble_Ke
  
  pure real(kind=dp) function int_field(fld)
    type(field), intent(in) :: fld
    integer :: i, j, gl

    gl = fld%glayer

    int_field = 0.0_dp

    do j = (1+gl), (fld%n+gl)
      do i = (1+gl), (fld%m+gl)
        int_field = int_field + fld%data(i,j)
      end do
    end do

  end function int_field
  
  pure real(kind=dp) function int_sq_field(fld)
    type(field), intent(in) :: fld
    integer :: i, j, gl

    gl = fld%glayer

    int_sq_field = 0.0_dp

    do j = (1+gl), (fld%n+gl)
      do i = (1+gl), (fld%m+gl)
        int_sq_field = int_sq_field + fld%data(i,j) ** 2
      end do
    end do

  end function int_sq_field

  pure real(kind=dp) function eval_field(fld, i, j)
    type(field), intent(in) :: fld
    integer, intent(in) :: i, j

    integer :: gl

    gl = fld%glayer

    ! Piecewise constant field
    eval_field = fld%data(i+gl,j+gl)

  end function eval_field

  subroutine indicator_field(fld, i, j)
    type(field), intent(inout) :: fld
    integer, intent(in) :: i, j

    integer :: gl

    gl = fld%glayer

    fld%data(i+gl,j+gl) = 1.0_dp

  end subroutine indicator_field

  pure integer function count_neg(fld)
    type(field), intent(in) :: fld
    integer :: i, j, gl

    gl = fld%glayer
    count_neg = 0

    do j = 1, fld%n
      do i = 1, fld%m
        if (fld%data(i+gl, j+gl) < 0.0_dp) count_neg = count_neg + 1
      end do
    end do

  end function count_neg

  pure real(kind=dp) function max_neg(fld)
    type(field), intent(in) :: fld
    integer :: i, j, gl

    gl = fld%glayer
    max_neg = 0.0_dp

    do j = 1, fld%n
      do i = 1, fld%m
        if (fld%data(i+gl,j+gl) < 0.0_dp) then
          max_neg = fld%data(i+gl,j+gl)
          !exit
        end if
      end do
    end do

    do j = 1, fld%n
      do i = 1, fld%m
        if ((fld%data(i+gl,j+gl) < 0.0_dp) .and. (fld%data(i+gl,j+gl) >= max_neg)) then
          max_neg = fld%data(i+gl,j+gl)
        end if
      end do
    end do

  end function max_neg

  pure real(kind=dp) function neg_int(fld)
    type(field), intent(in) :: fld
    integer :: i, j, gl

    gl = fld%glayer
    neg_int = 0.0_dp

    do j = 1, fld%n
      do i = 1, fld%m
        if (fld%data(i+gl,j+gl) < 0.0_dp) neg_int = neg_int + fld%data(i+gl,j+gl)
      end do
    end do

  end function neg_int

  subroutine print_info_field(fld)
    type(field), intent(in) :: fld

    write(6, "(a)") ""
    write(6, "(a,a)") " Display info: ", trim(fld%name)
    write(6, "(a,i0,a,i0,a,i0)") "|  m, n, dof: ", fld%m, ", ", fld%n, ", ", fld%m*fld%n
    write(6, "(a,i0)") "|  ghost layer: ", fld%glayer
    write(6, "(a,i0,a,i0)") "|  Negative dof = ", count_neg(fld), " out of ", fld%m*fld%n
    write(6, "(a,"//dp_chr//",a,"//dp_chr//")") "|  Max -ve dof = ", max_neg(fld), "; -ve_int ", neg_int(fld)
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
        if (abs(fld(i,j)) > 1e-13) then
          iszero = .false.
        end  if
      end do
    end do

   end function iszero

end module advdiff_field
