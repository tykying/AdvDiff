#define FFTW3 1
module advdiff_field
  use advdiff_precision
  use advdiff_debug

  implicit none

  private

  public :: field, field_ptr, allocate, deallocate, zeros, &
    & scale, set, addto, addto_product, &
    & dx_field, dy_field, manual_field, &
    & psi_to_uv, int_field, int_sq_field, &
    & eval_field, iszero
    
  public :: indicator_field
  public :: uv_signed, K_flux
  public :: print_info, diff_int
  public :: testcase_fields, fld_x, assign_q_periodic
  public :: q_stat, S_rhs_nonstat, S_rhs_stat
  public :: fourier_intpl, cosine_intpl, sine_intpl, bilinear_intpl
  public :: neg_int
  public :: sine_filter
  public :: del2
  public :: dt_CFL
  public :: sine_basis
  
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
  subroutine allocate_field(fld, m, n, name, glayer, type_id)
    type(field), intent(out) :: fld
    integer, intent(in) :: m
    integer, intent(in) :: n
    character(len = *), intent(in) :: name
    integer, optional, intent(in) :: glayer
    integer, optional, intent(in) :: type_id

    integer :: gl = 0
    integer :: tid = 2  ! Default = 2: centre-based


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

    if (fld%type_id .eq. 2) then
      ! Defined at cell centre
      allocate(fld%data(fld%m+2*gl, fld%n+2*gl))
    elseif (fld%type_id .eq. 1) then
      ! Defined at corner
      allocate(fld%data(fld%m+1+2*gl, fld%n+1+2*gl))
    elseif (fld%type_id .lt. 0) then
      ! DOF type: no ghost layer
      if (gl .ne. 0) &
        call abort_handle("Ghost layer for DOF should be zero", __FILE__, __LINE__)
      allocate(fld%data(fld%m, fld%n))
    elseif (fld%type_id .eq. 0) then
       call abort_handle("fld%type_id = 0: should correct", __FILE__, __LINE__)
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
  
  subroutine sine_basis(inout, wavenumber, amplitude)
#if FFTW3 == 1
    ! inout: no boundaries => size(in) = [2^m+1-2, 2^n+1-2]
    use, intrinsic :: iso_c_binding 
    include 'fftw3.f03'
    
    real(kind=dp), dimension(:, :), intent(inout) :: inout
    integer :: wavenumber
    real(kind=dp), intent(in) :: amplitude
    
    type(C_PTR) :: plan
    real(kind=dp), dimension(:, :), allocatable :: Bpq
    
    integer :: i, j
    
    
    allocate(Bpq(size(inout,1), size(inout,2)))
    Bpq = 0.0_dp
    
    Bpq(wavenumber, wavenumber) = amplitude
    
    plan = fftw_plan_r2r_2d(size(inout,2),size(inout,1), Bpq, &
                          inout, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE)
    call fftw_execute_r2r(plan, Bpq, inout)
    call fftw_destroy_plan(plan)
    
!     inout = inout/real(4*size(inout,1)*size(inout,2), kind=dp)
    inout = inout/4.0_dp
    
    deallocate(Bpq)
    call fftw_cleanup()
#else
    real(kind=dp), dimension(:, :), intent(inout) :: inout
    integer :: wavenumber
    real(kind=dp), intent(in) :: amplitude
    
    call abort_handle("No FFTW3 library", __FILE__, __LINE__)
#endif
  end subroutine sine_basis
  
  ! Fourier interpolations
  subroutine fourier_intpl(out, in, Mx, My)
#if FFTW3 == 1
    use, intrinsic :: iso_c_binding 
    include 'fftw3.f03'
    
    real(kind=dp), dimension(:, :), intent(out) :: out
    real(kind=dp), dimension(:, :), intent(in) :: in
    integer, intent(in) :: Mx, My  ! Factor of scaling up
    
    type(C_PTR) :: plan
    complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: Zk, Zk_pad
    real(kind=dp), dimension(:, :), allocatable :: in_copy
    
    integer :: i, j
    
    if (size(out, 1) .ne. size(in, 1) * Mx) then
      write(6, *) "size(out) = ", size(out, 1), size(out, 2), "size(in) = ", size(in, 1), size(in, 2)
      call abort_handle("fourier_intpl: size mismatch", __FILE__, __LINE__)
    end if
    if (size(out, 2) .ne. size(in, 2) * My) then
      write(6, *) "size(out) = ", size(out, 1), size(out, 2), "size(in) = ", size(in, 1), size(in, 2)
      call abort_handle("fourier_intpl: size mismatch", __FILE__, __LINE__)
    end if
    
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
    call fftw_cleanup()
#else
    real(kind=dp), dimension(:, :), intent(out) :: out
    real(kind=dp), dimension(:, :), intent(in) :: in
    integer, intent(in) :: Mx, My  ! Factor of scaling up
    
    call abort_handle("No FFTW3 library", __FILE__, __LINE__)
#endif
  end subroutine fourier_intpl

  subroutine cosine_intpl(out, in, Mx, My)
#if FFTW3 == 1
    ! NOT IN USE
    use, intrinsic :: iso_c_binding 
    include 'fftw3.f03'

    real(kind=dp), dimension(:, :), intent(out) :: out
    real(kind=dp), dimension(:, :), intent(in) :: in
    integer, intent(in) :: Mx, My  ! Factor of scaling up
  
    type(C_PTR) :: plan
    real(kind=dp), dimension(:, :), allocatable :: Bpq, Bpq_pad
    real(kind=dp), dimension(:, :), allocatable :: in_copy
    
    integer :: i, j
    
    if (size(out, 1) .ne. size(in, 1) * Mx) then
       write(6, *) "Mx = ", Mx
       write(6, *) "size(out) = ", size(out, 1), size(out, 2), "size(in) = ", size(in, 1), size(in, 2)
      call abort_handle("cosine_intpl: size mismatch", __FILE__, __LINE__)
    end if
    if (size(out, 2) .ne. size(in, 2) * My) then
      write(6, *) "My = ", My
      write(6, *) "size(out) = ", size(out, 1), size(out, 2), "size(in) = ", size(in, 1), size(in, 2)
      call abort_handle("cosine_intpl: size mismatch", __FILE__, __LINE__)
    end if
    
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

    out = out/real(4*size(in,1)*size(in,2), kind=dp)
    
    deallocate(Bpq)
    deallocate(Bpq_pad)
    deallocate(in_copy)
    call fftw_cleanup()
#else
    real(kind=dp), dimension(:, :), intent(out) :: out
    real(kind=dp), dimension(:, :), intent(in) :: in
    integer, intent(in) :: Mx, My  ! Factor of scaling up
    
    call abort_handle("No FFTW3 library", __FILE__, __LINE__)
#endif
  end subroutine cosine_intpl

  subroutine sine_intpl(out, in, Mx, My)
#if FFTW3 == 1
    ! in, out: no boundaries => size(in) = [2^m+1-2, 2^n+1-2]
    use, intrinsic :: iso_c_binding 
    include 'fftw3.f03'

    real(kind=dp), dimension(:, :), intent(out) :: out
    real(kind=dp), dimension(:, :), intent(in) :: in
    integer, intent(in) :: Mx, My  ! Factor of scaling up
  
    type(C_PTR) :: plan
    real(kind=dp), dimension(:, :), allocatable :: Bpq, Bpq_pad
    real(kind=dp), dimension(:, :), allocatable :: in_copy
    
    integer :: i, j
    
    if ((size(out, 1)+1) .ne. (size(in, 1)+1) * Mx) then
       write(6, *) "Mx = ", Mx
       write(6, *) "size(out) = ", size(out, 1), size(out, 2), "size(in) = ", size(in, 1), size(in, 2)
      call abort_handle("sine_intpl: size mismatch", __FILE__, __LINE__)
    end if
    if ((size(out, 2)+1) .ne. (size(in, 2)+1) * My) then
      write(6, *) "My = ", My
      write(6, *) "size(out) = ", size(out, 1), size(out, 2), "size(in) = ", size(in, 1), size(in, 2)
      call abort_handle("sine_intpl: size mismatch", __FILE__, __LINE__)
    end if
    
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
    Bpq_pad = 0.0_dp
    do j = 1, size(Bpq, 2)
      do i = 1, size(Bpq, 1)
        Bpq_pad(i, j) = Bpq(i, j)
      end do
    end do
    
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
    integer, intent(in) :: Mx, My  ! Factor of scaling up
    
    call abort_handle("No FFTW3 library", __FILE__, __LINE__)
#endif
  end subroutine sine_intpl
  
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
#endif
  end subroutine sine_filter
  
  subroutine bilinear_intpl(out, in, Mx, My)
    real(kind=dp), dimension(:, :), intent(out) :: out
    real(kind=dp), dimension(:, :), intent(in) :: in
    integer, intent(in) :: Mx, My  ! Factor of scaling up
    
    integer :: i, j, idiv, jdiv, imod, jmod
    real(kind=dp) :: ialpha, jalpha
    
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
    do j = 1, size(out, 2)
      do i = 1, size(out, 1)
        idiv = (i-1)/Mx+1
        jdiv = (j-1)/My+1
        imod = mod(i-1, Mx)
        jmod = mod(j-1, My)
        ialpha = real(imod, kind=dp)/real(Mx, kind=dp)
        jalpha = real(jmod, kind=dp)/real(My, kind=dp)

        if ((imod .eq. 0) .and. (jmod .eq. 0)) then
          ! Refined and coarse grid point
          ! No interpolations
          out(i,j) = in(idiv, jdiv)
        elseif (imod .eq. 0) then
          ! Grid point with x-coordinate aligned
          ! 1D interpolation
          out(i,j) = in(idiv, jdiv)*(1.0_dp-jalpha) &
                          + in(idiv, jdiv+1)*jalpha
        elseif (jmod .eq. 0) then
          ! Grid point with y-coordinate aligned
          ! 1D interpolation
          out(i,j) = in(idiv, jdiv)*(1.0_dp-ialpha) &
                          + in(idiv+1, jdiv)*ialpha
        else
          ! Bilinear interpolation
          out(i,j) = in(idiv, jdiv)*(1.0_dp-ialpha)*(1.0_dp-jalpha) &
                     + in(idiv+1, jdiv)*ialpha*(1.0_dp-jalpha) &
                     + in(idiv, jdiv+1)*(1.0_dp-ialpha)*jalpha &
                     + in(idiv+1, jdiv+1)*ialpha*jalpha
        end if
      end do
    end do

  end subroutine bilinear_intpl

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
        if (dabs(uv_fld%up(i,j)) > 1e-14) then
          uv_fld%u_sgn(i,j) = int(sign(1.0_dp, uv_fld%up(i,j)))
        else
          uv_fld%u_sgn(i,j) = 0
        end if

        uv_fld%um(i,j) = min(uv_fld%up(i, j), 0.0_dp)
        uv_fld%up(i,j) = max(uv_fld%up(i, j), 0.0_dp)
        uv_fld%u_abs(i,j) = uv_fld%up(i,j) - uv_fld%um(i,j)
      end do
    end do

    do j = 1, size(uv_fld%vp, 2)
      do i = 1, size(uv_fld%vp, 1)
        if (dabs(uv_fld%vp(i,j)) > 1e-14) then
          uv_fld%v_sgn(i,j) = int(sign(1.0_dp, uv_fld%vp(i,j)))
        else
          uv_fld%v_sgn(i,j) = 0
        end if

        uv_fld%vm(i,j) = min(uv_fld%vp(i, j), 0.0_dp)
        uv_fld%vp(i,j) = max(uv_fld%vp(i, j), 0.0_dp)
        uv_fld%v_abs(i,j) = uv_fld%vp(i,j) - uv_fld%vm(i,j)
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

    K11e => K_e%K11e
    K12e => K_e%K12e
    K21e => K_e%K21e
    K22e => K_e%K22e

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
        
        if (dabs(u(i, j)) < 1e-14) then
          u(i, j) = 0.0_dp
        end if
      end do
    end do

    do j = 1, psi_fld%n+1
      do i = 1, psi_fld%m
        v(i, j) = psi_fld%data(i+1, j) - psi_fld%data(i, j)
        
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
      if (psi%glayer .ne. 0) &
        call abort_handle("There should not be ghost layer in psi!", __FILE__, __LINE__)
      
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
      if (K11%glayer .ne. 0) &
        call abort_handle("There should not be ghost layer in K11!", __FILE__, __LINE__)
      
      if (K11%type_id .ne. 1) &
        call abort_handle("K11 not defined at corners", __FILE__, __LINE__)
      
      do j = 1, size(K11%data, 2)
        do i = 1, size(K11%data, 1)
          ! Analytic formula: Defined on [0,1]^2 -> Need a prefactor to scale to [0,m] x [0,n]
          x = fld_x(i, m, K11%type_id)
          y = fld_x(j, n, K11%type_id)
          K11%data(i, j) = m*m/((16.0*Pi)**2)* 4.0_dp*dsin(Pi*x)*dsin(Pi*y)
          K22%data(i, j) = n*n/((16.0*Pi)**2)* 2.0_dp*dsin(Pi*x)*dsin(Pi*y)
          K12%data(i, j) = m*n/((16.0*Pi)**2)* dsin(Pi*x)*dsin(Pi*y)
          
          if (dabs(K11%data(i, j)) .le. 1D-15) K11%data(i, j) = 0.0_dp
          if (dabs(K22%data(i, j)) .le. 1D-15) K22%data(i, j) = 0.0_dp
          if (dabs(K12%data(i, j)) .le. 1D-15) K12%data(i, j) = 0.0_dp
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
  
  subroutine del2(out, in)
    ! Define on [0, 1] times [0, 1]
    type(field), intent(inout) :: out
    type(field), intent(in) :: in
    
    integer :: m, n
    
    if ((size(out%data,1) .ne. size(in%data,1)) .or. &
          (size(out%data,2) .ne. size(in%data,2)) ) then
      call abort_handle("Inconsistent del2", __FILE__, __LINE__)
    end if
    
    m = size(in%data,1)
    n = size(in%data,2)
    
    out%data = 0.0_dp
    out%data(2:m-1,2:n-1) = in%data(3:m,2:n-1) + in%data(1:m-2,2:n-1) &
                         + in%data(2:m-1,3:n) + in%data(2:m-1,1:n-2) &
                         - 4.0_dp*in%data(2:m-1,2:n-1)
    
  end subroutine del2

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
    write(6, "(a,i0)") "|  type_id: ", fld%type_id
    write(6, "(a,i0,a,i0)") "|  size(fld%data, 1), size(fld%data, 2): ", size(fld%data, 1), ", ", size(fld%data, 2)
    write(6, "(a,i0,a,i0,a,i0)") "|  m, n, dof: ", fld%m, ", ", fld%n, ", ", fld%m*fld%n
    write(6, "(a,i0)") "|  ghost layer: ", fld%glayer
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
        if (abs(fld(i,j)) > 1e-13) then
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
