module advdiff_field
  use advdiff_precision
  use advdiff_debug

  implicit none

  private

  public :: field, allocate, deallocate, zeros, &
    & scale, set, addto, addto_product, &
    & dx_field, dy_field, &
    & psi_to_uv, int_field, int_sq_field, &
    & iszero, reset_Flux, &
    & dx_field_HoriEdge, dy_field_VertEdge
    
  public :: indicator_field
  public :: uv_signed, K_grid, FluxGrid
  public :: print_info, diff_int, fld_x
  public :: bilinear_intpl, pwc_intpl, coarsen_filter, linear_intpl
  public :: neg_int
  public :: dt_CFL
  
  type field
    integer :: m, n
    real(kind = dp), dimension(:, :), pointer :: data  ! It is an array, not a pointer
    character(len = field_name_len) :: name
    integer :: type_id
  end type field

  type uv_signed
    real(kind = dp), dimension(:, :), pointer :: up, um, vp, vm
    real(kind = dp), dimension(:, :), pointer :: u_abs, v_abs
    integer, dimension(:, :), pointer :: u_sgn, v_sgn
  end type uv_signed

  type K_grid
    real(kind = dp), dimension(:, :), pointer :: K11e, K21e, K12e, K22e
    real(kind = dp), dimension(:, :), pointer :: K11c, K12c, K22c
  end type K_grid

  type FluxGrid
    real(kind = dp), dimension(:, :), pointer :: F, G
    real(kind = dp), dimension(:, :), pointer :: dqx_F, dqy_F, dqx_G, dqy_G
  end type FluxGrid

  
  interface allocate
    module procedure allocate_field, allocate_uv, &
                     allocate_Kgrid, allocate_FluxGrid
  end interface allocate

  interface deallocate
    module procedure deallocate_field, deallocate_uv, & 
                     deallocate_Kgrid, deallocate_FluxGrid
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

  subroutine coarsen_filter(out, in, mp1, np1)
    ! in, out: defined at corners
    ! size(in) = [2^m+1, 2^n+1]
    real(kind=dp), dimension(:, :), intent(inout) :: out
    real(kind=dp), dimension(:, :), intent(in) :: in
    real(kind=dp), dimension(:, :), allocatable :: tmp
    integer, intent(in) :: mp1, np1
    
    integer :: i, j, Mx, My
    
    Mx = (size(in, 1)-1)/(mp1-1)
    My = (size(in, 2)-1)/(np1-1)
    
    allocate(tmp(mp1, np1))
    
    do j = 1, size(tmp, 2)
      do i = 1, size(tmp, 2)
        ! Pick up nodal values
        tmp(i, j) = in(Mx*(i-1)+1, My*(j-1)+1)
      end do
    end do
    
    Mx = (size(out, 1)-1)/(mp1-1)
    My = (size(out, 2)-1)/(np1-1)
    call bilinear_intpl(out, tmp, Mx, My)
    
    deallocate(tmp)
    
  end subroutine coarsen_filter

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
        uv_fld%u_sgn(i,j) = int(sign(1.0_dp, uv_fld%up(i,j)))
        uv_fld%u_abs(i,j) = dabs(uv_fld%up(i,j))
        uv_fld%um(i,j) = min(uv_fld%up(i, j), 0.0_dp)
        uv_fld%up(i,j) = max(uv_fld%up(i, j), 0.0_dp)
      end do
    end do

    do j = 1, size(uv_fld%vp, 2)
      do i = 1, size(uv_fld%vp, 1)
        uv_fld%v_sgn(i,j) = int(sign(1.0_dp, uv_fld%vp(i,j))) 
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

  subroutine allocate_Kgrid(K_g, K11, K22, K12)
    type(K_grid), intent(inout) :: K_g
    type(field), intent(in) :: K11, K22, K12

    integer :: m, n, i, j

    real(kind=dp), dimension(:,:), pointer :: K11e, K12e, K21e, K22e

    m = K11%m
    n = K11%n

    if (K11%type_id .ne. 1) then
      call abort_handle("Wrong K11%type_id", __FILE__, __LINE__)
    end if 
    
    allocate(K_g%K11e(m+1, n))  ! Define at vertical edges
    allocate(K_g%K12e(m+1, n))  ! Define at vertical edges
    allocate(K_g%K21e(m, n+1))  ! Define at horizonal edges
    allocate(K_g%K22e(m, n+1))  ! Define at horizonal edges
    
    allocate(K_g%K11c(m+1, n+1))  ! Define at corners
    allocate(K_g%K22c(m+1, n+1))  ! Define at corners
    allocate(K_g%K12c(m+1, n+1))  ! Define at corners
    
    K11e => K_g%K11e
    K12e => K_g%K12e
    K21e => K_g%K21e
    K22e => K_g%K22e

    K_g%K11c = K11%data
    K_g%K22c = K22%data
    K_g%K12c = K12%data
    
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
    
  end subroutine allocate_Kgrid

  subroutine deallocate_Kgrid(K_g)
    type(K_grid), intent(inout) :: K_g

    deallocate(K_g%K11e)
    deallocate(K_g%K21e)
    deallocate(K_g%K12e)
    deallocate(K_g%K22e)
    
    deallocate(K_g%K11c)
    deallocate(K_g%K22c)
    deallocate(K_g%K12c)
    
  end subroutine deallocate_Kgrid
  
  subroutine allocate_FluxGrid(Flux, m, n)
    type(FluxGrid), intent(inout) :: Flux
    integer, intent(in) :: m, n
    
    ! Flux terms
    allocate(Flux%F(m+1, n))
    allocate(Flux%G(m, n+1))
    
    ! Gradient terms
    allocate(Flux%dqx_F(m+1, n))
    allocate(Flux%dqy_F(m+1, n))
    allocate(Flux%dqx_G(m, n+1))
    allocate(Flux%dqy_G(m, n+1))
    
    call reset_Flux(Flux)
    
  end subroutine allocate_FluxGrid

  subroutine deallocate_FluxGrid(Flux)
    type(FluxGrid), intent(inout) :: Flux

    deallocate(Flux%F)
    deallocate(Flux%G)
   
    deallocate(Flux%dqx_F)
    deallocate(Flux%dqy_F)
    deallocate(Flux%dqx_G)
    deallocate(Flux%dqy_G)
    
  end subroutine deallocate_FluxGrid
  
  subroutine reset_Flux(Flux)
    type(FluxGrid), intent(inout) :: Flux

    Flux%F = 0.0_dp
    Flux%G = 0.0_dp
    Flux%dqx_F = 0.0_dp
    Flux%dqy_F = 0.0_dp
    Flux%dqx_G = 0.0_dp
    Flux%dqy_G = 0.0_dp
  end subroutine reset_Flux

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
    real(kind=dp), dimension(:,:), intent(inout) :: d_fld
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
    real(kind=dp), dimension(:,:), intent(inout) :: d_fld
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
    real(kind=dp), dimension(:,:), intent(inout) :: d_fld
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
    real(kind=dp), dimension(:,:), intent(inout) :: d_fld
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
