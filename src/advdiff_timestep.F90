module advdiff_timestep

  use advdiff_precision
  use advdiff_timing
  use advdiff_field
  use advdiff_debug

  implicit none

  private

  public :: print_info
  public :: imposeBC
  public :: advdiff_q

  public :: timestep_heun_Kflux, timestep_fe_Kflux
  public :: timestep_fe_Kflux_Src_stat, timestep_fe_Kflux_Src_nonstat
  public :: timestep_heun_Kflux_Src_stat, timestep_heun_Kflux_Src_nonstat
  public :: timestep_LMT_1DS, timestep_LaxWendoff
  public :: timestep_LMT_2D

  interface imposeBC
    module procedure imposeBC_q, imposeBC_FG, imposeBC_Kflux
  end interface imposeBC

  interface print_info
    module procedure print_info_timestepping
  end interface print_info

contains
  pure real(kind=dp) function CFL(u_max, dt, dx)
    real(kind=dp), intent(in) :: u_max, dt, dx

    ! CFL = u*dt/dx
    CFL = u_max*dt/dx

  end function CFL

  subroutine assemble_FG_DCU(F,G, q, uv_fld)
    real(kind=dp), dimension(:,:), intent(out) :: F,G
    type(field), intent(in) :: q
    type(uv_signed), intent(in) :: uv_fld

    integer :: gl
    real(kind = dp), dimension(:, :), pointer :: lq, lup, lum, lvp, lvm

    lq => q%data
    lup => uv_fld%up
    lum => uv_fld%um
    lvp => uv_fld%vp
    lvm => uv_fld%vm

    if (q%glayer .ne. 1) then
      call abort_handle("ghost layer not equal 1", __FILE__, __LINE__)
    else
      gl = q%glayer
    end if

    !#TODO: embarrassingly parallel between F and G
    call assemble_F_DCU(F, q, lup, lum)
    call assemble_G_DCU(G, q, lvp, lvm)

  end subroutine assemble_FG_DCU
  
  subroutine assemble_F_DCU(F, q, lup, lum)
    real(kind=dp), dimension(:,:), intent(out) :: F
    type(field), intent(in) :: q
    real(kind = dp), dimension(:, :), pointer, intent(in) :: lup, lum

    integer :: i, j, gl
    real(kind = dp), dimension(:, :), pointer :: lq

    lq => q%data
    gl = q%glayer

    do j = 1, size(F, 2)
      do i = 1, size(F, 1)
        F(i,j) = lum(i,j)*lq(i+gl, j+gl) + lup(i,j)*lq(i+gl-1, j+gl)
      end do
    end do

  end subroutine assemble_F_DCU
  
  subroutine assemble_G_DCU(G, q, lvp, lvm)
    real(kind=dp), dimension(:,:), intent(out) :: G
    type(field), intent(in) :: q
    real(kind = dp), dimension(:, :), pointer, intent(in) :: lvp, lvm

    integer :: i, j, gl
    real(kind = dp), dimension(:, :), pointer :: lq

    lq => q%data
    gl = q%glayer

    do j = 1, size(G, 2)
      do i = 1, size(G, 1)
        G(i,j) = lvm(i,j)*lq(i+gl, j+gl) + lvp(i,j)*lq(i+gl, j+gl-1)
      end do
    end do

  end subroutine assemble_G_DCU

  ! Corner transport upwind
  subroutine assemble_FG_CTU(F,G, dqx, dqy, uv_fld, dt)
    real(kind=dp), dimension(:,:), intent(inout) :: F,G
    real(kind=dp), dimension(:,:), intent(in) :: dqx, dqy
    type(uv_signed), intent(in) :: uv_fld
    real(kind=dp), intent(in) :: dt

    integer :: i, j
    real(kind = dp), dimension(:, :), pointer :: lup, lum, lvp, lvm
    
    lup => uv_fld%up
    lum => uv_fld%um
    lvp => uv_fld%vp
    lvm => uv_fld%vm

    if ( (size(F,1) .ne. size(dqx, 1)) .or. (size(F,2) .ne. size(dqx, 2)) ) then
      call abort_handle("F, dqx: Incompatible size", __FILE__, __LINE__)
    end if

    ! LeVeque p454 (Fig 20.4), eqn 20.24
    !! u(i,j) aligned with F(i,j); v(i,j) aligned with G(i,j)
    !! ---+---  Upper Left, Upper Right
    !!    |
    !! ---+---  Lower Left, Lower Right(i,j)
    do j = 1, size(F, 2)
      do i = 2, size(F, 1)-1
        G(i-1,j  ) = G(i-1,j  ) - 0.5_dp*dt *lvm(i-1,j  )*lum(i,j) *dqx(i,j) ! LL
        G(i-1,j+1) = G(i-1,j+1) - 0.5_dp*dt *lvp(i-1,j+1)*lum(i,j) *dqx(i,j) ! UL
        G(i  ,j  ) = G(i  ,j  ) - 0.5_dp*dt *lvm(i  ,j  )*lup(i,j) *dqx(i,j) ! LR
        G(i  ,j+1) = G(i  ,j+1) - 0.5_dp*dt *lvp(i  ,j+1)*lup(i,j) *dqx(i,j) ! UR
      end do
    end do

    ! Left boundary
      i = 1
      do j = 1, size(F, 2)
        G(i  ,j  ) = G(i  ,j  ) - 0.5_dp*dt *lvm(i  ,j  )*lup(i,j) *dqx(i,j) ! LR
        G(i  ,j+1) = G(i  ,j+1) - 0.5_dp*dt *lvp(i  ,j+1)*lup(i,j) *dqx(i,j) ! UR
      end do

    ! Right boundary
      i = size(F, 1)
      do j = 1, size(F, 2)
        G(i-1,j  ) = G(i-1,j  ) - 0.5_dp*dt *lvm(i-1,j  )*lum(i,j) *dqx(i,j) ! LL
        G(i-1,j+1) = G(i-1,j+1) - 0.5_dp*dt *lvp(i-1,j+1)*lum(i,j) *dqx(i,j) ! UL
      end do

    !! |     |  Upper Left(i,j), Upper Right
    !! +--=--+  (i,j)
    !! |     |  Lower Left, Lower Right
    ! Assign F
    do j = 2, size(G, 2)-1
      do i = 1, size(G, 1)
        F(i  ,j-1) = F(i  ,j-1) - 0.5_dp*dt *lum(i  ,j-1)*lvm(i,j) *dqy(i,j) ! LL
        F(i+1,j-1) = F(i+1,j-1) - 0.5_dp*dt *lup(i+1,j-1)*lvm(i,j) *dqy(i,j) ! LR
        F(i  ,j  ) = F(i  ,j  ) - 0.5_dp*dt *lum(i  ,j  )*lvp(i,j) *dqy(i,j) ! UL
        F(i+1,j  ) = F(i+1,j  ) - 0.5_dp*dt *lup(i+1,j  )*lvp(i,j) *dqy(i,j) ! UR
      end do
    end do

    ! Bottom boundary
      j = 1
    do i = 1, size(G, 1)
        F(i  ,j  ) = F(i  ,j  ) - 0.5_dp*dt *lum(i  ,j  )*lvp(i,j) *dqy(i,j) ! UL
        F(i+1,j  ) = F(i+1,j  ) - 0.5_dp*dt *lup(i+1,j  )*lvp(i,j) *dqy(i,j) ! UR
    end do

    ! Top boundary
      j = size(G, 2)
    do i = 1, size(G, 1)
        F(i  ,j-1) = F(i  ,j-1) - 0.5_dp*dt *lum(i  ,j-1)*lvm(i,j) *dqy(i,j) ! LL
        F(i+1,j-1) = F(i+1,j-1) - 0.5_dp*dt *lup(i+1,j-1)*lvm(i,j) *dqy(i,j) ! LR
    end do

  end subroutine assemble_FG_CTU

  pure real(kind=dp) function philim(theta)
    real(kind=dp), intent(in) :: theta

    ! minmod
    !philim = max(0.0_dp, min(1.0_dp, theta))

    ! Superbee
    !philim = max(0.0_dp, min(1.0_dp, 2.0_dp*theta), min(2.0_dp, theta))

    ! MC
    philim = max(0.0_dp, min(0.5_dp*(1.0_dp+theta), 2.0_dp, 2.0_dp*theta))

    ! Sweby
    !philim = max(0.0_dp, min(1.0_dp, 1.5_dp*theta), min(1.5_dp, theta))
    
    ! Lax-Wendoff
    !philim = 1.0_dp
    
    ! Beam-Warming
    !philim = theta
    
    ! Upwind scheme
    !philim = 0.0_dp

  end function philim

  subroutine assemble_F_limiter(F, dqx, lu_sgn, lu_abs , dt)
    real(kind=dp), dimension(:,:), intent(inout) :: F
    real(kind=dp), dimension(:,:), intent(in) :: dqx
    integer, dimension(:, :), pointer, intent(in) :: lu_sgn
    real(kind = dp), dimension(:, :), pointer, intent(in) :: lu_abs
    real(kind=dp), intent(in) :: dt
    
    integer :: i, j
    real(kind=dp) :: tmp
      
    ! LeVeque p119
    !! u(i,j) aligned with F(i,j) & dqx(i,j); v(i,j) aligned with G(i,j) & dqy(i,j)
    do j = 1, size(F, 2)
      do i = 2, size(F, 1)-1
!         if (dabs(dqx(i,j) .gt. 1D-16) then   !! Strange behavior with valgrind: with this line on it thinks there is a leak
          F(i,j) = F(i,j) + 0.5_dp *lu_abs(i,j)*(1.0_dp-dt*lu_abs(i,j)) &
                   * philim(dqx(i-lu_sgn(i,j), j)/dqx(i,j)) *dqx(i,j)
!         end if
      end do
    end do

  end subroutine assemble_F_limiter
  
  subroutine assemble_G_limiter(G, dqy, lv_sgn, lv_abs, dt)
    real(kind=dp), dimension(:,:), intent(inout) :: G
    real(kind=dp), dimension(:,:), intent(in) :: dqy
    integer, dimension(:, :), pointer, intent(in) :: lv_sgn
    real(kind = dp), dimension(:, :), pointer, intent(in) :: lv_abs
    real(kind=dp), intent(in) :: dt
    
    integer :: i, j
    
    ! LeVeque p119
    !! u(i,j) aligned with F(i,j) & dqx(i,j); v(i,j) aligned with G(i,j) & dqy(i,j)
    do j = 2, size(G, 2)-1
      do i = 1, size(G, 1)
!         if (dabs(dqy(i,j)) .gt. 1D-16) then
          G(i,j) = G(i,j) + 0.5_dp *lv_abs(i,j)*(1.0_dp-dt*lv_abs(i,j)) &
                * philim(dqy(i, j-lv_sgn(i,j))/dqy(i,j)) *dqy(i,j)
!         end if
      end do
    end do

  end subroutine assemble_G_limiter
  
  subroutine assemble_F_LW(F, q, uv_fld, dt)
    real(kind=dp), dimension(:,:), intent(inout) :: F
    type(field), intent(in) :: q
    type(uv_signed), intent(in) :: uv_fld
    real(kind=dp), intent(in) :: dt

    integer :: i, j, gl

    integer, dimension(:, :), pointer :: lu_sgn
    real(kind = dp), dimension(:, :), pointer :: lum, lup, lu_abs, lq


    lq => q%data
    gl = q%glayer
    
    lum => uv_fld%um
    lup => uv_fld%up
    lu_sgn => uv_fld%u_sgn
    lu_abs => uv_fld%u_abs

    ! LeVeque p119
    !! u(i,j) aligned with F(i,j) & dqx(i,j); v(i,j) aligned with G(i,j) & dqy(i,j)
    do j = 1, size(F, 2)
      do i = 2, size(F, 1)-1
        F(i,j) = lum(i,j)*lq(i+gl, j+gl) + lup(i,j)*lq(i+gl-1, j+gl) &  ! Upwind
                + 0.5_dp *lu_abs(i,j)*(1.0_dp-dt*lu_abs(i,j)) *(lq(i+gl, j+gl)-lq(i+gl-1, j+gl))
                ! Correction
      end do
    end do

  end subroutine assemble_F_LW
  
  subroutine assemble_G_LW(G, q, uv_fld, dt)
    real(kind=dp), dimension(:,:), intent(inout) :: G
    type(field), intent(in) :: q
    type(uv_signed), intent(in) :: uv_fld
    real(kind=dp), intent(in) :: dt

    integer :: i, j, gl

    integer, dimension(:, :), pointer :: lv_sgn
    real(kind = dp), dimension(:, :), pointer :: lvm, lvp, lv_abs, lq
    
    lq => q%data
    gl = q%glayer
    
    lvm => uv_fld%vm
    lvp => uv_fld%vp
    lv_sgn => uv_fld%v_sgn
    lv_abs => uv_fld%v_abs

    ! LeVeque p119
    !! u(i,j) aligned with F(i,j) & dqx(i,j); v(i,j) aligned with G(i,j) & dqy(i,j)
    do j = 2, size(G, 2)-1
      do i = 1, size(G, 1)
        G(i,j) = lvm(i,j)*lq(i+gl, j+gl) + lvp(i,j)*lq(i+gl, j+gl-1) & ! Upwind
                + 0.5_dp *lv_abs(i,j)*(1.0_dp-dt*lv_abs(i,j)) *(lq(i+gl, j+gl)-lq(i+gl, j+gl-1))
                ! Correction
      end do
    end do

  end subroutine assemble_G_LW
  
  subroutine imposeBC_Kflux(K_e)
    type(K_flux), intent(inout) :: K_e
    integer :: m, n, i, j

    m = size(K_e%K22e, 1)
    n = size(K_e%K11e, 2)

    ! Set diffusivity at domain boundaries to be zero to enforce zero flux
    do j = 1, n
      K_e%K11e(1,j) = 0.0_dp
      K_e%K12e(1,j) = 0.0_dp
      K_e%K11e(m+1,j) = 0.0_dp
      K_e%K12e(m+1,j) = 0.0_dp
    end do

    do i = 1, m
      K_e%K21e(i,1) = 0.0_dp
      K_e%K22e(i,1) = 0.0_dp
      K_e%K21e(i,n+1) = 0.0_dp
      K_e%K22e(i,n+1) = 0.0_dp
    end do

  end subroutine imposeBC_Kflux

  subroutine imposeBC_FG(F, G)
    real(kind=dp), dimension(:,:), intent(inout) :: F, G

    ! Set flux at domain boundaries to be zero
    call imposeBC_F(F)
    call imposeBC_G(G)
    
  end subroutine imposeBC_FG

  subroutine imposeBC_F(F)
    real(kind=dp), dimension(:,:), intent(inout) :: F
    integer :: j

    ! Set flux at domain boundaries to be zero
    do j = 1, size(F, 2)
      F(1,j) = 0.0_dp
      F(size(F, 1),j) = 0.0_dp
    end do

  end subroutine imposeBC_F

  subroutine imposeBC_G(G)
    real(kind=dp), dimension(:,:), intent(inout) :: G
    integer :: i

    ! Set flux at domain boundaries to be zero
    do i = 1, size(G, 1)
      G(i,1) = 0.0_dp
      G(i,size(G, 2)) = 0.0_dp
    end do

  end subroutine imposeBC_G
  
  subroutine imposeBC_q(q)
    type(field), intent(inout) :: q
    integer :: i, j, gl, m, n

    m = q%m
    n = q%n

    ! Set q in ghost layer to be the nearest value in the domain
    ! Along y
    do gl = 1, q%glayer  ! Inner to outer
      do j = 1+q%glayer, n+q%glayer  ! Exclude horizontal ghost layers
        q%data(1+q%glayer-gl,j) = q%data(1+q%glayer,j)
        q%data(m+q%glayer+gl,j) = q%data(m+q%glayer,j)
      end do
    end do

    ! Along x
    do gl = 1, q%glayer
      do i = 1, m+2*q%glayer    ! Now include the missing horizontal ghost layers
        q%data(i,1+q%glayer-gl) = q%data(i,1+q%glayer)
        q%data(i,n+q%glayer+gl) = q%data(i,n+q%glayer)
      end do
    end do

  end subroutine imposeBC_q
  
    ! Consider the ghost point
  subroutine eval_dqdt_Kfull(dqdt, q, K11e, K12e, K21e, K22e)
    ! 5-point stencil: ignore off diagonal terms
    ! Wrong BC??
    type(field), intent(inout) :: dqdt
    type(field), intent(in) :: q

    real(kind=dp), dimension(:,:) :: K11e, K12e, K21e, K22e

    integer :: i, j
    real(kind = dp), dimension(:, :), pointer :: ldqdt, lq

    integer :: gl   ! ghost layer

    lq => q%data
    ldqdt => dqdt%data

    gl = q%glayer

    ! 9points_Laplacian_FVM.mw
    do j = 1, q%n
      do i = 1, q%m
        ldqdt(i+gl,j+gl) = ldqdt(i+gl,j+gl) -( K11e(i+1,j) + K11e(i,j) + K22e(i,j+1) + K22e(i,j) ) *lq(i+gl, j+gl) &
                    + ( K11e(i+1,j) + 0.25_dp*(K21e(i,j+1) - K21e(i,j)) ) *lq(i+1+gl,j+gl) &
                    + ( K11e(i  ,j) - 0.25_dp*(K21e(i,j+1) - K21e(i,j)) ) *lq(i-1+gl,j+gl) &
                    + ( K22e(i,j+1) + 0.25_dp*(K12e(i+1,j) - K12e(i,j)) ) *lq(i+gl,j+1+gl) &
                    + ( K22e(i,j  ) - 0.25_dp*(K12e(i+1,j) - K12e(i,j)) ) *lq(i+gl,j-1+gl) &
                    + ( K12e(i+1,j) + K21e(i,j+1) )* 0.25_dp *lq(i+1+gl,j+1+gl) &
                    - ( K12e(i+1,j) + K21e(i,j  ) )* 0.25_dp *lq(i+1+gl,j-1+gl) &
                    - ( K12e(i  ,j) + K21e(i,j+1) )* 0.25_dp *lq(i-1+gl,j+1+gl) &
                    + ( K12e(i  ,j) + K21e(i,j  ) )* 0.25_dp *lq(i-1+gl,j-1+gl)
      end do
    end do
    
  end subroutine eval_dqdt_Kfull
  
    ! Consider the ghost point
  subroutine eval_dqdt_FG(dqdt, q, F, G)
    ! 5-point stencil: ignore off diagonal terms
    ! Assumed zero flux BC
    type(field), intent(out) :: dqdt
    type(field), intent(in) :: q

    real(kind=dp), dimension(:,:) :: F, G

    integer :: i, j
    real(kind = dp), dimension(:, :), pointer :: ldqdt, lq

    integer :: gl   ! ghost layer

    lq => q%data
    ldqdt => dqdt%data
    gl = q%glayer

    do j = 1, q%n
      do i = 1, q%m
        ldqdt(i+gl,j+gl) = -F(i+1,j) +F(i,j) -G(i,j+1) +G(i,j)
      end do
    end do
    
    ! Ensure no floating point issues
    call imposeBC(dqdt)
    
  end subroutine eval_dqdt_FG
  
  subroutine eval_dqdt_F(dqdt, q, F)
    ! 5-point stencil: ignore off diagonal terms
    ! Assumed zero flux BC
    type(field), intent(out) :: dqdt
    type(field), intent(in) :: q
    real(kind=dp), dimension(:,:), intent(in) :: F

    integer :: i, j
    real(kind = dp), dimension(:, :), pointer :: ldqdt, lq

    integer :: gl   ! ghost layer

    lq => q%data
    ldqdt => dqdt%data
    gl = q%glayer

    do j = 1, q%n
      do i = 1, q%m
        ldqdt(i+gl,j+gl) = -F(i+1,j) +F(i,j)
      end do
    end do
    
    ! Ensure no floating point issues
    call imposeBC(dqdt)
    
  end subroutine eval_dqdt_F
  
  subroutine eval_dqdt_G(dqdt, q, G)
    ! 5-point stencil: ignore off diagonal terms
    ! Assumed zero flux BC
    type(field), intent(out) :: dqdt
    type(field), intent(in) :: q
    real(kind=dp), dimension(:,:), intent(in) :: G

    integer :: i, j
    real(kind = dp), dimension(:, :), pointer :: ldqdt, lq

    integer :: gl   ! ghost layer

    lq => q%data
    ldqdt => dqdt%data
    gl = q%glayer

    do j = 1, q%n
      do i = 1, q%m
        ldqdt(i+gl,j+gl) = -G(i,j+1) +G(i,j)
      end do
    end do
    
    ! Ensure no floating point issues
    call imposeBC(dqdt)
    
  end subroutine eval_dqdt_G
  
  ! Consider the ghost point
  subroutine eval_dqdt_scr_nonstat(dqdt, q, t)
    ! 5-point stencil: ignore off diagonal terms
    ! Assumed zero flux BC
    type(field), intent(out) :: dqdt
    type(field), intent(in) :: q
    real(kind = dp), intent(in) :: t

    integer :: i, j
    real(kind = dp), dimension(:, :), pointer :: ldqdt, lq
    
    integer :: gl   ! ghost layer
    real(kind = dp) :: x, y

    lq => q%data
    ldqdt => dqdt%data
    gl = q%glayer
    
    do j = 1, q%n
      do i = 1, q%m
        x = fld_x(i, q%m, q%type_id)
        y = fld_x(j, q%n, q%type_id)
        ldqdt(i+gl,j+gl) = S_rhs_nonstat(x, y, t)
      end do
    end do
    
    ! Ensure no floating point issues
    call imposeBC(dqdt)

  end subroutine eval_dqdt_scr_nonstat
  
    ! Consider the ghost point
  subroutine eval_dqdt_scr_stat(dqdt, q, lambda)
    ! 5-point stencil: ignore off diagonal terms
    ! Assumed zero flux BC
    type(field), intent(out) :: dqdt
    type(field), intent(in) :: q
    real(kind = dp), intent(in) :: lambda
    
    integer :: i, j
    real(kind = dp), dimension(:, :), pointer :: ldqdt, lq
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    
    integer :: gl   ! ghost layer
    real(kind = dp) :: x, y, qxy

    lq => q%data
    ldqdt => dqdt%data
    gl = q%glayer
    
    do j = 1, q%n
      do i = 1, q%m
        x = fld_x(i, q%m, q%type_id)
        y = fld_x(j, q%n, q%type_id)
        qxy = lq(i+gl,j+gl)
        ldqdt(i+gl,j+gl) = S_rhs_stat(x,y, qxy, lambda)
      end do
    end do
    
    ! Ensure no floating point issues
    call imposeBC(dqdt)

  end subroutine eval_dqdt_scr_stat

  subroutine timestep_fe(q, dqdt, dt)
    type(field), intent(inout) :: q
    type(field), intent(in) :: dqdt
    real(kind = dp), intent(in) :: dt

    integer :: i, j
    real(kind = dp), dimension(:, :), pointer :: lq
    real(kind = dp), dimension(:, :), pointer :: ldqdt

    lq => q%data
    ldqdt => dqdt%data

    do j = 1, size(lq, 2)
      do i = 1, size(lq, 1)
        lq(i, j) = lq(i, j) + dt * ldqdt(i, j)
      end do
    end do

  end subroutine timestep_fe
  
  subroutine timestep_heun_Kflux(q, dt, K_e)
    ! Heun's method
    type(field), intent(inout) :: q
    type(K_flux), intent(in) :: K_e

    real(kind = dp), intent(in) :: dt

    real(kind = dp), dimension(:,:), pointer :: K11e, K12e, K21e, K22e

    type(field) :: dqdt
    type(field) :: qc

    call allocate(dqdt, q%m, q%n, 'dqdt', glayer=q%glayer, type_id=q%type_id)
    call allocate(qc, q%m, q%n, 'qc', glayer=q%glayer, type_id=q%type_id)

    K11e => K_e%K11e
    K12e => K_e%K12e
    K21e => K_e%K21e
    K22e => K_e%K22e

    ! Evaluate dqdt and copy q to qc
    call zeros(dqdt)
    call eval_dqdt_Kfull(dqdt, q, K11e, K12e, K21e, K22e)
    call set(qc, q)

    ! Prediction step
    call timestep_fe(qc, dqdt, dt)
    call imposeBC(qc)

    ! Correction step (add to dqdt, not rewriting it)
    call eval_dqdt_Kfull(dqdt, qc, K11e, K12e, K21e, K22e)
    
    ! Full timestepping
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    call imposeBC(q)

    call deallocate(dqdt)
    call deallocate(qc)

  end subroutine timestep_heun_Kflux
  
  subroutine timestep_fe_Kflux(q, dt, K_e)
    ! Heun's method
    type(field), intent(inout) :: q
    type(K_flux), intent(in) :: K_e

    real(kind = dp), intent(in) :: dt

    real(kind = dp), dimension(:,:), pointer :: K11e, K12e, K21e, K22e

    type(field) :: dqdt

    call allocate(dqdt, q%m, q%n, 'dqdt', glayer=q%glayer, type_id=q%type_id)

    K11e => K_e%K11e
    K12e => K_e%K12e
    K21e => K_e%K21e
    K22e => K_e%K22e

    ! Evaluate dqdt and copy q to qc
    call zeros(dqdt)
    call eval_dqdt_Kfull(dqdt, q, K11e, K12e, K21e, K22e)
    call timestep_fe(q, dqdt, dt)
    call imposeBC(q)

    call deallocate(dqdt)

  end subroutine timestep_fe_Kflux
  
  subroutine eval_dqdt_Kfull_scr_stat(dqdt, q, K11e, K12e, K21e, K22e, lambda)
    type(field), intent(inout) :: dqdt, q
    real(kind = dp), dimension(:,:), intent(in) :: K11e, K12e, K21e, K22e

    real(kind = dp), intent(in) :: lambda

    type(field) :: dqdt_S
    
    call allocate(dqdt_S, q%m, q%n, 'dqdt_S', glayer=q%glayer, type_id=q%type_id)

    call zeros(dqdt)
    call eval_dqdt_Kfull(dqdt, q, K11e, K12e, K21e, K22e) ! From diffusion

    call zeros(dqdt_S)
    call eval_dqdt_scr_stat(dqdt_S, q, lambda) ! From source term

    call addto(dqdt, dqdt_S, 1.0_dp)
    
    call deallocate(dqdt_S)
    
  end subroutine eval_dqdt_Kfull_scr_stat
  
  subroutine eval_dqdt_Kfull_scr_nonstat(dqdt, q, t, K11e, K12e, K21e, K22e)
    type(field), intent(inout) :: dqdt, q
    real(kind = dp), dimension(:,:), intent(in) :: K11e, K12e, K21e, K22e

    real(kind = dp), intent(in) :: t

    type(field) :: dqdt_S
    
    call allocate(dqdt_S, q%m, q%n, 'dqdt_S', glayer=q%glayer, type_id=q%type_id)

    call zeros(dqdt)
    call eval_dqdt_Kfull(dqdt, q, K11e, K12e, K21e, K22e) ! From diffusion

    call zeros(dqdt_S)
    call eval_dqdt_scr_nonstat(dqdt_S, q, t) ! From source term

    call addto(dqdt, dqdt_S, 1.0_dp)
    
    call deallocate(dqdt_S)
    
  end subroutine eval_dqdt_Kfull_scr_nonstat
  
  subroutine timestep_fe_Kflux_Src_stat(q, dt, K_e, lambda)
    type(field), intent(inout) :: q
    type(K_flux), intent(in) :: K_e

    real(kind = dp), intent(in) :: dt
    real(kind = dp), intent(in) :: lambda

    real(kind = dp), dimension(:,:), pointer :: K11e, K12e, K21e, K22e

    type(field) :: dqdt

    call allocate(dqdt, q%m, q%n, 'dqdt', glayer=q%glayer, type_id=q%type_id)

    K11e => K_e%K11e
    K12e => K_e%K12e
    K21e => K_e%K21e
    K22e => K_e%K22e

    ! Evaluate dqdt
    call eval_dqdt_Kfull_scr_stat(dqdt, q, K11e, K12e, K21e, K22e, lambda)
    
    call timestep_fe(q, dqdt, dt)
    call imposeBC(q)

    call deallocate(dqdt)
    
  end subroutine timestep_fe_Kflux_Src_stat
  
  subroutine timestep_heun_Kflux_Src_stat(q, dt, K_e, lambda)
    ! Heun's method
    type(field), intent(inout) :: q
    type(K_flux), intent(in) :: K_e

    real(kind = dp), intent(in) :: dt, lambda

    real(kind = dp), dimension(:,:), pointer :: K11e, K12e, K21e, K22e

    type(field) :: dqdt
    type(field) :: qc
    type(field) :: dqdtc

    call allocate(dqdt, q%m, q%n, 'dqdt', glayer=q%glayer, type_id=q%type_id)
    call allocate(qc, q%m, q%n, 'qc', glayer=q%glayer, type_id=q%type_id)
    call allocate(dqdtc, q%m, q%n, 'dqdtc', glayer=q%glayer, type_id=q%type_id)

    K11e => K_e%K11e
    K12e => K_e%K12e
    K21e => K_e%K21e
    K22e => K_e%K22e

    ! Evaluate dqdt and copy q to qc
    call zeros(dqdt)
    call eval_dqdt_Kfull_scr_stat(dqdt, q, K11e, K12e, K21e, K22e, lambda)
    call set(qc, q)

    ! Prediction step
    call timestep_fe(qc, dqdt, dt)
    call imposeBC(qc)

    ! Correction step
    call zeros(dqdtc)
    call eval_dqdt_Kfull_scr_stat(dqdt, q, K11e, K12e, K21e, K22e, lambda)
    call addto(dqdt, dqdtc, 1.0_dp)

    ! Full timestepping
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    call imposeBC(q)

    call deallocate(dqdt)
    call deallocate(qc)
    call deallocate(dqdtc)

  end subroutine timestep_heun_Kflux_Src_stat
  
  subroutine timestep_fe_Kflux_Src_nonstat(q, dt, t, K_e)
    ! Heun's method
    type(field), intent(inout) :: q
    type(K_flux), intent(in) :: K_e

    real(kind = dp), intent(in) :: dt, t

    real(kind = dp), dimension(:,:), pointer :: K11e, K12e, K21e, K22e

    type(field) :: dqdt

    call allocate(dqdt, q%m, q%n, 'dqdt', glayer=q%glayer, type_id=q%type_id)

    K11e => K_e%K11e
    K12e => K_e%K12e
    K21e => K_e%K21e
    K22e => K_e%K22e

    ! Evaluate dqdt
    call eval_dqdt_Kfull_scr_nonstat(dqdt, q, t, K11e, K12e, K21e, K22e)
    
    call timestep_fe(q, dqdt, dt)
    call imposeBC(q)

    call deallocate(dqdt)
    
  end subroutine timestep_fe_Kflux_Src_nonstat
  
  subroutine timestep_heun_Kflux_Src_nonstat(q, dt, t, K_e)
    ! Heun's method
    type(field), intent(inout) :: q
    type(K_flux), intent(in) :: K_e

    real(kind = dp), intent(in) :: dt, t

    real(kind = dp), dimension(:,:), pointer :: K11e, K12e, K21e, K22e

    type(field) :: dqdt
    type(field) :: qc
    type(field) :: dqdtc

    call allocate(dqdt, q%m, q%n, 'dqdt', glayer=q%glayer, type_id=q%type_id)
    call allocate(qc, q%m, q%n, 'qc', glayer=q%glayer, type_id=q%type_id)
    call allocate(dqdtc, q%m, q%n, 'dqdtc', glayer=q%glayer, type_id=q%type_id)

    K11e => K_e%K11e
    K12e => K_e%K12e
    K21e => K_e%K21e
    K22e => K_e%K22e

    ! Evaluate dqdt and copy q to qc
    call set(qc, q)
    call zeros(dqdt)
    call eval_dqdt_Kfull_scr_nonstat(dqdt, q, t, K11e, K12e, K21e, K22e)

    ! Prediction step
    call timestep_fe(qc, dqdt, dt)
    call imposeBC(qc)

    ! Correction step
    call zeros(dqdtc)
    call eval_dqdt_Kfull_scr_nonstat(dqdtc, qc, t+dt, K11e, K12e, K21e, K22e)
    call addto(dqdt, dqdtc, 1.0_dp)

    ! Full timestepping
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    call imposeBC(q)
    
    call deallocate(dqdt)
    call deallocate(qc)
    call deallocate(dqdtc)

  end subroutine timestep_heun_Kflux_Src_nonstat
  
  subroutine timestep_LMT_1DS(q, dt, uv_fld)
    type(field), intent(inout) :: q
    type(uv_signed), intent(in) :: uv_fld

    real(kind = dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), allocatable :: F, G, dqx, dqy

    type(field) :: dqdt
    
    call allocate(dqdt, q%m, q%n, 'dqdt', glayer=q%glayer, type_id=q%type_id)

    allocate(F(q%m+1, q%n))
    allocate(G(q%m, q%n+1))
    allocate(dqx(q%m+1, q%n))
    allocate(dqy(q%m, q%n+1))

    ! F, G at boundaries will not be touched
    F = 0.0_dp
    G = 0.0_dp

    ! Strang splitting: x, y direction - XYX
    ! X
    call assemble_F_DCU(F, q, uv_fld%up, uv_fld%um)
    call dx_field(dqx, q)
    call assemble_F_limiter(F, dqx, uv_fld%u_sgn, uv_fld%u_abs, 0.5_dp*dt)  ! Apply limiter
    
    call eval_dqdt_F(dqdt, q, F)
    call timestep_fe(q, dqdt, 0.5_dp*dt)

    ! Y
    call assemble_G_DCU(G, q, uv_fld%vp, uv_fld%vm)
    call dy_field(dqy, q)
    call assemble_G_limiter(G, dqy, uv_fld%v_sgn, uv_fld%v_abs, dt)  ! Apply limiter

    call eval_dqdt_G(dqdt, q, G)
    call timestep_fe(q, dqdt, dt)
    
    ! X
    call assemble_F_DCU(F, q, uv_fld%up, uv_fld%um)
    call dx_field(dqx, q)
    call assemble_F_limiter(F, dqx, uv_fld%u_sgn, uv_fld%u_abs, 0.5_dp*dt)  ! Apply limiter
    
    call eval_dqdt_F(dqdt, q, F)
    call timestep_fe(q, dqdt, 0.5_dp*dt)

    ! Deallocate
    call deallocate(dqdt)

    deallocate(dqx)
    deallocate(dqy)
    deallocate(F)
    deallocate(G)

  end subroutine timestep_LMT_1DS
  
  subroutine timestep_LaxWendoff(q, dt, uv_fld)
    type(field), intent(inout) :: q
    type(uv_signed), intent(in) :: uv_fld

    real(kind = dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), allocatable :: F, G

    type(field) :: dqdt
    
    call allocate(dqdt, q%m, q%n, 'dqdt', glayer=q%glayer, type_id=q%type_id)

    allocate(F(q%m+1, q%n))
    allocate(G(q%m, q%n+1))

    ! F, G at boundaries will not be touched
    F = 0.0_dp
    G = 0.0_dp

    ! Strang splitting: x, y direction - XYX
    ! X
    call assemble_F_LW(F, q, uv_fld, 0.5_dp*dt)  ! LaxWendoff
    
    call eval_dqdt_F(dqdt, q, F)
    call timestep_fe(q, dqdt, 0.5_dp*dt)

    ! Y
    call assemble_G_LW(G, q, uv_fld, dt)  ! LaxWendoff

    call eval_dqdt_G(dqdt, q, G)
    call timestep_fe(q, dqdt, dt)
    
    ! X
    call assemble_F_LW(F, q, uv_fld, 0.5_dp*dt)  ! LaxWendoff
    
    call eval_dqdt_F(dqdt, q, F)
    call timestep_fe(q, dqdt, 0.5_dp*dt)

    ! Deallocate
    call deallocate(dqdt)

    deallocate(F)
    deallocate(G)

  end subroutine timestep_LaxWendoff
  
  subroutine timestep_LMT_2D(q, dt, uv_fld)
    type(field), intent(inout) :: q
    type(uv_signed), intent(in) :: uv_fld

    real(kind = dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), allocatable :: F, G, dqx, dqy

    type(field) :: dqdt
    
    call allocate(dqdt, q%m, q%n, 'dqdt', glayer=q%glayer, type_id=q%type_id)

    allocate(F(q%m+1, q%n))
    allocate(G(q%m, q%n+1))
    allocate(dqx(q%m+1, q%n))
    allocate(dqy(q%m, q%n+1))

    ! F, G at boundaries will not be touched
    F = 0.0_dp
    G = 0.0_dp
    dqx = 0.0_dp
    dqy = 0.0_dp
    
    call assemble_F_DCU(F, q, uv_fld%up, uv_fld%um)
    call dx_field(dqx, q)
    call assemble_F_limiter(F, dqx, uv_fld%u_sgn, uv_fld%u_abs, 0.5_dp*dt)  ! Apply limiter

    call assemble_G_DCU(G, q, uv_fld%vp, uv_fld%vm)
    call dy_field(dqy, q)
    call assemble_G_limiter(G, dqy, uv_fld%v_sgn, uv_fld%v_abs, 0.5_dp*dt)  ! Apply limiter

    call eval_dqdt_FG(dqdt, q, F, G)
    call timestep_fe(q, dqdt, dt)
    call imposeBC(q)

    ! Deallocate
    call deallocate(dqdt)

    deallocate(dqx)
    deallocate(dqy)
    deallocate(F)
    deallocate(G)

  end subroutine timestep_LMT_2D
  
  subroutine print_info_timestepping(psi, dt)
    type(field), intent(in) :: psi
    real(kind=dp), intent(in) :: dt
    
    real(kind=dp) :: u_max, v_max

    real(kind=dp), dimension(:, :), allocatable :: u, v
    
    allocate(u(psi%m+1, psi%n))
    allocate(v(psi%m, psi%n+1))

    call psi_to_uv(u, v, psi)
    
    u_max = maxval(dabs(u))
    v_max = maxval(dabs(v))
    
    deallocate(u)
    deallocate(v)
    
    write(6, "(a)") ""
    write(6, "(a, a)") " Time stepping: ", trim(psi%name)
    write(6, "(a,"//dp_chr//",a,i0,a,"//dp_chr//")") "|  uv_max = ", max(u_max, v_max), "  dx = ", 1, " and dt = ", dt
    write(6, "(a,"//dp_chr//")") "|  Courant number in x-direction = ", CFL(u_max, dt, 1.0_dp)
    write(6, "(a,"//dp_chr//")") "|  Courant number in y-direction = ", CFL(v_max, dt, 1.0_dp)

  end subroutine print_info_timestepping

  subroutine advdiff_q(q, psi, K11, K22, K12, T, nts)
    type(field), intent(inout) :: q
    type(field), intent(in) :: psi, K11, K22, K12
    real(kind=dp), intent(in) :: T
    integer, intent(in) :: nts

    integer :: i
    type(uv_signed) :: uv_fld
    type(K_flux) :: K_e
    real(kind=dp) :: dt

    call imposeBC(q)
    
    call allocate(uv_fld, psi)
    call allocate(K_e, K11, K22, K12)
    call imposeBC(K_e)
    
    dt = T/nts
    do i = 1, nts
      call timestep_heun_Kflux(q, dt, K_e)
      call timestep_LMT_2D(q, dt, uv_fld)
!       call timestep_LMT_1DS(q, dt, uv_fld)
    end do

    call deallocate(uv_fld)
    call deallocate(K_e)

  end subroutine advdiff_q

end module advdiff_timestep
