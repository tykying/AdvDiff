module advdiff_timestep

  use advdiff_precision
  use advdiff_timing
  use advdiff_field
  use advdiff_debug

  implicit none

  private

  public :: timestep_heun_K, timestep_LMT, timestep_fe_Kflux
  public :: imposeBC, eval_dqdt
  public :: advdiff_q
  public :: timestep_heun_T_ops

  type(timer), save :: timestep_timer, assemble_timer
  public :: reset_timestep_timers, print_timestep_timers

  public :: maxmax, print_info
  public :: timestep_fe_Kflux_Src_stat, timestep_fe_Kflux_Src_nonstat
  public :: timestep_heun_Kflux_Src_stat, timestep_heun_Kflux_Src_nonstat
  public :: timestep_LMT_1DS, timestep_LMT_1DS2
  
  
  interface eval_dqdt
    module procedure eval_dqdt_Kfull, eval_dqdt_FG
  end interface eval_dqdt

  interface imposeBC
    module procedure imposeBC_q, imposeBC_FG, imposeBC_Kflux
  end interface imposeBC

  interface timestep_heun_K
    module procedure timestep_heun_Kflux
  end interface timestep_heun_K

  interface print_info
    module procedure print_info_timestepPing
  end interface print_info

  interface advdiff_q
     module procedure advdiff_qC
!    module procedure advdiff_qMPFA
  end interface advdiff_q

contains
  pure real(kind=dp) function maxmax(u_abs, v_abs)
    real(kind=dp), dimension(:,:), intent(in) :: u_abs, v_abs

    ! Assume dx = dy = 1
    ! CFL = u*dt/dx
    maxmax = max(maxval(u_abs), maxval(v_abs))

  end function maxmax

  pure real(kind=dp) function CFL(uv_max, dt, dx)
    real(kind=dp), intent(in) :: uv_max, dt, dx

    ! Assume dx = dy = 1
    ! CFL = u*dt/dx
    CFL = uv_max*dt/dx

  end function CFL

  subroutine assemble_FG_DCU(F,G, q, uv_fld)
    real(kind=dp), dimension(:,:), intent(out) :: F,G
    type(field), intent(in) :: q
    type(uv_signed), intent(in) :: uv_fld

    integer :: i, j, gl
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
    call assemble_F_DCU(F, q, uv_fld)
    call assemble_G_DCU(G, q, uv_fld)

  end subroutine assemble_FG_DCU
  
  subroutine assemble_F_DCU(F, q, uv_fld)
    real(kind=dp), dimension(:,:), intent(out) :: F
    type(field), intent(in) :: q
    type(uv_signed), intent(in) :: uv_fld

    integer :: i, j, gl
    real(kind = dp), dimension(:, :), pointer :: lq, lup, lum, lvp, lvm

    lq => q%data
    lup => uv_fld%up
    lum => uv_fld%um

    if (q%glayer .ne. 1) then
      call abort_handle("ghost layer not equal 1", __FILE__, __LINE__)
    else
      gl = q%glayer
    end if

    do j = 1, size(F, 2)
      do i = 1, size(F, 1)
        F(i,j) = lum(i,j)*lq(i+gl, j+gl) + lup(i,j)*lq(i+gl-1, j+gl)
      end do
    end do

  end subroutine assemble_F_DCU

  subroutine assemble_F_DCU2(F, q, lup, lum)
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

  end subroutine assemble_F_DCU2
  
  subroutine assemble_G_DCU(G, q, uv_fld)
    real(kind=dp), dimension(:,:), intent(out) :: G
    type(field), intent(in) :: q
    type(uv_signed), intent(in) :: uv_fld

    integer :: i, j, gl
    real(kind = dp), dimension(:, :), pointer :: lq, lup, lum, lvp, lvm

    lq => q%data
    lvp => uv_fld%vp
    lvm => uv_fld%vm

    if (q%glayer .ne. 1) then
      call abort_handle("ghost layer not equal 1", __FILE__, __LINE__)
    else
      gl = q%glayer
    end if

    do j = 1, size(G, 2)
      do i = 1, size(G, 1)
        G(i,j) = lvm(i,j)*lq(i+gl, j+gl) + lvp(i,j)*lq(i+gl, j+gl-1)
      end do
    end do

  end subroutine assemble_G_DCU
  
  subroutine assemble_G_DCU2(G, q, lvp, lvm)
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

  end subroutine assemble_G_DCU2

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

  pure real(kind=dp) function phi_limiter(theta)
    real(kind=dp), intent(in) :: theta

    ! Superbee
    !phi_limiter = max(0.0_dp, min(1.0_dp, 2.0_dp*theta), min(2.0_dp, theta))

    ! MC
    !phi_limiter = max(0.0_dp, min(0.5_dp*(1.0_dp+theta), 2.0_dp, 2.0_dp*theta))

    ! Sweby
    !phi_limiter = max(0.0_dp, min(1.0_dp, 1.5_dp*theta), min(1.5_dp, theta))
    !phi_limiter = 0.0_dp
    
    ! Lax-Wendoff
    phi_limiter = 1.0_dp
  end function phi_limiter
  
  subroutine assemble_phi_x(phi_x, dqx, lu_sgn)
    real(kind=dp), dimension(:,:), intent(out) :: phi_x
    real(kind=dp), dimension(:,:), intent(in) :: dqx
    integer, dimension(:, :), pointer, intent(in) :: lu_sgn
    
    integer :: i, j

    real(kind = dp) :: theta_x

    ! LeVeque p119
    !! phi_x(i,j) aligned with F(i,j) & dqx(i,j) & u(i,j); v(i,j) aligned with G(i,j) & dqy(i,j)
    do j = 1, size(phi_x, 2)
      do i = 2, size(phi_x, 1)-1
        if (dabs(dqx(i,j)) < 1.0e-13) then
          theta_x = sign(1.0_dp, dqx(i,j)) * 1.0e14
        else
          theta_x = dqx(i-lu_sgn(i,j), j)/dqx(i,j)
        end if
        
        phi_x(i,j)  = phi_limiter(theta_x)
      end do
    end do
    
  end subroutine assemble_phi_x
  
  subroutine assemble_phi_y(phi_y, dqy, lv_sgn)
    real(kind=dp), dimension(:,:), intent(out) :: phi_y
    real(kind=dp), dimension(:,:), intent(in) :: dqy
    integer, dimension(:, :), pointer, intent(in) :: lv_sgn

    integer :: i, j

    real(kind = dp) :: theta_y

    ! LeVeque p119
    !! phi_x(i,j) aligned with F(i,j) & dqx(i,j) & u(i,j); v(i,j) aligned with G(i,j) & dqy(i,j)
    do j = 2, size(phi_y, 2)-1
      do i = 1, size(phi_y, 1)
        if (dabs(dqy(i,j)) < 1.0e-13) then
          theta_y = sign(1.0_dp, dqy(i,j)) * 1.0e14
        else
          theta_y = dqy(i, j-lv_sgn(i,j))/dqy(i,j)
        end if
        
        phi_y(i,j)  = phi_limiter(theta_y)
      end do
    end do
    
  end subroutine assemble_phi_y

  subroutine assemble_FG_limiter(F,G, dqx, dqy, uv_fld, dt)
    real(kind=dp), dimension(:,:), intent(inout) :: F,G
    real(kind=dp), dimension(:,:), intent(in) :: dqx, dqy
    type(uv_signed), intent(in) :: uv_fld
    real(kind=dp), intent(in) :: dt

    call assemble_F_limiter(F, dqx, uv_fld, dt)
    call assemble_G_limiter(G, dqy, uv_fld, dt)

  end subroutine assemble_FG_limiter

  subroutine assemble_F_limiter(F, dqx, uv_fld, dt)
    real(kind=dp), dimension(:,:), intent(inout) :: F
    real(kind=dp), dimension(:,:), intent(in) :: dqx
    type(uv_signed), intent(in) :: uv_fld
    real(kind=dp), intent(in) :: dt

    integer :: i, j

    integer, dimension(:, :), pointer :: lu_sgn
    real(kind = dp), dimension(:, :), pointer :: lu_abs

    real(kind = dp) :: theta_x, phi_x

    lu_sgn => uv_fld%u_sgn
    lu_abs => uv_fld%u_abs

    ! LeVeque p119
    !! u(i,j) aligned with F(i,j) & dqx(i,j); v(i,j) aligned with G(i,j) & dqy(i,j)
    do j = 1, size(F, 2)
      do i = 2, size(F, 1)-1
        if (dabs(dqx(i,j)) < 1.0e-13) then
          theta_x = sign(1.0_dp, dqx(i,j)) * 1.0e14
        else
          theta_x = dqx(i-lu_sgn(i,j), j)/dqx(i,j)
        end if
        phi_x  = phi_limiter(theta_x)

        F(i,j) = F(i,j) + 0.5_dp *lu_abs(i,j)*(1.0_dp-dt*lu_abs(i,j)) *phi_x *dqx(i,j)

!         if (1.0_dp-dt*lu_abs(i,j) < 0.0_dp) then
!           write(6, *) "Check limiter CFL conditions: "
!           write(6, "(a,"//dp_chr//",a,"//dp_chr//",a,"//dp_chr//")") "dt = ", dt, "; lu_abs(i,j) = ", lu_abs(i,j), "; dt*lu_abs(i,j) = ", dt*lu_abs(i,j)
!         end if
      end do
    end do

  end subroutine assemble_F_limiter
  
  subroutine assemble_G_limiter(G, dqy, uv_fld, dt)
    real(kind=dp), dimension(:,:), intent(inout) :: G
    real(kind=dp), dimension(:,:), intent(in) :: dqy
    type(uv_signed), intent(in) :: uv_fld
    real(kind=dp), intent(in) :: dt

    integer :: i, j
    real(kind = dp), dimension(:, :), pointer :: lvp, lvm

    integer, dimension(:, :), pointer :: lv_sgn
    real(kind = dp), dimension(:, :), pointer :: lv_abs

    real(kind = dp) :: theta_y, phi_y

    lvp => uv_fld%vp
    lvm => uv_fld%vm

    lv_sgn => uv_fld%v_sgn
    lv_abs => uv_fld%v_abs

    ! LeVeque p119
    !! u(i,j) aligned with F(i,j) & dqx(i,j); v(i,j) aligned with G(i,j) & dqy(i,j)
    do j = 2, size(G, 2)-1
      do i = 1, size(G, 1)
        if (dabs(dqy(i,j)) < 1.0e-13) then
          theta_y = sign(1.0_dp, dqy(i,j)) * 1.0e14
        else
          theta_y = dqy(i, j-lv_sgn(i,j))/dqy(i,j)
        end if
        phi_y  = phi_limiter(theta_y)

        G(i,j) = G(i,j) + 0.5_dp *lv_abs(i,j)*(1.0_dp-dt*lv_abs(i,j)) *phi_y *dqy(i,j)

!         if (1.0_dp-dt*lv_abs(i,j) < 0.0_dp) then
!           write(6, *) "Check limiter CFL conditions: "
!           write(6, "(a,"//dp_chr//",a,"//dp_chr//",a,"//dp_chr//")") "dt = ", dt, "; lv_abs(i,j) = ", lv_abs(i,j), "; dt*lv_abs(i,j) = ", dt*lv_abs(i,j)
!         end if
      end do
    end do

  end subroutine assemble_G_limiter
  
  
  subroutine assemble_F_limiter2(F, dqx, lu_abs, phi_x, dt)
    real(kind=dp), dimension(:,:), intent(inout) :: F
    real(kind=dp), dimension(:,:), intent(in) :: dqx, lu_abs, phi_x
    real(kind=dp), intent(in) :: dt

    integer :: i, j

    ! LeVeque p119
    !! u(i,j) aligned with F(i,j) & dqx(i,j); v(i,j) aligned with G(i,j) & dqy(i,j)
    do j = 1, size(F, 2)
      do i = 2, size(F, 1)-1
        F(i,j) = F(i,j) + 0.5_dp*lu_abs(i,j)*(1.0_dp-dt*lu_abs(i,j))*phi_x(i,j)*dqx(i,j)
      end do
    end do

  end subroutine assemble_F_limiter2
  
!   subroutine assemble_F_limiter3(F, q, uv_fld, dt)
!     real(kind=dp), dimension(:,:), intent(inout) :: F
!     type(field), intent(in) :: q
!     type(uv_signed), intent(in) :: uv_fld
! 
!     real(kind=dp), intent(in) :: dt
! 
!     integer :: i, j
!     
!     lu_sgn => uv_fld%u_sgn
!     lu_abs => uv_fld%u_abs
! 
!     ! LeVeque p119
!     !! u(i,j) aligned with F(i,j) & dqx(i,j); v(i,j) aligned with G(i,j) & dqy(i,j)
!     do j = 1, size(F, 2)
!       do i = 2, size(F, 1)-1
!         F(i,j) = F(i,j) + 0.5_dp*lu_abs(i,j)*(1.0_dp-dt*lu_abs(i,j))*phi_x(i,j)*dqx(i,j)
!       end do
!     end do
! 
!   end subroutine assemble_F_limiter3
  
  subroutine assemble_G_limiter2(G, dqy, lv_abs, phi_y, dt)
    real(kind=dp), dimension(:,:), intent(inout) :: G
    real(kind=dp), dimension(:,:), intent(in) :: dqy, lv_abs, phi_y
    real(kind=dp), intent(in) :: dt

    integer :: i, j

    ! LeVeque p119
    !! u(i,j) aligned with F(i,j) & dqx(i,j); v(i,j) aligned with G(i,j) & dqy(i,j)
    do j = 2, size(G, 2)-1
      do i = 1, size(G, 1)
        G(i,j) = G(i,j) + 0.5_dp*lv_abs(i,j)*(1.0_dp-dt*lv_abs(i,j))*phi_y(i,j)*dqy(i,j)
      end do
    end do

  end subroutine assemble_G_limiter2
  
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
    integer :: i, j

    ! Set flux at domain boundaries to be zero
    do j = 1, size(F, 2)
      F(1,j) = 0.0_dp
      F(size(F, 1),j) = 0.0_dp
    end do

    do i = 1, size(G, 1)
      G(i,1) = 0.0_dp
      G(i,size(G, 2)) = 0.0_dp
    end do

  end subroutine imposeBC_FG

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
    ! Assumed zero flux BC
    type(field), intent(out) :: dqdt
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
        ldqdt(i+gl,j+gl) = -( K11e(i+1,j) + K11e(i,j) + K22e(i,j+1) + K22e(i,j) ) *lq(i+gl, j+gl) &
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
    
    ! Ensure no floating point issues
    call imposeBC(dqdt)

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
    real(kind = dp) :: x, y

    lq => q%data
    ldqdt => dqdt%data
    gl = q%glayer
    
    do j = 1, q%n
      do i = 1, q%m
        x = fld_x(i, q%m, q%type_id)
        y = fld_x(j, q%n, q%type_id)
        ldqdt(i+gl,j+gl) = lambda*(q_stat(x,y)-lq(i+gl,j+gl)) + S_rhs_stat(x,y)
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

  subroutine timestep_heun_T_ops(q, dt, T_ops)
    ! Heun's method
    type(field), intent(inout) :: q
    type(T_operator), dimension(:, :), pointer, intent(in) :: T_ops

    real(kind = dp), intent(in) :: dt

    real(kind = dp), dimension(:,:), allocatable :: F, G

    type(field) :: dqdt
    type(field) :: qc
    type(field) :: dqdtc

    call allocate(dqdt, q%m, q%n, 'dqdt', q%glayer)
    call allocate(qc, q%m, q%n, 'qc', q%glayer)
    call allocate(dqdtc, q%m, q%n, 'dqdtc', q%glayer)
    allocate(F(q%m+1, q%n))
    allocate(G(q%m, q%n+1))

    ! Evaluate dqdt and copy q to qc
    call zeros(dqdt)
    call assemble_Flux_MPFA(F, G, q, T_ops)
    call imposeBC(F, G)

    call eval_dqdt(dqdt, q, F, G)
    call set(qc, q)

    ! Prediction step
    call timestep_fe(qc, dqdt, dt)
    call imposeBC(qc)

    ! Correction step
    call zeros(dqdtc)
    call assemble_Flux_MPFA(F, G, qc, T_ops)
    call imposeBC(F, G)

    call eval_dqdt(dqdtc, qc, F, G)
    call addto(dqdt, dqdtc, 1.0_dp)

    ! Full timestepPing
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    call imposeBC(q)

    deallocate(F)
    deallocate(G)

    call deallocate(dqdt)
    call deallocate(qc)
    call deallocate(dqdtc)

  end subroutine timestep_heun_T_ops

  subroutine timestep_heun_Kflux(q, dt, K_e)
    ! Heun's method
    type(field), intent(inout) :: q
    type(K_flux), intent(in) :: K_e

    real(kind = dp), intent(in) :: dt

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
    call eval_dqdt(dqdt, q, K11e, K12e, K21e, K22e)
    call set(qc, q)

    ! Prediction step
    call timestep_fe(qc, dqdt, dt)
    call imposeBC(qc)

    ! Correction step
    call zeros(dqdtc)
    call eval_dqdt(dqdtc, qc, K11e, K12e, K21e, K22e)
    call addto(dqdt, dqdtc, 1.0_dp)

    ! Full timestepPing
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    call imposeBC(q)

    call deallocate(dqdt)
    call deallocate(qc)
    call deallocate(dqdtc)

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
    call eval_dqdt(dqdt, q, K11e, K12e, K21e, K22e)
    call timestep_fe(q, dqdt, dt)
    call imposeBC(q)

    call deallocate(dqdt)

  end subroutine timestep_fe_Kflux
  
  
  subroutine eval_dqdt_Kfull_scr_stat(dqdt, q, K11e, K12e, K21e, K22e, lambda)
    ! Heun's method
    type(field), intent(inout) :: dqdt, q
    real(kind = dp), dimension(:,:), intent(in) :: K11e, K12e, K21e, K22e

    real(kind = dp), intent(in) :: lambda

    type(field) :: dqdt_S
    
    call allocate(dqdt_S, q%m, q%n, 'dqdt_S', glayer=q%glayer, type_id=q%type_id)

    call zeros(dqdt)
    call eval_dqdt(dqdt, q, K11e, K12e, K21e, K22e) ! From diffusion

    call zeros(dqdt_S)
    call eval_dqdt_scr_stat(dqdt_S, q, lambda) ! From source term

    call addto(dqdt, dqdt_S, 1.0_dp)
    
    call deallocate(dqdt_S)
    
  end subroutine eval_dqdt_Kfull_scr_stat
  
  subroutine eval_dqdt_Kfull_scr_nonstat(dqdt, q, t, K11e, K12e, K21e, K22e)
    ! Heun's method
    type(field), intent(inout) :: dqdt, q
    real(kind = dp), dimension(:,:), intent(in) :: K11e, K12e, K21e, K22e

    real(kind = dp), intent(in) :: t

    type(field) :: dqdt_S
    
    call allocate(dqdt_S, q%m, q%n, 'dqdt_S', glayer=q%glayer, type_id=q%type_id)

    call zeros(dqdt)
    call eval_dqdt(dqdt, q, K11e, K12e, K21e, K22e) ! From diffusion

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

    ! Full timestepPing
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
    call zeros(dqdt)
    call eval_dqdt_Kfull_scr_nonstat(dqdt, q, t, K11e, K12e, K21e, K22e)
    call set(qc, q)

    ! Prediction step
    call timestep_fe(qc, dqdt, dt)
    call imposeBC(qc)

    ! Correction step
    call zeros(dqdtc)
    call eval_dqdt_Kfull_scr_nonstat(dqdt, q, t, K11e, K12e, K21e, K22e)
    call addto(dqdt, dqdtc, 1.0_dp)

    ! Full timestepPing
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    call imposeBC(q)
    
    call deallocate(dqdt)
    call deallocate(qc)
    call deallocate(dqdtc)

  end subroutine timestep_heun_Kflux_Src_nonstat
  
  
  subroutine timestep_LMT(q, dt, uv_fld)
    ! Heun's method
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

    ! DCU
    F = 0.0_dp
    G = 0.0_dp
    call assemble_FG_DCU(F,G, q, uv_fld)

    ! CTU
    call dx_field(dqx, q)
    call dy_field(dqy, q)
    call assemble_FG_CTU(F,G, dqx, dqy, uv_fld, dt)
    call assemble_FG_limiter(F,G, dqx, dqy, uv_fld, dt)  ! Apply limiter

    call imposeBC(F, G)

    call eval_dqdt(dqdt, q, F, G)

    ! Prediction step
    call timestep_fe(q, dqdt, dt)
    call imposeBC(q)

    ! Deallocate
    call deallocate(dqdt)

    deallocate(dqx)
    deallocate(dqy)
    deallocate(F)
    deallocate(G)

  end subroutine timestep_LMT
  
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
    
    ! Strang splitting: x, y direction
    ! XYYX
    G = 0.0_dp
    call dx_field(dqx, q)
    call assemble_F_DCU(F, q, uv_fld)
    call assemble_F_limiter(F, dqx, uv_fld, 0.5_dp*dt)  ! Apply limiter
    
    call imposeBC(F, G)
    call eval_dqdt(dqdt, q, F, G)
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    call imposeBC(q)

    F = 0.0_dp
    call dy_field(dqy, q)
    call assemble_G_DCU(G, q, uv_fld)
    call assemble_G_limiter(G, dqy, uv_fld, 0.5_dp*dt)  ! Apply limiter
    
    call imposeBC(F, G)
    call eval_dqdt(dqdt, q, F, G)
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    call imposeBC(q)

    F = 0.0_dp
    call dy_field(dqy, q)
    call assemble_G_DCU(G, q, uv_fld)
    call assemble_G_limiter(G, dqy, uv_fld, 0.5_dp*dt)  ! Apply limiter
    
    call imposeBC(F, G)
    call eval_dqdt(dqdt, q, F, G)
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    call imposeBC(q)

    G = 0.0_dp
    call dx_field(dqx, q)
    call assemble_F_DCU(F, q, uv_fld)
    call assemble_F_limiter(F, dqx, uv_fld, 0.5_dp*dt)  ! Apply limiter
    
    call imposeBC(F, G)
    call eval_dqdt(dqdt, q, F, G)
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    call imposeBC(q)

    ! Deallocate
    call deallocate(dqdt)

    deallocate(dqx)
    deallocate(dqy)
    deallocate(F)
    deallocate(G)

  end subroutine timestep_LMT_1DS
  
  subroutine timestep_LMT_1DS2(q, dt, uv_fld)
    type(field), intent(inout) :: q
    type(uv_signed), intent(in) :: uv_fld

    real(kind = dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), allocatable :: F, G, dqx, dqy, phi_x, phi_y

    type(field) :: dqdt
    
    if (q%glayer .ne. 1) then
      call abort_handle("ghost layer not equal 1", __FILE__, __LINE__)
    end if
    
    call allocate(dqdt, q%m, q%n, 'dqdt', glayer=q%glayer, type_id=q%type_id)

    allocate(F(q%m+1, q%n))
    allocate(G(q%m, q%n+1))
    allocate(dqx(q%m+1, q%n))
    allocate(dqy(q%m, q%n+1))
!     allocate(phi_x(q%m+1, q%n))
!     allocate(phi_y(q%m, q%n+1))

    ! F, G at boundaries will not be touched
    F = 0.0_dp
    G = 0.0_dp

    ! Strang splitting: x, y direction - XYYX
    ! X
    call assemble_F_DCU2(F, q, uv_fld%up, uv_fld%um)
    call dx_field(dqx, q)
!     call assemble_phi_x(phi_x, dqx, uv_fld%u_sgn)
!     call assemble_F_limiter2(F, dqx, uv_fld%u_abs, phi_x, 0.5_dp*dt)  ! Apply limiter
    call assemble_F_limiter(F, dqx, uv_fld, 0.5_dp*dt)  ! Apply limiter
    
    call eval_dqdt_F(dqdt, q, F)
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    call imposeBC(q)

    ! Y
    call assemble_G_DCU2(G, q, uv_fld%vp, uv_fld%vm)
    call dy_field(dqy, q)
!     call assemble_phi_y(phi_y, dqy, uv_fld%v_sgn)
!     call assemble_G_limiter2(G, dqy, uv_fld%v_abs, phi_y, 0.5_dp*dt)  ! Apply limiter 
    call assemble_G_limiter(G, dqy, uv_fld, 0.5_dp*dt)  ! Apply limiter

    call eval_dqdt_G(dqdt, q, G)
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    call imposeBC(q)

    ! Y
    call assemble_G_DCU2(G, q, uv_fld%vp, uv_fld%vm)
    call dy_field(dqy, q)
!     call assemble_phi_y(phi_y, dqy, uv_fld%v_sgn)
!     call assemble_G_limiter2(G, dqy, uv_fld%v_abs, phi_y, 0.5_dp*dt)  ! Apply limiter 
    call assemble_G_limiter(G, dqy, uv_fld, 0.5_dp*dt)  ! Apply limiter

    call eval_dqdt_G(dqdt, q, G)
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    call imposeBC(q)
    
    ! X
    call assemble_F_DCU2(F, q, uv_fld%up, uv_fld%um)
    call dx_field(dqx, q)
!     call assemble_phi_x(phi_x, dqx, uv_fld%u_sgn)
!     call assemble_F_limiter2(F, dqx, uv_fld%u_abs, phi_x, 0.5_dp*dt)  ! Apply limiter
    call assemble_F_limiter(F, dqx, uv_fld, 0.5_dp*dt)  ! Apply limiter
    
    call eval_dqdt_F(dqdt, q, F)
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    call imposeBC(q)

    ! Deallocate
    call deallocate(dqdt)

    deallocate(dqx)
    deallocate(dqy)
!     deallocate(phi_x)
!     deallocate(phi_y)
    deallocate(F)
    deallocate(G)

  end subroutine timestep_LMT_1DS2

  subroutine print_info_timestepPing(psi, dt)
    type(field), intent(in) :: psi
    real(kind=dp), intent(in) :: dt
    
    real(kind=dp) :: uv_max

    real(kind=dp), dimension(:, :), allocatable :: u, v
    
    allocate(u(psi%m+1, psi%n))
    allocate(v(psi%m, psi%n+1))

    call psi_to_uv(u, v, psi)
    uv_max = maxmax(dabs(u), dabs(v))
    
    deallocate(u)
    deallocate(v)
    
    write(6, "(a)") ""
    write(6, "(a, a)") " Time stepPing: ", trim(psi%name)
    write(6, "(a,"//dp_chr//",a,i0,a,"//dp_chr//")") "|  uv_max = ", uv_max, "  dx = ", 1, " and dt = ", dt
    write(6, "(a,"//dp_chr//")") "|  Courant number = ", CFL(uv_max, dt, 1.0_dp)

  end subroutine print_info_timestepPing

  subroutine reset_timestep_timers()
    call reset(timestep_timer)
    call reset(assemble_timer)
  end subroutine reset_timestep_timers

  subroutine print_timestep_timers()
    write(6, "(a)") "Timestep timers:"
    call print(timestep_timer, "Timestep increment", prefix = "  ")
    call print(assemble_timer, "Assemble increment", prefix = "  ")
  end subroutine print_timestep_timers

  subroutine advdiff_qC(q, psi, K11, K22, K12, T, nts)
    type(field), intent(inout) :: q
    type(field), intent(in) :: psi, K11, K22, K12
    real(kind=dp), intent(in) :: T
    integer, intent(in) :: nts

    integer :: i
    type(uv_signed) :: uv_fld
    type(K_flux) :: K_e
    real(kind=dp) :: dt

    call start(assemble_timer)

    call allocate(uv_fld, psi)
    call allocate(K_e, K11, K22, K12)
    call imposeBC(K_e)

    call stop(assemble_timer)

    call start(timestep_timer)

    dt = T/nts
    do i = 1, nts
      ! Strang splitting
      call timestep_fe_Kflux(q, 0.5_dp*dt, K_e)
      call timestep_LMT(q, dt, uv_fld)
      call timestep_fe_Kflux(q, 0.5_dp*dt, K_e)
    end do

    call stop(timestep_timer)

    call start(assemble_timer)

    call deallocate(uv_fld)
    call deallocate(K_e)

    call stop(assemble_timer)

  end subroutine advdiff_qC

  subroutine advdiff_qMPFA(q, psi, K11, K22, K12, T, nts)
    type(field), intent(inout) :: q
    type(field), intent(in) :: psi, K11, K22, K12
    real(kind=dp), intent(in) :: T
    integer, intent(in) :: nts

    integer :: i
    type(T_operator), dimension(:,:), pointer :: T_ops
    type(uv_signed) :: uv_fld
    real(kind=dp) :: dt

    call start(assemble_timer)

    call allocate(uv_fld, psi)
    call allocate(T_ops, K11%m, K11%n)

    call assemble_T_MPFA(T_ops, K11, K22, K12)

    call stop(assemble_timer)

    call start(timestep_timer)

    dt = T/nts
    do i = 1, nts
      call timestep_heun_T_ops(q, dt, T_ops)
      call timestep_LMT(q, dt, uv_fld)
    end do

    call stop(timestep_timer)

    call start(assemble_timer)

    call deallocate(uv_fld)
    call deallocate(T_ops)

    call stop(assemble_timer)

  end subroutine advdiff_qMPFA

end module advdiff_timestep
