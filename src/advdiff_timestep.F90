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

  public :: timestep_heun_K, timestep_LMT_2D
  public :: timestep_fe, eval_dqdt_K, eval_dqdt_FG

  interface imposeBC
    module procedure imposeBC_FG, imposeBC_Grad
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
  
  subroutine assemble_F_DCU(F, q, lup, lum)
    real(kind=dp), dimension(:,:), intent(out) :: F
    type(field), intent(in) :: q
    real(kind = dp), dimension(:, :), pointer, intent(in) :: lup, lum

    integer :: i, j
    real(kind = dp), dimension(:, :), pointer :: lq

    lq => q%data

    ! Interior
    do j = 1, size(F, 2)
      do i = 2, (size(F, 1)-1)
        F(i,j) = lum(i,j)*lq(i, j) + lup(i,j)*lq(i-1, j)
      end do
    end do
    
    ! Boundary
    do j = 1, size(F, 2)
      F(1,j) = 0.0_dp
      F(size(F, 1), j) = 0.0_dp
    end do
    
  end subroutine assemble_F_DCU
  
  subroutine assemble_G_DCU(G, q, lvp, lvm)
    real(kind=dp), dimension(:,:), intent(out) :: G
    type(field), intent(in) :: q
    real(kind = dp), dimension(:, :), pointer, intent(in) :: lvp, lvm

    integer :: i, j
    real(kind = dp), dimension(:, :), pointer :: lq

    lq => q%data

    ! Interior
    do j = 2, (size(G, 2)-1)
      do i = 1, size(G, 1)
        G(i,j) = lvm(i,j)*lq(i, j) + lvp(i,j)*lq(i, j-1)
      end do
    end do
    
    ! Boundary
    do i = 1, size(G, 1)
      G(i, 1) = 0.0_dp
      G(i, size(G, 2)) = 0.0_dp
    end do

  end subroutine assemble_G_DCU
  
  subroutine assemble_FG_CTU(F, G, dqx, dqy, lup, lum, lvp, lvm, dt)
    real(kind=dp), dimension(:,:), intent(inout) :: F,G
    real(kind=dp), dimension(:,:), intent(in) :: dqx, dqy
    real(kind = dp), dimension(:, :), pointer, intent(in) :: lup, lum, lvp, lvm
    real(kind=dp), intent(in) :: dt

    integer :: i, j

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

#include "advdiff_configuration.h"
#if MC0LW1 == 0
    ! MC
    philim = max(0.0_dp, min(0.5_dp*(1.0_dp+theta), 2.0_dp, 2.0_dp*theta))
#elif MC0LW1 == 1
    ! Lax-Wendroff
    philim = 1.0_dp
#endif

  end function philim

  subroutine assemble_F_limiter(F, dqx, lu_sgn, lu_abs , dt)
    real(kind=dp), dimension(:,:), intent(inout) :: F
    real(kind=dp), dimension(:,:), intent(in) :: dqx
    integer, dimension(:, :), pointer, intent(in) :: lu_sgn
    real(kind = dp), dimension(:, :), pointer, intent(in) :: lu_abs
    real(kind=dp), intent(in) :: dt
    
    integer :: i, j
      
    ! LeVeque p119
    !! u(i,j) aligned with F(i,j) & dqx(i,j); v(i,j) aligned with G(i,j) & dqy(i,j)
    do j = 1, size(F, 2)
      do i = 2, size(F, 1)-1
        !! Strange behavior with valgrind: with the `if' on it thinks there is a leak
        if (dabs(dqx(i,j)) .gt. 1D-16) then
          F(i,j) = F(i,j) + 0.5_dp *lu_abs(i,j)*(1.0_dp-dt*lu_abs(i,j)) &
                   * philim(dqx(i-lu_sgn(i,j), j)/dqx(i,j)) *dqx(i,j)
        end if
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
        if (dabs(dqy(i,j)) .gt. 1D-16) then
          G(i,j) = G(i,j) + 0.5_dp *lv_abs(i,j)*(1.0_dp-dt*lv_abs(i,j)) &
                * philim(dqy(i, j-lv_sgn(i,j))/dqy(i,j)) *dqy(i,j)
        end if
      end do
    end do

  end subroutine assemble_G_limiter

  subroutine imposeBC_FG(F, G)
    real(kind=dp), dimension(:,:), intent(inout) :: F, G
    integer :: j, i

    ! Set flux at domain boundaries to be zero
    do j = 1, size(F, 2)
      F(1,j) = 0.0_dp
      F(size(F,1),j) = 0.0_dp
    end do

    do i = 1, size(G, 1)
      G(i,1) = 0.0_dp
      G(i,size(G,2)) = 0.0_dp
    end do
    
  end subroutine imposeBC_FG
  
  subroutine imposeBC_Grad(Flux, KGrid)
    type(FluxGrid), intent(inout) :: Flux
    type(K_grid), intent(in) :: KGrid
    integer :: i, j, m, n

    real(kind = dp), dimension(:, :), pointer :: dqx_F, dqx_G, dqy_F, dqy_G
    real(kind = dp), dimension(:, :), pointer :: K11, K22, K12

    dqx_F => Flux%dqx_F
    dqx_G => Flux%dqx_G
    dqy_F => Flux%dqy_F
    dqy_G => Flux%dqy_G

    K11 => KGrid%K11c
    K22 => KGrid%K22c
    K12 => KGrid%K12c
    
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
    
    ! Bottom & Top, Near Boundary; Use BC
    do i = 2, m
      dqy_F(i, 1) = 0.5_dp*( &
        - K12(i,1)/K22(i,1)*(1.5_dp*dqx_F(i,1) - 0.5_dp*dqx_F(i,2)) &
        + 0.5_dp* dqy_G(i,2) + 0.5_dp* dqy_G(i-1,2) )
      dqy_F(i, n) = 0.5_dp*( &
        - K12(i,n+1)/K22(i,n+1)*(1.5_dp*dqx_F(i,n) - 0.5_dp*dqx_F(i,n-1)) &
        + 0.5_dp* dqy_G(i,n) + 0.5_dp* dqy_G(i-1,n) )
    end do
    
  end subroutine imposeBC_Grad
  
  subroutine eval_dqdt_K(dqdt, Flux, q, KGrid)
    type(field), intent(inout) :: dqdt
    type(FluxGrid), intent(inout) :: Flux
    type(field), intent(in) :: q
    type(K_grid), intent(in) :: KGrid

    ! Compute gradients
    call dx_field(Flux%dqx_F, q)
    call dy_field_VertEdge(Flux%dqy_F, q)
    
    call dx_field_HoriEdge(Flux%dqx_G, q)
    call dy_field(Flux%dqy_G, q)
    
    ! Impose BC for the gradient terms
    call imposeBC(Flux, KGrid)
    
    Flux%F = -(KGrid%K11e * Flux%dqx_F + KGrid%K12e * Flux%dqy_F)
    Flux%G = -(KGrid%K21e * Flux%dqx_G + KGrid%K22e * Flux%dqy_G)

    ! Impose BC for the flux terms
    call imposeBC(Flux%F, Flux%G)
    
    ! Evaluate dqdt
    call eval_dqdt_FG(dqdt, q, Flux%F, Flux%G)
    
  end subroutine eval_dqdt_K
  
  subroutine eval_dqdt_FG(dqdt, q, F, G)
    type(field), intent(inout) :: dqdt
    type(field), intent(in) :: q
    
    real(kind=dp), dimension(:,:) :: F, G
    integer :: i, j

    do j = 1, q%n
      do i = 1, q%m
        dqdt%data(i,j) = dqdt%data(i,j) -F(i+1,j) +F(i,j) -G(i,j+1) +G(i,j)
      end do
    end do
    
  end subroutine eval_dqdt_FG

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
  
  subroutine timestep_heun_K(q, Flux, dt, KGrid)
    ! Heun's method
    type(field), intent(inout) :: q
    type(FluxGrid), intent(inout) :: Flux
    type(K_grid), intent(in) :: KGrid
    real(kind = dp), intent(in) :: dt

    type(field) :: dqdt
    type(field) :: qc

    call allocate(dqdt, q%m, q%n, 'dqdt', type_id=q%type_id)
    call allocate(qc, q%m, q%n, 'qc', type_id=q%type_id)

    ! Evaluate dqdt and copy q to qc
    call zeros(dqdt)
    call eval_dqdt_K(dqdt, Flux, q, KGrid)
    call set(qc, q)

    ! Prediction step
    call timestep_fe(qc, dqdt, dt)

    ! Correction step (add to dqdt)
    call eval_dqdt_K(dqdt, Flux, qc, KGrid)
    
    ! Full timestepping
    call timestep_fe(q, dqdt, 0.5_dp*dt)

    call deallocate(qc)
    call deallocate(dqdt)

  end subroutine timestep_heun_K
  
  subroutine timestep_LMT_2D(q, Flux, dt, uv_fld)
    type(field), intent(inout) :: q
    type(FluxGrid), intent(inout) :: Flux
    type(uv_signed), intent(in) :: uv_fld
    real(kind = dp), intent(in) :: dt

    type(field) :: dqdt
    
    call allocate(dqdt, q%m, q%n, 'dqdt', type_id=q%type_id)
    
    ! DCU + 1D limiter
    call assemble_F_DCU(Flux%F, q, uv_fld%up, uv_fld%um)
    call dx_field(Flux%dqx_F, q)
    call assemble_F_limiter(Flux%F, Flux%dqx_F, uv_fld%u_sgn, uv_fld%u_abs, dt)  ! Apply limiter

    call assemble_G_DCU(Flux%G, q, uv_fld%vp, uv_fld%vm)
    call dy_field(Flux%dqy_G, q)
    call assemble_G_limiter(Flux%G, Flux%dqy_G, uv_fld%v_sgn, uv_fld%v_abs, dt)  ! Apply limiter
    
    ! CTU ( for 2nd order )
    call assemble_FG_CTU(Flux%F, Flux%G, Flux%dqx_F, Flux%dqy_G, uv_fld%up, uv_fld%um, uv_fld%vp, uv_fld%vm, dt)

    call imposeBC(Flux%F, Flux%G)
    
    call zeros(dqdt)
    call eval_dqdt_FG(dqdt, q, Flux%F, Flux%G)
    call timestep_fe(q, dqdt, dt)

    ! Deallocate
    call deallocate(dqdt)

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
    type(K_grid) :: KGrid
    type(FluxGrid) :: Flux
    real(kind=dp) :: dt
    
    call allocate(uv_fld, psi)
    call allocate(KGrid, K11, K22, K12)
    call allocate(Flux, q%m, q%n)
    
    dt = T/nts
    do i = 1, nts
      call timestep_heun_K(q, Flux, 0.5_dp*dt, KGrid)
      call timestep_LMT_2D(q, Flux, dt, uv_fld)
      call timestep_heun_K(q, Flux, 0.5_dp*dt, KGrid)
    end do

    call deallocate(uv_fld)
    call deallocate(KGrid)
    call deallocate(Flux)

  end subroutine advdiff_q

end module advdiff_timestep
