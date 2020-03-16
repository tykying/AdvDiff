module advdiff_unittest
  use advdiff_precision
  use advdiff_debug
  use advdiff_timing
  use advdiff_io
  use advdiff_field
  use advdiff_timestep
  use advdiff_trajdata
  use advdiff_complib
  use advdiff_inference

  implicit none
  
  public :: unittest_comptools, unittest_trajdata, unittest_solver_properties
  public :: unittest_operator_order, unittest_BCflux_order, unittest_solver_convergence
  public :: unittest_IO, unittest_timer
  public :: unittest_FPsolver_INPUT, unittest_interpolation
  public :: try_evaluate_logPost, evaluate_numerical_diffusion
  
  integer, dimension(4) :: SEED = (/ 314, 159, 26, 535 /)

contains
#include "advdiff_configuration.h"
  subroutine unittest_comptools()
    real(kind=dp) :: mean, stdv, mu, sigma, gamma, sample
    real(kind=dp), dimension(:), allocatable :: stdnRV, UniRV
    real(kind=dp), dimension(2,2) :: K, Krt
    integer :: i, n, m, r
    real(kind=dp) :: sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy
    
    real(kind=dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    real(kind=dp), parameter :: eps = 1.0D-14
    
    ! Mod and division
    n = 11
    i = 3
    
    m = n/i
    r = mod(n, i)
    write(6, *) "n = ", n, "i = ",  i
    write(6, *) "n/i = ", n/i
    write(6, *) "mod(n, i) = ", r
    
    n = 9
    i = 11
    
    m = n/i
    r = mod(n, i)
    write(6, *) "n = ", n, "i = ",  i
    write(6, *) "n/i = ", n/i
    write(6, *) "mod(n, i) = ", r

    n = 11
    i = 11
    
    m = n/i
    r = mod(n, i)
    write(6, *) "n = ", n, "i = ",  i
    write(6, *) "n/i = ", n/i
    write(6, *) "mod(n, i) = ", r
    
    
    ! randn generator
    mean = 2.0_dp
    stdv = 5.0_dp
    n = 1000000

    ! Test randn()
    mu = 0.0_dp
    sigma = 0.0_dp
    gamma = 0.0_dp
    
    allocate(stdnRV(n))
    call randn(stdnRV, SEED)
    
    do i = 1, n
      sample = mean + stdv * stdnRV(i)
!      sample = 0.0_dp
      
      mu = mu + sample
      sigma = sigma + sample*sample
      gamma = gamma + sample*sample*sample
    end do

    mu = mu/n
    sigma = dsqrt(sigma/n - mu*mu)
    gamma = (gamma/n - 3.0_dp*mu*sigma*sigma - mu*mu*mu)/(sigma*sigma*sigma)

    deallocate(stdnRV)
    
    write(6, "(a, i0)") " Testing standard normal random number generator: n = ", n
    write(6, "(a,"//dp_chr//","//dp_chr//","//dp_chr//")") "| Truth = ", mean, stdv, 0.0_dp
    write(6, "(a,"//dp_chr//","//dp_chr//","//dp_chr//")") "| Computed", mu, sigma, gamma

    if ( (dabs((mean - mu)/mu) > 0.02) .or. &
       (dabs((stdv - sigma)/sigma) > 0.02) .or. &
       (dabs((0.0_dp - gamma)) > 0.01) ) then
      call abort_handle("  !! E: Random number generator", __FILE__, __LINE__)
    end if
    
    ! randu
    mean = 0.5_dp
    stdv = dsqrt(1.0_dp/12.0_dp)
    n = 1000000

    ! Test randn()
    mu = 0.0_dp
    sigma = 0.0_dp
    gamma = 0.0_dp
    
    allocate(UniRV(n))
    call randu(UniRV, SEED)
    
    do i = 1, n
      sample = UniRV(i)
!      sample = 0.0_dp
      
      mu = mu + sample
      sigma = sigma + sample*sample
      gamma = gamma + sample*sample*sample
    end do

    mu = mu/n
    sigma = dsqrt(sigma/n - mu*mu)
    gamma = (gamma/n - 3.0_dp*mu*sigma*sigma - mu*mu*mu)/(sigma*sigma*sigma)

    deallocate(UniRV)
    
    write(6, "(a, i0)") " Testing uniform(0, 1) random number generator: n = ", n
    write(6, "(a,"//dp_chr//","//dp_chr//","//dp_chr//")") "| Truth = ", mean, stdv, 0.0_dp
    write(6, "(a,"//dp_chr//","//dp_chr//","//dp_chr//")") "| Computed", mu, sigma, gamma

    if ( (dabs((mean - mu)/mu) > 0.02) .or. &
       (dabs((stdv - sigma)/sigma) > 0.02) .or. &
       (dabs((0.0_dp - gamma)) > 0.01) ) then
      call abort_handle("  !! E: Random number generator", __FILE__, __LINE__)
    end if
    
    
    write(6, "(a)") "  -- P: Passed simple random number generator test"

    ! sqrt of a matrix
    K = reshape((/ 5.0_dp, -5.0_dp, -5.0_dp, 10.0_dp /), shape(K))

    call sqrt_matrix(Krt, K)

   write(6, "(a)") ""
   write(6, "(a)") " Testing sqrt of a matrix: "
   write(6, "(a, "//dp_chr//","//dp_chr//","//dp_chr//")") "|  K = ", K(1,1), K(1,2), K(2,2)
   write(6, "(a, "//dp_chr//","//dp_chr//","//dp_chr//")") "|  Krt = ", Krt(1,1), Krt(1,2), Krt(2,2)

    if (dabs(Krt(1,1) - 2.0_dp) > eps) then
      call abort_handle("Krt(1,1):", __FILE__, __LINE__)
    end if
    if (dabs(Krt(1,2) + 1.0_dp) > eps) then
      call abort_handle("Krt(1,2):", __FILE__, __LINE__)
    end if
    if (dabs(Krt(2,2) - 3.0_dp) > eps) then
      call abort_handle("Krt(2,2):", __FILE__, __LINE__)
    end if
    write(6, "(a)") "  -- P: Passed sqrt(Matrix) test"

    ! Test case - a
    sigma1_sq = 2.0_dp
    sigma2_sq = 6.0_dp
    phi_K = (4.0_dp*atan(1.0_dp))/6.0_dp

    call KCanon_to_KCarte(Kxx, Kyy, Kxy, sigma1_sq, sigma2_sq, phi_K)
    if ((abs(Kxx-3.0_dp) > eps) .or. (abs(Kyy-5.0_dp) > eps) .or. &
         (abs(Kxy+sqrt(3.0_dp)) > eps) ) then
       call abort_handle("KCanon_to_KCarte:", __FILE__, __LINE__)
    end if

    ! Test case - b
    Kxx = 9.0_dp
    Kyy = 9.0_dp
    Kxy = 6.0_dp
    call KCarte_to_KCanon(sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy)
    if ((abs(sigma1_sq-15.0_dp) > eps) .or. (abs(sigma2_sq-3.0_dp) > eps) .or. &
         (abs(phi_K-datan(1.0_dp)) > eps) ) then
       call abort_handle("KCarte_to_KCanon:", __FILE__, __LINE__)
    end if
    write(6, "(a)") "  -- P: Passed KCarte <-> KCanon conversion test"
    
    write(6, "(a)") &
      "! ----------- Passed unittest_comptools ----------- !"
  end subroutine unittest_comptools

  subroutine unittest_trajdata()
    type(jumpsdat), dimension(:), allocatable :: jumps
    integer :: m, n
    type(meshdat) :: mesh

    type(field) :: q
    type(trajdat) :: traj
    integer :: nparticles, nts0, k, m_part, k_i, jump
    real(kind=dp), dimension(:,:), allocatable :: lik, lik_exact

    ! Part 2
    m = 64
    n = m
    call allocate(mesh, m, n)
    call allocate(q, mesh%m, mesh%n, 'q', type_id=2)

    ! Assign particle
    m_part = 16
    nparticles = 1
    
    nts0 = m_part*m_part+1
    call allocate(traj, nparticles, nts0)
    allocate(lik_exact(nparticles, nts0-1))
    allocate(lik(nparticles, nts0-1))
    
    do k = 2, nts0
      traj%x(1, k) = particle_x(k-1, m_part)
      traj%y(1, k) = particle_y(k-1, m_part)
      traj%t(:, k) = real(k-1, kind=dp)
    end do
    traj%x(:, 1) = traj%x(:, 2)
    traj%y(:, 1) = traj%y(:, 2)
    traj%t(:, 1) = 0.0_dp
    
    do k = 2, nts0
      lik_exact(1, k-1) = testcase_traj(traj%x(1, k), traj%y(1, k))
    end do
    
!     call print_array(traj%x, "traj%x")
!     call print_array(traj%y, "traj%y")
    
    call allocate(jumps, mesh)
    ! Convert trajectory into jumps
    call traj2jumps(jumps, traj, mesh)
    
    call assign_q_traj(q)
    
    k = 0
    do k_i = 1, size(jumps, 1)
      do jump = 1, jumps(k_i)%njumps
        k = k + 1
        lik(1, k) = eval_fldpt(q, mesh, jumps(k_i)%k_f(jump), jumps(k_i)%alpha_f(:,jump))
      end do
    end do   
    
    if (.not. iszero(lik-lik_exact)) then
      call print_array(lik, "lik")
      call print_array(lik_exact, "lik_exact")
      call abort_handle("eval_fldpt may be wrong! Note it is meant to be so near boundaries.", __FILE__, __LINE__)
    end if
    
    deallocate(lik)
    deallocate(lik_exact)
    call deallocate(q)
    call deallocate(traj)
    call deallocate(jumps)
    
    write(6, "(a)") &
      "! ----------- Passed unittest_trajdata ----------- !"
  end subroutine unittest_trajdata
  
  pure real(kind=dp) function particle_x(k,m)
    ! Output: on [0, 1]
    integer, intent(in) :: k, m
    integer :: i
    
    i = k2i(k,m)
    particle_x = (real(i, kind=dp)-0.5_dp)/real(m, kind=dp)
    
  end function particle_x
  
  pure real(kind=dp) function particle_y(k,m)
    ! Output: on [0, 1]
    integer, intent(in) :: k, m
    integer :: j
    
    j = k2j(k,m)
    particle_y = (real(j, kind=dp)-0.5_dp)/real(m, kind=dp)
    
  end function particle_y
  
  pure real(kind=dp) function testcase_traj(x, y)
    real(kind=dp), intent(in) :: x, y
    
    testcase_traj = 2.0_dp * x - 3.0_dp * y + 5.0_dp * x * y
  end function testcase_traj
  
  subroutine assign_q_traj(q)
    type(field), intent(inout) :: q
    integer :: i, j
    real(kind=dp) :: x, y
    
    do j = 1, q%n
      do i = 1, q%m
        x = (real(i, kind=dp)-0.5_dp)/real(q%m, kind=dp)
        y = (real(j, kind=dp)-0.5_dp)/real(q%n, kind=dp)
        
        q%data(i, j) = testcase_traj(x, y)
      end do
    end do
    
  end subroutine assign_q_traj

  subroutine unittest_operator_order()
    type(field) :: q
    integer :: m, n, type_id
    
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    integer :: i, j, h, h_range
    real(kind = dp) :: x, y, p_F, p_G
    
    real(kind=dp), dimension(:,:), allocatable :: dqx_F, dqy_F, dqx_G, dqy_G
    real(kind=dp), dimension(:,:), allocatable :: dqx_G_exact, dqy_F_exact
    real(kind=dp), dimension(:), allocatable :: h_arr, err_arr_F, err_arr_G
    
    h_range = 10
    allocate(err_arr_F(h_range))
    allocate(err_arr_G(h_range))
    allocate(h_arr(h_range))

    do h = 1, h_range
      m = 2**(h+2)
      n = m
      h_arr(h) = 1.0_dp/real(m, kind=dp)
      
      type_id = 2  ! Defined at cell centre; =1: defined at corners
      call allocate(q, m, n, 'q', type_id=type_id)
      
      allocate(dqx_F(q%m+1, q%n))
      allocate(dqy_F(q%m+1, q%n))
      allocate(dqx_G(q%m, q%n+1))
      allocate(dqy_G(q%m, q%n+1))
      allocate(dqx_G_exact(q%m, q%n+1))
      allocate(dqy_F_exact(q%m+1, q%n))
        
      !call zeros(q)
      do j = 1, q%n
        do i = 1, q%m
          x = fld_x(i, q%m, q%type_id)
          y = fld_x(j, q%n, q%type_id)
          q%data(i, j) = dcos(0.5_dp*Pi*x) * dcos(Pi*y)
        end do
      end do
          
      dqx_F = 999.0_dp
      dqy_F = 999.0_dp
      dqx_G = 999.0_dp
      dqy_G = 999.0_dp
      
      call dx_field(dqx_F, q)
      call dx_field_HoriEdge(dqx_G, q)
      
      call dy_field_VertEdge(dqy_F, q)
      call dy_field(dqy_G, q)
      
      ! Compare with exact solution
      do j = 1, size(dqx_G, 2)
        do i = 1, size(dqx_G, 1)
          x = fld_x(i, q%m, type_id=2)  ! type_id=2: Centre
          y = fld_x(j, q%n, type_id=1)  ! type_id=1: Corner
          dqx_G_exact(i, j) = -0.5_dp*Pi*dsin(0.5_dp*Pi*x) * dcos(Pi*y)
        end do
      end do
      
      do j = 1, size(dqy_F, 2)
        do i = 1, size(dqy_F, 1)
          x = fld_x(i, q%m, type_id=1)  ! type_id=1: Corner
          y = fld_x(j, q%n, type_id=2)  ! type_id=2: Centre
          dqy_F_exact(i, j) = -Pi*dcos(0.5_dp*Pi*x) * dsin(Pi*y)
        end do
      end do
      
      dqx_G = dqx_G*real(m, kind=dp)
      dqy_F = dqy_F*real(n, kind=dp)
      
      err_arr_G(h) = max_diff(dqx_G, dqx_G_exact)
      err_arr_F(h) = max_diff(dqy_F, dqy_F_exact)
      
      deallocate(dqy_F_exact)
      deallocate(dqx_G_exact)
      deallocate(dqy_G)
      deallocate(dqx_G)
      deallocate(dqy_F)
      deallocate(dqx_F)
      call deallocate(q)
    end do
    
    p_G = estimate_order(h_arr, err_arr_G)
    write(6, *) "  G order = ", p_G
    if (dabs(p_G-2.0_dp) > 1D-3) then
      call abort_handle("G: Not 2nd order gradient operator", __FILE__, __LINE__)
    end if
    
    
    p_F = estimate_order(h_arr, err_arr_F)
    write(6, *) "  F order = ", p_F
    if (dabs(p_F-2.0_dp) > 1D-3) then
      call abort_handle("F: Not 2nd order gradient operator", __FILE__, __LINE__)
    end if
    
    deallocate(h_arr)
    deallocate(err_arr_F)
    deallocate(err_arr_G)
    
    write(6, "(a)") &
      "! ----------- Passed unittest_operator_order ----------- !"
  end subroutine unittest_operator_order
  
  subroutine unittest_BCflux_order()
    type(field) :: q, K11, K22, K12, psi
    integer :: m, n
    
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    integer :: i, j, h, h_range
    real(kind = dp) :: x, y, p_F, p_G
    real(kind = dp), parameter ::  t = 0.0_dp
    
    real(kind=dp), dimension(:,:), allocatable :: G_exact, F_exact
    real(kind=dp), dimension(:), allocatable :: h_arr, err_arr_F, err_arr_G
    
    type(uv_signed) :: uv_fld
    type(K_grid) :: KGrid
    type(FluxGrid) :: Flux
    
    h_range = 10
    allocate(err_arr_F(h_range))
    allocate(err_arr_G(h_range))
    allocate(h_arr(h_range))

    do h = 1, h_range
      m = 2**(h+2)
      n = m
      h_arr(h) = 1.0_dp/real(m, kind=dp)
      
      ! type_id = 2 Defined at cell centre; =1: defined at corners
      call allocate(q, m, n, 'q', type_id=2)

      call allocate(psi, m, n, 'psi', type_id=1)
      call allocate(K11, m, n, 'K11', type_id=1)
      call allocate(K22, m, n, 'K22', type_id=1)
      call allocate(K12, m, n, 'K12', type_id=1)

      ! Assign fields (q = q(0))
      call testcase_fields(psi, K11, K22, K12, q)
      call assign_q(q, t)
      
      call allocate(uv_fld, psi)
      call allocate(KGrid, K11, K22, K12)
      call allocate(Flux, q%m, q%n)
      
      call dx_field(Flux%dqx_F, q)
      call dy_field_VertEdge(Flux%dqy_F, q)
      
      call dx_field_HoriEdge(Flux%dqx_G, q)
      call dy_field(Flux%dqy_G, q)
      
      ! Impose BC for the gradient terms
      call imposeBC(Flux, KGrid)
      call updateBC_S(Flux, KGrid, t)
      
      Flux%F = -(KGrid%K11e * Flux%dqx_F + KGrid%K12e * Flux%dqy_F)
      Flux%G = -(KGrid%K21e * Flux%dqx_G + KGrid%K22e * Flux%dqy_G)
      
      ! Update boundary flux
      call correctBFlux_S(Flux, t)
      
      ! Compare with exact solution
      allocate(G_exact(q%m, q%n+1))
      allocate(F_exact(q%m+1, q%n))
      
      do j = 1, size(Flux%G, 2)
        do i = 1, size(Flux%G, 1)
          x = fld_x(i, q%m, type_id=2)  ! type_id=2: Centre
          y = fld_x(j, q%n, type_id=1)  ! type_id=1: Corner
          G_exact(i, j) = -G_BC(x,y,t)
        end do
      end do
      
      do j = 1, size(Flux%F, 2)
        do i = 1, size(Flux%F, 1)
          x = fld_x(i, q%m, type_id=1)  ! type_id=1: Corner
          y = fld_x(j, q%n, type_id=2)  ! type_id=2: Centre
          F_exact(i, j) = -F_BC(x,y,t)
        end do
      end do
      
      Flux%G = Flux%G/real(m, kind=dp)  ! or n?
      Flux%F = Flux%F/real(n, kind=dp)
      
      err_arr_G(h) = max_diff(G_exact, Flux%G)
      err_arr_F(h) = max_diff(F_exact, Flux%F)
      
      call deallocate(q)

      call deallocate(psi)
      call deallocate(K11)
      call deallocate(K22)
      call deallocate(K12)
      
      call deallocate(uv_fld)
      call deallocate(KGrid)
      call deallocate(Flux)
      
      deallocate(F_exact)
      deallocate(G_exact)
    end do
    
    p_G = estimate_order(h_arr, err_arr_G)
    write(6, *) "  G order = ", p_G
    if (dabs(p_G-2.0_dp) > 1D-3) then
      call print_array(err_arr_G, "err_arr_G")
      call print_array(err_arr_F, "err_arr_F")
      call abort_handle("G: Not 2nd order flux approximation", __FILE__, __LINE__)
    end if
    
    p_F = estimate_order(h_arr, err_arr_F)
    write(6, *) "  F order = ", p_F
    if (dabs(p_F-2.0_dp) > 1D-3) then
      call abort_handle("F: Not 2nd order flux approximation", __FILE__, __LINE__)
    end if
    
    deallocate(h_arr)
    deallocate(err_arr_F)
    deallocate(err_arr_G)
    
    write(6, "(a)") &
      "! ----------- Passed unittest_BCflux_order ----------- !"
  end subroutine unittest_BCflux_order

  
  pure real(kind=dp) function max_diff(arr1, arr2)
    real(kind=dp), dimension(:,:), intent(in) :: arr1, arr2

    integer :: i, j
    real(kind=dp) :: max_err, diff_arr
    
    max_err = 0.0_dp
    
    do j = 1, size(arr1, 2)
      do i = 1, size(arr1, 1)
        diff_arr = dabs(arr1(i,j)-arr2(i,j))
        if (diff_arr .gt. max_err) then
          max_err = diff_arr
        end if
      end do
    end do
    
    max_diff = max_err

  end function max_diff 
  
  pure real(kind=dp) function estimate_order(h_arr, err_arr)
    real(kind = dp), dimension(:), intent(in) :: h_arr, err_arr
    real(kind = dp) :: p
    integer :: N, i
    
    N = size(h_arr, 1)
    p = 0.0_dp
    ! Use average of last three to estimate the order
    do i = 1, 4
      p = p + dlog(err_arr(N-i+1)/err_arr(N-i))/dlog(h_arr(N-i+1)/h_arr(N-i))
    end do
    
    estimate_order = p/4.0_dp
    
  end function estimate_order
  
  subroutine unittest_solver_properties()
    type(field) :: q, K11, K22, K12, psi
    integer :: m, n
    integer :: i, j

    real(kind=dp) :: t
    real(kind=dp) :: dt
    integer :: nts

    type(field) :: q0
    real(kind=dp) :: q_int_0, x, y
    real(kind=dp), parameter :: eps = 1D-14
    
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    
    m = 64+1
    n = 64+1
    call allocate(q, m, n, 'q', type_id=2)
    call allocate(q0, m, n, 'q0', type_id=2)
    call allocate(psi, m, n, 'psi', type_id=1)
    call allocate(K11, m, n, 'K11', type_id=1)
    call allocate(K22, m, n, 'K22', type_id=1)
    call allocate(K12, m, n, 'K12', type_id=1)

    ! Sin rotation for psi
    do j = 1, size(psi%data,2)
      do i = 1, size(psi%data,1)
        x = fld_x(i, psi%m, psi%type_id)  ! in [0, 1]
        y = fld_x(j, psi%n, psi%type_id)  ! in [0, 1]
        psi%data(i, j) = psi%m*psi%n/(Pi**2) *dsin(Pi*x)*dsin(Pi*y)
      end do
    end do
    call set(K11, 0.1_dp)
    call set(K22, 0.1_dp)
    call set(K12, 0.0_dp)
    
    ! Reset q
    call zeros(q)
    do i = -2, 2
      do j = -2, 2
        call indicator_field(q, (q%m+1)/2+i, (q%n+1)/2+j)
      end do
    end do
    call set(q0, q)
    q_int_0 = int_field(q)
    
    t = 1.0_dp * Pi
    nts = 1000
    call print_info(psi, t/nts)
    
    if (rotsym(q) .ge. eps) then
      call abort_handle(" ! FAIL @ Rotational symmetry at initialisation", __FILE__, __LINE__)
    end if
    
    call advdiff_q(q, psi, K11, K22, K12, t, nts)
    
    ! Test rotational symmetry
    if (rotsym(q) .le. eps) then
      write(6, *) " -- P: Rotational symmetry preserved --"
    else
      call print_array(q%data, "Final q")
      write(6, *) "Maximum discrepency = ", rotsym(q)
      call abort_handle(" ! FAIL @ Rotational symmetry", __FILE__, __LINE__)
    end if

    ! Test mass conservation
    if (diff_int(q, q0)/q_int_0 .le. 10.0_dp*eps) then
      write(6, *) " -- P: Total mass conserved --"
    else
      write(6, "(a,"//dp_chr//")") "Difference in mass = ", diff_int(q, q0)
      call abort_handle(" ! FAIL  @ Mass conversation", __FILE__, __LINE__)
    end if
    
    call deallocate(q)
    call deallocate(q0)
    call deallocate(psi)
    call deallocate(K11)
    call deallocate(K22)
    call deallocate(K12)
    
    ! Test Divergence consistent: constant initial condition
    m = 64
    n = 64
    call allocate(q, m, n, 'q', type_id=2)
    call allocate(q0, m, n, 'q0', type_id=2)
    
    call allocate(psi, m, n, 'psi', type_id=1)
    call allocate(K11, m, n, 'K11', type_id=1)
    call allocate(K22, m, n, 'K22', type_id=1)
    call allocate(K12, m, n, 'K12', type_id=1)

    call set(q, 1.0_dp)
    call set(q0, q)

    ! Initialise psi
    do j = 1, size(psi%data,2)
      do i = 1, size(psi%data,1)
        x = fld_x(i, psi%m, psi%type_id)  ! in [0, 1]
        y = fld_x(j, psi%n, psi%type_id)  ! in [0, 1]
        psi%data(i, j) = (x + 2.0_dp*y)*dsin(Pi*x)*dsin(Pi*2.0_dp*y)
      end do
    end do

    ! Initialise K
    do j = 1, size(K11%data,2)
      do i = 1, size(K11%data,1)
        x = fld_x(i, K11%m, K11%type_id)  ! in [0, 1]
        y = fld_x(j, K11%n, K11%type_id)  ! in [0, 1]
        K11%data(i, j) = 0.02_dp*(1.0_dp + dsin(2.0_dp*Pi*x) * dsin(2.0_dp*Pi*y))**2
        K22%data(i, j) = 0.03_dp*(1.0_dp + dsin(2.0_dp*Pi*x) * dsin(2.0_dp*Pi*y))**2
        K12%data(i, j) = 0.01_dp*(1.0_dp + dsin(2.0_dp*Pi*x) * dsin(2.0_dp*Pi*y))**2
      end do
    end do

    call advdiff_q(q, psi, K11, K22, K12, t, nts)

    call addto(q0, q, -1.0_dp)

    if (iszero(q0%data)) then
      write(6, *) " -- P: Divergence consistent --"
    else
      call print_array(q%data, 'q-q0')
      call abort_handle(" ! FAIL  @ Divergence consistent", __FILE__, __LINE__)
    end if
    
    call deallocate(q)
    call deallocate(q0)
    call deallocate(psi)
    call deallocate(K11)
    call deallocate(K22)
    call deallocate(K12)
    
    ! Test positivity
    m = 16
    n = 16
    call allocate(q, m, n, 'q', type_id=2)
    call allocate(q0, m, n, 'q0', type_id=2)
    
    call allocate(psi, m, n, 'psi', type_id=1)
    call allocate(K11, m, n, 'K11', type_id=1)
    call allocate(K22, m, n, 'K22', type_id=1)
    call allocate(K12, m, n, 'K12', type_id=1)

    ! Initialise q
    do j = 1, q%n
      do i = 1, q%m
        x = fld_x(i, q%m, q%type_id)  ! in [0, 1]
        y = fld_x(j, q%n, q%type_id)  ! in [0, 1]
        q%data(i, j) = dexp(-(((x-0.3_dp)/0.05_dp)**2 + ((y-0.5_dp)/0.15_dp)**2))
      end do
    end do
!     call print_array(q%data, "q0")
    
    ! Initialise psi
    do j = 1, size(psi%data,2)
      do i = 1, size(psi%data,1)
        x = fld_x(i, psi%m, psi%type_id)  ! in [0, 1]
        y = fld_x(j, psi%n, psi%type_id)  ! in [0, 1]
        psi%data(i, j) = (x + 2.0_dp*y)*dsin(Pi*x)*dsin(Pi*2.0_dp*y)
      end do
    end do

    ! Initialise K
    do j = 1, size(K11%data,2)
      do i = 1, size(K11%data,1)
        x = fld_x(i, K11%m, K11%type_id)  ! in [0, 1]
        y = fld_x(j, K11%n, K11%type_id)  ! in [0, 1]
        K11%data(i, j) = 0.02_dp*(1.0_dp + dsin(2.0_dp*Pi*x) * dsin(2.0_dp*Pi*y))**2
        K22%data(i, j) = 0.03_dp*(1.0_dp + dsin(2.0_dp*Pi*x) * dsin(2.0_dp*Pi*y))**2
        K12%data(i, j) = 0.00_dp*(1.0_dp + dsin(2.0_dp*Pi*x) * dsin(2.0_dp*Pi*y))**2
      end do
    end do
    
    call set(q0, q)
    t = 10.0_dp*Pi
    nts = 1000
    dt = t/nts
    
    call print_info(psi, dt)
    
    do i = 1, nts
      call advdiff_q(q, psi, K11, K22, K12, dt, 1)
#if MC0LW1 == 0
      ! Check at every step
      if (dabs(neg_int(q)) .ge. eps) then
        write(6, *) i, "-th step"
        call print_array(q%data, 'q')
        call print_info(q)
        call abort_handle(" ! FAIL  @ Preserved positivity", __FILE__, __LINE__)
      end if
#endif
    end do

    if (dabs(neg_int(q)) .lt. eps) then
      write(6, *) " -- P: Preserved positivity --"
    end if
    
    
    call deallocate(q)
    call deallocate(q0)
    call deallocate(psi)
    call deallocate(K11)
    call deallocate(K22)
    call deallocate(K12)

    
    write(6, "(a)") &
      "! ----------- Passed unittest_solver_properties ----------- !"
  end subroutine unittest_solver_properties
  
  pure real(kind=dp) function rotsym(fld)
    type(field), intent(in) :: fld
    integer :: m, n, i, j
    real(kind=dp) :: max_err, err

    m = fld%m
    n = fld%n

    max_err = 0.0_dp

    ! Test rotational symmetry
    if ((m .ne. n) .or. (mod(m,2) .ne. 1)) then
      rotsym = 1D15
    else
      do j = 1, n
        do i = 1, m
          ! Rotate by pi: x->-x; y->-y
          err = dabs(fld%data(i, j)-fld%data(m-i+1, n-j+1))
          if (err .gt. max_err) then
            max_err = err
          end if
          ! Rotate by 3pi/2: x->y;  y->-x
          err = dabs(fld%data(i, j)-fld%data(j, m-i+1))
          if (err .gt. max_err) then
            max_err = err
          end if
          ! Rotate by pi/2: x->-y;  y->x
          err = dabs(fld%data(i, j)-fld%data(n-j+1, i))
          if (err .gt. max_err) then
            max_err = err
          end if
        end do
      end do
      
      rotsym = max_err
    end if
  end function rotsym
  
  subroutine unittest_solver_convergence()
    integer :: i, m
    integer :: nts, C, testcase
    character(len = 256) :: scr_lne, filename
    
    ! Test temporal accuarcy (periodic solution)
    C = 16
#if TESTCASE == 0
    testcase = 0
#elif TESTCASE == 1
    testcase = 1
#else
    stop
#endif
    
#if MC0LW1 == 0
    write(filename, "(a, i0, a)") "./output/MC_Strang_case",testcase,".txt"
#elif MC0LW1 == 1
    write(filename, "(a, i0, a)") "./output/LW_Strang_case",testcase,".txt"
#else
    stop
#endif

    open(unit = output_unit, file = trim(filename), &
      & status = "replace", action = "write")
    write(output_unit, "(a)") "% Timestepping to periodic state"
    write(output_unit, "(a)") "% Scale dt = 1/C dx"
    write(output_unit, "(a,i0,a)") "% nts =", C,"*m"
    do i = 1, 8
      m = 2 ** (i+3)
      nts = C* m
      write(scr_lne, "(a,i0,a,i6,a,i0,a,i6,a,i0,a,"//dp_chr//")") &
          "m(", i,") =", m,"; nts(",i,") =", nts,"; error(", i,") =", compute_temporal_err(m, nts)
      write(output_unit, *) trim(scr_lne)
      
      ! Print to screen
      write(6, *) trim(scr_lne)
      flush(6)
    end do
    close(output_unit)
    
    write(6, "(a)") &
      "! ----------- Passed unittest_solver_convergence ----------- !"
  end subroutine unittest_solver_convergence
  
  real(kind=dp) function compute_temporal_err(m, nts)
    integer, intent(in) :: m, nts

    type(field) :: q, K11, K22, K12, psi
    integer :: n, i

    type(uv_signed) :: uv_fld
    type(K_grid) :: KGrid
    type(FluxGrid) :: Flux

    real(kind=dp), parameter :: Pi = 4.0_dp* datan (1.0_dp)
!    real(kind=dp), parameter :: T = 0.2_dp * 2.0_dp
    real(kind=dp), parameter :: T = 1.0_dp
   
    type(field) :: q_exact, q_err
    real(kind=dp) :: t_i, err_int, err_sup, dt
    
    n = m
    call allocate(q, m, n, 'q', type_id=2)
    call allocate(q_exact, m, n, 'q_exact', type_id=2)
    call allocate(q_err, m, n, 'q_err', type_id=2)
    call allocate(psi, m, n, 'psi', type_id=1)

    call allocate(K11, m, n, 'K11', type_id=1)
    call allocate(K22, m, n, 'K22', type_id=1)
    call allocate(K12, m, n, 'K12', type_id=1)

    ! test_id = 1: q(t) 1-periodic
    call testcase_fields(psi, K11, K22, K12, q_exact)
    
    call allocate(uv_fld, psi)
    call allocate(KGrid, K11, K22, K12)
    call allocate(Flux, q%m, q%n)
    
    ! Set inital condition to zero
    call assign_q(q, 0.0_dp)
    
    err_sup = 0.0_dp
    
    dt = T/real(nts,kind=dp)
    do i = 1, nts
      t_i = real(i-1, kind=dp)*dt
      
      ! Strang splitting ABA
      call timestep_heun_K_S(q, Flux, 0.5_dp*dt, t_i, KGrid)
      call timestep_LMT_2D(q, Flux, dt, uv_fld)
      call timestep_heun_K_S(q, Flux, 0.5_dp*dt, t_i+0.5_dp*dt, KGrid)

      call set(q_err, q)
      ! Compare q with exact stat state
      call assign_q(q_exact, t_i+dt)
      call addto(q_err, q_exact, -1.0_dp)
      ! Method 1: L2-norm
      err_int = dsqrt(int_sq_field(q_err)/(m * n))
!       ! Method 2: inf-norm
!       err_int = max_diff(q%data, q_exact%data)

      err_sup = max(err_int, err_sup)
    end do
        
    call deallocate(q)
    call deallocate(q_exact)
    call deallocate(q_err)
    call deallocate(psi)

    call deallocate(K11)
    call deallocate(K22)
    call deallocate(K12)

    call deallocate(uv_fld)
    call deallocate(KGrid)
    call deallocate(Flux)

    compute_temporal_err = err_sup
  end function compute_temporal_err
  
  subroutine testcase_fields(psi, K11, K22, K12, q)
    type(field), intent(inout) :: psi, K11, K22, K12, q
    integer :: i, j, m, n
    
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    real(kind = dp) :: x, y
    
    m = q%m
    n = q%n

    ! test case = 2-Period q
    ! psi
    if ((psi%type_id .ne. 1) .or. (K11%type_id .ne. 1))&
      call abort_handle("psi/K11 not defined at corners", __FILE__, __LINE__)
    
      
    do j = 1, size(psi%data, 2)
      do i = 1, size(psi%data, 1)
        ! Analytic formula: Defined on [0,1]^2 -> Need a prefactor to scale to [0,m] x [0,n]
        x = fld_x(i, m, psi%type_id)  ! in [0, 1]
        y = fld_x(j, n, psi%type_id)  ! in [0, 1]
#if TESTCASE == 0
        psi%data(i, j) = m*n/((2.0_dp*Pi)**2) *dsin(2.0_dp*Pi*x)*dsin(Pi*y)
        K11%data(i, j) = m*m/(1000.0_dp)* 4.0_dp*dexp(-0.5_dp*((x-0.5_dp)**2+(y-0.75_dp)**2))
        K22%data(i, j) = n*n/(1000.0_dp)* 2.0_dp*dexp(-0.5_dp*((x-0.5_dp)**2+(y-0.75_dp)**2))
        K12%data(i, j) = m*n/(1000.0_dp)* dexp(-0.5_dp*((x-0.5_dp)**2+(y-0.75_dp)**2))
#elif TESTCASE == 1
        psi%data(i, j) = 0.0_dp
        K11%data(i, j) = m*m/(1000.0_dp)* ((-y**2-2.0_dp*y+8.0_dp)*(dcos(3.0_dp*Pi*(x-y)))**2 + (y+1.0_dp)**2)
        K22%data(i, j) = n*n/(1000.0_dp)* (9.0_dp + (y**2+2.0_dp*y-8.0_dp)*(dcos(3.0_dp*Pi*(x-y)))**2)
        K12%data(i, j) = m*n/(1000.0_dp)* (-dcos(3.0_dp*Pi*(x-y))*dsin(3.0_dp*Pi*(x-y))*(y+4.0_dp)*(y-2.0_dp))
#else
        psi%data(i, j) = 0.0_dp
        K11%data(i, j) = m*m/(1000.0_dp)* 5.0_dp
        K22%data(i, j) = n*n/(1000.0_dp)* 5.0_dp
        K12%data(i, j) = m*n/(1000.0_dp)* 4.0_dp
#endif
        if (dabs(psi%data(i, j)) .le. 1D-15) psi%data(i, j) = 0.0_dp
      end do
    end do
    
    ! q0
    call assign_q(q, 0.0_dp)
    
  end subroutine testcase_fields
  
  subroutine assign_q(q, t)
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
    
  end subroutine assign_q
  
  pure real(kind=dp) function q_periodic(x, y, t)
    ! Define on [0, 1] times [0, 1]
    real(kind = dp), intent(in) :: x, y, t
    real(kind=dp), parameter :: Pi =  4.0_dp* datan (1.0_dp)
#if TESTCASE == 0
    q_periodic = 1.0_dp + dcos(Pi*t) * Pi**2 * (4.0_dp*x*(1.0_dp-y))**2*(4.0_dp*y*(1.0_dp-x))**3
#elif TESTCASE == 1
    q_periodic = 1.0_dp + dcos(Pi*t) * Pi**2 * (x**2 + 2.0_dp*y)
#else
    q_periodic = 1.0_dp + dcos(Pi*t) * Pi**2 * (x**3 + 2.0_dp*y)
#endif
  end function q_periodic
  
  pure real(kind=dp) function S_rhs(x, y, t)
    ! Define on [0, m] times [0, n]
    real(kind = dp), intent(in) :: x, y, t
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    
    real(kind = dp) :: dqt, udqx, vdqy, divKgradq
#if TESTCASE == 0
    dqt = -1024.0_dp*Pi**3*x**2*(1.0_dp-y)**2*y**3*(1.0_dp-x)**3*dsin(Pi*t)
    udqx = 1280.0_dp*dcos(Pi*y)*Pi*y**3*dsin(2.0_dp*Pi*x)*(-2.0_dp/5.0_dp+x)*(-1.0_dp+y)**2*dcos(Pi*t)*x*(-1.0_dp+x)**2
    vdqy = -512.0_dp*Pi*dcos(2.0_dp*Pi*x)*dsin(Pi*y)*x**2*y**2*(-1.0_dp+x)**3*dcos(Pi*t)*(5.0_dp*y**2-8.0_dp*y+3)
    divKgradq = 0.01_dp* ( (512.0_dp*(-1.0_dp+x))*exp(-(0.5_dp)*x**2+(0.5_dp)*x-13.0_dp/32.0_dp-(0.5_dp)*y**2+0.75_dp*y)*Pi**2 &
                * cos(Pi*t)*y*((y**3-8.0_dp/5.0_dp*(y**2)+3.0_dp/5.0_dp*y)*x**5 &
                + (-12.0_dp/5.0_dp+6.0_dp*y**4-76.0_dp/5.0_dp*(y**3)+18.0_dp/5.0_dp*(y**2)+36.0_dp/5.0_dp*y)*x**4 &
                + (24.0_dp/5.0_dp+y**5-287.0_dp/20.0_dp*(y**4)+191.0_dp/10.0_dp*(y**3)+53.0_dp/4.0_dp*(y**2)-111.0_dp/5.0_dp*y)*x**3 &
                +(-7.0_dp/5.0_dp*(y**5)-23.0_dp/4.0_dp*(y**4)+57.0_dp/2.0_dp*(y**3)-731.0_dp/20.0_dp*(y**2)+84.0_dp/5.0_dp*y-12.0_dp/5.0_dp)*x**2 &
                +(0.2_dp)*(2.0_dp*(-1.0_dp+y))*(y**3+113.0_dp/4.0_dp*(y**2)-157.0_dp/4.0_dp*y+6.0_dp)*y*x-8.0_dp*y**2*(-1.0_dp+y)**2*(0.2_dp)) )
#elif TESTCASE == 1
    dqt = -Pi**3*(x**2+2.0_dp*y)*dsin(Pi*t)
    udqx = 0.0_dp
    vdqy = 0.0_dp
    divKgradq = (1.0_dp/1000.0_dp) * dcos(Pi*t) * Pi**2 * &
                ( (12.0_dp*(y+4.0_dp)*(y-2.0_dp)*(x-1.0_dp)*Pi &
                + 20.0_dp-2.0_dp*y**2)*(dcos(3.0_dp*Pi*(x-y)))**2  &
                + (12.0_dp*(y+4.0_dp)*(y-2.0_dp)*(x+1.0_dp)*Pi-4.0_dp*x*(y+1.0_dp)) &
                * dsin(3.0_dp*Pi*(x-y))*dcos(3.0_dp*Pi*(x-y)) &
                - 6.0_dp*(y+4.0_dp)*(y-2.0_dp)*(x-1.0_dp)*Pi + 2.0_dp*(y+1.0_dp)**2 )
#else
    dqt = -Pi**3*(x**3+2.0_dp*y)*dsin(Pi*t)
    udqx = 0.0_dp
    vdqy = 0.0_dp
    divKgradq = (1.0_dp/1000.0_dp) * 30.0_dp*PI**2*x*dcos(Pi*t)
#endif
    S_rhs = dqt + udqx + vdqy - divKgradq
    
  end function S_rhs
  
  pure real(kind=dp) function F_BC(x, y, t)
    ! Define on [0, m] times [0, n]??
    real(kind = dp), intent(in) :: x, y, t
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    
#if TESTCASE == 0
    F_BC = 0.0_dp
#elif TESTCASE == 1
    F_BC = (1.0_dp/1000.0_dp) * &
        ( -2.0_dp*dcos(Pi*t)*(x*(y+4.0_dp)*(y-2.0_dp)*(dcos(3.0_dp*Pi*(x-y)))**2 &
          + dcos(3.0_dp*Pi*(x-y))*dsin(3.0_dp*Pi*(x-y))*(y+4.0_dp)*(y-2.0_dp) &
          - x*(y+1.0_dp)**2)*Pi**2 )
#else
    F_BC = (1.0_dp/1000.0_dp) * Pi**2*dcos(Pi*t)*(15.0_dp*x**2+8.0_dp)
#endif
  end function F_BC
  
  pure real(kind=dp) function G_BC(x, y, t)
    ! Define on [0, m] times [0, n]??
    real(kind = dp), intent(in) :: x, y, t
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    
#if TESTCASE == 0
    G_BC = 0.0_dp
#elif TESTCASE == 1
    G_BC = (1.0_dp/1000.0_dp) * &
        ( -2.0_dp*dcos(Pi*t)*(-9.0_dp+(-y**2-2.0_dp*y+8.0_dp)*(dcos(3.0_dp*Pi*(x-y)))**2 &
          +x*dsin(3.0_dp*Pi*(x-y))*(y+4.0_dp)*(y-2.0_dp)*dcos(3.0_dp*Pi*(x-y)))*Pi**2 )
#else
    G_BC = (1.0_dp/1000.0_dp) * (12.0_dp*x**2+10.0_dp)*Pi**2*cos(Pi*t)
#endif
  end function G_BC
  
  subroutine updateBC_S(Flux, KGrid, t)
    type(FluxGrid), intent(inout) :: Flux
    type(K_grid), intent(in):: KGrid
    real(kind = dp), intent(in) :: t

    integer :: i, j, m, n
    real(kind = dp) :: x, y
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
      y = fld_x(j, n, 1)  ! Corners
      dqx_G(1,j) = dqx_G(1,j) + 0.5_dp*(m*F_BC(0.0_dp,y,t)/K11(1,j)  )
      dqx_G(m,j) = dqx_G(m,j) + 0.5_dp*(m*F_BC(1.0_dp,y,t)/K11(m+1,j))
    end do
    
!     !! Ad-hoc: Test sensitivity
!     dqx_G(m,n) = dqx_G(m-1,n)
    
    ! Bottom & Top, Near Boundary; Use BC
    do i = 2, m
      x = fld_x(i, m, 1)  ! Corners
      dqy_F(i,1) = dqy_F(i,1) + 0.5_dp*(n*G_BC(x,0.0_dp,t)/K22(i,1)  )
      dqy_F(i,n) = dqy_F(i,n) + 0.5_dp*(n*G_BC(x,1.0_dp,t)/K22(i,n+1))
    end do
    
  end subroutine updateBC_S

  subroutine correctBFlux_S(Flux, t)
    type(FluxGrid), intent(inout) :: Flux
    real(kind = dp), intent(in) :: t

    integer :: i, j, m, n
    real(kind = dp) :: x, y
    real(kind = dp), dimension(:, :), pointer :: dqx_F, dqx_G, dqy_F, dqy_G

    dqx_F => Flux%dqx_F
    dqx_G => Flux%dqx_G
    dqy_F => Flux%dqy_F
    dqy_G => Flux%dqy_G
    
    m = size(dqy_G, 1)
    n = size(dqx_F, 2)
    
    do j = 1, n
      y = fld_x(j, n, 2)  ! centre
      Flux%F(1, j) = -m*F_BC(0.0_dp,y,t)
      Flux%F(m+1, j) = -m*F_BC(1.0_dp,y,t)
    end do

    do i = 1, m
      x = fld_x(i, m, 2)  ! centre
      Flux%G(i, 1) = -n*G_BC(x,0.0_dp,t)
      Flux%G(i, n+1) = -n*G_BC(x,1.0_dp,t)
    end do
    
  end subroutine correctBFlux_S

  
  subroutine eval_dqdt_K_S(dqdt, Flux, q, KGrid, t)
    type(field), intent(inout) :: dqdt
    type(FluxGrid), intent(inout) :: Flux
    type(field), intent(in) :: q
    type(K_grid), intent(in):: KGrid
    real(kind = dp), intent(in) :: t

    integer :: i, j
    real(kind = dp) :: x, y
    
    call dx_field(Flux%dqx_F, q)
    call dy_field_VertEdge(Flux%dqy_F, q)
    
    call dx_field_HoriEdge(Flux%dqx_G, q)
    call dy_field(Flux%dqy_G, q)
    
    ! Impose BC for the gradient terms
    call imposeBC(Flux, KGrid)
    call updateBC_S(Flux, KGrid, t)
    
    Flux%F = -(KGrid%K11e * Flux%dqx_F + KGrid%K12e * Flux%dqy_F)
    Flux%G = -(KGrid%K21e * Flux%dqx_G + KGrid%K22e * Flux%dqy_G)
    
    ! Update boundary flux
    call correctBFlux_S(Flux, t)
    
    ! Updated Diffusive flux
    call eval_dqdt_FG(dqdt, q, Flux%F, Flux%G)
    
    ! Source term
    do j = 1, q%n
      do i = 1, q%m
        x = fld_x(i, q%m, q%type_id)
        y = fld_x(j, q%n, q%type_id)
        dqdt%data(i,j) = dqdt%data(i,j) + S_rhs(x, y, t)
      end do
    end do

  end subroutine eval_dqdt_K_S
  
  subroutine timestep_heun_K_S(q, Flux, dt, t, KGrid)
    ! Heun's method
    type(field), intent(inout) :: q
    type(FluxGrid), intent(inout) :: Flux
    real(kind = dp), intent(in) :: dt, t
    type(K_grid), intent(in):: KGrid

    type(field) :: dqdt
    type(field) :: qc
    
    call allocate(dqdt, q%m, q%n, 'dqdt', type_id=q%type_id)
    call allocate(qc, q%m, q%n, 'qc', type_id=q%type_id)
    
    ! Evaluate dqdt and copy q to qc
    call zeros(dqdt)
    call eval_dqdt_K_S(dqdt, Flux, q, KGrid, t) ! diffusion + source term
    call set(qc, q)

    ! Prediction step
    call timestep_fe(qc, dqdt, dt)
    
    ! Correction step
    call eval_dqdt_K_S(dqdt, Flux, qc, KGrid, t+dt) ! diffusion + source term

    ! Full timestepping
    call timestep_fe(q, dqdt, 0.5_dp*dt)
    
    call deallocate(qc)
    call deallocate(dqdt)
  end subroutine timestep_heun_K_S
  
  subroutine unittest_IO()
    type(trajdat) :: traj
    real(kind=dp) :: t
    
    integer :: m, n, reflv
    type(meshdat) :: mesh

    type(dofdat) :: dof

    real(kind=dp), parameter :: L = 3840.0_dp*1000.0_dp
    real(kind=dp), parameter :: kappa_scale = 10000.0_dp
    real(kind=dp), parameter :: psi_scale = 100.0_dp*1000.0_dp
    real(kind=dp) :: sc
    
    ! Trajectories
    call read_header(traj, "./unittest/test", t)
    call read(traj, "./unittest/test")

    if ( (dabs(traj%x(3,2) - 1.0_dp) > 1D-14) .or. &
         (dabs(traj%y(1,2) - (-2.0_dp)) > 1D-14) .or. &
         (dabs(traj%t(3,4) - 3.0_dp) > 1D-14) ) then
      call print_array(traj%x)
      call print_array(traj%y)
      call print_array(traj%t)
      call abort_handle("!!! --F: testcase IO", __FILE__, __LINE__)
    end if
    
    call deallocate(traj)
!     
    ! m, n: DOF/control
    m = 4
    n = m
    ! reflv: Solver grid = (m*reflv, n*reflv)
    reflv = 2
    
    call allocate(mesh, m*reflv, n*reflv)
    call allocate(dof, m-1, n-1, m, n)

    ! Initialise fields
    sc = real(mesh%m,kind=dp)/L  ! mesh%m = solver grid, not DOF
    call zeros(dof%psi)

    call set(dof%K11, 0.5_dp*kappa_scale*sc)
    call set(dof%K22, 0.5_dp*kappa_scale*sc)
    call set(dof%K12, 0.1_dp*sc)

    call print_info(mesh)
    call print_info(dof)
    
    call print_array(dof%K22%data)
    call print_array(dof%K12%data)
    
    call write_theta(dof, "./unittest/test_theta", sc, 1)
    call print_info(dof)
    
    call deallocate(dof)
    
    call read_theta(dof, "./unittest/test_theta", sc, 1)
    
    call print_info(dof)
    call print_array(dof%K22%data)
    call print_array(dof%K12%data)

    call deallocate(dof)

    write(6, "(a)") &
      "! ----------- Passed unittest_IO ----------- !"
  end subroutine unittest_IO
  
  subroutine unittest_timer()
    integer :: i, j
    real(kind=dp) :: t, t2
    real(kind=dp), dimension(200000) :: dummy_t
    ! Timer
    type(timer), save :: total_timer
    
    call start(total_timer)
    
    t = 2.0D0
    !$OMP PARALLEL DO PRIVATE(j) num_threads(32)
    do i = 1, size(dummy_t)
      t2 = (1.0D0 + dlog(i*t)*3.0D0)/dlog(t)
      
      do j = 1, 200000
        t2 = t2 * (1.0D0 + dlog(i*t)*3.0D0)/dlog(t)
      end do
      
      dummy_t(i) = t2 * (1.0D0 + dlog(i*t)*3.0D0)/dlog(t)
    end do
    !$OMP END PARALLEL DO
    
    call stop(total_timer)

    write(6, *) sum(dummy_t)

    
    call print(total_timer, "total_timer")
    
    write(6, "(a)") &
      "! ----------- Passed unittest_timer ----------- !"
  end subroutine unittest_timer
  
  subroutine unittest_FPsolver_INPUT()
      use omp_lib
    type(jumpsdat), dimension(:), allocatable :: jumps
    type(trajdat) :: traj

    integer :: m, n, m_solver, m_Ind, n_Ind
    type(meshdat) :: mesh

    type(dofdat) :: dof, dof_old
    type(IndFndat) :: IndFn

    real(kind=dp), parameter :: L = 3840.0_dp*1000.0_dp
    real(kind=dp), parameter :: kappa_scale = 10000.0_dp
    real(kind=dp), parameter :: psi_scale = 100.0_dp*1000.0_dp
    real(kind=dp) :: sc, h, dt
    integer :: nts, dt_rf

    integer, parameter :: Td = 32
    character(len = *), parameter :: RunProfile = "QGM2_L1_NPART2704"
    character(len = 256) :: output_fld, Td_char, resol_param
    character(len = 256) :: output_q

    type(dofdat) :: dof_solver, dof_solver_inv
    type(field) :: q
    integer, dimension(:), allocatable :: klist
    integer :: INDk, i, k
    
    real(kind=dp), dimension(:), allocatable :: theta, theta_old
    
    integer, dimension(6) :: INDk_list
    
    ! m, n: DOF/control
    m = 16
    n = 16
    ! m_solver: solver grid
    m_solver = 64
    ! m_Ind, n_Ind: indicator functions
    m_Ind = 16
    n_Ind = m_Ind
    
    write(Td_char, "(a,i0,a)") "h", Td, "d"
    write(resol_param, "(a,i0,a,i0,a,i0)") "N",m_solver*m_solver,"_D", m*n, "_I", m_Ind*n_Ind
    write(output_fld, "(a,a,a,a,a,a,a)") "./output/", trim(resol_param), "/", trim(RunProfile), "/", trim(Td_char), "/"
    write(6, "(a, a)") "Output path: ", trim(output_fld)
    
    call allocate(mesh, m_solver, m_solver)  ! Needs to be identical in both directions
    call allocate(IndFn, m_Ind, n_Ind)

    ! Initialise fields
    sc = real(mesh%m,kind=dp)/L

    call read_theta(dof, trim(output_fld)//"theta_sine", sc, 40)
    
    call read(traj, "./trajdat/"//RunProfile//"/"//trim(Td_char))
    ! Ensure the positions are normalised
    if (.not. check_normalised_pos(traj)) then
      call abort_handle("E: particle out of range [0, 1]", __FILE__, __LINE__)
    end if

    ! Allocate array of jumps with initial positions
    call allocate(jumps, mesh)

    ! Convert trajectory into jumps
    call traj2jumps(jumps, traj, mesh)
    call deallocate(traj)

    !call print_info(jumps, mesh)
    
    ! Determine h: Assumed h = uniform; assumed mesh%m = mesh%n
    h = read_uniform_h(jumps) *real(mesh%m, kind=dp)   ! *m to rescale it to [0, mesh%m]

    call allocate(dof_solver, mesh)
    call intpl_dof_solver(dof_solver, dof, mesh)
    
    ! Test intpl_dof_solver
    call print_array(dof%psi%data)
    call print_array(dof_solver%psi%data)
    stop
    ! END Test intpl_dof_solver
    
    call allocate(dof_solver_inv, mesh)
    call reverse_dof_to_dof_inv(dof_solver_inv, dof_solver)
    
    INDk = m_Ind*m_Ind/2+2
    call INDk_to_klist(klist, INDk, IndFn, mesh)
    call allocate(q, mesh%m, mesh%n, 'q', type_id=2)
    
    !! dt = 3 hours
    dt_rf = 1
    dt = 3.0_dp*3600.0_dp/real(dt_rf, kind=dp)  ! 1.5 hours, in seconds
    nts = int((h/sc)/dt)    ! 2 hours time step
    
    !! Output video
    INDk_list = (/m_Ind*m_Ind/2+2, m_Ind*m_Ind/2+m_Ind+2, m_Ind*m_Ind/2-m_Ind+2, &
               &  m_Ind*m_Ind/2+1, m_Ind*m_Ind/2+m_Ind+1, m_Ind*m_Ind/2-m_Ind+1 /)
    
    write(6, *) INDk_list
    do k = 1, size(INDk_list, 1)
      INDk = INDk_list(k)
      call INDk_to_klist(klist, INDk, IndFn, mesh)
      
      write(6, *) "Solving indicator function ", INDk, " with h = ", h/sc/(3600.0_dp*24.0_dp), " days"
      call initialise_q(q, klist, mesh)
    
      write(output_q, "(a, i0)") "./unittest/q64/QGM2_L1/q", INDk
      call write(q, trim(output_q), 0.0_dp, 0)

      do i = 1, nts
        call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h/nts, 1)
        call write(q, trim(output_q), (h/sc)/nts*i, i)
      end do
      
      write(6, *) "Reverse - Solving indicator function ", INDk, " with h = ", h/sc/(3600.0_dp*24.0_dp), " days"
      call initialise_q(q, klist, mesh)
    
      write(output_q, "(a, i0)") "./unittest/q64/QGM2_L1/rev_q", INDk
      call write(q, trim(output_q), 0.0_dp, 0)

      do i = 1, nts
        call advdiff_q(q, dof_solver_inv%psi, dof_solver_inv%K11, dof_solver_inv%K22, dof_solver_inv%K12, h/nts, 1)
        call write(q, trim(output_q), (h/sc)/nts*i, i)
      end do
      
      flush(6)
    end do
    
    ! Test canon conversion
    call allocate(dof_old, dof%m_psi, dof%n_psi, dof%m_K, dof%n_K)
    call set(dof_old, dof)
    
    call allocate(theta, dof)
    call allocate(theta_old, dof)
    
    call convert_dof_to_canon(theta, dof)
    call convert_canon_to_dof(dof, theta)
    
    if (iszero(dof_old%psi%data-dof%psi%data)) then
      write(6, *) "P: Consistent psi"
    else
      call abort_handle("F: Inconsistent canonical conversion - psi", __FILE__, __LINE__)
    end if
    
    if (iszero(dof_old%K11%data-dof%K11%data)) then
      write(6, *) "P: Consistent K11"
    else
      call abort_handle("F: Inconsistent canonical conversion - K11", __FILE__, __LINE__)
    end if
    
    if (iszero(dof_old%K22%data-dof%K22%data)) then
      write(6, *) "P: Consistent K22"
    else
      call abort_handle("F: Inconsistent canonical conversion - K22", __FILE__, __LINE__)
    end if
    
    if (iszero(dof_old%K12%data-dof%K12%data)) then
      write(6, *) "P: Consistent K12"
    else
      call abort_handle("F: Inconsistent canonical conversion - K12", __FILE__, __LINE__)
    end if
    
    call deallocate(dof_old)
    deallocate(theta_old)
    deallocate(theta)
    call deallocate(jumps)
    call deallocate(q)
    deallocate(klist)
    call deallocate(dof_solver)
    
    call deallocate(dof)
    call deallocate(IndFn)

  end subroutine unittest_FPsolver_INPUT
  
  subroutine unittest_interpolation()
    real(kind=dp), dimension(:,:), allocatable :: in, out, exact
    integer :: m, n, Mx, My
    
    integer :: i, j, type_id
    real(kind=dp) :: x, y
    
    type_id = 1  ! Corner
    m = 4
    n = 8
    Mx = 3
    My = 5
    
    allocate(in(m+1, n+1))
    allocate(out(m*Mx+1, n*My+1))
    allocate(exact(m*Mx+1, n*My+1))
    
    ! Evaluate input
    do j = 1, size(in, 2)
      do i = 1, size(in, 1)
        x = fld_x(i, m, type_id)  ! in [0, 1]
        y = fld_x(j, n, type_id)  ! in [0, 1]
        in(i, j) = testcase_intpl(x, y)
      end do
    end do
    
    ! Perform interpolation
    call bilinear_intpl(out, in, Mx, My)

    ! Evaluate exact field
    do j = 1, size(exact, 2)
      do i = 1, size(exact, 1)
        x = fld_x(i, size(exact,1)-1, type_id)  ! in [0, 1]
        y = fld_x(j, size(exact,2)-1, type_id)  ! in [0, 1]
        exact(i, j) = testcase_intpl(x, y)
      end do
    end do
    
    if (.not. iszero(out-exact)) then
      call abort_handle("Bilinear interpolator broke.", __FILE__, __LINE__)
    end if
    
    deallocate(exact)
    deallocate(out)
    deallocate(in)    
    
    write(6, "(a)") &
      "! ----------- Passed unittest_interpolation ----------- !"
  end subroutine unittest_interpolation
  
  pure real(kind=dp) function testcase_intpl(x, y)
    real(kind=dp), intent(in) :: x, y
    
    testcase_intpl = 2.0_dp * x + 3.0_dp * y
  end function testcase_intpl
  
  subroutine try_evaluate_logPost(Td, layer, NPart, Phase, Seed_ID)
    use mpi
    integer, intent(in) :: Td, layer, NPart, Phase, Seed_ID
    
    type(jumpsdat), dimension(:), allocatable :: jumps
    type(trajdat) :: traj

    integer :: m, n, m_Ind, n_Ind
    integer, parameter :: m_solver = 64  ! solver grid
    type(meshdat) :: mesh

    type(dofdat) :: dof
    type(IndFndat) :: IndFn
    
    real(kind=dp) :: sc, h, dt
    integer :: nts
    
    real(kind=dp), parameter :: L = 3840.0_dp*1000.0_dp

    character(len = 128) :: RunProfile
    character(len = 256) :: Td_char, resol_param
    character(len = 256) :: input_fld, fld_tmp
    
#if EM_MEAN == 1
    ! Use Eulerian time-average
    type(field) :: psi_EM
    real(kind=dp) :: t_avg
#endif
    
    ! Prior
    type(priordat) :: prior_param
    
    integer :: restart_ind
    integer :: output_dn

    ! MPI
    integer :: my_id, num_procs
    integer :: ierr !, PROVIDED, !! PROVIDED: For MPI-OpenMP
    real(kind=dp) :: SlogPost_local, SlogPost_global
    
    m = 16
    n = 16
    m_Ind = 16
    n_Ind = m_Ind
   
    write(Td_char, "(a,i0,a)") "h", Td, "d"
    write(resol_param, "(a,i0,a,i0,a,i0)") "N",m_solver*m_solver,"_D", m*n, "_I", m_Ind*n_Ind
    
    output_dn = 0
    restart_ind = -1
    write(RunProfile, "(a,i0,a,i0)") "QGM2_L", layer, "_NPART", NPart
    write(fld_tmp, "(a,a,a,a,a,a,a)") "./eddie/", trim(resol_param), "/", trim(RunProfile), "/", trim(Td_char), "/"
    if (Phase .eq. 1) then
      write(input_fld, "(a,a,i0,a)") trim(fld_tmp), "Seed", Seed_ID, "/SpinUp/"
    elseif (Phase .eq. 2) then
      write(input_fld, "(a,a,i0,a)") trim(fld_tmp), "Seed", Seed_ID, "/Tuned/"
    elseif (Phase .eq. 3) then
      write(input_fld, "(a,a,i0,a)") trim(fld_tmp), "Seed", Seed_ID, "/Data/"
    end if
    
    call allocate(mesh, m_solver, m_solver)  ! Needs to be identical in both directions
    ! N.B. assumed zeros at boundaries of psi-> hence only need (m-1)*(n-1) instead of (m+1)*(n+1)

    call MPI_INIT( ierr )
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
    if (ierr .ne. MPI_SUCCESS) then 
      call abort_handle("MPI Failed", __FILE__, __LINE__) 
    end if
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    if (ierr .ne. MPI_SUCCESS) then 
      call abort_handle("MPI Failed", __FILE__, __LINE__) 
    end if
    
    call allocate(IndFn, m_Ind, n_Ind, my_id, num_procs)
    
    ! Initialise fields
    sc = real(mesh%m,kind=dp)/L     ! Needs mesh to be identical in both directions

    if (restart_ind .ne. 0) then
      call read_theta(dof, trim(input_fld)//"theta", sc, restart_ind)
    else
      call allocate(dof, m-1, n-1, m+1, n+1)
      call init_dof(dof, sc)
    end if

    call allocate(prior_param, sc)
    
    ! Read trajectory
    call read(traj, "./trajdat/"//trim(RunProfile)//"/"//trim(Td_char))

    ! Ensure the positions are normalised
    if (.not. check_normalised_pos(traj)) then
      call abort_handle("E: particle out of range [0, 1]", __FILE__, __LINE__)
    end if

    ! Allocate array of jumps with initial positions
    call allocate(jumps, mesh)
    ! Convert trajectory into jumps
    call traj2jumps(jumps, traj, mesh)
    call deallocate(traj)
    
    ! Determine h: Assumed h = uniform; assumed mesh%m = mesh%n
    h = read_uniform_h(jumps) *real(mesh%m, kind=dp)   ! *m to rescale it to [0, mesh%m]

    ! Calculate nts = number of timesteps needed in solving FP
    dt = 6.0_dp*3600.0_dp  ! 6 hours, in seconds
    nts = int((h/sc)/dt)
    
    SlogPost_local = evaluate_logPost_OMP(prior_param, jumps, IndFn, mesh, dof, h, nts)
    
    SlogPost_global = 0.0_dp
    call MPI_Reduce(SlogPost_local, SlogPost_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
             (num_procs-1), MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) then 
      call abort_handle("MPI Failed", __FILE__, __LINE__) 
    end if
    
    dof%SlogPost = SlogPost_global
    
    if (my_id .eq. (num_procs-1)) then
      write(6, "(a, "//dp_chr//")")  "logPost = ", dof%SlogPost
    end if
    
    ! Release memory
    call deallocate(prior_param)
    
    call deallocate(dof)
    call deallocate(IndFn)
    call deallocate(jumps)

    call MPI_FINALIZE(ierr)

  end subroutine try_evaluate_logPost
  
  subroutine evaluate_numerical_diffusion()
    type(field) :: q, K11, K22, K12, psi
    integer :: m, n
    integer :: i, j, k, k10

    real(kind=dp) :: t
    real(kind=dp) :: dt
    integer :: nts

    type(field) :: q0
    real(kind=dp) :: q_int_0, x, y
    real(kind=dp), parameter :: eps = 1D-14
    
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    type(meshdat) :: mesh

    real(kind=dp), dimension(3) :: Gauss_param    !=(mean_x, mean_y, sigma)
    real(kind=dp), parameter :: Gauss_sigma = 1.0_dp/32.0_dp
    real(kind=dp) :: q_sq_int, dq_sq_int, q_sq_int_old, K_inst, K_avg
    real(kind=dp), dimension(:), allocatable :: K_small_arr, K_avg_arr
    integer :: K_avg_counter, K_nscale
    
    m = 64+1
    n = 64+1
!     m = 8+1
!     n = 8+1
    call allocate(mesh, m, n)  ! Needs to be identical in both directions
    call allocate(q, m, n, 'q', type_id=2)
    call allocate(q0, m, n, 'q0', type_id=2)
    call allocate(psi, m, n, 'psi', type_id=1)
    call allocate(K11, m, n, 'K11', type_id=1)
    call allocate(K22, m, n, 'K22', type_id=1)
    call allocate(K12, m, n, 'K12', type_id=1)

    ! Rigid body rotation for psi
    do j = 1, size(psi%data,2)
      do i = 1, size(psi%data,1)
        x = fld_x(i, psi%m, psi%type_id)  ! in [0, 1]
        y = fld_x(j, psi%n, psi%type_id)  ! in [0, 1]
        psi%data(i, j) = psi%m*psi%n/(Pi**2) *dsin(Pi*x)*dsin(Pi*y)
      end do
    end do
    
    K_nscale = 10
    allocate(K_small_arr(K_nscale))
    allocate(K_avg_arr(K_nscale))
    
    do k10 = 1, K_nscale
      K_small_arr(k10) = 10.0_dp**(-k10+4)
      
      call set(K11, K_small_arr(k10))
      call set(K22, K_small_arr(k10))
      call set(K12, 0.0_dp)
      
      ! Reset q
      call zeros(q)
      
      Gauss_param = (/ 0.5_dp, 0.5_dp, Gauss_sigma /)
      call initialise_q_Gauss(q, Gauss_param, mesh)

      call set(q0, q)
      q_int_0 = int_field(q)
      
!       t = 0.001_dp * Pi
      t = 0.0000001_dp
      nts = 10000
      dt = t/nts
!       call print_info(psi, t/nts)
  !     call print_array(q%data, "q%data")
      
      K_avg = 0.0
      K_avg_counter = 0
      
      q_sq_int_old = int_grad_sq_field(q)
      do k = 1, nts
        call advdiff_q(q, psi, K11, K22, K12, dt, 1)
        dq_sq_int = int_grad_sq_field(q)
        q_sq_int = int_sq_field(q)
        
        K_inst = -(q_sq_int-q_sq_int_old)/(2.0_dp*dt*dq_sq_int)
        if (k .ge. nts/2) then
          K_avg_counter = K_avg_counter + 1
          K_avg  = K_avg + K_inst
        end if 
        
        q_sq_int_old = q_sq_int
      end do
      
      K_avg_arr(k10) = K_avg/K_avg_counter
      !write(6, *) K_avg_arr(k10), " ", K_inst
      
      ! Test mass conservation
      if (diff_int(q, q0)/q_int_0 .le. 100.0_dp*eps) then
        write(6, *) " -- P: Total mass conserved --"
      else
        write(6, "(a,"//dp_chr//")") "Difference in mass = ", diff_int(q, q0)
        call abort_handle(" ! FAIL  @ Mass conversation", __FILE__, __LINE__)
      end if
    end do
    
    call print_array(K_small_arr, "K")
    call print_array(K_avg_arr, "K_avg")

    
    deallocate(K_small_arr)
    deallocate(K_avg_arr)
    
    call deallocate(q)
    call deallocate(q0)
    call deallocate(psi)
    call deallocate(K11)
    call deallocate(K22)
    call deallocate(K12)
  
  end subroutine evaluate_numerical_diffusion
  
  pure real(kind=dp) function int_grad_sq_field(fld)
    type(field), intent(in) :: fld
    integer :: i, j

    int_grad_sq_field = 0.0_dp

    do j = 2, (fld%n-1)
      do i = 2, (fld%m-1)
        int_grad_sq_field = int_grad_sq_field + (fld%data(i+1,j)-fld%data(i-1,j))**2 &
                      + (fld%data(i,j+1)-fld%data(i,j-1))**2
      end do
    end do

  end function int_grad_sq_field
  
end module advdiff_unittest

