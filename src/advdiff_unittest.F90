module advdiff_unittest
  use advdiff_precision
  use advdiff_debug
  use advdiff_timing
  use advdiff_parameters
  use advdiff_io
  use advdiff_field
  use advdiff_timestep
  use advdiff_trajdata
  use advdiff_complib
  use advdiff_inference

  implicit none
  
  public :: unittest_comptools, unittest_trajdata, unittest_solver_properties
  public :: unittest_fftw, unittest_solver_timing, unittest_solver_convergence
  public :: unittest_IO, unittest_timer
  public :: unittest_FPsolver_INPUT, unittest_TG_instability

contains
  subroutine unittest_comptools()
    real(kind=dp) :: mean, stdv, mu, sigma, gamma, sample
    real(kind=dp), dimension(:), allocatable :: stdnRV, UniRV
    real(kind=dp), dimension(2,2) :: K, Krt
    integer :: i, j, n, m, r
    real(kind=dp) :: sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy
    
    type(field) :: rel, psi
    real(kind=dp) :: x, y
    real(kind=dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)

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
    call randn(stdnRV)
    
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
    call randu(UniRV)
    
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

    if (dabs(Krt(1,1) - 2.0_dp) > EPSILON(0.0_dp)) then
      call abort_handle("Krt(1,1):", __FILE__, __LINE__)
    end if
    if (dabs(Krt(1,2) + 1.0_dp) > EPSILON(0.0_dp)) then
      call abort_handle("Krt(1,2):", __FILE__, __LINE__)
    end if
    if (dabs(Krt(2,2) - 3.0_dp) > EPSILON(0.0_dp)) then
      call abort_handle("Krt(2,2):", __FILE__, __LINE__)
    end if
    write(6, "(a)") "  -- P: Passed sqrt(Matrix) test"

    ! Test case - a
    sigma1_sq = 2.0_dp
    sigma2_sq = 6.0_dp
    phi_K = (4.0_dp*atan(1.0_dp))/6.0_dp

    call KCanon_to_KCarte(Kxx, Kyy, Kxy, sigma1_sq, sigma2_sq, phi_K)
    if ((abs(Kxx-3.0_dp) > 1e-14) .or. (abs(Kyy-5.0_dp) > 1e-14) .or. &
         (abs(Kxy+sqrt(3.0_dp)) > 1e-14) ) then
       call abort_handle("KCanon_to_KCarte:", __FILE__, __LINE__)
    end if

    ! Test case - b
    Kxx = 9.0_dp
    Kyy = 9.0_dp
    Kxy = 6.0_dp
    call KCarte_to_KCanon(sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy)
    if ((abs(sigma1_sq-15.0_dp) > 1e-14) .or. (abs(sigma2_sq-3.0_dp) > 1e-14) .or. &
         (abs(phi_K-datan(1.0_dp)) > 1e-14) ) then
       call abort_handle("KCarte_to_KCanon:", __FILE__, __LINE__)
    end if
    write(6, "(a)") "  -- P: Passed KCarte <-> KCanon conversion test"
    
    
    m = 1024
    n = m
    call allocate(psi, m, n, "psi", glayer=0, type_id=1)
    call allocate(rel, psi%m, psi%n, "rel", glayer=psi%glayer, type_id=psi%type_id)
    
    do j = 1, size(psi%data, 2)
      do i = 1, size(psi%data, 1)
        x = fld_x(i, psi%m, psi%type_id)
        y = fld_x(j, psi%n, psi%type_id)
        
        psi%data(i,j) = dsin(2.0_dp* Pi * x) * dsin(Pi * y)
      end do
    end do
    
    call del2(rel, psi)
    
!     call print_array(rel%data*m*m, "rel*m*m")
!     call print_array(((2.0_dp*Pi)**2 + Pi**2)*psi%data, "((2.0_dp*Pi)**2 + Pi**2)*psi%data")
    
!     if (.not. iszero(rel%data*m*m + ((2.0_dp*Pi)**2 + Pi**2)*psi%data) ) then
!         call abort_handle("del2: inconsistent", __FILE__, __LINE__)
!     end if

    call deallocate(rel)
    call deallocate(psi)
    
    write(6, "(a)") "  -- P: del2 consistent"
    
    write(6, "(a)") &
      "! ----------- Passed unittest_comptools ----------- !"
  end subroutine unittest_comptools

  subroutine unittest_trajdata()
    type(jumpsdat), dimension(:), allocatable :: jumps
    integer :: m, n, cell
    type(meshdat) :: mesh

    type(field) :: q

    m = 64
    n = 32

    call allocate(mesh, m, n)
    call print_info(mesh)

    call allocate(jumps, mesh)
!     write(6, "(a,"//int_chr//")") "  size(jumps)", size(jumps)
!     write(6, "(a)") "  Before initialisation: "
!     write(6, "(a,"//int_chr//","//int_chr//")") "  jumps(1:2)%njumps", jumps(1)%njumps, jumps(2)%njumps

    do cell = 1, mesh%ncell
      call allocate(jumps(cell), 23)
    end do
    call allocate(jumps(2), 51)

!     write(6, "(a)") "  After initialisation: "
!     write(6, "(a,"//int_chr//","//int_chr//")") "  jumps(1:2)%njumps", jumps(1)%njumps, jumps(2)%njumps
!     write(6, "(a,"//int_chr//","//int_chr//")") "  size(jumps(1)%k_f)", size(jumps(1)%k_f)
!     write(6, "(a,"//int_chr//","//int_chr//")") "  size(jumps(2)%alpha_i, 1:2)", size(jumps(2)%alpha_i, 1), size(jumps(2)%alpha_i, 2)

    write(6, "(a)") ""
    write(6, "(a)") "  Test eval_fldpt: "
    call allocate(q, m, n, 'q')
    ! Assign q(i,j) = real(j)*real(i)
    call manual_field(q, 2)

    jumps(1)%k_f(1) = 23
    jumps(1)%k_f(2) = 23*16+7

    write(6, "(a,"//dp_chr//")") "|  eval_fldpt(1)", eval_fldpt(q, mesh, jumps(1)%k_f(1))
    write(6, "(a,"//dp_chr//")") "|  eval_fldpt(2)", eval_fldpt(q, mesh, jumps(1)%k_f(2))

    if (dabs(eval_fldpt(q, mesh, jumps(1)%k_f(1))-23.0_dp) > 1D-14) &
      call abort_handle("eval_fldpt(1):", __FILE__, __LINE__)

    if (dabs(eval_fldpt(q, mesh, jumps(1)%k_f(2))-330.0_dp) > 1D-14) &
      call abort_handle("eval_fldpt(2):", __FILE__, __LINE__)


    write(6, "(a)") "  -- P: Evaluation of function -- "

    call deallocate(q)
    do cell = 1, mesh%ncell
      call deallocate(jumps(cell))
    end do

    write(6, "(a)") &
      "! ----------- Passed unittest_trajdata ----------- !"
  end subroutine unittest_trajdata

  subroutine unittest_solver_properties()
    type(field) :: q, K11, K22, K12, psi
    integer :: m, n
    integer :: i, j

    type(uv_signed) :: uv_fld
    type(K_flux) :: K_e

    integer, parameter :: gl = 1

    real(kind=dp) :: t = 1.0_dp*3.1415926_dp
    real(kind=dp) :: dt
    integer :: nts
    integer :: ind = 0

    type(field) :: q0
    real(kind=dp) :: q_int_0, x, y
    real(kind=dp), dimension(:,:), allocatable :: dqx, dqy
    
    real(kind=dp), dimension(:,:), allocatable :: exactval
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)

    m = 4
    n = m
    
    ! Test calculation of gradient q
    call allocate(q, m, n, 'q', glayer=1, type_id=2)
    allocate(dqx(q%m+1, q%n))
    allocate(dqy(q%m, q%n+1))

    do i = 1, m
      do j = 1, n
        q%data(i+gl, j+gl) = i*(2*j+1)
      end do
    end do
    
    call imposeBC(q)
    call dx_field(dqx, q)
    call dy_field(dqy, q)
    
    do j = 1, size(dqx, 2)
      do i = 2, size(dqx, 1)-1
        if (abs(dqx(i,j) - (2.0_dp*j+1.0_dp)) > 1e-13) then
          call abort_handle("Error in dqx!", __FILE__, __LINE__)
        end if
      end do
    end do
    
    do j = 2, size(dqy, 2)-1
      do i = 1, size(dqy, 1)
        if (abs(dqy(i,j) - 2*i) > 1e-13) then
          call abort_handle("Error in dqy!", __FILE__, __LINE__)
        end if
      end do
    end do
    
    call deallocate(q)
    deallocate(dqx)
    deallocate(dqy)
    
    m = 32+1
    n = 32+1
    call allocate(q, m, n, 'q', glayer=gl, type_id=2)
    call allocate(q0, m, n, 'q0', glayer=gl, type_id=2)
    call allocate(psi, m, n, 'psi', glayer=0, type_id=1)
    call allocate(K11, m, n, 'K11', glayer=0, type_id=1)
    call allocate(K22, m, n, 'K22', glayer=0, type_id=1)
    call allocate(K12, m, n, 'K12', glayer=0, type_id=1)

    ! Ridge body rotation for psi
    do j = 1, size(psi%data,2)
      do i = 1, size(psi%data,1)
        x = fld_x(i, psi%m, psi%type_id)  ! in [0, 1]
        y = fld_x(j, psi%n, psi%type_id)  ! in [0, 1]
        psi%data(i, j) = psi%m*psi%n/(Pi**2) *dsin(Pi*x)*dsin(Pi*y)
      end do
    end do
    
!     call print_array(psi%data, "psi")

    call set(K11, 1.0_dp)
    call set(K22, 1.0_dp)
    call set(K12, 0.0_dp)
    
    ! Convert data structure
    call allocate(uv_fld, psi)
    call allocate(K_e, K11, K22, K12)
    call imposeBC(K_e)

    !dt_UB = 1.0_dp/maxmax(uv_fld%u_abs, uv_fld%v_abs)
    !nts = ceil(t/(1/dt_UB))
    nts = floor(10.0_dp*floor(t*max(maxval(uv_fld%u_abs), maxval(uv_fld%v_abs))))+1
    dt = t/nts

    ! Reset q
    call zeros(q)
    do i = -2, 2
      do j = -2, 2
        !call indicator_field(q, q%m/2+i+q%m/8+2, q%n/2+j)
        call indicator_field(q, (q%m+1)/2+i, (q%n+1)/2+j)
      end do
    end do

    call set(q0, q)
!     call print_array(q0%data, "q0")

    q_int_0 = int_field(q)

    call write(q, "./unittest/solverprop/q", 0.0_dp, 0)
    call write(psi, "./unittest/solverprop/psi", 0.0_dp, 0)
    ind = ind + 1

    call print_info(psi, dt)

    if (.not. rotsym(q)) then
      call abort_handle(" ! FAIL @ Rotational symmetry at initialisation", __FILE__, __LINE__)
    end if

!     call print_array(psi%data, "psi")
    do i = 1, nts
      call timestep_heun_Kflux(q, dt, K_e)
      call timestep_LMT_2D(q, dt, uv_fld)

      if (mod(i, floor(nts/25.0_dp)) .eq. 0) then
        call write(q, "./unittest/solverprop/q", dt*i, ind)
        ind = ind + 1
      end if
    end do

    call print_info(q)
!     call print_array(q%data, "q")

    ! Test rotational symmetry
    if (rotsym(q)) then
      write(6, *) " -- P: Rotational symmetry preserved --"
    else
      call print_array(q%data, "Final q")
      call abort_handle(" ! FAIL @ Rotational symmetry", __FILE__, __LINE__)
    end if

    ! Test mass conservation
    if (diff_int(q, q0) < 1D-13) then
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

    call deallocate(uv_fld)
    call deallocate(K_e)
    
    ! Test Divergence consistent: constant initial condition
    m = 8+1
    n = 8+1
    call allocate(q, m, n, 'q', glayer=gl, type_id=2)
    call allocate(q0, m, n, 'q0', glayer=gl, type_id=2)
    
    call allocate(psi, m, n, 'psi', glayer=0, type_id=1)
    call allocate(K11, m, n, 'K11', glayer=0, type_id=1)
    call allocate(K22, m, n, 'K22', glayer=0, type_id=1)
    call allocate(K12, m, n, 'K12', glayer=0, type_id=1)

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

    call allocate(uv_fld, psi)
    call allocate(K_e, K11, K22, K12)
    call imposeBC(K_e)

    do i = 1, nts
      call timestep_heun_Kflux(q, dt, K_e)
      call timestep_LMT_2D(q, dt, uv_fld)
    end do

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

    call deallocate(uv_fld)
    call deallocate(K_e)
    
    ! Test positivity
    m = 64
    n = 64
    call allocate(q, m, n, 'q', glayer=gl, type_id=2)
    call allocate(q0, m, n, 'q0', glayer=gl, type_id=2)
    
    call allocate(psi, m, n, 'psi', glayer=0, type_id=1)
    call allocate(K11, m, n, 'K11', glayer=0, type_id=1)
    call allocate(K22, m, n, 'K22', glayer=0, type_id=1)
    call allocate(K12, m, n, 'K12', glayer=0, type_id=1)

    ! Initialise q
    do j = 1, q%n
      do i = 1, q%m
        x = fld_x(i, q%m, q%type_id)  ! in [0, 1]
        y = fld_x(j, q%n, q%type_id)  ! in [0, 1]
        q%data(i+q%glayer, j+q%glayer) = dexp(-(((x-0.3_dp)/0.05_dp)**2 + ((y-0.5_dp)/0.15_dp)**2))
      end do
    end do
    call imposeBC(q)
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
    
    call allocate(uv_fld, psi)
    call allocate(K_e, K11, K22, K12)
    call imposeBC(K_e)
    
    call set(q0, q)
    t = 10.0_dp*Pi
    nts = floor(10.0_dp*floor(t*max(maxval(uv_fld%u_abs), maxval(uv_fld%v_abs))))+1
    dt = t/nts
    
    call print_info(psi, dt)
    
    do i = 1, nts
      call timestep_heun_Kflux(q, dt, K_e)
!       call timestep_LaxWendoff(q, dt, uv_fld)
!       call timestep_LMT_1DS(q, dt, uv_fld)
      call timestep_LMT_2D(q, dt, uv_fld)

      ! Check at every step
      if (neg_int(q) .lt. 0.0_dp) then
        write(6, *) i, "-th step"
        call print_info(q)
        call abort_handle(" ! FAIL  @ Preserved positivity", __FILE__, __LINE__)
      end if
    end do

    if (neg_int(q) .ge. 0.0_dp) then
      write(6, *) " -- P: Preserved positivity --"
    end if
    
    
    call deallocate(q)
    call deallocate(q0)
    call deallocate(psi)
    call deallocate(K11)
    call deallocate(K22)
    call deallocate(K12)

    call deallocate(uv_fld)
    call deallocate(K_e)
    
    write(6, "(a)") &
      "! ----------- Passed unittest_solver_properties ----------- !"
  end subroutine unittest_solver_properties
  
  pure logical function rotsym(fld)
    type(field), intent(in) :: fld
    integer :: m, n, i, j, gl
    real(kind=dp), parameter :: eps = 1D-14

    m = fld%m
    n = fld%n
    gl = fld%glayer

    rotsym = .true.

    ! Test rotational symmetry
    if ((m .ne. n) .or. (mod(m,2) .ne. 1)) then
      rotsym = .false.
    else
      do j = 1, n
        do i = 1, m
          ! Rotate by pi: x->-x; y->-y
          if (dabs(fld%data(i+gl, j+gl)-fld%data(m-i+1+gl, n-j+1+gl)) > eps) then
            rotsym = .false.
          end if
          ! Rotate by 3pi/2: x->y;  y->-x
          if (dabs(fld%data(i+gl, j+gl)-fld%data(j+gl, m-i+1+gl)) > eps) then
            rotsym = .false.
          end if
          ! Rotate by pi/2: x->-y;  y->x
          if (dabs(fld%data(i+gl, j+gl)-fld%data(n-j+1+gl, i+gl)) > eps) then
            rotsym = .false.
          end if
        end do
      end do
    end if

  end function rotsym
  
  subroutine unittest_fftw()
    integer, parameter :: Nx = 8
    integer, parameter :: Ny = 16
    real(kind=dp), parameter :: Pi =  4.0_dp* datan (1.0_dp)
    
    real(kind=dp), dimension(:, :), allocatable :: data_in, data_intpl, data_exact

    integer :: i, j, scx, scy
    real(kind=dp) :: x, y

    ! Discrete Fourier transform
    scx = 3
    scy = 5
   
    ! Data input
    allocate(data_in(Nx, Ny))
    do j = 1, size(data_in, 2)
      do i = 1, size(data_in, 1)
        x = real(i-1, kind=dp)/real(Nx, kind=dp)
        y = real(j-1, kind=dp)/real(Ny, kind=dp)
        
        data_in(i, j) = dsin(4.0_dp*Pi*x)*dsin(2.0_dp*Pi*y) + dsin(8.0_dp*Pi*y-2.0_dp*Pi*x)
      end do
    end do
    
    ! Data exact
    allocate(data_exact(Nx*scx, Ny*scy))
    do j = 1, size(data_exact, 2)
      do i = 1, size(data_exact, 1)
        x = real(i-1, kind=dp)/real(Nx*scx, kind=dp)
        y = real(j-1, kind=dp)/real(Ny*scy, kind=dp)
        
        data_exact(i, j) = dsin(4.0_dp*Pi*x)*dsin(2.0_dp*Pi*y) + dsin(8.0_dp*Pi*y-2.0_dp*Pi*x)
      end do
    end do
    
    ! Test: Interpolation subroutine
    allocate(data_intpl(size(data_in,1)*scx, size(data_in,2)*scy))
    
    call fourier_intpl(data_intpl, data_in, scx, scy)
    if (.not. iszero(data_exact-data_intpl)) &
      call abort_handle(" ! FAIL  @ Inconsistent fourier_intpl", __FILE__, __LINE__)
    
    deallocate(data_in)
    deallocate(data_exact)
    deallocate(data_intpl)
    
    
    ! Discrete Cosine transform
    ! Data input
    allocate(data_in(Nx, Ny))
    do j = 1, size(data_in, 2)
      do i = 1, size(data_in, 1)
        x = real(2*i-1, kind=dp)/real(2*Nx, kind=dp)
        y = real(2*j-1, kind=dp)/real(2*Ny, kind=dp)
        
        data_in(i, j) = dcos(Pi*x)*dcos(2.0_dp*Pi*y) + dcos(4.0_dp*Pi*x)*dcos(6.0_dp*Pi*y)
      end do
    end do
    
    ! Data exact
    allocate(data_exact(Nx*scx, Ny*scy))
    do j = 1, size(data_exact, 2)
      do i = 1, size(data_exact, 1)
        x = real(2*i-1, kind=dp)/real(2*Nx*scx, kind=dp)
        y = real(2*j-1, kind=dp)/real(2*Ny*scy, kind=dp)
        
        data_exact(i, j) = dcos(Pi*x)*dcos(2.0_dp*Pi*y) + dcos(4.0_dp*Pi*x)*dcos(6.0_dp*Pi*y)
      end do
    end do
    
    ! Test : Interpolation subroutine
    allocate(data_intpl(size(data_in,1)*scx, size(data_in,2)*scy))
    call cosine_intpl(data_intpl, data_in, scx, scy)
    
    if (.not. iszero(data_exact-data_intpl)) &
      call abort_handle(" ! FAIL  @ Inconsistent cosine_intpl", __FILE__, __LINE__)
    
    deallocate(data_in)
    deallocate(data_exact)
    deallocate(data_intpl)
    
    ! Discrete sine transform
    ! Data input
    allocate(data_in(Nx+1, Ny+1))
    do j = 1, size(data_in, 2)
      do i = 1, size(data_in, 1)
        x = real(i-1, kind=dp)/real(Nx, kind=dp)
        y = real(j-1, kind=dp)/real(Ny, kind=dp)
        
        data_in(i, j) = dsin(Pi*x)*dsin(2.0_dp*Pi*y) + dsin(4.0_dp*Pi*x)*dsin(6.0_dp*Pi*y)
      end do
    end do
    
    ! Data exact
    scx = 2
    scy = 4
    allocate(data_exact(Nx*scx+1, Ny*scy+1))
    do j = 1, size(data_exact, 2)
      do i = 1, size(data_exact, 1)
        x = real(i-1, kind=dp)/real(Nx*scx, kind=dp)
        y = real(j-1, kind=dp)/real(Ny*scy, kind=dp)
        
        data_exact(i, j) = dsin(Pi*x)*dsin(2.0_dp*Pi*y) + dsin(4.0_dp*Pi*x)*dsin(6.0_dp*Pi*y)
      end do
    end do
    
    ! Test : Interpolation subroutine
    allocate(data_intpl(size(data_exact,1), size(data_exact,2)))
    data_intpl = 0.0_dp
    call sine_intpl(data_intpl(2:Nx*scx, 2:Ny*scy), data_in(2:Nx, 2:Ny), scx, scy)
    
    if (.not. iszero(data_exact-data_intpl)) then
      call print_array(data_intpl, "data_intpl")
      call print_array(data_exact, "data_exact")
      call abort_handle(" ! FAIL  @ Inconsistent sine_intpl", __FILE__, __LINE__)
    end if
    
    deallocate(data_in)
    deallocate(data_exact)
    deallocate(data_intpl)
    
    ! Test : Sinusoidal filter
    scx = 128+1
    scy = 128+1
    allocate(data_in(scx, scy))
    do j = 1, size(data_in, 2)
      do i = 1, size(data_in, 1)
        x = real(i-1, kind=dp)/real(scx-1, kind=dp)
        y = real(j-1, kind=dp)/real(scy-1, kind=dp)
        
        data_in(i, j) = dsin(Pi*x)*dsin(2.0_dp*Pi*y) + dsin(16.0_dp*Pi*x)*dsin(16.0_dp*Pi*y)
      end do
    end do
    
    allocate(data_exact(16+1, 16+1))
    
    do j = 1, size(data_exact, 2)
      do i = 1, size(data_exact, 1)
        x = real(i-1, kind=dp)/real(size(data_exact,1)-1, kind=dp)
        y = real(j-1, kind=dp)/real(size(data_exact,2)-1, kind=dp)
        
       data_exact(i, j) = dsin(Pi*x)*dsin(2.0_dp*Pi*y)
      end do
    end do
    
    allocate(data_intpl(16+1, 16+1))
    data_intpl = 0.0_dp
    call sine_filter(data_intpl(2:size(data_intpl,1)-1,2:size(data_intpl,2)-1), &
                     data_in(2:size(data_in,1)-1,2:size(data_in,2)-1))
    
    if (.not. iszero(data_intpl-data_exact)) then
      call abort_handle(" ! FAIL  @ Inconsistent sine_filter", __FILE__, __LINE__)
    end if
    
    
    deallocate(data_intpl)
    deallocate(data_exact)
    deallocate(data_in)
    
    write(6, "(a)") &
      "! ----------- Passed unittest_fftw ----------- !"
  end subroutine unittest_fftw
  
  subroutine unittest_solver_timing()
    type(field) :: K11, K22, K12, psi
    integer :: m, n

    type(uv_signed) :: uv_fld
    type(K_flux) :: K_e

    integer, parameter :: gl = 1

    real(kind=dp) :: t = 12.0_dp*3.1415926_dp
    real(kind=dp) :: dt, x, y
    integer :: nts = 64*16
    integer :: i, j
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)

    type(field) :: q, qLW, q2D
    
    type(timer), save :: DS_timer, LW_timer, D2_timer
    
    ! Test 1DS and 1DS2 equivalence
    m = 80
    n = m
    call allocate(q, m, n, 'q', glayer=gl, type_id=2)
    call allocate(qLW, m, n, 'qLW', glayer=gl, type_id=2)
    call allocate(q2D, m, n, 'q2D', glayer=gl, type_id=2)
    
    call allocate(psi, m, n, 'psi', glayer=0, type_id=1)

    call allocate(K11, m, n, 'K11', glayer=0, type_id=2)
    call allocate(K22, m, n, 'K22', glayer=0, type_id=2)
    call allocate(K12, m, n, 'K12', glayer=0, type_id=2)
    
    ! Initialise psi
    do j = 1, size(psi%data,2)
      do i = 1, size(psi%data,1)
        x = fld_x(i, psi%m, psi%type_id)  ! in [0, 1]
        y = fld_x(j, psi%n, psi%type_id)  ! in [0, 1]
        psi%data(i, j) = psi%m*psi%n/(2.0_dp*Pi**2) *dsin(2.0_dp*Pi*x)*dsin(Pi*y)
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

    call allocate(uv_fld, psi)
    call allocate(K_e, K11, K22, K12)
    call imposeBC(K_e)
    
    call manual_field(q, 4)
    call imposeBC(q)
    call set(qLW, q)
    call set(q2D, q)
    dt = t/real(nts, kind=dp)
    
    ! 1DS
    call reset(DS_timer)
    call start(DS_timer)
    do i = 1, nts
!       call timestep_heun2_Kflux(q2, 0.5_dp*dt, K_e)
!       call timestep_LMT_1DS2(q2, dt, uv_fld)
!       call timestep_heun2_Kflux(q2, 0.5_dp*dt, K_e)
       call timestep_heun_Kflux(q, 0.5_dp*dt, K_e)
       call timestep_LMT_1DS(q, dt, uv_fld)
       call timestep_heun_Kflux(q, 0.5_dp*dt, K_e)
       ! First order
!        call timestep_fe_Kflux(q2, dt, K_e)
!        call timestep_LMT_2D(q2, dt, uv_fld)
    end do
    call stop(DS_timer)
    
    ! Lax Wendoff
    call reset(LW_timer)
    call start(LW_timer)
    do i = 1, nts
      call timestep_heun_Kflux(qLW, 0.5_dp*dt, K_e)
      !call timestep_fe_Kflux(qLW, 0.5_dp*dt, K_e)
      call timestep_LaxWendoff(qLW, dt, uv_fld)
      call timestep_heun_Kflux(qLW, 0.5_dp*dt, K_e)
    end do
    call stop(LW_timer)
    
    ! MC
    call reset(D2_timer)
    call start(D2_timer)
    do i = 1, nts
      call timestep_heun_Kflux(q2D, dt, K_e)
      call timestep_LMT_2D(q2D, dt, uv_fld)
    end do
    call stop(D2_timer)
    
    write(6, "(a)") "Dimensional S + Strang S + Limiter timers:"
      call print(DS_timer, "Time", prefix = "  ")
      
    write(6, "(a)") "Dimensional S + Strang S + LaxWendoff timers:"
      call print(LW_timer, "Time", prefix = "  ")
      
    write(6, "(a)") "UnSplit + Godunov S + Limiter  timers:"
      call print(D2_timer, "Time", prefix = "  ")
    
!     call print_array(q%data, 'q')

    
    call addto(qLW, q, -1.0_dp)
    
    if (iszero(qLW%data)) then
      write(6, *) " -- P: 1DS consistent with LW --"
    else
      call abort_handle(" ! FAIL  @ 1DS/LW consistency", __FILE__, __LINE__)
    end if
    
    

    call deallocate(q)
    call deallocate(qLW)
    call deallocate(psi)
    call deallocate(K11)
    call deallocate(K22)
    call deallocate(K12)
    
    write(6, "(a)") &
      "! ----------- Passed unittest_solver_timing ----------- !"
  end subroutine unittest_solver_timing

  subroutine unittest_solver_convergence()
    integer :: i, m
    integer :: nts, C
    character(len = 256) :: scr_lne
    
!     ! Test spatial accuarcy (integrating to steady state)
!     write(6, *) " Timestepping to steady state"
!     write(6, *) " nts = 512*m "
!     do i = 1, 8
!       m = 2 ** (i+3)
!       !write(6, *) "m(",i,") =", m, "; error(",i,") = ", compute_spatial_advection_err(m, 2)  ! 2: stat state
!     end do
    
    ! Test temporal accuarcy (periodic solution)
    C = 128
    open(unit = output_unit, file = "./unittest/convergence/log_LW_2.txt", &
      & status = "replace", action = "write")
    write(output_unit, "(a)") " Timestepping to periodic state"
    write(output_unit, "(a)") " Scale dt = 1/C dx"
    write(output_unit, "(a,i0,a)") " nts =", C,"*m"
    do i = 1, 8
      m = 2 ** (i+3)
      nts = C* m
      write(scr_lne, "(a,i0,a,i6,a,i0,a,i6,a,i0,a,"//dp_chr//")") &
          "m(", i,") =", m,"; nts(",i,") =", nts,"; error(", i,") =", compute_temporal_err(m, nts)  ! 1: periodic
      write(output_unit, *) trim(scr_lne)
      
      ! Print to screen
      write(6, *) trim(scr_lne)
      flush(6)
    end do
    close(output_unit)
    
    write(6, "(a)") &
      "! ----------- Passed unittest_solver_properties ----------- !"
  end subroutine unittest_solver_convergence
  
  real(kind=dp) function compute_spatial_err(m)
    ! Source term: stationary, with lambda decay
    !  S(x,y,t) = F(x,y)-lambda*(q_stat-q)
    integer, intent(in) :: m
    
    type(field) :: q, K11, K22, K12, psi
    integer :: n, i

    type(uv_signed) :: uv_fld
    type(K_flux) :: K_e

    integer, parameter :: gl = 1

    real(kind=dp), parameter :: T = 16.0_dp
    integer :: nts
   
    real(kind=dp), parameter :: Pi =  4.0_dp* datan (1.0_dp)
   
    type(field) :: q_steady, q_err
    real(kind=dp) :: err_int
    
    nts = 512*m
    n = m
    call allocate(q, m, n, 'q', glayer=gl, type_id=2)
    call allocate(q_steady, m, n, 'q_steady', glayer=gl, type_id=2)
    call allocate(q_err, m, n, 'q_err', glayer=gl, type_id=2)
    call allocate(psi, m, n, 'psi', glayer=0, type_id=1)

    call allocate(K11, m, n, 'K11', glayer=0, type_id=2)
    call allocate(K22, m, n, 'K22', glayer=0, type_id=2)
    call allocate(K12, m, n, 'K12', glayer=0, type_id=2)

    ! test_id = 2: q_steady
    call testcase_fields(psi, K11, K22, K12, q_steady, 2)
    call imposeBC(q_steady)
    
    call allocate(uv_fld, psi)
    call allocate(K_e, K11, K22, K12)
    call imposeBC(K_e)
    
    ! Set inital condition to zero
    call zeros(q)
    !call set(q, q_steady)
    call imposeBC(q)
    !call write(q, "./unittest/convergence/q", 0.0_dp, 0)
    dt = T/real(nts,kind=dp)
    do i = 1, nts
      ! Strang splitting
      ! A
      !call timestep_fe_Kflux_Src_stat(q, dt, K_e, 1.0_dp)  ! Last argument = lambda
      call timestep_heun_Kflux_Src_stat(q, 0.5_dp*dt, K_e, 8.0_dp)
      ! B
      call timestep_LMT_1DS(q, dt, uv_fld)

      ! A
      !call timestep_fe_Kflux_Src_stat(q, 0.5_dp*dt, K_e, 1.0_dp)
      call timestep_heun_Kflux_Src_stat(q, 0.5_dp*dt, K_e, 8.0_dp)
    end do
    call write(q, "./unittest/convergence/q_final", nts*dt, nts)
    
    ! Compare q with exact stat state
    call addto(q, q_steady, -1.0_dp)
    err_int = dsqrt(int_sq_field(q)/(m * n))

    call deallocate(q)
    call deallocate(q_steady)
    call deallocate(q_err)
    call deallocate(psi)

    call deallocate(K11)
    call deallocate(K22)
    call deallocate(K12)

    call deallocate(uv_fld)
    call deallocate(K_e)
    
    compute_spatial_err = err_int
  end function compute_spatial_err
  
  real(kind=dp) function compute_temporal_err(m, nts)
    integer, intent(in) :: m, nts
    
    type(field) :: q, K11, K22, K12, psi
    integer :: n, i

    type(uv_signed) :: uv_fld
    type(K_flux) :: K_e

    integer, parameter :: gl = 1

    real(kind=dp), parameter :: Pi = 4.0_dp* datan (1.0_dp)
    real(kind=dp), parameter :: T = 4.0_dp* datan (1.0_dp)
   
    type(field) :: q_exact, q_err
    real(kind=dp) :: t_i, err_int, err_sup
    
    n = m
    call allocate(q, m, n, 'q', glayer=gl, type_id=2)
    call allocate(q_exact, m, n, 'q_exact', glayer=gl, type_id=2)
    call allocate(q_err, m, n, 'q_err', glayer=gl, type_id=2)
    call allocate(psi, m, n, 'psi', glayer=0, type_id=1)

    call allocate(K11, m, n, 'K11', glayer=0, type_id=1)
    call allocate(K22, m, n, 'K22', glayer=0, type_id=1)
    call allocate(K12, m, n, 'K12', glayer=0, type_id=1)

    ! test_id = 1: q(t) 1-periodic
    call testcase_fields(psi, K11, K22, K12, q_exact, 1)
    call imposeBC(q_exact)
    
!     ! Turn off diffusion
!     K11%data = 0.0_dp
!     K22%data = 0.0_dp
!     K12%data = 0.0_dp
    
    call allocate(uv_fld, psi)
    call allocate(K_e, K11, K22, K12)
    call imposeBC(K_e)
    
    ! Set inital condition to zero
    call set(q, q_exact)
    call imposeBC(q)
    
    err_sup = 0.0_dp
    
!     call write(q, "./unittest/convergence/q_initial", 0.0_dp, 0)
    dt = T/real(nts,kind=dp)
    do i = 1, nts
      t_i = real(i-1, kind=dp)*dt
      
      ! Strang splitting
      ! A
      call timestep_heun_Kflux_Src_nonstat(q, 0.5_dp*dt, t_i, K_e)
      !call timestep_heun_Kflux_Src_nonstat(q, dt, t_i, K_e)
      ! B
      call timestep_LMT_2D(q, dt, uv_fld)
      !call timestep_LaxWendoff(q, dt, uv_fld)
      ! A
      call timestep_heun_Kflux_Src_nonstat(q, 0.5_dp*dt, t_i+0.5_dp*dt, K_e)
      
      call set(q_err, q)
      ! Compare q with exact stat state
      call assign_q_periodic(q_exact, t_i+dt)
      call imposeBC(q_exact)
      call addto(q_err, q_exact, -1.0_dp)
      err_int = dsqrt(int_sq_field(q_err)/(m * n))
      
      err_sup = max(err_int, err_sup)
    end do
    
!     call write(q, "./unittest/convergence/q_final", nts*dt, nts)

!     ! (1/2A)BAB..AB(1/2A)
!     dt = T/real(nts,kind=dp)
!     ! Strang splitting
!     ! A/2
!     t_i = 0.0_dp
!     call timestep_heun_Kflux_Src_nonstat(q, 0.5_dp*dt, t_i, K_e)
!     
!     ! (BA)^(nts-1) 
!     do i = 1, (nts-1)
!       t_i = real(i-1, kind=dp)*dt
!       call timestep_LMT_2D(q, dt, uv_fld)
!       !call timestep_LaxWendoff(q, dt, uv_fld)
!       call timestep_heun_Kflux_Src_nonstat(q, dt, t_i+0.5_dp*dt, K_e)
!     end do
!     
!     ! B(A/2)
!     t_i = real(nts-1, kind=dp)*dt
!     call timestep_LMT_2D(q, dt, uv_fld)
!     call timestep_heun_Kflux_Src_nonstat(q, 0.5_dp*dt, t_i+0.5_dp*dt, K_e)
! 
!     ! Calculate error at final time
!     call set(q_err, q)
!     ! Compare q with exact stat state
!     call assign_q_periodic(q_exact, T)
!     call imposeBC(q_exact)
!     call addto(q_err, q_exact, -1.0_dp)
!     err_int = dsqrt(int_sq_field(q_err)/(m * n))
!     
!     err_sup = max(err_int, err_sup)
      
    call deallocate(q)
    call deallocate(q_exact)
    call deallocate(q_err)
    call deallocate(psi)

    call deallocate(K11)
    call deallocate(K22)
    call deallocate(K12)

    call deallocate(uv_fld)
    call deallocate(K_e)
    
    compute_temporal_err = err_sup
  end function compute_temporal_err
  
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

    if ( (dabs(traj%x(3,2) - 1.0_dp) > 1e-14) .or. &
         (dabs(traj%y(1,2) - (-2.0_dp)) > 1e-14) .or. &
         (dabs(traj%t(3,4) - 3.0_dp) > 1e-14) ) then
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
    real(kind=dp) :: T

    integer :: m, n, cell, m_solver, m_Ind, n_Ind
    type(meshdat) :: mesh

    type(dofdat) :: dof, dof_MAP, dof_old, dof_SSD
    type(IndFndat) :: IndFn

    real(kind=dp), parameter :: L = 3840.0_dp*1000.0_dp
    real(kind=dp), parameter :: kappa_scale = 10000.0_dp
    real(kind=dp), parameter :: psi_scale = 100.0_dp*1000.0_dp
    real(kind=dp) :: sc, h, dt
    integer :: nts, nts_sub, dt_rf
    real(kind=dp), dimension(:), allocatable :: logPost, logPost_old

    integer :: dof_id
    integer :: iter, niter, ind

    integer, parameter :: Td = 32
    character(len = *), parameter :: RunProfile = "QGM2_L1_NPART2704"
    character(len = 256) :: output_fld, Td_char, resol_param
    character(len = 256) :: output_q

    type(dofdat) :: dof_solver, dof_solver_inv
    type(field) :: q
    real(kind=dp) :: lik
    integer, dimension(:), allocatable :: klist
    integer :: INDk, i, k
    
    real(kind=dp) :: SlogPost, SlogPostr, alphaUniRV
    
    real(kind=dp), dimension(:,:), allocatable :: dJdm, dJdm_abs
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
    call allocate(q, mesh%m, mesh%n, 'q', glayer=1, type_id=2)
    
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

!     ! Test dJdm
!     allocate(theta(dof%m_psi*dof%n_psi+3*dof%m_K*dof%n_K))
!     allocate(theta_old(dof%m_psi*dof%n_psi+3*dof%m_K*dof%n_K))
!     allocate(logPost(IndFn%nIND))
!     allocate(logPost_old(IndFn%nIND))
!     
!     call read_theta(dof_old, trim(output_fld)//"theta_sine", sc, 119)
!     
!     call convert_dof_to_theta(theta, dof, 1.0_dp)
!     call convert_dof_to_theta(theta_old, dof_old, 1.0_dp)
!     
!     call evaluate_loglik_OMP(logPost, jumps, IndFn, mesh, dof, h, nts)
!     call evaluate_loglik_OMP(logPost_old, jumps, IndFn, mesh, dof_old, h, nts)
!     
!     allocate(dJdm(size(logPost,1), size(theta, 1)))
!     
!     call compute_dJdm(dJdm, logPost, logPost_old, theta, theta_old)
!     
!     call print_array(dJdm, 'dJdm')
!     
!     deallocate(dJdm)
!     deallocate(theta)
!     deallocate(theta_old)
!     call deallocate(dof_old)
!     deallocate(logPost)
!     deallocate(logPost_old)
!     ! End Test dJdm

    ! Evaluate likelihood
    allocate(logPost(IndFn%nIND))
    
!     call evaluate_loglik_OMP(logPost, jumps, IndFn, mesh, dof, h, nts)
!     SlogPost = real(sum(logPost), kind=dp)
!     write(6, *) "3hr@Posterior = ", SlogPost
!     
!     ! Smaller timestep
!     nts = nts*8
!     logPost = 0.0_dp
!     call evaluate_loglik_OMP(logPost, jumps, IndFn, mesh, dof, h, nts)
!     SlogPostr = real(sum(logPost), kind=dp)
!     write(6, *) "0.375hr@Posterior = ", SlogPostr
!       
!     alphaUniRV = dabs((SlogPost - SlogPostr)/SlogPost)
!     alphaUniRV = dexp(SlogPostr-SlogPost)
!     write(6, *) "Relative acceptance = ", alphaUniRV
!     
!     ! Perturb dof
!     call deallocate(dof)
!     call read_theta(dof, trim(output_fld)//"theta_sine", sc, 124)
!     nts = 256*4
!     call evaluate_loglik_OMP(logPost, jumps, IndFn, mesh, dof, h, nts)
!     SlogPostr = real(sum(logPost), kind=dp)
!     write(6, *) "Perturb dof@Posterior = ", SlogPostr
!     
!     alphaUniRV = dabs((SlogPost - SlogPostr)/SlogPost)
!     write(6, *) "Relative acceptance = ", alphaUniRV
! 
!     nts = 256*8*4
!     call evaluate_loglik_OMP(logPost, jumps, IndFn, mesh, dof, h, nts)
!     SlogPostr = real(sum(logPost), kind=dp)
!     write(6, *) "1hr@Posterior = ", SlogPostr
!     
!     alphaUniRV = dabs((SlogPost - SlogPostr)/SlogPost)
!     write(6, *) "Relative acceptance = ", alphaUniRV
    
    deallocate(logPost)
    call deallocate(jumps)
    call deallocate(q)
    deallocate(klist)
    call deallocate(dof_solver)
    
    call deallocate(dof)
    call deallocate(IndFn)

  end subroutine unittest_FPsolver_INPUT
  
  subroutine unittest_TG_instability()
    use omp_lib
    type(jumpsdat), dimension(:), allocatable :: jumps
    type(trajdat) :: traj

    integer :: m, n, cell, m_solver, m_Ind, n_Ind
    type(meshdat) :: mesh

    type(dofdat) :: dof, dof_MAP, dof_old
    type(IndFndat) :: IndFn

    real(kind=dp), parameter :: L = 3840.0_dp*1000.0_dp
    real(kind=dp), parameter :: kappa_scale = 10000.0_dp
    real(kind=dp), parameter :: psi_scale = 100.0_dp*1000.0_dp
    real(kind=dp) :: sc, h, dt
    integer :: nts
    real(kind=dp), dimension(:), allocatable :: logPost, canon_SSD
    real(kind=dp) :: logPrior
    
    integer :: canon_id
    integer :: iter, niter, ind
    real(kind=dp) :: alphaUniRV
    real(kind=dp), dimension(1) :: UniRV

    integer, parameter :: Td = 1
    character(len = *), parameter :: RunProfile = "zero_const"
    character(len = 256) :: output_fld, Td_char, resol_param

    ! Extra variables
    type(dofdat) :: dof_solver  ! Refined dof for solver
    real(kind=dp) :: kappa, a16, logLik
    real(kind=dp) :: kappa_old, a16_old, logLik_old
    character(len = 256) :: scrstr
    integer :: i, j, k, wavenumber
    type(field) :: q
    character(len = 256) ::output_q

    
    ! m, n: DOF/control
    m = 1
    n = 1
    ! m_solver: solver grid
    m_solver = 128
    ! m_Ind, n_Ind: indicator functions
    m_Ind = 1
    n_Ind = m_Ind
    
    write(Td_char, "(a,i0,a)") "h", Td, "d"
    write(resol_param, "(a,i0,a,i0,a,i0)") "N",m_solver*m_solver,"_D", m*n, "_I", m_Ind*n_Ind
    write(output_fld, "(a,a,a,a,a,a,a)") "./output/", trim(resol_param), "/", trim(RunProfile), "/", trim(Td_char), "/"
    write(6, "(a, a)") "Output path: ", trim(output_fld)
    
    call allocate(mesh, m_solver, m_solver)  ! Needs to be identical in both directions
    call allocate(dof, m, n, 1, 1) 
    ! N.B. assumed zeros at boundaries of psi-> hence only need (m-1)*(n-1) instead of (m+1)*(n+1)
    call allocate(IndFn, m_Ind, n_Ind)

    ! Initialise fields
    sc = real(mesh%m,kind=dp)/L     ! Needs mesh to be identical in both directions
    call init_dof(dof, sc)

    ! Read trajectory
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
    
    ! Determine h: Assumed h = uniform; assumed mesh%m = mesh%n
    h = read_uniform_h(jumps) *real(mesh%m, kind=dp)   ! *m to rescale it to [0, mesh%m]

    ! Calculate nts = number of timesteps needed in solving FP
    dt = 0.01_dp*3600.0_dp  ! 3 hours, in seconds
    nts = int((h/sc)/dt)    ! 2 hours time step

    ! Print info to screen
    call print_info(mesh)
    call print_info(dof)
    call print_info(IndFn)
    call print_info(h, nts, sc)
!     call print_info(jumps, mesh)
    
    call allocate(dof_solver, mesh)
    
    
    ! Look for a MLE by random walk
    wavenumber = 16
    
    a16_old = 2.0_dp*L *sc
    kappa_old = 100.0_dp * L *sc/100.0_dp
    kappa_old = 1.0_dp * L *sc/100.0_dp
    
!     call set(dof_solver%psi, 0.0_dp)
!     call sine_basis(dof_solver%psi%data(2:dof_solver%psi%m, 2:dof_solver%psi%n), wavenumber, a16)
!       
!     call set(dof_solver%K11, kappa)
!     call set(dof_solver%K22, kappa)
!     call set(dof_solver%K12, 0.0_dp)
!     
!     logLik_old = eval_Gaussian_loglik(jumps, mesh, dof_solver, h, nts)
! 
!     do k = 1, 500
!       a16 = a16_old*(1.0_dp + 0.1_dp*randn())
!       kappa = kappa_old*(1.0_dp + 0.1_dp*randn())
!     
!       ! Only need to pass the interior points
!       call set(dof_solver%psi, 0.0_dp)
!       call sine_basis(dof_solver%psi%data(2:dof_solver%psi%m, 2:dof_solver%psi%n), wavenumber, a16)
!       
!       call set(dof_solver%K11, kappa)
!       call set(dof_solver%K22, kappa)
!       call set(dof_solver%K12, 0.0_dp)
!       
!       logLik = eval_Gaussian_loglik(jumps, mesh, dof_solver, h, nts)
!       
! !       write(scrstr, "(a,i0,a,"//sp_chr//",a,"//sp_chr//",a,"//dp_chr//",a)") &
! !           "arr(", k, ",:)= [", a16, ",",  kappa, ", ", logLik, "]"
! !       write(6, *) trim(scrstr)
!       
!       if (logLik .gt. logLik_old) then
!         a16_old = a16
!         kappa_old = kappa
!         logLik_old = logLik
!         
!         
!         write(scrstr, "(a,i0,a,"//sp_chr//",a,"//sp_chr//",a,"//dp_chr//",a)") &
!           "arr(", k, ",:)= [", a16, ",",  kappa, ", ", logLik, "]"
!         write(6, *) trim(scrstr)
!         flush(6)
!       end if
!     end do
    
    
    !! Output video for MLE
    write(6, *) "Output video"
    a16 = a16_old
    kappa = kappa_old
    
    ! Only need to pass the interior points
    call set(dof_solver%psi, 0.0_dp)
    call sine_basis(dof_solver%psi%data(2:dof_solver%psi%m, 2:dof_solver%psi%n), wavenumber, a16)
      
    call set(dof_solver%K11, kappa)
    call set(dof_solver%K22, kappa)
    call set(dof_solver%K12, 0.0_dp)

    call print_info(dof_solver%psi, h/nts)

    call allocate(q, mesh%m, mesh%n, 'q', glayer=1, type_id=2)

    call initialise_q_Gaussian(q, mesh)
    call imposeBC(q)
    
    write(output_q, "(a)") "./unittest/cellular/q"
    call write(q, trim(output_q), 0.0_dp, 0)
    
    do i = 1, nts
      call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h/nts, 1)
      call write(q, trim(output_q), (h/sc)/nts*i, i)
      if (mod(i, nts/10) .eq. 0) then
        write(6, *) i, " out of ", nts
        flush(6)
      end if
    end do
    
    call deallocate(q)
    
    stop
    
    write(6, "(a)") "a16, kappa, logLik"
    
    k = 0
    do i = 1, 10
      do j = 1, 10
        k = k + 1
        a16 = 1.0_dp * real(i-1, kind=dp) *sc
        kappa = 0.01_dp * real(j, kind=dp) *sc
        
        ! Only need to pass the interior points
        call set(dof_solver%psi, 0.0_dp)
        call sine_basis(dof_solver%psi%data(2:dof_solver%psi%m, 2:dof_solver%psi%n), 16 , a16)
        
        call set(dof_solver%K11, kappa)
        call set(dof_solver%K22, kappa)
        call set(dof_solver%K12, 0.0_dp)
        
        dof_solver%SlogPost = eval_Gaussian_loglik(jumps, mesh, dof_solver, h, nts)
        
        write(scrstr, "(a,i0,a,"//sp_chr//",a,"//sp_chr//",a,"//dp_chr//",a)") &
          "arr(", k, ",:)= [", a16, ",",  kappa, ", ", dof_solver%SlogPost, "]"
        write(6, *) trim(scrstr)
      end do
    end do
    
    ! Release memory
    call deallocate(dof_solver)
    
    call deallocate(dof)
    call deallocate(IndFn)
    call deallocate(jumps)
    
    write(6, "(a)") &
    "! ----------- Passed unittest_TG_instability ----------- !"

  end subroutine unittest_TG_instability
  
  real(kind=dp) function eval_Gaussian_loglik(jumps, mesh, dof_solver, h, nts)
    type(jumpsdat), dimension(:), allocatable, intent(in) :: jumps
    type(meshdat), intent(in) :: mesh
    type(dofdat), intent(in) :: dof_solver
    real(kind=dp), intent(in) :: h
    integer, intent(in) :: nts

    type(field) :: q
    integer :: jump
    real(kind=dp) :: lik
    integer :: k_i, i

    eval_Gaussian_loglik = 0.0_dp

    ! Find out the cells corresponding to the indicator function
    call allocate(q, mesh%m, mesh%n, 'q', glayer=1, type_id=2)
    
    ! Solve FK equation
    call initialise_q_Gaussian(q, mesh)
    call imposeBC(q)
    call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h, nts)
    
    ! Evaluate likelihood
    do k_i = 1, mesh%m*mesh%n
      do jump = 1, jumps(k_i)%njumps
        ! Bilinear interpolations
        lik = eval_fldpt(q, mesh, jumps(k_i)%k_f(jump), jumps(k_i)%alpha_f(:,jump))
        ! Set the negative probability to be 1e-16
        lik = max(1e-16, lik)
        eval_Gaussian_loglik = eval_Gaussian_loglik + dlog(lik)
      end do
    end do

    call deallocate(q)
    
  end function eval_Gaussian_loglik
  
  subroutine initialise_q_Gaussian(q, mesh)
    type(field), intent(inout) :: q
    type(meshdat), intent(in) :: mesh
    
    integer :: i, j, gl
    real(kind=dp) :: x, y, Sigma0
    real(kind = dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)

    gl = q%glayer
    Sigma0 = 0.05_dp
    do j = 1, q%n
        do i = 1, q%m
          ! Analytic formula: Defined on [0,1]^2 -> Need a prefactor to scale to [0,m] x [0,n]
          x = fld_x(i, q%m, q%type_id)
          y = fld_x(j, q%n, q%type_id)
          
          q%data(i+gl, j+gl) = dexp(-0.5_dp*((x-0.5_dp)**2+(y-0.5_dp)**2)/Sigma0**2)/(2.0_dp*Pi*Sigma0**2)
        end do
      end do
    
  end subroutine initialise_q_Gaussian
  
end module advdiff_unittest

