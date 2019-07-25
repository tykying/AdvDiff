program advdiff
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

  ! Main program
  call unittest()

contains

  subroutine unittest()
    write(6, *) "!!! ----------- Unittests begin ----------- !!!"
    !call unittest_comptools()
    !call unittest_trajdata()
    !call unittest_solver_properties()
    !call unittest_fftw()
    !call unittest_solver_timing()
    !call unittest_IO()
    !call unittest_2DLik()
    !call unittest_FPsolver_INPUT()
    !call unittest_solver_convergence()
    call unittest_inference_OMP()
    write(6, *) "!!! ----------- All unittests passed ----------- !!!"
  end subroutine unittest

  subroutine unittest_comptools()
    real(kind = dp) :: mean, stdv, mu, sigma, gamma, sample
    real(kind = dp), dimension(2,2) :: K, Krt
    integer :: i, n
    real(kind = dp) :: sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy

    ! randn generator
    mean = 2.0_dp
    stdv = 5.0_dp
    n = 10000

    call init_random_seed();
    mu = 0.0_dp
    sigma = 0.0_dp

    do i = 1, n
      sample = mean + stdv * randn()

      mu = mu + sample
      sigma = sigma + sample*sample
      gamma = gamma + sample*sample*sample
    end do

    mu = mu/n
    sigma = dsqrt(sigma/n - mu*mu)
    gamma = (gamma/n - 3.0_dp*mu*sigma*sigma - mu*mu*mu)/(sigma*sigma*sigma)

    write(6, "(a, i0)") " Testing random number generator: n = ", n
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
    
    write(6, "(a)") &
      "! ----------- Passed unittest_comptools ----------- !"
  end subroutine unittest_comptools

  subroutine unittest_trajdata()
    type(jumpsdat), dimension(:), pointer :: jumps
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

    real(kind=dp) :: t = 12.0_dp*3.1415926_dp
    real(kind=dp) :: dt
    integer :: nts = 10000
    integer :: ind = 0

    type(field) :: q0
    real(kind=dp) :: q_int_0
    real(kind=dp), dimension(:,:), allocatable :: dqx, dqy
    
    real(kind=dp), dimension(:,:), allocatable :: exactval

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
    
    ! Test calculation of psi and uv_fld
    call allocate(psi, m, n, 'q', glayer=0, type_id=1)
    call manual_field(psi, 6)  ! Rigid body rotation
    call allocate(uv_fld, psi)
    

    call print_array(uv_fld%up, "up")
    call print_array(uv_fld%um, "um")
    call print_array(real(uv_fld%u_sgn, kind=dp), "u_sgn")
    call print_array(uv_fld%u_abs, "u_abs")
    
    call print_array(uv_fld%vp, "vp")
    call print_array(uv_fld%vm, "vm")
    !call print_array(uv_fld%v_sgn, "v_sgn")
    !call print_array(uv_fld%v_abs, "v_abs")
    
    call deallocate(psi)
    call deallocate(uv_fld)
    
    call allocate(K11, m, n, 'K11', glayer=0, type_id=2)
    call allocate(K22, m, n, 'K22', glayer=0, type_id=2)
    call allocate(K12, m, n, 'K12', glayer=0, type_id=2)
    
    call assemble_K(K11, 11)
    call assemble_K(K22, 22)
    call assemble_K(K12, 12)

    call allocate(K_e, K11, K22, K12)
    call imposeBC(K_e)

    ! K11
    allocate(exactval(m+1, n))
    call assemble_Ke(exactval, 11)
    exactval(1,:) = 0.0_dp
    exactval(m+1,:) = 0.0_dp
    if (iszero(K_e%K11e-exactval)) then
      write(6, *) " -- P: Interpolation K11 --"
    else
      call print_array(K_e%K11e, "K11e")
      call print_array(exactval, "Exact K11e")
      call abort_handle(" ! FAIL  @ Interpolation K11", __FILE__, __LINE__)
    end if
    deallocate(exactval)
    
    ! K12
    allocate(exactval(m+1, n))
    call assemble_Ke(exactval, 12)
    exactval(1,:) = 0.0_dp
    exactval(m+1,:) = 0.0_dp
    if (iszero(K_e%K12e-exactval)) then
      write(6, *) " -- P: Interpolation K12 --"
    else
      call print_array(K_e%K12e, "K12e")
      call print_array(exactval, "Exact K12e")
      call abort_handle(" ! FAIL  @ Interpolation K12", __FILE__, __LINE__)
    end if
    deallocate(exactval)

    ! K21
    allocate(exactval(m, n+1))
    call assemble_Ke(exactval, 12)
    exactval(:,1) = 0.0_dp
    exactval(:,n+1) = 0.0_dp
    if (iszero(K_e%K21e-exactval)) then
      write(6, *) " -- P: Interpolation K21 --"
    else
      call print_array(K_e%K21e, "K21e")
      call print_array(exactval, "Exact K21e")
      call abort_handle(" ! FAIL  @ Interpolation K21", __FILE__, __LINE__)
    end if
    deallocate(exactval)
    
    ! K22
    allocate(exactval(m, n+1))
    call assemble_Ke(exactval, 22)
    exactval(:,1) = 0.0_dp
    exactval(:,n+1) = 0.0_dp
    if (iszero(K_e%K22e-exactval)) then
      write(6, *) " -- P: Interpolation K22 --"
    else
      call print_array(K_e%K22e, "K22e")
      call print_array(exactval, "Exact K22e")
      call abort_handle(" ! FAIL  @ Interpolation K22", __FILE__, __LINE__)
    end if
    deallocate(exactval)

    call deallocate(K_e)
    call deallocate(K11)
    call deallocate(K22)
    call deallocate(K12)
    

    m = 8+1
    n = 8+1
    call allocate(q, m, n, 'q', glayer=gl, type_id=2)
    call allocate(q0, m, n, 'q0', glayer=gl, type_id=2)
    call allocate(psi, m, n, 'psi', glayer=0, type_id=1)

    call allocate(K11, m, n, 'K11', glayer=0, type_id=2)
    call allocate(K22, m, n, 'K22', glayer=0, type_id=2)
    call allocate(K12, m, n, 'K12', glayer=0, type_id=2)


    ! Assign fields
    call manual_field(psi, 6)  ! Rigid body rotation
    call scale(psi, 2.0_dp)

    call set(K11, 100.0_dp)
    call set(K22, 100.0_dp)
    call set(K12, 0.0_dp)

    call scale(K11, 0.01_dp)
    call scale(K22, 0.01_dp)
    call scale(K12, 0.0_dp)

    ! Convert data structure
    call allocate(uv_fld, psi)
    call allocate(K_e, K11, K22, K12)
    call imposeBC(K_e)

    !dt_UB = 1.0_dp/maxmax(uv_fld%u_abs, uv_fld%v_abs)
    !nts = ceil(t/(1/dt_UB))
    nts = floor(2.0_dp*floor(t*maxmax(uv_fld%u_abs, uv_fld%v_abs)))+1
    dt = t/nts

    ! Reset q
    call zeros(q)
    do i = -4, 4
      do j = -4, 4
        !call indicator_field(q, q%m/2+i+q%m/8+2, q%n/2+j)
        call indicator_field(q, (q%m+1)/2+i, (q%n+1)/2+j)
      end do
    end do

    call set(q0, q)

    q_int_0 = int_field(q)

    call write(q, "./unittest/solverprop/q", 0.0_dp, 0)
    call write(psi, "./unittest/solverprop/psi", 0.0_dp, 0)
    ind = ind + 1

    call print_info(psi, dt)

    if (.not. rotsym(q)) then
      call abort_handle(" ! FAIL @ Rotational symmetry at initialisation", __FILE__, __LINE__)
    end if

    call print_array(psi%data, "psi")
    do i = 1, nts
      call timestep_heun_Kflux(q, dt, K_e)
      call timestep_LMT_2D(q, dt, uv_fld)

      if (mod(i, floor(nts/25.0_dp)) .eq. 0) then
        call write(q, "./unittest/solverprop/q", dt*i, ind)
        ind = ind + 1
      end if
    end do

    call print_info(q)

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


    ! Test Divergence consistent
    call deallocate(uv_fld)

    call set(q, 1.0_dp)
    call set(q0, q)

    call assemble_psi(psi, 2)  ! Ensure zero flux condition holds
    call allocate(uv_fld, psi)

    call assemble_K(K11, 2)  ! Ensure zero flux condition holds
    call assemble_K(K22, 2)  ! Ensure zero flux condition holds
    call assemble_K(K12, 2)  ! Ensure zero flux condition holds
    call scale(K11, 0.02_dp)
    call scale(K22, 0.03_dp)
    call scale(K12, 0.01_dp)
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
    
    write(6, "(a)") &
      "! ----------- Passed unittest_solver_properties ----------- !"
  end subroutine unittest_solver_properties
  
  subroutine unittest_fftw()
    use, intrinsic :: iso_c_binding 
    include 'fftw3.f03'
  
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
        
!        data_in(i, j) = dcos(Pi*(x+y))*dcos(2.0_dp*Pi*(x-y));
        data_in(i, j) = dcos(Pi*x)*dcos(2.0_dp*Pi*y) + dcos(4.0_dp*Pi*x)*dcos(6.0_dp*Pi*y)
        !data_in(i, j) = (x + y**2)*dsin(Pi*x)*dcos(2.0_dp*Pi*(x-y))
        !data_in(i, j) = dsin(4.0_dp*Pi*x)*dsin(2.0_dp*Pi*y) + dsin(8.0_dp*Pi*y-2.0_dp*Pi*x)
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
    
  end subroutine unittest_fftw
  
  subroutine unittest_solver_timing()
    type(field) :: K11, K22, K12, psi
    integer :: m, n
    integer :: i

    type(uv_signed) :: uv_fld
    type(K_flux) :: K_e

    integer, parameter :: gl = 1

    real(kind=dp) :: t = 12.0_dp*3.1415926_dp
    real(kind=dp) :: dt
    integer :: nts = 64*16

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
    
    call manual_field(psi, 7)  ! TG

    call set(K11, 2.0_dp)
    call set(K22, 1.0_dp)
    call set(K12, 0.5_dp)

    call allocate(uv_fld, psi)
    call allocate(K_e, K11, K22, K12)
    call imposeBC(K_e)
    
    call manual_field(q, 4)
    call imposeBC(q)
    call set(qLW, q)
    call set(q2D, q)
    dt = t / real(nts, kind=dp)
    
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
    integer :: nts 
      
!     ! Test spatial accuarcy (integrating to steady state)
!     write(6, *) " Timestepping to steady state"
!     write(6, *) " nts = 512*m "
!     do i = 1, 8
!       m = 2 ** (i+3)
!       !write(6, *) "m(",i,") =", m, "; error(",i,") = ", compute_spatial_advection_err(m, 2)  ! 2: stat state
!     end do
    
    ! Test temporal accuarcy (periodic solution)
    write(6, *) " Timestepping to periodic state"
    write(6, *) " Scale dt = C dx"
    write(6, *) " nts = 128*m "
    do i = 1, 8
      m = 2 ** (i+3)
      nts = 128* m
      write(6, *) "nts(",i,") =", nts, "; error(",i,") = ", compute_temporal_err(m, nts)  ! 1: periodic
    end do
    
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

    real(kind=dp), parameter :: T = 2.0_dp
   
    real(kind=dp), parameter :: Pi =  4.0_dp* datan (1.0_dp)
   
    type(field) :: q_exact, q_err
    real(kind=dp) :: t_i, err_int, err_sup
    
    n = m
    call allocate(q, m, n, 'q', glayer=gl, type_id=2)
    call allocate(q_exact, m, n, 'q_exact', glayer=gl, type_id=2)
    call allocate(q_err, m, n, 'q_err', glayer=gl, type_id=2)
    call allocate(psi, m, n, 'psi', glayer=0, type_id=1)

    call allocate(K11, m, n, 'K11', glayer=0, type_id=2)
    call allocate(K22, m, n, 'K22', glayer=0, type_id=2)
    call allocate(K12, m, n, 'K12', glayer=0, type_id=2)

    ! test_id = 1: q(t) 1-periodic
    call testcase_fields(psi, K11, K22, K12, q_exact, 1)
    call imposeBC(q_exact)
    
    ! TODO
!     call set(K11, 0.0_dp)
!     call set(K22, 0.0_dp)
!     call set(K12, 0.0_dp)    
    
    call allocate(uv_fld, psi)
    call allocate(K_e, K11, K22, K12)
    call imposeBC(K_e)
    
    ! Set inital condition to zero
    call set(q, q_exact)
    call imposeBC(q)
    
    err_sup = 0.0_dp
    
    call write(q, "./unittest/convergence/q_initial", 0.0_dp, 0)
    dt = T/real(nts,kind=dp)
    do i = 1, nts
      t_i = real(i-1, kind=dp)*dt
      
      ! Strang splitting
      ! A
      call timestep_heun_Kflux_Src_nonstat(q, 0.5_dp*dt, t_i, K_e)
      !call timestep_heun_Kflux_Src_nonstat(q, dt, t_i, K_e)
      ! B
      call timestep_LMT_2D(q, dt, uv_fld)
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
    call write(q, "./unittest/convergence/q_final", nts*dt, nts)


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
  
  pure logical function rotsym(fld)
    type(field), intent(in) :: fld
    integer :: m, n, i, j, gl

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
          if (dabs(fld%data(i+gl, j+gl)-fld%data(m-i+1+gl, n-j+1+gl)) > 1D-14) then
            rotsym = .false.
          end if
          ! Rotate by 3pi/2: x->y;  y->-x
          if (dabs(fld%data(i+gl, j+gl)-fld%data(j+gl, m-i+1+gl)) > 1D-14) then
            rotsym = .false.
          end if
          ! Rotate by pi/2: x->-y;  y->x
          if (dabs(fld%data(i+gl, j+gl)-fld%data(n-j+1+gl, i+gl)) > 1D-14) then
            rotsym = .false.
          end if
        end do
      end do
    end if

  end function rotsym
  
  subroutine unittest_IO()
    type(trajdat) :: traj
    real(kind=dp) :: t
    
    integer :: m, n, reflv
    type(meshdat) :: mesh

    type(dofdat) :: dof

    real(kind=dp), parameter :: L = 3840.0_dp*1000.0_dp
    real(kind=dp), parameter :: kappa_scale = 10000.0_dp
!     real(kind=dp), parameter :: psi_scale = 0.05_dp*L
    real(kind=dp), parameter :: psi_scale = 0.1_dp*L
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
    call set(dof%K12, 0.0_dp*sc)

    call print_info(mesh)
    call print_info(dof)
    
    call write_theta(dof, "./unittest/test_theta", sc, 1)

    call print_info(dof)

    write(6, "(a)") &
      "! ----------- Passed unittest_IO ----------- !"
  end subroutine unittest_IO
  
  subroutine unittest_inference_OMP()
    use omp_lib
    type(jumpsdat), dimension(:), pointer :: jumps
    type(trajdat) :: traj
    real(kind=dp) :: T

    integer :: m, n, cell, m_solver, n_solver, m_Ind, n_Ind
    type(meshdat) :: mesh

    type(dofdat) :: dof, dof_MAP, dof_old, dof_SSD
    type(IndFndat) :: IndFn

    real(kind=dp), parameter :: L = 3840.0_dp*1000.0_dp
    real(kind=dp), parameter :: kappa_scale = 10000.0_dp
!     real(kind=dp), parameter :: psi_scale = 0.05_dp*L
    real(kind=dp), parameter :: psi_scale = 0.1_dp*L
    real(kind=dp) :: sc, h, dt
    integer :: nts
    real(kind=dp), dimension(:), allocatable :: logPost, logPost_old
    real(kind=dp) :: SlogPost_old, SlogPost, SlogPost_MAP

    integer :: dof_id
    integer :: iter, niter, ind
    real(kind=dp) :: alphaUniRV
    real(kind=dp), dimension(1) :: UniRV

    integer, parameter :: Td = 30
!    character(len = *), parameter :: RunProfile = "K_TG_iso"
!    character(len = *), parameter :: RunProfile = "K_sinusoidal"
!    character(len = *), parameter :: RunProfile = "TTG_const"
    character(len = *), parameter :: RunProfile = "TTG_sinusoidal"
!     character(len = *), parameter :: RunProfile = "K_const"
!    character(len = *), parameter :: RunProfile = "TTG100_sinusoidal"
!     character(len = *), parameter :: RunProfile = "QGM2_L1_NPART676"
!    character(len = *), parameter :: RunProfile = "QGM2_L1_NPART2704"
    character(len = 256) :: output_fld, Td_char, resol_param

    ! Timer
    type(timer), save :: total_timer, commun_timer, propose_timer

    ! Timer
    call reset(total_timer)
    call reset(commun_timer)
    call reset(propose_timer)
    call start(total_timer)

    ! m, n: DOF/control
    m = 4
    n = 8
    ! m_solver: solver grid
    m_solver = 16
    ! m_Ind, n_Ind: indicator functions
    m_Ind = 8
    n_Ind = m_Ind
    
    write(Td_char, "(a,i0,a)") "h", Td, "d"
    write(resol_param, "(a,i0,a,i0,a,i0)") "N",m_solver*m_solver,"_D", m*n, "_I", m_Ind*n_Ind
    write(output_fld, "(a,a,a,a,a,a,a)") "./output/", trim(resol_param), "/", trim(RunProfile), "/", trim(Td_char), "/"
    write(6, "(a, a)") "Output path: ", trim(output_fld)
    
    ! AD-HOC
    write(output_fld, "(a)") "./output/test/"
    ! END AD-HOC

    call allocate(mesh, m_solver, m_solver)  ! Needs to be identical in both directions
    call allocate(dof, m-1, n-1, m, n) 
    ! N.B. assumed zeros at boundaries of psi-> hence only need (m-1)*(n-1) instead of (m+1)*(n+1)
    call allocate(IndFn, m_Ind, n_Ind)

    ! Initialise fields
    sc = real(mesh%m,kind=dp)/L
    call zeros(dof%psi)

    call set(dof%K11, 0.5_dp*kappa_scale*sc)
    call set(dof%K22, 0.5_dp*kappa_scale*sc)
    call set(dof%K12, 0.0_dp*sc)

    ! Read trajectory
    call read_header(traj, "./trajdat/"//RunProfile//"/"//trim(Td_char), T)
    call read(traj, "./trajdat/"//RunProfile//"/"//trim(Td_char))

    ! Ensure the positions are normalised
    if (.not. check_normalised_pos(traj)) then
      call abort_handle("E: particle out of range [0, 1]", __FILE__, __LINE__)
    end if

    ! Allocate array of jumps with initial positions
    call allocate(jumps, mesh)

    ! Convert trajectory into jumps
    call traj2jumps(jumps, traj, mesh)
    !call print_info(jumps, mesh)
    
    ! Determine h: Assumed h = uniform; assumed mesh%m = mesh%n
    h = read_uniform_h(jumps) *real(mesh%m, kind=dp)   ! *m to rescale it to [0, mesh%m]

    ! Calculate nts = number of timesteps needed in solving FP
    dt = 1.5_dp*3600.0_dp  ! 1.5 hours, in seconds
    dt = 3.0_dp*3600.0_dp  ! 1.5 hours, in seconds
    nts = int((h/sc)/dt)    ! 2 hours time step

    call print_info(mesh)
    call print_info(dof)
    call print_info(IndFn)
    call print_info(h, nts, sc)
    

!   likelihood for each indicator function
    allocate(logPost(IndFn%nIND))
    allocate(logPost_old(IndFn%nIND))

    ! MCMC
    call allocate(dof_old, dof%m_psi, dof%n_psi, dof%m_K, dof%n_K)
    call allocate(dof_MAP, dof%m_psi, dof%n_psi, dof%m_K, dof%n_K)
    call allocate(dof_SSD, dof%m_psi, dof%n_psi, dof%m_K, dof%n_K)

    call init_random_seed();
    
    ! Set up stepsize for random sampling
    call set(dof_SSD%psi, 0.05_dp*psi_scale*sc)

    call set(dof_SSD%K11, 0.0125_dp*kappa_scale*sc)
    call set(dof_SSD%K22, 0.0125_dp*kappa_scale*sc)
    call set(dof_SSD%K12, 0.050_dp)
    
    ! Initialise
    call evaluate_loglik_OMP(logPost, jumps, IndFn, mesh, dof, h, nts)
    SlogPost = sum(logPost)
    
    write(6, *) "SlogPost = ", SlogPost
    
    SlogPost_MAP = SlogPost
    SlogPost_old = SlogPost
    logPost_old = logPost

    call set(dof_old, dof)
    call set(dof_MAP, dof)

    ind = 0
    write(6, "(i0, a, "//dp_chr//")") 0, "-th step: logPost = ", SlogPost_old
    FLUSH(6)

    call write_theta(dof_old, trim(output_fld)//"theta", sc, ind)
    call reset_timestep_timers()
    call reset_inference_timers()

    ! TODO: delete Adhoc
    ! write(6, *) "propose_dof_K: isotropic diffusivity"
    ! TODO: delete Adhoc
    
    
    niter = max(500, m*m)
    niter = 8
    do iter = 1, niter
      ! In each iteration loop over all cells
      do dof_id = 1, 1
        ! Initialise K_prop
        call set(dof, dof_old)
        logPost = logPost_old

        ! Proposal
        call start(propose_timer)
        !call propose_dof(dof, dof_id, dof_SSD)   ! TODO: to implement
        call propose_dof_all(dof, iter, dof_SSD)
        call stop(propose_timer)

!         call print_array(dof%psi%data, "psi")
!         call print_array(dof%K11%data, "K11")
!         call print_array(dof%K22%data, "K22")
!         call print_array(dof%K12%data, "K12")
          
        ! Evaluate likelihood
        call evaluate_loglik_OMP(logPost, jumps, IndFn, mesh, dof, h, nts)
        SlogPost = sum(logPost)
        
        ! Metropolis-Hastings
        alphaUniRV = dexp(SlogPost - SlogPost_old);
        call RANDOM_NUMBER(UniRV)

        if (UniRV(1) < alphaUniRV) then
          logPost_old = logPost
          SlogPost_old = SlogPost
          call set(dof_old, dof)
        end if

        ! Record MAP
        if (SlogPost > SlogPost_MAP) then
          SlogPost_MAP = SlogPost
          call set(dof_MAP, dof)
        end if
      end do

      ! I/O
      if (mod(iter, 4) .eq. 0) then
        write(6, "(i0, a, "//dp_chr//")") iter, "-th step: logPost = ", SlogPost_old
        FLUSH(6)
        
        ind = ind + 1
        call write_theta(dof_old, trim(output_fld)//"theta", sc, ind)
      end if
    end do

    ! Write MAP
    call write_theta(dof_MAP, trim(output_fld)//"theta", sc, -1)

    ! Release memory
    call deallocate(traj)
    call deallocate(dof)
    call deallocate(IndFn)
    call deallocate(jumps)
    
    ! MH
    deallocate(logPost_old)
    deallocate(logPost)
    call deallocate(dof_old)
    call deallocate(dof_MAP)
    call deallocate(dof_SSD)

    ! Timer
    call stop(total_timer)
    call print_timestep_timers()
    call print_inference_timers()

    write(6, "(a)") "Proposal timers:"
    call print(propose_timer, "Proposal on node 0", prefix = "  ")

    write(6, "(a)") "Communication timers:"
    call print(commun_timer, "communication on node 0", prefix = "  ")

    write(6, "(a)") "Total timers:"
    call print(total_timer, "Total time on a node", prefix = "  ")

    write(6, "(a)") &
    "! ----------- Passed unittest_inference_OMP ----------- !"

  end subroutine unittest_inference_OMP

!   subroutine unittest_FPsolver_INPUT()
!     use mpi
!     type(jumpsdat), dimension(:), pointer :: jumps
!     type(trajdat) :: traj
!     real(kind=dp) :: T
! 
!     integer :: m, n, cell, reflv
!     type(meshdat) :: mesh
! 
!     type(dofdat) :: dof
!     type(IndFndat) :: IndFn
! 
!     real(kind=dp), parameter :: L = 3840.0_dp*1000.0_dp
!     real(kind=dp), parameter :: kappa_scale = 10000.0_dp
!     real(kind=dp), parameter :: psi_scale = 0.05_dp*L
!     real(kind=dp) :: sc, h, dt
!     integer :: nts
!     real(kind=dp), dimension(:), allocatable :: logPost
! 
!     integer, parameter :: Td = 30
!     character(len = *), parameter :: RunProfile = "TTG_sinusoidal"
!     character(len = 256) :: output_fld, Td_char, q_output_fld
! 
!     ! MPI
!     integer :: ierr, my_id, num_procs
! 
!     ! FP
!     integer, dimension(:), allocatable :: klist
!     type(dofdat) :: rdof
!     type(field) :: q
!     integer :: INDk, IND_read, write_ind
!     
!     ! Formulate dense matrix
!     integer :: densemat_switch, dof_id
!     real(kind=dp), dimension(:,:), allocatable :: QK
!     real(kind=dp), dimension(:), allocatable :: q_unfold
!     integer :: i0, j0, cmp, k
! 
!     call MPI_INIT( ierr )
!     ! find out MY process ID, and how many processes were started.
!     call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
!     call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
! 
!     ! m, n: For DOF
!     m = 16
!     n = m
!     reflv = 4
!     
!     write(output_fld, "(a,i0,a,a,a)") "./output/output", m, "/", RunProfile, "/"
!     write(q_output_fld, "(a,i0,a,a,a)") "./unittest/q", m*reflv, "/", RunProfile, "/"
!     if (my_id .eq. (num_procs-1)) then
!       write(6, "(a, a)") "Reading from Output path: ", trim(output_fld)
!     end if
! 
!     call allocate(mesh, m*reflv, n*reflv)
!     call allocate(IndFn, 8, 8)
!  
!     ! Initialise fields
!     sc = real(mesh%m,kind=dp)/L
!     
!     ! Read dof
!     IND_read = 50
!     call read(dof, output_fld, sc, IND_read)
!     !call set(dof%K12, 0.0_dp)
!     
!     ! Read trajectory
!     write(Td_char, "(a,i0,a)") "h", Td, "d"
!     call read_header(traj, "./trajdat/"//RunProfile//"/"//trim(Td_char), T)
!     call read(traj, "./trajdat/"//RunProfile//"/"//trim(Td_char))
! 
!     ! Ensure the positions are normalised
!     if (.not. check_normalised_pos(traj)) then
!       call abort_handle("E: particle out of range [0, 1]", __FILE__, __LINE__)
!     end if
! 
!     ! Allocate array of jumps with initial positions
!     call allocate(jumps, mesh)
! 
!     ! Convert trajectory
!     call traj2jumps(jumps, traj, mesh)
! 
!     ! Determine h: Assumed h = uniform; assumed mesh%m = mesh%n
!     ! *m to rescale it from [0, 1] to [0, mesh%m]
!     h = read_uniform_h(jumps) *real(mesh%m, kind=dp)
!     
!     ! Calculate nts = number of timesteps needed in solving FP
!     dt = 1.5_dp*3600.0_dp  ! 6 hours, in seconds
!     nts = int((h/sc)/dt)
! 
!     if (my_id .eq. (num_procs-1)) then
!       call print_info(dof)
!       call print_info(mesh)
!       call print_info(h, nts, sc)
!     end if
! 
!     !   likelihood for each indicator function
!     allocate(logPost(IndFn%nIND))
!     
!     ! Choose the tracer
!     INDk = 37
!     
!     call allocate(rdof, dof, mesh%m_reflv, mesh%n_reflv)
!     call allocate(q, rdof%m, rdof%n, 'q', glayer=1,type_id=2)
!     
!     ! Solve FK equation
!     call INDk_to_klist(klist, INDk, IndFn, mesh)
!     
!     call initialise_q(q, klist, mesh)
!     call write(q, trim(q_output_fld)//"q", 0.0_dp, 0)
! 
!     do write_ind = 1, nts
!       call advdiff_q(q, rdof%psi, rdof%K11, rdof%K22, rdof%K12, h/real(nts, kind=dp), 1)
!   !     call advdiff_q(q, rdof%psi, rdof%K11, rdof%K22, rdof%K12, h, nts)
!       call write(q, trim(q_output_fld)//"q", real(write_ind, kind=dp), write_ind)
! 
!       if (my_id .eq. (num_procs-1)) then
!         call print_info(q)
!         call print_info(rdof%psi, h/real(nts, kind=dp))
!       end if
!     end do
!     
!     call deallocate(rdof)
!     call deallocate(q)
!     deallocate(klist)
! 
!     ! Release memory
!     call deallocate(traj)
!     call deallocate(dof)
! 
! 
!     do cell = 1, mesh%ncell
!       call deallocate(jumps(cell))
!     end do
! 
!     deallocate(logPost)
!     
!     ! Assemble Qk
!     densemat_switch = 1
!     if (densemat_switch .eq. 1) then
!       write(6, *) "Now assemble Qk"
!       m = 16
!       n = m
!       reflv = 4
!       call allocate(mesh, m*reflv, n*reflv)
!       call allocate(dof, m, n, m, n)
!       
!       ! To determine size of Qk
!       call allocate(rdof, dof, mesh%m_reflv, mesh%n_reflv)
!       call allocate(q, rdof%m, rdof%n, 'q', glayer=1,type_id=2)
! 
!       allocate(Qk(size(q%data,1)*size(q%data,2), dof%ndof*3))
!       allocate(q_unfold(size(q%data,1)*size(q%data, 2)))
!       
!       call deallocate(q)
!       call deallocate(rdof)
!       
!       INDk = 37
!       call INDk_to_klist(klist, INDk, IndFn, mesh)
!       write(6, *) klist
!       
!       do dof_id = 1, dof%ndof
!         i0 = k2i(dof_id, dof)
!         j0 = k2j(dof_id, dof)
!         call zeros(dof%psi)
!           
!         do cmp = 1, 3
!           call zeros(dof%K11)
!           call zeros(dof%K22)
!           call zeros(dof%K12)
!           if (cmp .eq. 1) dof%K11%data(i0, j0) = 1.0_dp
!           if (cmp .eq. 2) dof%K22%data(i0, j0) = 1.0_dp
!           if (cmp .eq. 3) dof%K12%data(i0, j0) = 1.0_dp
!           
!           call allocate(rdof, dof, mesh%m_reflv, mesh%n_reflv)
!           call allocate(q, rdof%m, rdof%n, 'q', glayer=1,type_id=2)
!           call initialise_q(q, klist, mesh)
!           !call scale(q, 1000000.0_dp)
!           call advdiff_q(q, rdof%psi, rdof%K11, rdof%K22, rdof%K12, h, nts)
!           call print_info(h, nts, sc)
!           call print_info(q)
! 
!           call unfold_matrix(q_unfold, q%data)
!           do k = 1, size(q_unfold, 1)
!             Qk(k, 3*(dof_id-1)+cmp) = q_unfold(k)
!           end do
!           
!           call deallocate(q)
!           call deallocate(rdof)
!         end do
!       end do
!       
!       
!       !call print_array(Qk)
!       if (my_id .eq. (num_procs-1)) then
!         write(output_fld, "(a,i0,a)") "./unittest/Qk", m, "/"
!         write(6, *) "Now writing text file"
!         call write(Qk, trim(output_fld)//"Qk")
!       end if
!       
!       deallocate(q_unfold)
!       deallocate(Qk)
!     end if
!     
!     call deallocate(dof)
!     call deallocate(IndFn)
!     deallocate(klist)
!     
!     
!     call MPI_FINALIZE( ierr )
! 
!   end subroutine unittest_FPsolver_INPUT
!   
end program advdiff
