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
    call unittest_FPsolver_INPUT()
    !call unittest_solver_convergence()
    !call unittest_inference_OMP()
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
    
    ! Discrete sine transform
    ! Data input
    allocate(data_in(Nx+1, Ny+1))
    do j = 1, size(data_in, 2)
      do i = 1, size(data_in, 1)
        x = real(i-1, kind=dp)/real(Nx, kind=dp)
        y = real(j-1, kind=dp)/real(Ny, kind=dp)
        
!        data_in(i, j) = dcos(Pi*(x+y))*dcos(2.0_dp*Pi*(x-y));
        data_in(i, j) = dsin(Pi*x)*dsin(2.0_dp*Pi*y) + dsin(4.0_dp*Pi*x)*dsin(6.0_dp*Pi*y)
        !data_in(i, j) = (x + y**2)*dsin(Pi*x)*dcos(2.0_dp*Pi*(x-y))
        !data_in(i, j) = dsin(4.0_dp*Pi*x)*dsin(2.0_dp*Pi*y) + dsin(8.0_dp*Pi*y-2.0_dp*Pi*x)
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
    
!     call print_array(data_intpl, "data_intpl")
!     call print_array(data_exact, "data_exact")
    if (.not. iszero(data_exact-data_intpl)) then
      call abort_handle(" ! FAIL  @ Inconsistent sine_intpl", __FILE__, __LINE__)
    end if
    
    deallocate(data_in)
    deallocate(data_exact)
    deallocate(data_intpl)
    
    write(6, "(a)") &
      "! ----------- Passed unittest_fftw ----------- !"
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
    
    call print_array(dof%K22%data)
    
     call write_theta(dof, "./unittest/test_theta", sc, 1)

    call print_info(dof)
    
    call deallocate(dof)

    call read_theta(dof, "./unittest/test_theta", sc, 1)
    
    call print_info(dof)
    call print_array(dof%K22%data)
    
    call deallocate(dof)

    write(6, "(a)") &
      "! ----------- Passed unittest_IO ----------- !"
  end subroutine unittest_IO
  
  subroutine unittest_inference_OMP()
    use omp_lib
    type(jumpsdat), dimension(:), pointer :: jumps
    type(jumpsdat), dimension(:), pointer :: jumps_inv
    type(trajdat) :: traj
!     real(kind=dp) :: T

    integer :: m, n, cell, m_solver, m_Ind, n_Ind
    type(meshdat) :: mesh

    type(dofdat) :: dof, dof_MAP, dof_old, dof_SSD
    type(dofdat) :: dof_inv
    type(IndFndat) :: IndFn

    real(kind=dp), parameter :: L = 3840.0_dp*1000.0_dp
    real(kind=dp), parameter :: kappa_scale = 10000.0_dp
!     real(kind=dp), parameter :: psi_scale = 0.05_dp*L
    real(kind=dp), parameter :: psi_scale = 0.1_dp*L
    real(kind=dp) :: sc, h, dt
    integer :: nts
    real(kind=dp), dimension(:), allocatable :: logPost, logPost_old
    real(kind=dp), dimension(:), allocatable :: canon, canon_old

    integer :: canon_id
    integer :: iter, niter, ind
    real(kind=dp) :: alphaUniRV
    real(kind=dp), dimension(1) :: UniRV

    integer, parameter :: Td = 32
!    character(len = *), parameter :: RunProfile = "TTG_sinusoidal"
!     character(len = *), parameter :: RunProfile = "QGM2_L1_NPART676"
    character(len = *), parameter :: RunProfile = "QGM2_L1_NPART2704"
    character(len = 256) :: output_fld, Td_char, resol_param

    ! Timer
    type(timer), save :: total_timer, commun_timer, propose_timer

    ! dJdm
    real(kind=dp), dimension(:, :), allocatable :: dJdm, dJdm_sqsum, criteria
    
    ! Timer
    call reset(total_timer)
    call reset(commun_timer)
    call reset(propose_timer)
    call start(total_timer)

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
    call allocate(dof, m-1, n-1, 8+1, 8+1) 
!     call allocate(dof, 1, 1, 8+1, 8+1)   !! Ad-hoc No advection
    ! N.B. assumed zeros at boundaries of psi-> hence only need (m-1)*(n-1) instead of (m+1)*(n+1)
    call allocate(IndFn, m_Ind, n_Ind)

    ! Initialise fields
    sc = real(mesh%m,kind=dp)/L
    call zeros(dof%psi)
    
    call set(dof%K11, 5.0_dp*kappa_scale*sc)
    call set(dof%K22, 5.0_dp*kappa_scale*sc)
    call set(dof%K12, 0.0_dp*sc)
    

    ! Read trajectory
!     call read_header(traj, "./trajdat/"//RunProfile//"/"//trim(Td_char), T)
    call read(traj, "./trajdat/"//RunProfile//"/"//trim(Td_char))

    ! Ensure the positions are normalised
    if (.not. check_normalised_pos(traj)) then
      call abort_handle("E: particle out of range [0, 1]", __FILE__, __LINE__)
    end if

    ! Allocate array of jumps with initial positions
    call allocate(jumps, mesh)
    call allocate(jumps_inv, mesh)

    ! Convert trajectory into jumps
    call traj2jumps(jumps, traj, mesh)
    call traj2jumps_inv(jumps_inv, traj, mesh)
    call deallocate(traj)

    !call print_info(jumps, mesh)
    
    ! Determine h: Assumed h = uniform; assumed mesh%m = mesh%n
    h = read_uniform_h(jumps) *real(mesh%m, kind=dp)   ! *m to rescale it to [0, mesh%m]

    ! Calculate nts = number of timesteps needed in solving FP
    dt = 1.5_dp*3600.0_dp  ! 1.5 hours, in seconds
    dt = 3.0_dp*3600.0_dp  ! 1.5 hours, in seconds
!     dt = 6.0_dp*3600.0_dp  ! 1.5 hours, in seconds
    nts = int((h/sc)/dt)    ! 2 hours time step

    call print_info(mesh)
    call print_info(dof)
    call print_info(IndFn)
    call print_info(h, nts, sc)

!   likelihood for each indicator function
    allocate(logPost(IndFn%nIND))
    allocate(logPost_old(IndFn%nIND))
    call allocate(canon, dof)
    call allocate(canon_old, dof)
    
    allocate(dJdm(size(logPost,1), size(canon,1)))
    allocate(dJdm_sqsum(size(logPost,1), size(canon,1)))
    allocate(criteria(size(logPost,1), size(canon,1)))
    dJdm = 0.0_dp
    dJdm_sqsum = 0.0_dp
    criteria = 0.0_dp 
    
    ! MCMC
    call allocate(dof_inv, dof%m_psi, dof%n_psi, dof%m_K, dof%n_K)
    call allocate(dof_old, dof%m_psi, dof%n_psi, dof%m_K, dof%n_K)
    call allocate(dof_MAP, dof%m_psi, dof%n_psi, dof%m_K, dof%n_K)
    call allocate(dof_SSD, dof%m_psi, dof%n_psi, dof%m_K, dof%n_K)

    call init_random_seed();
    
    ! Set up stepsize for random sampling
    ! dof_SSD: By default zeros
    ! Obtain canonical
    do canon_id = 1, dof%ndof
      call propose_dof_canon(dof, canon, canon_id, dof_SSD)
    end do
    
    ! Tune proposal distribution: SSD
    call set(dof_SSD%psi, 0.05_dp*psi_scale*sc)
    call set(dof_SSD%K11, 0.1_dp*kappa_scale*sc)
    call set(dof_SSD%K22, 0.1_dp*kappa_scale*sc)
    call set(dof_SSD%K12, 0.1_dp)
    
    ! Set reverse psi
    call reverse_dof_to_dof_inv(dof_inv, dof)

    ! Initialise
    call evaluate_loglik_OMP(logPost, jumps, IndFn, mesh, dof, h, nts)
    dof%SlogPost = real(sum(logPost), kind=dp)
!     call evaluate_loglik_OMP(logPost, jumps_inv, IndFn, mesh, dof_inv, h, nts)
!     dof%SlogPost = dof%SlogPost + real(sum(logPost), kind=dp)
    
    call set(dof_old, dof)
    logPost_old = logPost
    canon_old = canon
    call set(dof_MAP, dof)

    ind = 0
    write(6, "(i0, a, "//dp_chr//")") 0, "-th step: logPost = ", dof_old%SlogPost
    FLUSH(6)

    call write_theta(dof_old, trim(output_fld)//"theta_canon", sc, ind)
    
    call reset_timestep_timers()
    call reset_inference_timers()

    ! TODO: delete Adhoc
    ! write(6, *) "propose_dof_K: isotropic diffusivity"
    ! TODO: delete Adhoc
    
    niter = max(2500, m*m)
    do iter = 1, niter
      ! In each iteration loop over all cells
      do canon_id = 1, dof%ndof
        ! Initialise K_prop
        call set(dof, dof_old)
        logPost = logPost_old
        canon = canon_old
        
        ! Proposal
        call propose_dof_canon(dof, canon, canon_id, dof_SSD)
!         call propose_dof(dof, canon_id, dof_SSD)
        call reverse_dof_to_dof_inv(dof_inv, dof)
!         !! Ad-hoc No advection
!         dof%psi%data(1,1) = 0.0_dp  ! Ad-hoc No advection
!         !! Ad-hoc No advection
        !call propose_dof_all(dof, iter, dof_SSD)
        
        ! Evaluate likelihood
!        call evaluate_loglik_OMP(logPost, jumps, IndFn, mesh, dof, h, nts)
        call evaluate_loglik_guided(logPost, jumps, IndFn, mesh, dof, h, nts, canon_id, criteria)
        dof%SlogPost = real(sum(logPost), kind=dp)
!         call evaluate_loglik_OMP(logPost, jumps_inv, IndFn, mesh, dof_inv, h, nts)
!         dof%SlogPost = dof%SlogPost + real(sum(logPost), kind=dp)
        
        ! Evaluate dJdm
        call compute_dJdm_Cid(dJdm, logPost, logPost_old, canon, canon_old, canon_id)
        
        ! Metropolis-Hastings
        alphaUniRV = dexp(dof%SlogPost - dof_old%SlogPost)
        call RANDOM_NUMBER(UniRV)

        if (UniRV(1) < alphaUniRV) then
          call set(dof_old, dof)
          logPost_old = logPost
          canon_old = canon
          dof_old%occurrence = 1
        else
          dof%occurrence = dof%occurrence + 1
        end if

        ! Record MAP
        if (dof%SlogPost > dof_MAP%SlogPost) then
          call set(dof_MAP, dof)
        end if
      end do
      
      ! dJdm
      dJdm_sqsum = dJdm_sqsum + dJdm*dJdm
      dJdm = 0.0_dp
      
      ! I/O
      if (mod(iter, 10) .eq. 0) then
        write(6, "(i0, a, "//dp_chr//")") iter, "-th step: logPost = ", dof_old%SlogPost
        FLUSH(6)
        
        ind = ind + 1
        call write_theta(dof_old, trim(output_fld)//"theta_canon", sc, ind)
        
        if (iter .eq. 50) then
          criteria = dJdm_sqsum
          dJdm_sqsum = 0.0
          call write(criteria, trim(output_fld)//"criteria")
        end if
      end if
    end do

    ! Write MAP
    call write_theta(dof_MAP, trim(output_fld)//"theta_canon", sc, -1)

    ! Release memory
    ! MH
    deallocate(dJdm)
    deallocate(dJdm_sqsum)
    deallocate(criteria)
    deallocate(canon)
    deallocate(canon_old)
    
    deallocate(logPost_old)
    deallocate(logPost)
    call deallocate(dof_old)
    call deallocate(dof_MAP)
    call deallocate(dof_SSD)

    call deallocate(dof)
    call deallocate(dof_inv)
    call deallocate(IndFn)
    call deallocate(jumps)
    call deallocate(jumps_inv)
    
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

  subroutine unittest_FPsolver_INPUT()
      use omp_lib
    type(jumpsdat), dimension(:), pointer :: jumps
    type(trajdat) :: traj
    real(kind=dp) :: T

    integer :: m, n, cell, m_solver, m_Ind, n_Ind
    type(meshdat) :: mesh

    type(dofdat) :: dof, dof_MAP, dof_old, dof_SSD
    type(IndFndat) :: IndFn

    real(kind=dp), parameter :: L = 3840.0_dp*1000.0_dp
    real(kind=dp), parameter :: kappa_scale = 10000.0_dp
!     real(kind=dp), parameter :: psi_scale = 0.05_dp*L
    real(kind=dp), parameter :: psi_scale = 0.1_dp*L
    real(kind=dp) :: sc, h, dt
    integer :: nts, nts_sub, dt_rf
    real(kind=dp), dimension(:), allocatable :: logPost, logPost_old

    integer :: dof_id
    integer :: iter, niter, ind

    integer, parameter :: Td = 32
    character(len = *), parameter :: RunProfile = "QGM2_L1_NPART2704"
    character(len = 256) :: output_fld, Td_char, resol_param
    character(len = 256) :: output_q

    type(dofdat) :: dof_solver
    type(field) :: q
    real(kind=dp) :: lik
    integer, dimension(:), allocatable :: klist
    integer :: INDk, i
    
    real(kind=dp) :: SlogPost, SlogPostr, alphaUniRV
    
    real(kind=dp), dimension(:,:), allocatable :: dJdm, dJdm_abs
    real(kind=dp), dimension(:), allocatable :: theta, theta_old
    
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

    call read_theta(dof, trim(output_fld)//"theta_canon", sc, 38)
    
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
    
    INDk = m_Ind*m_Ind/2+2
    call INDk_to_klist(klist, INDk, IndFn, mesh)
    call allocate(q, mesh%m, mesh%n, 'q', glayer=1, type_id=2)
    
    !! dt = 3 hours
    dt_rf = 1
    dt = 3.0_dp*3600.0_dp/real(dt_rf, kind=dp)  ! 1.5 hours, in seconds
    nts = int((h/sc)/dt)    ! 2 hours time step
    call initialise_q(q, klist, mesh)

    call write(q, "./unittest/q64/QGM2_L1/q32_sine1", 0.0_dp, 0)

    call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h*0.75_dp, 3*nts/4)
    call write(q, "./unittest/q64/QGM2_L1/q32_sine1", h*0.75_dp, 1)
!      3hr@Posterior =   -1498342.09318808     
!  0.375hr@Posterior =   -1498396.61044673

    ! Output video
    INDk = m_Ind*m_Ind/2+2
    call INDk_to_klist(klist, INDk, IndFn, mesh)
    write(6, *) INDk
    write(6, *) h/sc/(3600.0_dp*24.0_dp)
    call initialise_q(q, klist, mesh)
    
    write(output_q, "(a, i0)") "./unittest/q64/QGM2_L1/q", INDk
    call write(q, trim(output_q), 0.0_dp, 0)

    do i = 1, nts
      call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h/nts, 1)
      call write(q, trim(output_q), (h/sc)/nts*i, i)
    end do
    
    INDk = m_Ind*m_Ind/2+m_Ind+2
    call INDk_to_klist(klist, INDk, IndFn, mesh)
    write(6, *) INDk
    write(6, *) h/sc/(3600.0_dp*24.0_dp)
    call initialise_q(q, klist, mesh)
    
    write(output_q, "(a, i0)") "./unittest/q64/QGM2_L1/q", INDk
    call write(q, trim(output_q), 0.0_dp, 0)

    do i = 1, nts
      call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h/nts, 1)
      call write(q, trim(output_q), (h/sc)/nts*i, i)
    end do
    
    INDk = m_Ind*m_Ind/2+m_Ind+2
    call INDk_to_klist(klist, INDk, IndFn, mesh)
    write(6, *) INDk
    write(6, *) h/sc/(3600.0_dp*24.0_dp)
    call initialise_q(q, klist, mesh)
    
    write(output_q, "(a, i0)") "./unittest/q64/QGM2_L1/qr", INDk
    call write(q, trim(output_q), 0.0_dp, 0)

    do i = 1, nts
      call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h/nts, 32)
      call write(q, trim(output_q), (h/sc)/nts*i, i)
    end do

    INDk = m_Ind*m_Ind/2-m_Ind+2
    call INDk_to_klist(klist, INDk, IndFn, mesh)
    write(6, *) INDk
    write(6, *) h/sc/(3600.0_dp*24.0_dp)
    call initialise_q(q, klist, mesh)
    
    write(output_q, "(a, i0)") "./unittest/q64/QGM2_L1/q", INDk
    call write(q, trim(output_q), 0.0_dp, 0)

    do i = 1, nts
      call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h/nts, 1)
      call write(q, trim(output_q), (h/sc)/nts*i, i)
    end do

    INDk = m_Ind*m_Ind/2-m_Ind+1
    call INDk_to_klist(klist, INDk, IndFn, mesh)
    write(6, *) INDk
    write(6, *) h/sc/(3600.0_dp*24.0_dp)
    call initialise_q(q, klist, mesh)
    
    write(output_q, "(a, i0)") "./unittest/q64/QGM2_L1/q", INDk
    call write(q, trim(output_q), 0.0_dp, 0)

    do i = 1, nts
      call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h/nts, 1)
      call write(q, trim(output_q), (h/sc)/nts*i, i)
    end do
    
    INDk = m_Ind*m_Ind/2+1
    call INDk_to_klist(klist, INDk, IndFn, mesh)
    write(6, *) INDk
    write(6, *) h/sc/(3600.0_dp*24.0_dp)
    call initialise_q(q, klist, mesh)
    
    write(output_q, "(a, i0)") "./unittest/q64/QGM2_L1/q", INDk
    call write(q, trim(output_q), 0.0_dp, 0)

    do i = 1, nts
      call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h/nts, 1)
      call write(q, trim(output_q), (h/sc)/nts*i, i)
    end do
    
    INDk = m_Ind*m_Ind/2+m_Ind+1
    call INDk_to_klist(klist, INDk, IndFn, mesh)
    write(6, *) INDk
    write(6, *) h/sc/(3600.0_dp*24.0_dp)
    call initialise_q(q, klist, mesh)
    
    write(output_q, "(a, i0)") "./unittest/q64/QGM2_L1/q", INDk
    call write(q, trim(output_q), 0.0_dp, 0)

    do i = 1, nts
      call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h/nts, 1)
      call write(q, trim(output_q), (h/sc)/nts*i, i)
    end do
    
    INDk = m_Ind*m_Ind/2+m_Ind+1
    call INDk_to_klist(klist, INDk, IndFn, mesh)
    write(6, *) INDk
    write(6, *) h/sc/(3600.0_dp*24.0_dp)
    call initialise_q(q, klist, mesh)
    
    write(output_q, "(a, i0)") "./unittest/q64/QGM2_L1/qr", INDk
    call write(q, trim(output_q), 0.0_dp, 0)

    do i = 1, nts
      call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h/nts, 256)
      call write(q, trim(output_q), (h/sc)/nts*i, i)
    end do
    
    INDk = m_Ind*m_Ind/2+m_Ind/2
    call INDk_to_klist(klist, INDk, IndFn, mesh)
    write(6, *) INDk
    write(6, *) h/sc/(3600.0_dp*24.0_dp)
    call initialise_q(q, klist, mesh)
    
    write(output_q, "(a, i0)") "./unittest/q64/QGM2_L1/q", INDk
    call write(q, trim(output_q), 0.0_dp, 0)

    do i = 1, nts
      call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h/nts, 1)
      call write(q, trim(output_q), (h/sc)/nts*i, i)
    end do
    
    call print_info(dof_solver%psi, h/real(nts,kind=dp))
    
    !! dt  = 0.375 hours
    dt_rf = 8
    dt = 3.0_dp*3600.0_dp/real(dt_rf, kind=dp)  ! 1.5 hours, in seconds
    nts = int((h/sc)/dt)    ! 2 hours time step
    call initialise_q(q, klist, mesh)

    call write(q, "./unittest/q64/QGM2_L1/q32_sine8", 0.0_dp, 0)

    call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h, nts)
    call write(q, "./unittest/q64/QGM2_L1/q32_sine8", h, 1)
    
    call print_info(dof_solver%psi, h/real(nts,kind=dp))
    

    
    ! Test canon conversion
    call allocate(dof_old, dof%m_psi, dof%n_psi, dof%m_K, dof%n_K)
    call set(dof_old, dof)
    
    call allocate(theta, dof)
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
  
end program advdiff
