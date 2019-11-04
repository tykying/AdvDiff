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
  use advdiff_unittest

  implicit none

  character(len=32) :: arg
  integer :: Td, layer, NPart, Phase, Seed_ID
  character(len=128) :: path_arg
  
  ! Main program
  if (command_argument_count() > 0) then
!     write(6, "(a)") &
!       "! ----------- Reading Command Line Arguments ----------- !"      
    call get_command_argument(1, arg)
    READ(arg,*) Td     ! cast string into integer
  
    call get_command_argument(2, arg)
    READ(arg,*) layer

    call get_command_argument(3, arg)
    READ(arg,*) NPart
    
    call get_command_argument(4, arg)
    READ(arg,*) Phase
    
    call get_command_argument(5, arg)
    READ(arg,*) Seed_ID
    
    call get_command_argument(6, path_arg)
    
!     write(6, "(a, i0, a, i0, a, i0, a, i0, a, i0)") " SampInt = ", Td, &
!               & "d; Layer = ", layer, "; NPart = ", NPart, "; Phase = ", Phase, & 
!               & "; Seed_ID = ", Seed_ID
  else
    call abort_handle("Invalid input to the programme", __FILE__, __LINE__)
  end if
  
  call inference(Td, layer, NPart, Phase, Seed_ID, path_arg)
  !call unittest()

contains

  subroutine unittest()
    write(6, *) "!!! ----------- Unittests begin ----------- !!!"
    call unittest_comptools()
    !call unittest_trajdata()
    !call unittest_solver_properties()
    !call unittest_fftw()
    !call unittest_solver_timing()
    !call unittest_IO()
    !call unittest_timer()
    !call unittest_FPsolver_INPUT()
    !call unittest_solver_convergence()
    !call unittest_optimisation_OMP()
    !call unittest_TG_instability()
    write(6, *) "!!! ----------- All unittests passed ----------- !!!"
  end subroutine unittest

  subroutine validation()
    integer :: Td_i = 32
    integer :: layer = 2
    
    do Td_i = 1, 3
      call validate_inference(Td=Td_i*32, layer=layer)
    end do
    write(6, *) "!!! ----------- Validation finished ----------- !!!"
  end subroutine validation
  
#define OMP0MPI1 1
#define IBACKWARD 0
#define EM_MEAN 0

  subroutine inference(Td, layer, NPart, Phase, Seed_ID, path_arg)
#if OMP0MPI1 == 0
    use omp_lib
#elif OMP0MPI1 == 1
    use mpi
#endif
    integer, intent(in) :: Td, layer, NPart, Phase, Seed_ID
    character(len=128), intent(in) :: path_arg
    
    type(jumpsdat), dimension(:), allocatable :: jumps
    type(trajdat) :: traj

    integer :: m, n, m_solver, m_Ind, n_Ind
    type(meshdat) :: mesh

    type(dofdat) :: dof, dof_MAP, dof_old
    type(IndFndat) :: IndFn

#if IBACKWARD == 1
    type(jumpsdat), dimension(:), allocatable :: jumps_inv
    type(dofdat) :: dof_inv
#endif
    
    real(kind=dp) :: sc, h, dt
    integer :: nts
    
    real(kind=dp), parameter :: L = 3840.0_dp*1000.0_dp
    integer, dimension(:), allocatable :: accept_counter
    real(kind=dp), dimension(:), allocatable :: canon_SSD
    
    integer :: niter = 500
    integer :: iter, canon_id, ind
    real(kind=dp) :: alphaUniRV

    character(len = 128) :: RunProfile
    character(len = 256) :: output_fld, Td_char, resol_param
    character(len = 256) :: input_fld, fld_tmp
    
    ! Timer
    type(timer), save :: total_timer
    
#if EM_MEAN == 1
    ! Use Eulerian time-average
    type(field) :: psi_EM
    real(kind=dp) :: t_avg
#endif
    
    ! Prior
    type(priordat) :: prior_param
    
    ! Random number generator
    real(kind=dp), dimension(:), allocatable :: stdnRV, UniRV
    integer :: RV_ind, accepted

    integer :: restart_ind = 0
    integer :: output_dn = 100

    ! MPI
    integer :: my_id, num_procs
#if OMP0MPI1 == 1
    integer :: ierr !, PROVIDED, !! PROVIDED: For MPI-OpenMP
    real(kind=dp) :: SlogPost_local
    real(kind=dp) :: SlogPost_global = 0.0_dp
#endif
    
    ! m, n: DOF/control
    if (index(path_arg, "QGM2") .gt. 0) then
      write(RunProfile, "(a,i0,a,i0)") "QGM2_L", layer, "_NPART", NPart
      m = 16
      n = 16
      ! m_solver: solver grid
      m_solver = 64
      ! m_Ind, n_Ind: indicator functions
      m_Ind = 16
      n_Ind = m_Ind
    elseif (index(path_arg, "TTG_sinusoidal") .gt. 0) then
      write(RunProfile, "(a)") "TTG_sinusoidal"
      m = 8
      n = 8
      ! m_solver: solver grid
      m_solver = 64
      ! m_Ind, n_Ind: indicator functions
      m_Ind = 8
      n_Ind = m_Ind
    else
      ! Debug mode
      m = 2
      n = 2
      ! m_solver: solver grid
      m_solver = 8
      ! m_Ind, n_Ind: indicator functions
      m_Ind = 2
      n_Ind = m_Ind
      
      call abort_handle("Invalid path", __FILE__, __LINE__)
    end if

    write(Td_char, "(a,i0,a)") "h", Td, "d"
    write(resol_param, "(a,i0,a,i0,a,i0)") "N",m_solver*m_solver,"_D", m*n, "_I", m_Ind*n_Ind
    !write(fld_tmp, "(a,a,a,a,a,a)") "./output/", trim(resol_param), "/", trim(RunProfile), "/", trim(Td_char)
    
    write(fld_tmp, "(a,a)") trim(path_arg), "/"
!     write(6, "(a, a)") "New path: ", trim(fld_tmp)
    
    if (Phase .eq. 1) then
      niter = 500
      output_dn = 100
      restart_ind = 0
      write(input_fld, "(a,a)") trim(fld_tmp), "" 
      write(output_fld, "(a,a)") trim(fld_tmp), "SpinUp/"
    elseif (Phase .eq. 2) then
      niter = 500
      output_dn = 100
      restart_ind = -1
      write(input_fld, "(a,a)") trim(fld_tmp), "SpinUp/" 
      write(output_fld, "(a,a)") trim(fld_tmp), "Tuned/"
    elseif (Phase .eq. 3) then
      niter = 800
      output_dn = 1
      restart_ind = -1
      write(input_fld, "(a,a)") trim(fld_tmp), "Tuned/" 
      write(output_fld, "(a,a)") trim(fld_tmp), "Data/"
    end if
    
    call allocate(mesh, m_solver, m_solver)  ! Needs to be identical in both directions
    ! N.B. assumed zeros at boundaries of psi-> hence only need (m-1)*(n-1) instead of (m+1)*(n+1)
#if OMP0MPI1 == 0
    my_id = 0
    num_procs = 1
    
    !call omp_set_num_threads(16);
    call allocate(IndFn, m_Ind, n_Ind)
#elif OMP0MPI1 == 1
    call MPI_INIT( ierr )
    !call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, PROVIDED, ierr)
    !if (PROVIDED .ne. MPI_THREAD_FUNNELED) call abort_handle( &
    !    & "Failed to initialise MPI multi-threads", __FILE__, __LINE__)
    !call omp_set_num_threads(16/num_procs);
    !write(6, "(a,i0,a,i0)") " #run = ", num_procs, "; each with #threads = ", 16/num_procs
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
    if (ierr .ne. MPI_SUCCESS) then 
      call abort_handle("MPI Failed", __FILE__, __LINE__) 
    end if
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    if (ierr .ne. MPI_SUCCESS) then 
      call abort_handle("MPI Failed", __FILE__, __LINE__) 
    end if
    
    call allocate(IndFn, m_Ind, n_Ind, my_id, num_procs)
#endif
    
    ! Initialise fields
    sc = real(mesh%m,kind=dp)/L     ! Needs mesh to be identical in both directions

    if (restart_ind .ne. 0) then
      call read_theta(dof, trim(input_fld)//"theta", sc, restart_ind)
    else
      call allocate(dof, m-1, n-1, m+1, n+1)
      call init_dof(dof, sc)
    end if
    
#if EM_MEAN == 1
    !! Use Eulerian time-average
    call read_QGfield(psi_EM, t_avg, "./meanflow/psi_int_final", layer)
    
    ! Time-average + rescale wrt mesh
    psi_EM%data = (psi_EM%data - psi_EM%data(1,1))/t_avg
    psi_EM%data = psi_EM%data/(100.0_dp*real(psi_EM%m, kind=dp)) * real(mesh%m, kind=dp)  ! 100: from cm to m
    t_avg =  t_avg/(psi_EM%m/(L*100.0_dp))/(3600.0_dp*24.0_dp*365.25_dp)
    write(6, *) "Using Eulerian time-averaged mean flow of ", t_avg, " years"
    
    ! Filter out
    call sine_filter(dof%psi%data, &
      psi_EM%data(2:size(psi_EM%data,1)-1, 2:size(psi_EM%data,1)-1) )
    dof%ndof = 3*(dof%m_K*dof%n_K)
    call deallocate(psi_EM)
    !! End: Use Eulerian time-average
#endif

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
    
    ! Determine h: Assumed h = uniform; assumed mesh%m = mesh%n
    h = read_uniform_h(jumps) *real(mesh%m, kind=dp)   ! *m to rescale it to [0, mesh%m]

    ! Calculate nts = number of timesteps needed in solving FP
    dt = 6.0_dp*3600.0_dp  ! 6 hours, in seconds
    nts = int((h/sc)/dt)
    
!   likelihood for each indicator function
    allocate(accept_counter(dof%ndof))
    call allocate(canon_SSD, dof)
    
    ! MCMC
    call allocate(dof_old, dof%m_psi, dof%n_psi, dof%m_K, dof%n_K)
    call allocate(dof_MAP, dof%m_psi, dof%n_psi, dof%m_K, dof%n_K)
#if IBACKWARD == 1
    call allocate(jumps_inv, mesh)
    call traj2jumps_inv(jumps_inv, traj, mesh)
    call allocate(dof_inv, dof%m_psi, dof%n_psi, dof%m_K, dof%n_K)
    call reverse_dof_to_dof_inv(dof_inv, dof)   ! Set reverse psi
#endif

    call deallocate(traj)
    
    ! Tune proposal distribution: SSD
    accept_counter = 0
    
    if (Phase .eq. 1) then
      call init_canon_SSD(canon_SSD, dof, sc)
      ind = restart_ind
    else
      ! Read SSD
      call read(canon_SSD, trim(input_fld)//"canon_SSD")
      if (restart_ind .eq. -1) then
        ind = 0
      else
        ind = restart_ind
      end if
    end if
    
#if OMP0MPI1 == 0
    dof%SlogPost = evaluate_logPost_OMP(prior_param, jumps, IndFn, mesh, dof, h, nts)
#elif OMP0MPI1 == 1
    SlogPost_local = evaluate_logPost_OMP(prior_param, jumps, IndFn, mesh, dof, h, nts)
    
    SlogPost_global = 0.0_dp
    call MPI_AllReduce(SlogPost_local, SlogPost_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
             MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) then 
      call abort_handle("MPI Failed", __FILE__, __LINE__) 
    end if
    
    dof%SlogPost = SlogPost_global
#endif

    call set(dof_old, dof)
    call set(dof_MAP, dof)
    
    ! Generate random numbers
    allocate(stdnRV(niter*dof%ndof))
    allocate(UniRV(niter*dof%ndof))
    ! Seed must lie in (0, 4095), last entry must be odd
    call randn(stdnRV, (/ mod(1989*Seed_ID, 4096), mod(28*Seed_ID, 4096), &
                        & mod(11*m_solver, 4096), mod(Td*6*Phase, 4096)+1 /))
    call randu(UniRV, (/ mod(nts*(4+Phase), 4096), mod(26*n, 4096), &
                        & mod(8*niter, 4096), mod(8*Seed_ID*1991, 4096)+1 /))
    
    ! I/O
    if (my_id .eq. (num_procs-1)) then
      write(6, "(a, a)") "Output path: ", trim(output_fld)
      
      call print_info(mesh)
      call print_info(dof)
      call print_info(IndFn)
      call print_info(h, nts, sc)
!     call print_info(jumps, mesh)
      write(6, "(a)") " Normal RV samples: "
      call print_info(stdnRV)
      write(6, "(a)") " Uniform RV samples: "
      call print_info(UniRV)
      
      write(6, "(i0, a, "//dp_chr//")") ind, "-th step: logPost = ", dof_old%SlogPost
      FLUSH(6)
      
      call write_theta(dof_old, trim(output_fld)//"theta", sc, ind)
    end if
    
    ! Timer
    call reset(total_timer)
    call start(total_timer)
    
    do iter = 1, niter
      ! In each iteration loop over all cells
      do canon_id = 1, dof%ndof
        RV_ind = (iter-1)*dof%ndof + canon_id
        
        call set(dof, dof_old)
        
        ! Proposal
#if EM_MEAN == 1
        call propose_dof_EM(dof, canon_SSD, canon_id, stdnRV(RV_ind))
#else
        call propose_dof_canon(dof, canon_SSD, canon_id, stdnRV(RV_ind))
#endif
        
#if IBACKWARD == 1
        call reverse_dof_to_dof_inv(dof_inv, dof)
#endif
        
        ! Evaluate likelihood
#if OMP0MPI1 == 0
        dof%SlogPost = evaluate_logPost_OMP(prior_param, jumps, IndFn, mesh, dof, h, nts)
#elif OMP0MPI1 == 1
        SlogPost_local = evaluate_logPost_OMP(prior_param, jumps, IndFn, mesh, dof, h, nts)
        
        SlogPost_global = 0.0_dp
        call MPI_Reduce(SlogPost_local, SlogPost_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                 (num_procs-1), MPI_COMM_WORLD, ierr)
        if (ierr .ne. MPI_SUCCESS) then 
          call abort_handle("MPI Failed", __FILE__, __LINE__) 
        end if
        
        dof%SlogPost = SlogPost_global
#endif
        
        ! Metropolis-Hastings
        if (my_id .eq. (num_procs-1)) then
          alphaUniRV = dexp(dof%SlogPost - dof_old%SlogPost)
          
          if (UniRV(RV_ind) .lt. alphaUniRV) then
            accepted = 1
          else
            accepted = 0
          end if
        end if 

#if OMP0MPI1 == 1
        call MPI_Bcast(accepted, 1, MPI_INT, (num_procs-1), MPI_COMM_WORLD,  ierr)
        if (ierr .ne. MPI_SUCCESS) then 
          call abort_handle("MPI Failed", __FILE__, __LINE__) 
        end if
#endif
        
        if (accepted .eq. 1) then
          ! Accepted
          call set(dof_old, dof)
          dof_old%occurrence = 1
          accept_counter(canon_id) = accept_counter(canon_id) + 1
        else
          ! Rejected
          dof_old%occurrence = dof_old%occurrence + 1
        end if
        
        ! Record MAP
        if (dof%SlogPost .gt. dof_MAP%SlogPost) then
          call set(dof_MAP, dof)
        end if
      end do
      
      ind = ind + 1
      
      ! Tune
      if ((Phase .eq. 2) .and. (mod(iter, output_dn) .eq. 0)) then
        call tune_SSD(canon_SSD, real(accept_counter,kind=dp)/real(output_dn,kind=dp))
        if (iter .lt. niter) then
          accept_counter = 0
        end if
      end if
      
      ! I/O
      if (my_id .eq. (num_procs-1)) then
        if (mod(iter, 5) .eq. 0) then
          write(6, "(i0, a, "//dp_chr//")") ind, "-th step: logPost = ", dof_old%SlogPost
          FLUSH(6)
        end if

        if (mod(iter, output_dn) .eq. 0) then
          call write_theta(dof_old, trim(output_fld)//"theta", sc, ind)
        end if
      end if
      
    end do
      
    ! I/O
    if ((my_id .eq. (num_procs-1))) then
      ! Write MAP
      call write_theta(dof_MAP, trim(output_fld)//"theta", sc, -1)
      ! Write SSD
      call write(canon_SSD, trim(output_fld)//"canon_SSD")
      
      ! Write TXT for easy verification
      call write_txt(real(accept_counter,kind=dp), trim(output_fld)//"accept_counter")
      call write_txt(canon_SSD, trim(output_fld)//"canon_SSD")
    end if
    
    ! Timer
    call stop(total_timer)
    
    ! Release memory
    ! MH
    deallocate(stdnRV)
    deallocate(UniRV)
    call deallocate(prior_param)
    deallocate(accept_counter)
    deallocate(canon_SSD)
    
    call deallocate(dof_old)
    call deallocate(dof_MAP)

    call deallocate(dof)
    call deallocate(IndFn)
    call deallocate(jumps)
    
#if IBACKWARD == 1
    call deallocate(dof_inv)
    call deallocate(jumps_inv)
#endif

#if OMP0MPI1 == 1
    call MPI_FINALIZE(ierr)
#endif

    if ((my_id .eq. (num_procs-1))) then
      write(6, "(a, "//dp_chr//")")  "Final step: logPost = ", dof_old%SlogPost
      write(6, "(a)") "Total timers:"
      call print(total_timer, "Total time on a node", prefix = "  ")
      
      write(6, "(a)") &
      "! ----------- Passed inference ----------- !"
    end if

  end subroutine inference
  
  subroutine pwc_intpl(rfld, fld)
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
  
  subroutine validate_inference(Td, layer)
    use omp_lib
    integer, intent(in) :: Td, layer

    type(jumpsdat), dimension(:), allocatable :: jumps
    type(trajdat) :: traj

    integer :: m, n, m_solver, m_Ind, n_Ind
    type(meshdat) :: mesh

    type(dofdat) :: dof, dof_Davis
    type(IndFndat) :: IndFn

    real(kind=dp) :: sc, h, dt
    integer :: nts
    
    real(kind=dp), parameter :: L = 3840.0_dp*1000.0_dp
    
    character(len=256) :: RunProfile
    character(len=256) :: output_fld, Td_char, resol_param

    ! Timer
    type(timer), save :: total_timer
    
    ! Use Eulerian time-average
    type(field) :: psi_EM
    real(kind=dp) :: t_avg
    integer :: nEMx, nEMy 
    character(len=256) :: KDavis_fld
    ! Use Eulerian time-average
    
    integer, parameter :: restart_ind = -1
    
    integer :: INDk, nts_adapted

    type(dofdat) :: dof_solver  ! Refined dof for solver
    real(kind=dp) :: dt_min, dt_arg
    type(field) :: q
    character(len = 256) ::output_q, suffix
    real(kind=dp), dimension(3) :: Gauss_param    !=(mean_x, mean_y, sigma)
    real(kind=dp), parameter :: Gauss_sigma = 1.0_dp/32.0_dp
    integer :: i, j, nx, ny, DavisID
    
    write(RunProfile, "(a,i0,a)") "QGM2_L", layer, "_NPART676"
    
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
    
    m_solver = 64
    call allocate(mesh, m_solver, m_solver)  ! Needs to be identical in both directions
    ! N.B. assumed zeros at boundaries of psi-> hence only need (m-1)*(n-1) instead of (m+1)*(n+1)
    
    ! Initialise fields
    sc = real(mesh%m,kind=dp)/L     ! Needs mesh to be identical in both directions
    call read_theta(dof, trim(output_fld)//"theta_sigma", sc, restart_ind)
!     call allocate(dof, m-1, n-1, m/2+1, n/2+1)
!     call init_dof(dof, sc)
    
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
    dt = 12.0_dp*3600.0_dp  ! 1.5 hours, in seconds
    nts = int((h/sc)/dt)    ! 2 hours time step
    
    ! Alternative
    call allocate(IndFn, m_Ind, n_Ind)
    
    ! Print info to screen
    call print_info(mesh)
    call print_info(dof)
    call print_info(IndFn)
    call print_info(h, nts, sc)
!     call print_info(jumps, mesh)
    
    ! Timer
    call reset(total_timer)
    call start(total_timer)
    
    call allocate(dof_solver, mesh)
    call intpl_dof_solver(dof_solver, dof, mesh)
    
    dt_arg = h/nts
    dt_min = dt_CFL(dof_solver%psi, 0.1_dp)
    
    nts_adapted = max(int(h/dt_min + 0.5_dp), nts)
    write(6, *) "nts_adapted = ", nts_adapted
    
    ! MAP diffusivity
    ! Solve FK equation
    nx = 6
    ny = 6
    do j = 1, ny
      do i = 1, nx
        call allocate(q, mesh%m, mesh%n, 'q', glayer=1, type_id=2)
        
        INDk = i + (j-1)*ny
        Gauss_param = (/ fld_x(i, nx, 2), fld_x(j, ny, 2), Gauss_sigma /)
        
        call initialise_q_Gauss(q, Gauss_param, mesh)
        write(output_q, "(a,a,a,a,a,a,a)") "./unittest/validation/", trim(resol_param), "/", &
                                          trim(RunProfile), "/", trim(Td_char), "/"
        call write(q, trim(output_q)//"q0", 0.0_dp, INDk)
        
        call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h, 2*nts_adapted)
        
        write(output_q, "(a,a,a,a,a,a,a)") "./unittest/validation/", trim(resol_param), "/", &
                                          trim(RunProfile), "/", trim(Td_char), "/"
        call write(q, trim(output_q)//"q", h, INDk)
        
        call deallocate(q)
      end do
    end do

    ! Write theta_MAP
    call write_theta(dof_solver, trim(output_q)//"theta_MAP", sc)
    
    
    ! Davis diffusivity
    !! Use Eulerian time-average
    call read_QGfield(psi_EM, t_avg, "./meanflow/psi_int_final", layer)

    ! Time-average + rescale wrt mesh
    psi_EM%data = (psi_EM%data - psi_EM%data(1,1))/t_avg
    psi_EM%data = psi_EM%data/(100.0_dp*real(psi_EM%m, kind=dp)) * real(mesh%m, kind=dp)  ! 100: from cm to m
    t_avg =  t_avg/(psi_EM%m/(L*100.0_dp))/(3600.0_dp*24.0_dp*365.25_dp)
    ! call write(psi_EM, trim("./unittest/meanflow/psi_test"), t_avg, 0)
    write(6, *) "Using Eulerian time-averaged mean flow of ", t_avg, " years"
    
    !N.B.: psi_EM%m = 512
    
    ! Filter out
    nEMx = size(psi_EM%data, 1) - 1
    nEMy = size(psi_EM%data, 2) - 1

    dof_solver%psi%data = 0.0_dp
    call sine_filter(dof_solver%psi%data(2:m_solver, 2:m_solver), psi_EM%data(2:nEMx, 2:nEMy) )
    call deallocate(psi_EM)
    
    write(KDavis_fld, "(a,i0,a,i0,a)") "KDavis_L",layer,"_h",Td, "d"

    do DavisID = 1, 2
      if (DavisID .eq. 1) then
        write(suffix, "(a)") "_676"
        call read_theta(dof_Davis, "/home/s1046972/opt/qgm2_particle_diagnosis/&
              &production/DavisDiffusivity_ScarceData/postprocess_output/out/"&
              &//trim(KDavis_fld), sc)
      else
        write(suffix, "(a)") "_40000"
        call read_theta(dof_Davis, "./LocalInf/"//trim(KDavis_fld), sc)
      end if
      
!     call write_theta(dof_Davis, trim("./unittest/meanflow/theta_Davis_test"), sc, 0)
      call pwc_intpl(dof_solver%K11%data(1:m_solver,1:m_solver), dof_Davis%K11%data)
      call pwc_intpl(dof_solver%K22%data(1:m_solver,1:m_solver), dof_Davis%K22%data)
      call pwc_intpl(dof_solver%K12%data(1:m_solver,1:m_solver), dof_Davis%K12%data)
      
      call deallocate(dof_Davis)
      
      ! Ad-hoc way
      dof_solver%K11%data(m_solver+1,:) = dof_solver%K11%data(m_solver,:)
      dof_solver%K22%data(m_solver+1,:) = dof_solver%K22%data(m_solver,:)
      dof_solver%K12%data(m_solver+1,:) = dof_solver%K12%data(m_solver,:)
      dof_solver%K11%data(:,m_solver+1) = dof_solver%K11%data(:,m_solver)
      dof_solver%K22%data(:,m_solver+1) = dof_solver%K22%data(:,m_solver)
      dof_solver%K12%data(:,m_solver+1) = dof_solver%K12%data(:,m_solver)
      
      dt_arg = h/nts
      dt_min = dt_CFL(dof_solver%psi, 0.1_dp)
      
      nts_adapted = max(int(h/dt_min + 0.5_dp), nts)
      write(6, *) "DavisID = ", DavisID, "; nts_adapted = ", nts_adapted, "; dt_min =", dt_min
      
      write(output_q, "(a,a,a,a,a,a,a)") "./unittest/validation/", trim(resol_param), "/", &
                                          trim(RunProfile), "/", trim(Td_char), "/"
      
      do j = 1, ny
        do i = 1, nx
          call allocate(q, mesh%m, mesh%n, 'q', glayer=1, type_id=2)
          
          INDk = i + (j-1)*ny
          Gauss_param = (/ fld_x(i, nx, 2), fld_x(j, ny, 2), Gauss_sigma /)
          
          call initialise_q_Gauss(q, Gauss_param, mesh)
          call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h, nts_adapted)
          
          call write(q, trim(output_q)//"Davis"//trim(suffix), h, INDk)
          
          call deallocate(q)
        end do
      end do
      
      ! Write theta_MAP
      if (DavisID .eq. 1) then
        call write_theta(dof_solver, trim(output_q)//"theta_Davis", sc, 676)
      else
        call write_theta(dof_solver, trim(output_q)//"theta_Davis", sc, 40000)
      end if
      
      
    end do
!     !! End: Use Eulerian time-average
        
    ! Timer
    call stop(total_timer)
    write(6, "(a)") "Total timers:"
    call print(total_timer, "Total time on a node", prefix = "  ")
    
    ! Release memory
    call deallocate(dof_solver)
    call deallocate(dof)
    call deallocate(IndFn)
    call deallocate(jumps)
    
    write(6, "(a)") &
    "! ----------- Passed inference ----------- !"

  end subroutine validate_inference

end program advdiff
