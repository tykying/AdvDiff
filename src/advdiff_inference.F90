module advdiff_inference
  use advdiff_precision
  use advdiff_debug
  use advdiff_timing
  use advdiff_parameters
  use advdiff_field
  use advdiff_timestep
  use advdiff_trajdata
  use advdiff_complib

  implicit none

!   public :: initialise_loglik, eval_cell_loglik, update_logPost
  public :: dofdat, IndFndat
  public :: allocate, deallocate, set, reset_inference_timers, print_inference_timers
  public :: validate_IndFn
  public :: propose_dof
  public :: eval_INDk_loglik

  type dofdat
    integer :: m_psi, n_psi
    integer :: m_K, n_K
    integer :: ndof
    type(field) :: psi, K11, K22, K12
    integer :: occurrence
    real(kind=dp) :: SlogPost
  end type dofdat
  
  type IndFndat
    integer :: m_Ind, n_Ind
    integer :: NInd
  end type IndFndat

  type(timer), save :: loglik_timer, intpl_timer, reffld_timer

  interface allocate
    module procedure allocate_dof, allocate_IndFn, allocate_dof_solver
  end interface allocate

  interface deallocate
    module procedure deallocate_dof, deallocate_IndFn
  end interface deallocate

  interface set
    module procedure set_dof
  end interface set
  
  interface print_info
    module procedure print_info_inference, print_info_rescaled, print_info_IndFn
  end interface print_info

contains
  subroutine allocate_dof(dof, m_psi, n_psi, m_K, n_K)
    type(dofdat), intent(inout) :: dof
    integer, intent(in) :: m_psi, n_psi, m_K, n_K

    integer :: type_id ! 1= corner; 2= centres; x-1 = DOF
    
    dof%m_psi = m_psi
    dof%n_psi = n_psi
    dof%m_K = m_K
    dof%n_K = n_K
    dof%ndof = m_psi*n_psi + m_K*n_K  ! Truth: m_psi*n_psi + 3(m_K*n_K)

    ! Defined at corners + DOF -> -1
    type_id = -1  ! NOT including the boundaries (assumed zeros, hence not DOF)
    call allocate(dof%psi, m_psi, n_psi, 'psi', glayer=0, type_id=type_id)

    ! Defined at cell centre + DOF -> -2
    type_id = -2
    call allocate(dof%K11, m_K, n_K, 'K11', glayer=0, type_id=type_id)
    call allocate(dof%K22, m_K, n_K, 'K22', glayer=0, type_id=type_id)
    call allocate(dof%K12, m_K, n_K, 'K12', glayer=0, type_id=type_id)

    dof%occurrence = 0
    dof%SlogPost = 0.0_dp
    
  end subroutine allocate_dof

  subroutine deallocate_dof(dof)
    type(dofdat), intent(inout) :: dof

    call deallocate(dof%psi)
    call deallocate(dof%K11)
    call deallocate(dof%K22)
    call deallocate(dof%K12)

  end subroutine deallocate_dof
  
  subroutine allocate_IndFn(InfFn, m_Ind, n_Ind)
    type(IndFndat), intent(inout) :: InfFn
    integer, intent(in) :: m_Ind, n_Ind

    ! Assume no overlapping
    InfFn%m_Ind = m_Ind
    InfFn%n_Ind = n_Ind
    InfFn%NInd = InfFn%m_Ind * InfFn%n_Ind
    
  end subroutine allocate_IndFn
  
  subroutine deallocate_IndFn(InfFn)
    type(IndFndat), intent(inout) :: InfFn

  end subroutine deallocate_IndFn
  
  subroutine validate_IndFn(InfFn, mesh)
    type(IndFndat), intent(in) :: InfFn
    type(meshdat), intent(in) :: mesh
   
    if ((mod(mesh%m, InfFn%m_ind) .ne. 0) .or. (mod(mesh%n, InfFn%n_ind) .ne. 0)) &
      call abort_handle("dim(mesh) inconsistent with indicator functions m_ind, n_ind", __FILE__, __LINE__)
    
  end subroutine validate_IndFn

  subroutine set_dof(dof, dof_in)
    type(dofdat), intent(out) :: dof
    type(dofdat), intent(in) :: dof_in

    call set(dof%psi, dof_in%psi)
    call set(dof%K11, dof_in%K11)
    call set(dof%K22, dof_in%K22)
    call set(dof%K12, dof_in%K12)

  end subroutine set_dof
  
  subroutine evaluate_loglik_OMP(loglik, jumps, IndFn, mesh, dof, h, nts)
    real(kind=dp), dimension(:), intent(out) :: loglik
    type(jumpsdat), dimension(:), pointer, intent(in) :: jumps
    type(IndFndat), intent(in) :: IndFn
    type(meshdat), intent(in) :: mesh
    type(dofdat), intent(in) :: dof
    real(kind=dp), intent(in) :: h
    integer, intent(in) :: nts

    integer :: INDk

    type(dofdat) :: dof_solver  ! Refined dof for solver
    
    call allocate(dof_solver, mesh)
    call intpl_dof_solver(dof_solver, dof, mesh)

    !! Control number of thread = 32 by 
    !! "!$OMP PARALLEL DO num_threads(32)"
    
    !$OMP PARALLEL DO num_threads(4)
    do INDk = 1, IndFn%nIND
      loglik(INDk) = eval_INDk_loglik(INDk, jumps, IndFn, mesh, dof_solver, h, nts)
    end do
    !$OMP END PARALLEL DO
    
    call deallocate(dof_solver)
    
  end subroutine evaluate_loglik_OMP
  
  subroutine INDk_to_klist(klist, INDk, IndFn, mesh)
    integer, dimension(:), allocatable, intent(out) :: klist
    integer, intent(in) :: INDk
    type(IndFndat), intent(in) :: IndFn
    type(meshdat), intent(in) :: mesh
    
    integer :: INDi, INDj, di, dj, i0, j0, i1, j1
    integer :: klist_len, i, j
    
    ! Map k to (i,j)
    INDi = mod(INDk-1, IndFn%m_IND) + 1
    INDj = (INDk-1)/ IndFn%m_IND + 1
    
    ! Determine number of computational grid required
    di = mesh%m/IndFn%m_IND
    dj = mesh%n/IndFn%n_IND
    
    i0 = (INDi-1)*di+1
    i1 = INDi*di
    j0 = (INDj-1)*dj+1
    j1 = INDj*dj
    
    klist_len = (i1-i0+1)*(j1-j0+1)
    
    allocate(klist(klist_len))
    
    klist_len = 0
    do j = j0, j1
      do i = i0, i1
        klist_len = klist_len + 1
        klist(klist_len) = ij2k(i,j,mesh%m)
      end do
    end do
    
  end subroutine INDk_to_klist
  
  subroutine initialise_q(q, klist, mesh)
    type(field), intent(inout) :: q
    integer, dimension(:), intent(in) :: klist
    type(meshdat), intent(in) :: mesh

    integer :: k_i, i, nzc
    
    nzc = size(klist, 1)  ! Number of non-zero cells

    call zeros(q)
    do i = 1, nzc
      k_i = klist(i)
      call indicator_field(q, k2i(k_i, mesh%m), k2j(k_i, mesh%m))
    end do
    call scale(q, 1.0_dp/real(nzc, kind=dp))
  
  end subroutine initialise_q
  
  real(kind=dp) function eval_INDk_loglik(INDk, jumps, IndFn, mesh, dof_solver, h, nts)
    integer, intent(in) :: INDk
    type(jumpsdat), dimension(:), pointer, intent(in) :: jumps
    type(IndFndat), intent(in) :: IndFn
    type(meshdat), intent(in) :: mesh
    type(dofdat), intent(in) :: dof_solver
    real(kind=dp), intent(in) :: h
    integer, intent(in) :: nts
    integer, dimension(:), allocatable :: klist

    type(field) :: q
    integer :: jump
    real(kind=dp) :: lik
    integer :: k_i, i

    call start(loglik_timer)

    eval_INDk_loglik = 0.0_dp

    ! Find out the cells corresponding to the indicator function
    call INDk_to_klist(klist, INDk, IndFn, mesh)
    call allocate(q, mesh%m, mesh%n, 'q', glayer=1, type_id=2)
    
    ! Solve FK equation
    call initialise_q(q, klist, mesh)

!     call print_array(q%data, "q0%data")
    call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h, nts)
!       if (.not. is_nneg(q%data)) then
!         write(6, "(a,"//dp_chr//")") "Non-negative q!: ", minval(q%data)
!       end if

!     call print_array(dof_solver%psi%data, "dof_solver%psi")
!     call print_array(dof_solver%K11%data, "dof_solver%K11")
!     call print_array(dof_solver%K22%data, "dof_solver%K22")
!     call print_array(dof_solver%K12%data, "dof_solver%K12")
!     call print_array(q%data, "q%data")

    ! Evaluate likelihood
    call start(intpl_timer)
    do i = 1, size(klist, 1)
      k_i = klist(i)
      do jump = 1, jumps(k_i)%njumps
        lik = eval_fldpt(q, mesh, jumps(k_i)%k_f(jump), jumps(k_i)%alpha_f(:,jump))
        ! Set the negative probability to be 1e-16
        lik = max(1e-16, lik)
        eval_INDk_loglik = eval_INDk_loglik + dlog(lik)
      end do
    end do
    call stop(intpl_timer)

    call deallocate(q)

    deallocate(klist)
    
    call stop(loglik_timer)

  end function eval_INDk_loglik
  
  subroutine allocate_dof_solver(dof_solver, mesh)
    type(dofdat), intent(out) :: dof_solver
    type(meshdat), intent(in) :: mesh

    integer :: type_id ! 1= corner; 2= centres; *-1 = DOF
    
    dof_solver%m_psi = mesh%m
    dof_solver%n_psi = mesh%n
    dof_solver%m_K = mesh%m
    dof_solver%n_K = mesh%n
    dof_solver%ndof = 0

    ! Defined at corners
    type_id = 1
    call allocate(dof_solver%psi, mesh%m, mesh%n, 'psi', glayer=0, type_id=type_id)

    ! Defined at cell centre
    type_id = 2
    call allocate(dof_solver%K11, mesh%m, mesh%n, 'K11', glayer=0, type_id=type_id)
    call allocate(dof_solver%K22, mesh%m, mesh%n, 'K22', glayer=0, type_id=type_id)
    call allocate(dof_solver%K12, mesh%m, mesh%n, 'K12', glayer=0, type_id=type_id)

    dof_solver%occurrence = 0
    dof_solver%SlogPost = 0.0_dp
    
  end subroutine allocate_dof_solver
  
  subroutine intpl_dof_solver(dof_solver, dof, mesh)
    ! Interpolate from DOF to fields for the solver
    type(dofdat), intent(inout) :: dof_solver  ! Field to use for solver
    type(dofdat), intent(in) :: dof
    type(meshdat), intent(in) :: mesh
    
    integer :: m, n, Mx, My
    real(kind=dp), dimension(:, :), allocatable :: fldij
    
    ! fldij(i, j) = field(x_i, y_j)
    ! DOF: can be coarse field(x_i, y_j) or Fourier coefficients
    ! dof_to_fldij: an interface from DOF to fld(i,j)
    
    ! Streamfunction: psi:
    ! N.B. DOF%psi = non-bondary values, assumed boundary value = 0
    m = dof%psi%m
    n = dof%psi%n

    allocate(fldij(m+1, n+1))
    fldij = 0.0_dp  ! Boundary condition for stream function
    fldij(2:(m+1), 2:(n+1)) = dof%psi%data
    
    Mx = mesh%m/(m+1)
    My = mesh%n/(n+1)
    ! Assumed periodic(= 0) at first and final row/column
    call fourier_intpl(dof_solver%psi%data(1:mesh%m, 1:mesh%n), fldij, Mx, My)
    dof_solver%psi%data(mesh%m+1, :) = dof_solver%psi%data(1, :)  ! Impose periodicity
    dof_solver%psi%data(:, mesh%n+1) = dof_solver%psi%data(:, 1)  ! Impose periodicity

    deallocate(fldij)
    
    ! Diffusivity
    ! N.B. DOF%psi = Cartesian grid at cell centre
    Mx = mesh%m/dof%K11%m
    My = mesh%n/dof%K11%n
    
    call cosine_intpl(dof_solver%K11%data, dof%K11%data, Mx, My)
    call cosine_intpl(dof_solver%K22%data, dof%K22%data, Mx, My)
    call cosine_intpl(dof_solver%K12%data, dof%K12%data, Mx, My)
    
  end subroutine intpl_dof_solver

  pure logical function is_nneg(data)
    real(kind=dp), dimension(:, :), intent(in) :: data

    integer :: i, j

    is_nneg = .true.
    do j = 1, size(data, 2)
      do i = 1, size(data, 1)
        ! Not stictly negative but with a threshold -1e-14
        is_nneg = (is_nneg .and. (data(i,j) .ge. -1e-10))
      end do
    end do

  end function is_nneg

  subroutine propose_dof(dof, dof_id, dof_SSD)
    type(dofdat), intent(inout) :: dof
    integer, intent(in) :: dof_id
    type(dofdat), intent(in) :: dof_SSD
    
    integer :: K_id, psi_id

    call abort_handle("TODO: Non overlapping DOF in K and psi", __FILE__, __LINE__)
    
    if (dof_id .le. dof%m_psi*dof%n_psi) then
      K_id = dof_id
      psi_id = dof%ndof - dof_id + 1
    else
      K_id = 1
      psi_id = dof_id
    end if
    
    call propose_dof_K(dof, K_id, dof_SSD)
    call propose_dof_psi(dof, psi_id, dof_SSD)
    
  end subroutine propose_dof
  
  subroutine propose_dof_all(dof, iter, dof_SSD)
    type(dofdat), intent(inout) :: dof
    integer, intent(in) :: iter
    type(dofdat), intent(in) :: dof_SSD

    integer :: dof_id
    
    !write(6, *) "propose_dof: Proposing two changes at a time"
    if (mod(iter, 2) .eq. 0) then
      do dof_id = 1, dof%m_K*dof%n_K
        call propose_dof_K(dof, dof_id, dof_SSD)
      end do
    else
      do dof_id = 1, dof%m_psi*dof%n_psi
        call propose_dof_psi(dof, dof_id, dof_SSD)
      end do
    end if
  
  end subroutine propose_dof_all

  subroutine propose_dof_K(dof, dof_id, dof_SSD)
    type(dofdat), intent(inout) :: dof
    integer, intent(in) :: dof_id
    type(dofdat), intent(in) :: dof_SSD

    integer :: i0, j0, i, j
    real(kind = dp) :: sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy

    i0 = k2i(dof_id, dof%m_K)
    j0 = k2j(dof_id, dof%m_K)

    if ( (dof%K11%type_id .eq. -2) .and. &
      (dof%K22%type_id .eq. -2) .and. (dof%K12%type_id .eq. -2) ) then
      ! K: defined at cell-centres
      Kxx = dof%K11%data(i0,j0)
      Kyy = dof%K22%data(i0,j0)
      Kxy = dof%K12%data(i0,j0)

      call KCarte_to_KCanon(sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy)

      ! Perturbation
      sigma1_sq = dabs(sigma1_sq + dof_SSD%K11%data(i0,j0)*randn() )
      sigma2_sq = dabs(sigma2_sq + dof_SSD%K22%data(i0,j0)*randn() )
      phi_K = phi_K + dof_SSD%K12%data(i0,j0)*randn()

      !! TODO: delete Adhoc
      !phi_K = 0.0_dp
      !! TODO: delete Adhoc

      call KCanon_to_KCarte(Kxx, Kyy, Kxy, sigma1_sq, sigma2_sq, phi_K)

      dof%K11%data(i0,j0) = Kxx
      dof%K22%data(i0,j0) = Kyy
      dof%K12%data(i0,j0) = Kxy
    else
      call abort_handle("Invalid type_id for K", __FILE__, __LINE__)
    end if

  end subroutine propose_dof_K

  subroutine propose_dof_psi(dof, dof_id, dof_SSD)
    type(dofdat), intent(inout) :: dof
    integer, intent(in) :: dof_id
    type(dofdat), intent(in) :: dof_SSD

    integer :: i0, j0

    i0 = k2i(dof_id, dof%m_psi)
    j0 = k2j(dof_id, dof%m_psi)
    
    ! Fourier basis
    if (dof%psi%type_id .ne. -1) &
      call abort_handle("Invalid type_id for psi", __FILE__, __LINE__)
    
    dof%psi%data(i0,j0) = dof%psi%data(i0,j0) + dof_SSD%psi%data(i0,j0)*randn()
    
  end subroutine propose_dof_psi

  subroutine reset_inference_timers()
    call reset(loglik_timer)
    call reset(intpl_timer)
    call reset(reffld_timer)
  end subroutine reset_inference_timers
  
  subroutine print_inference_timers()
    write(6, "(a)") "Inference timers:"
    call print(loglik_timer, "Evaluate loglik", prefix = "  ")
    call print(intpl_timer, "  -> interpolations", prefix = "  ")
    call print(reffld_timer, "  -> refine fields", prefix = "  ")
  end subroutine print_inference_timers
  
  subroutine print_info_IndFn(IndFn)
    type(IndFndat), intent(in) :: IndFn
    
    write(6, "(a)") ""
    write(6, "(a)") " Indicator function info: "
    write(6, "(a,"//int_chr//",a,"//int_chr//")") &
      "|  IndFn%m_ind, n_ind", IndFn%m_ind, " , ", IndFn%n_ind
    write(6, "(a,"//int_chr//")") "|  IndFn%nIND", IndFn%nIND
 
   end subroutine print_info_IndFn
  
  subroutine print_info_inference(dof)
    type(dofdat), intent(in) :: dof

    write(6, "(a)") ""
    write(6, "(a)") " Inference info: "
    write(6, "(a,"//int_chr//",a,"//int_chr//")") & 
      "|  dof%m_psi, n_psi", dof%m_psi, " , ", dof%n_psi
    write(6, "(a,"//int_chr//",a,"//int_chr//")") & 
      "|  dof%m_K, n_K", dof%m_K, " , ", dof%n_K
    write(6, "(a,"//int_chr//")") "|  dof%ndof", dof%ndof
    write(6, "(a)") "| [type_id =2 -> centre; =1 -> corner]"
    write(6, "(a,"//int_chr//")") "|  psi%type_id", dof%psi%type_id
    write(6, "(a,"//int_chr//")") "|  K11%type_id", dof%K11%type_id
    write(6, "(a,"//int_chr//")") "|  K22%type_id", dof%K22%type_id
    write(6, "(a,"//int_chr//",a,"//int_chr//")") "| dof%occurrence", dof%occurrence
    write(6, "(a,"//sp_chr//")") "|  dof%SlogPost", dof%SlogPost
    
    FLUSH(6)
  end subroutine print_info_inference
  
  subroutine print_info_rescaled(h, nts, sc)
    real(kind = dp), intent(in) :: h, sc
    integer, intent(in) :: nts

    write(6, "(a)") ""
    write(6, "(a)") " Physical quantities info: "
    write(6, "(a,"//sp_chr//", a)") "|  h = ", (h/sc)/(24.0_dp*3600.0_dp), " days"
    write(6, "(a,"//int_chr//")") "|  nts = ", nts
    write(6, "(a,"//sp_chr//", a)") "|  dt = ", ((h/sc)/real(nts,kind=dp))/(3600.0_dp), " hours"

    FLUSH(6)
  end subroutine print_info_rescaled

end module advdiff_inference
