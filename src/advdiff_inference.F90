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
  public :: eval_INDk_loglik
  public :: reverse_dof_to_dof_inv
  public :: compute_dJdm_Cid
  public :: convert_dof_to_theta, convert_theta_to_dof
  public :: convert_dof_to_canon, convert_canon_to_dof
  public :: propose_dof_canon, propose_dof_EM

  type dofdat
    integer :: m_psi, n_psi
    integer :: m_K, n_K
    integer :: ndof
    type(field) :: psi, K11, K22, K12
    integer :: occurrence
    real(kind=dp) :: SlogPost
    real(kind=dp), dimension(:), pointer :: canon
  end type dofdat
  
  type IndFndat
    integer :: m_Ind, n_Ind
    integer :: NInd
  end type IndFndat

  type(timer), save :: loglik_timer, intpl_timer, reffld_timer

  interface allocate
    module procedure allocate_dof, allocate_IndFn, allocate_dof_solver, allocate_theta
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
    dof%ndof = m_psi*n_psi + 3*(m_K*n_K)

    ! Defined at corners + DOF -> -1
    type_id = -1  ! NOT including the boundaries (assumed zeros, hence not DOF)
    call allocate(dof%psi, m_psi, n_psi, 'psi', glayer=0, type_id=type_id)

    ! Defined at cell centre + DOF -> -2
    type_id = -1  ! Including the boundaries
    call allocate(dof%K11, m_K, n_K, 'K11', glayer=0, type_id=type_id)
    call allocate(dof%K22, m_K, n_K, 'K22', glayer=0, type_id=type_id)
    call allocate(dof%K12, m_K, n_K, 'K12', glayer=0, type_id=type_id)
    
    ! By default zeros
    call zeros(dof%psi)
    call zeros(dof%K11)
    call zeros(dof%K22)
    call zeros(dof%K12)

    dof%occurrence = 0
    dof%SlogPost = 0.0_dp
   
    ! canonical
    allocate(dof%canon(dof%ndof))

  end subroutine allocate_dof

  subroutine deallocate_dof(dof)
    type(dofdat), intent(inout) :: dof

    call deallocate(dof%psi)
    call deallocate(dof%K11)
    call deallocate(dof%K22)
    call deallocate(dof%K12)
    
    deallocate(dof%canon)

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
  
  subroutine init_dof(dof, sc)
    type(dofdat), intent(inout) :: dof
    real(kind=dp), optional, intent(in) :: sc
    
    real(kind=dp) :: ssc
    real(kind=dp), parameter :: L = 3840.0_dp*1000.0_dp
    real(kind=dp), parameter :: kappa_scale = 10000.0_dp
!    real(kind=dp), parameter :: kappa_scale = 0.2_dp*10000.0_dp   ! Layer 2
    real(kind=dp), parameter :: psi_scale = 100.0_dp*1000.0_dp  ! Irrelevant
    
    if (present(sc)) then
      ssc = sc
    else
      ssc = 1.0_dp
    end if
    
    ! Initialise dof
    call zeros(dof%psi)
    call set(dof%K11, 5.0_dp*kappa_scale*ssc)
    call set(dof%K22, 5.0_dp*kappa_scale*ssc)
    call set(dof%K12, 0.0_dp*kappa_scale*ssc)
    
    call imprint_canon(dof)
    
  end subroutine init_dof

  subroutine set_dof(dof, dof_in)
    type(dofdat), intent(out) :: dof
    type(dofdat), intent(in) :: dof_in

    call set(dof%psi, dof_in%psi)
    call set(dof%K11, dof_in%K11)
    call set(dof%K22, dof_in%K22)
    call set(dof%K12, dof_in%K12)

    dof%SlogPost = dof_in%SlogPost
    dof%occurrence = dof_in%occurrence
    dof%ndof = dof_in%ndof

    dof%canon = dof_in%canon
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
    
!     call print_info(dof_solver%psi, h/nts)

    !! Control number of thread = 32 by 
    !! "!$OMP PARALLEL DO num_threads(32)"
    
    !$OMP PARALLEL DO
    do INDk = 1, IndFn%nIND
      loglik(INDk) = eval_INDk_loglik(INDk, jumps, IndFn, mesh, dof_solver, h, nts)
    end do
    !$OMP END PARALLEL DO
    
    call deallocate(dof_solver)
    
  end subroutine evaluate_loglik_OMP
  
  subroutine obtain_logLik_guide(IndFn_List, canon_id, criteria)
    logical, dimension(:), intent(out) :: IndFn_List
    integer, intent(in) :: canon_id
    real(kind=dp), dimension(:, :), intent(in) :: criteria

    real(kind=dp) :: threshold
    integer :: INDk, IndFn_List_len, i
    
    real(kind=dp) :: TE  ! total energy
    
    IndFn_List = .FALSE.
    
    ! Select by threshold
    TE = sum(criteria(:, canon_id))
    threshold = 1E-8 * TE/real(size(criteria, 1), kind=dp)
!    threshold = 1E-10
    
    IndFn_List_len = 0
    do INDk = 1, size(criteria, 1)
      if (criteria(INDk, canon_id) .ge. threshold) then
        IndFn_List(INDk) = .TRUE.
      end if
    end do
    
  end subroutine obtain_logLik_guide
  
  subroutine evaluate_loglik_guided(loglik, jumps, IndFn, mesh, dof, h, nts, canon_id, dJdm_sqsum)
    real(kind=dp), dimension(:), intent(inout) :: loglik
    type(jumpsdat), dimension(:), pointer, intent(in) :: jumps
    type(IndFndat), intent(in) :: IndFn
    type(meshdat), intent(in) :: mesh
    type(dofdat), intent(in) :: dof
    real(kind=dp), intent(in) :: h
    integer, intent(in) :: nts
    integer, intent(in) :: canon_id
    real(kind=dp), dimension(:, :), intent(in) :: dJdm_sqsum
    
    integer :: INDk, i
    logical, dimension(:), allocatable :: IndFn_List

    type(dofdat) :: dof_solver  ! Refined dof for solver
    
    call allocate(dof_solver, mesh)
    call intpl_dof_solver(dof_solver, dof, mesh)
    
    allocate(IndFn_List(IndFn%nIND))
    
    call obtain_logLik_guide(IndFn_List, canon_id, dJdm_sqsum)
    
!     if (mod(canon_id, 40) .eq. 2) then
!       write(6, *) "canon_id = ", canon_id
!       write(6, *) IndFn_List
!     end if
    
    !$OMP PARALLEL DO
    do INDk = 1, IndFn%nIND
      if (IndFn_List(INDk)) then
        loglik(INDk) = eval_INDk_loglik(INDk, jumps, IndFn, mesh, dof_solver, h, nts)
      end if
    end do
    !$OMP END PARALLEL DO
    
    call deallocate(dof_solver)
    deallocate(IndFn_List)
    
  end subroutine evaluate_loglik_guided
  
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
    call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h, nts)
!       if (.not. is_nneg(q%data)) then
!         write(6, "(a,"//dp_chr//")") "Non-negative q!: ", minval(q%data)
!       end if
    
    ! Evaluate likelihood
    call start(intpl_timer)
    do i = 1, size(klist, 1)
      k_i = klist(i)
      do jump = 1, jumps(k_i)%njumps
        ! Bilinear interpolations
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
    
    dof_solver%m_psi = mesh%m
    dof_solver%n_psi = mesh%n
    dof_solver%m_K = mesh%m
    dof_solver%n_K = mesh%n
    dof_solver%ndof = 0

    ! Defined at corners (! 1= corner; 2= centres; *-1 = DOF)
    call allocate(dof_solver%psi, mesh%m, mesh%n, 'psi', glayer=0, type_id=1)

    call allocate(dof_solver%K11, mesh%m, mesh%n, 'K11', glayer=0, type_id=1)
    call allocate(dof_solver%K22, mesh%m, mesh%n, 'K22', glayer=0, type_id=1)
    call allocate(dof_solver%K12, mesh%m, mesh%n, 'K12', glayer=0, type_id=1)

    dof_solver%occurrence = 0
    dof_solver%SlogPost = 0.0_dp
    
    allocate(dof_solver%canon(1))
    dof_solver%canon(1) = 0.0_dp
    
  end subroutine allocate_dof_solver
  
  subroutine intpl_dof_solver(dof_solver, dof, mesh)
    ! Interpolate from DOF to fields for the solver
    type(dofdat), intent(inout) :: dof_solver  ! Field to use for solver
    type(dofdat), intent(in) :: dof
    type(meshdat), intent(in) :: mesh
    
    integer :: m, n, Mx, My
    real(kind=dp), dimension(:, :), allocatable :: fldij ! Bilinear intepolation
    
    ! Unstructured grid
    integer, dimension(:), allocatable :: xind, yind
    integer :: i, j, k
    ! END Unstructured grid
    
    ! fldij(i, j) = field(x_i, y_j)
    ! DOF: can be coarse field(x_i, y_j) or Fourier coefficients
    ! dof_to_fldij: an interface from DOF to fld(i,j)
    
    ! Streamfunction: psi
    ! N.B. DOF%psi = non-bondary values, assumed boundary value = 0
    m = dof%psi%m
    n = dof%psi%n

!     ! Fourier intepolation
!     allocate(fldij(m+1, n+1))
!     fldij = 0.0_dp  ! Boundary condition for stream function
!     fldij(2:(m+1), 2:(n+1)) = dof%psi%data
!     
!     Mx = mesh%m/(m+1)
!     My = mesh%n/(n+1)
!     ! Assumed periodic(= 0) at first and final row/column
!     call fourier_intpl(dof_solver%psi%data(1:mesh%m, 1:mesh%n), fldij, Mx, My)
!     dof_solver%psi%data(mesh%m+1, :) = dof_solver%psi%data(1, :)  ! Impose periodicity
!     dof_solver%psi%data(:, mesh%n+1) = dof_solver%psi%data(:, 1)  ! Impose periodicity
! 
!     deallocate(fldij)

! !     ! Bilinear intepolation
!     allocate(fldij(m+2, n+2))
!     fldij = 0.0_dp  ! Boundary condition for stream function
!     fldij(2:(m+1), 2:(n+1)) = dof%psi%data
!    
!     Mx = mesh%m/(m+1)
!     My = mesh%n/(n+1)
!     ! Assumed periodic(= 0) at first and final row/column
!     call bilinear_intpl(dof_solver%psi%data, fldij, Mx, My)
!     deallocate(fldij)

!     ! Sine intepolation
!     Mx = mesh%m/(m+1)
!     My = mesh%n/(n+1)
!     dof_solver%psi%data = 0.0_dp
!     call sine_intpl(dof_solver%psi%data(2:mesh%m, 2:mesh%n), dof%psi%data, Mx, My)

!     ! Unstructured bilinear intepolation: dof_solver%psi%data
!     if ((m .eq. 15) .and. (size(dof_solver%psi%data,1) .eq. 65)) then
!       allocate(fldij(m+2, n+2))
!       fldij = 0.0_dp  ! Boundary condition for streamfunction
!       fldij(2:(m+1), 2:(n+1)) = dof%psi%data
! 
!       dof_solver%psi%data = 0.0_dp
!       allocate(xind(17))
!       xind = (/0, 1, 2, 4, 6, 9, 12, 16, 20, 25, 30, 36, 42, 49, 55, 60, 64 /) ! 65: end
!       xind = xind + 1
!       allocate(yind(17))
! !       yind = (/0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64 /) ! 65: end
!       yind = (/0, 4, 10, 18, 24, 29, 33, 35, 37, 38, 39, 41, 43, 47, 52, 60, 64 /)
!       yind = yind + 1
!       
!       ! Interpolate along x
!       do j = 1, size(yind, 1)
!         do i = 1, (size(xind, 1)-1)
!           do k = xind(i), (xind(i+1)-1)
!             dof_solver%psi%data(k, yind(j)) = &
!               ( fldij(i, j)  *(xind(i+1)-k) + &
!                 fldij(i+1, j)*(k  -xind(i)) ) / &
!               (xind(i+1)-xind(i))
!           end do
!         end do
!         
!         ! Final point: size(xind, 1)
!         dof_solver%psi%data(xind(size(xind,1)), yind(j)) = 0.0_dp
!       end do
!       
!       ! Interpolate along y
!       do i = 1, size(dof_solver%psi%data, 2)
!         do j = 1, (size(yind, 1)-1)
!           do k = yind(j), (yind(j+1)-1)
!             dof_solver%psi%data(i, k) = &
!               ( dof_solver%psi%data(i, yind(j))  *(yind(j+1)-k) + &
!                 dof_solver%psi%data(i, yind(j+1))*(k  -yind(j)) ) / &
!               (yind(j+1)-yind(j))
!           end do
!         end do
!         
!         ! Final point: size(yind, 2)
!         dof_solver%psi%data(i, yind(size(yind,1))) = 0.0_dp
!       end do
!       
!       deallocate(fldij)
!       deallocate(xind)
!       deallocate(yind)
!     end if
!     ! END Unstructured grid
    
    !! Use Eulerian time-average
    dof_solver%psi%data = 0.0_dp
    dof_solver%psi%data(2:size(dof_solver%psi%data,1)-1, &
                        2:size(dof_solver%psi%data,2)-1) = dof%psi%data
    !! End: Use Eulerian time-average

    ! Diffusivity
    ! N.B. DOF%psi = Cartesian grid at cell centre
    Mx = mesh%m/(dof%K11%m-1)
    My = mesh%n/(dof%K11%n-1)
    
    call bilinear_intpl(dof_solver%K11%data, dof%K11%data, Mx, My)
    call bilinear_intpl(dof_solver%K22%data, dof%K22%data, Mx, My)
    call bilinear_intpl(dof_solver%K12%data, dof%K12%data, Mx, My)
    
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
  
  subroutine init_canon_SSD(canon_SSD, dof, sc)
    real(kind=dp), dimension(:), intent(inout) :: canon_SSD
    type(dofdat), intent(in) :: dof
    real(kind=dp), optional, intent(in) :: sc
    
    real(kind=dp) :: ssc
    real(kind=dp), parameter :: L = 3840.0_dp*1000.0_dp
    real(kind=dp), parameter :: kappa_scale = 10000.0_dp
    real(kind=dp), parameter :: psi_scale = 100.0_dp*1000.0_dp  ! Layer 1
!    real(kind=dp), parameter :: psi_scale = 10.0_dp*1000.0_dp  ! Layer 2
    
    integer :: k, list_i, list_f
    
    if (present(sc)) then
      ssc = sc
    else
      ssc = 1.0_dp
    end if
    
    ! psi
    list_i = 1
    list_f = dof%m_psi*dof%n_psi
    do k = list_i, list_f
      canon_SSD(k) = 0.1_dp*psi_scale*ssc
    end do
    
    ! sigma1_sq, sigma2_sq 
    list_i = dof%m_psi*dof%n_psi + 1
    list_f = dof%m_psi*dof%n_psi + 2*dof%m_K*dof%n_K
    do k = list_i, list_f
      canon_SSD(k) = 0.1_dp*kappa_scale*ssc
    end do
    
    ! phi
    list_i = dof%m_psi*dof%n_psi + 2*dof%m_K*dof%n_K + 1
    list_f = dof%m_psi*dof%n_psi + 3*dof%m_K*dof%n_K
    do k = list_i, list_f
      canon_SSD(k) = 0.1_dp
    end do

  end subroutine init_canon_SSD
  
  subroutine imprint_canon(dof, canon_id)
    type(dofdat), intent(inout) :: dof
    integer, optional, intent(in) :: canon_id
    
    integer :: loop_i, loop_f, cid, i0, j0
    integer :: psi_id, K_id, Kloc_id
    integer :: ncid, Kcmp1, Kcmp2, Kcmp3
    
    real(kind=dp) :: sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy
    
    ncid = dof%m_psi*dof%n_psi + dof%m_K*dof%n_K
    
    if (present(canon_id)) then
      loop_i = canon_id
      loop_f = canon_id
    else
      loop_i = 1
      loop_f = ncid
    end if
    
    do cid = loop_i, loop_f
      if (cid .le. dof%m_psi*dof%n_psi) then
        ! psi
        psi_id = cid
      
        i0 = k2i(psi_id, dof%m_psi)
        j0 = k2j(psi_id, dof%m_psi)
    
        dof%canon(cid) = dof%psi%data(i0,j0)
      else
        ! K
        K_id = cid - dof%m_psi*dof%n_psi
        
!         cmp = (K_id-1)/(dof%m_K*dof%n_K)
        Kloc_id = mod(K_id-1, dof%m_K*dof%n_K)+1
        
        i0 = k2i(Kloc_id, dof%m_K)
        j0 = k2j(Kloc_id, dof%m_K)
      
        ! K: defined at corners
        Kxx = dof%K11%data(i0,j0)
        Kyy = dof%K22%data(i0,j0)
        Kxy = dof%K12%data(i0,j0)
      
        call KCarte_to_KCanon(sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy)
      
        ! Update all 
        Kcmp1 = dof%m_psi*dof%n_psi
        Kcmp2 = dof%m_psi*dof%n_psi + dof%m_K*dof%n_K
        Kcmp3 = dof%m_psi*dof%n_psi + 2*dof%m_K*dof%n_K

        dof%canon(Kcmp1 + Kloc_id) = sigma1_sq
        dof%canon(Kcmp2 + Kloc_id) = sigma2_sq
        dof%canon(Kcmp3 + Kloc_id) = phi_K
      end if
      
!       write(6, *) "cid, i0, j0 = ", cid, i0, j0
    end do
    
  end subroutine imprint_canon
  
  subroutine convert_dof_to_canon(canon, dof, canon_id)
    real(kind = dp), dimension(:), intent(inout) :: canon
    type(dofdat), intent(in) :: dof
    integer, optional, intent(in) :: canon_id
    
    integer :: loop_i, loop_f, cid, i0, j0
    integer :: psi_id, K_id, Kloc_id
    integer :: ncid, Kcmp1, Kcmp2, Kcmp3
    
    real(kind=dp) :: sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy
    
    ncid = dof%m_psi*dof%n_psi + dof%m_K*dof%n_K
    
    if (present(canon_id)) then
      loop_i = canon_id
      loop_f = canon_id
    else
      loop_i = 1
      loop_f = ncid
    end if
    
    do cid = loop_i, loop_f
      if (cid .le. dof%m_psi*dof%n_psi) then
        ! psi
        psi_id = cid
      
        i0 = k2i(psi_id, dof%m_psi)
        j0 = k2j(psi_id, dof%m_psi)
    
        canon(cid) = dof%psi%data(i0,j0)
      else
        ! K
        K_id = cid - dof%m_psi*dof%n_psi
        
!         cmp = (K_id-1)/(dof%m_K*dof%n_K)
        Kloc_id = mod(K_id-1, dof%m_K*dof%n_K)+1
        
        i0 = k2i(Kloc_id, dof%m_K)
        j0 = k2j(Kloc_id, dof%m_K)
      
        ! K: defined at corners
        Kxx = dof%K11%data(i0,j0)
        Kyy = dof%K22%data(i0,j0)
        Kxy = dof%K12%data(i0,j0)
      
        call KCarte_to_KCanon(sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy)
      
        ! Update all 
        Kcmp1 = dof%m_psi*dof%n_psi
        Kcmp2 = dof%m_psi*dof%n_psi + dof%m_K*dof%n_K
        Kcmp3 = dof%m_psi*dof%n_psi + 2*dof%m_K*dof%n_K

        canon(Kcmp1 + Kloc_id) = sigma1_sq
        canon(Kcmp2 + Kloc_id) = sigma2_sq
        canon(Kcmp3 + Kloc_id) = phi_K
      end if
      
!       write(6, *) "cid, i0, j0 = ", cid, i0, j0
    end do
    
  end subroutine convert_dof_to_canon
  
  subroutine convert_canon_to_dof(dof, canon, canon_id)
    type(dofdat), intent(inout) :: dof
    real(kind = dp), dimension(:), intent(in) :: canon
    integer, optional, intent(in) :: canon_id
    
    integer :: loop_i, loop_f, cid, i0, j0
    integer :: psi_id, K_id, Kloc_id
    integer :: ncid, Kcmp1, Kcmp2, Kcmp3
    
    real(kind=dp) :: sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy
    
    ncid = dof%m_psi*dof%n_psi + dof%m_K*dof%n_K
    
    ! No given canon_id => all
    if (present(canon_id)) then
      loop_i = canon_id
      loop_f = canon_id
    else
      loop_i = 1
      loop_f = ncid
    end if
    
    do cid = loop_i, loop_f
      if (cid .le. dof%m_psi*dof%n_psi) then
        ! psi
        psi_id = cid
      
        i0 = k2i(psi_id, dof%m_psi)
        j0 = k2j(psi_id, dof%m_psi)
    
        dof%psi%data(i0,j0) = canon(cid)
      else
        ! K
        K_id = cid - dof%m_psi*dof%n_psi
        
        Kloc_id = mod(K_id-1, dof%m_K*dof%n_K)+1
        
        Kcmp1 = dof%m_psi*dof%n_psi
        Kcmp2 = dof%m_psi*dof%n_psi + dof%m_K*dof%n_K
        Kcmp3 = dof%m_psi*dof%n_psi + 2*dof%m_K*dof%n_K
        
        sigma1_sq = canon(Kcmp1 + Kloc_id)
        sigma2_sq = canon(Kcmp2 + Kloc_id)
        phi_K = canon(Kcmp3 + Kloc_id)
        call KCanon_to_KCarte(Kxx, Kyy, Kxy, sigma1_sq, sigma2_sq, phi_K)
        
        i0 = k2i(Kloc_id, dof%m_K)
        j0 = k2j(Kloc_id, dof%m_K)
        
        ! K: defined at corners
        dof%K11%data(i0,j0) = Kxx
        dof%K22%data(i0,j0) = Kyy
        dof%K12%data(i0,j0) = Kxy
      end if
      
!       write(6, *) "cid, i0, j0 = ", cid, i0, j0
    end do
    
  end subroutine convert_canon_to_dof
  
  subroutine evaluate_logPrior(logPrior, dof, sc)
    real(kind=dp), intent(out) :: logPrior
    type(dofdat), intent(in) :: dof
    real(kind=dp), intent(in) :: sc
    integer :: cid, list_i, list_f
    integer :: i0, j0
    type(meshdat) :: mesh
    type(dofdat) :: dof_solver
    real(kind=dp), dimension(:,:), allocatable :: u, v

    logical :: inrange
    
    call allocate(mesh, dof%m_psi+1, dof%m_psi+1)  ! Needs to be identical in both directions
    
    list_i = dof%m_psi*dof%n_psi+1
    list_f = dof%m_psi*dof%n_psi+2*dof%m_K*dof%n_K
    
    call allocate(dof_solver, mesh)
    call intpl_dof_solver(dof_solver, dof, mesh)
    
    inrange = .TRUE.
    
    ! u and v
    allocate(u(mesh%m, mesh%n-1))
    allocate(v(mesh%m-1, mesh%n))
    
    ! u = diff(psi)*sc/(1*sc)
    u = -(dof_solver%psi%data(:,2:mesh%n) - dof_solver%psi%data(:,1:mesh%n-1))
    v = (dof_solver%psi%data(2:mesh%m,:) - dof_solver%psi%data(1:mesh%m-1,:))
    
    ! L1: u = 1.5*sc/sc; v = 2.5 works
    do j0 = 1, size(u, 2)
      do i0 = 1, size(u, 1)
        if (dabs(u(i0, j0)) .ge. 10.0_dp) then
          !inrange = .FALSE.
        end if
      end do
    end do
    
    do j0 = 1, size(v, 2)
      do i0 = 1, size(v, 1)
        if (dabs(v(i0, j0)) .ge. 10.0_dp) then
          !inrange = .FALSE.
        end if
      end do
    end do
    
    deallocate(u)
    deallocate(v)
    call deallocate(dof_solver)
    
    ! Sigma 1 and Sigma 2
    do cid = list_i, list_f
      if ((dof%canon(cid) .le. 10.0_dp*sc) .or. (dof%canon(cid) .gt. 1E8*sc)) then
        inrange = .FALSE.
      end if
    end do
    
    if (inrange) then
      logPrior = 0.0_dp
    else
      logPrior = -1E15
    end if
    

  end subroutine evaluate_logPrior
  
  subroutine propose_dof_canon(dof, canon_SSD, canon_id)
    type(dofdat), intent(inout) :: dof
    real(kind = dp), dimension(:), intent(in) :: canon_SSD
    integer, intent(in) :: canon_id
    
    real(kind=dp), dimension(:), allocatable :: canon
        
    dof%canon(canon_id) = dof%canon(canon_id) + canon_SSD(canon_id)*randn()
    
    call convert_canon_to_dof(dof, dof%canon, canon_id)
    
  end subroutine propose_dof_canon
  
  subroutine propose_dof_EM(dof, canon_SSD, canon_id)
    type(dofdat), intent(inout) :: dof
    real(kind = dp), dimension(:), intent(in) :: canon_SSD
    integer, intent(in) :: canon_id
    
    real(kind=dp), dimension(:), allocatable :: canon
    integer :: canon_K
    
    ! No need to propose psi
    canon_K = canon_id + dof%m_psi*dof%n_psi
    
    dof%canon(canon_K) = dof%canon(canon_K) + canon_SSD(canon_K)*randn()
    
    call convert_canon_to_dof(dof, dof%canon, canon_K)
    
  end subroutine propose_dof_EM
  
  subroutine reverse_dof_to_dof_inv (dof_inv, dof)
    type(dofdat), intent(inout) :: dof_inv
    type(dofdat), intent(in) :: dof
    
    call set(dof_inv, dof)
    call scale(dof_inv%psi, -1.0_dp)
    
  end subroutine reverse_dof_to_dof_inv
  
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

  subroutine propose_dof_Kcmp(dof, dof_id, dof_SSD)
    type(dofdat), intent(inout) :: dof
    integer, intent(in) :: dof_id
    type(dofdat), intent(in) :: dof_SSD

    integer :: i0, j0, i, j
    real(kind = dp) :: sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy
    integer :: cmp, K_id

    ! Define at corners
    if (dof%K11%type_id .ne. -1) &
      call abort_handle("Invalid type_id for K11", __FILE__, __LINE__)
    
    if (dof_id .le. dof%m_K*dof%n_K) then
      K_id = dof_id
      cmp = 0       ! Kxx
    elseif (dof_id .le. 2*dof%m_K*dof%n_K) then
      K_id = dof_id - dof%m_K*dof%n_K
      cmp = 1       ! Kyy
    else 
      K_id = dof_id - 2*dof%m_K*dof%n_K
      cmp = 2       ! Kxy
    end if
    
    i0 = k2i(K_id, dof%m_K)
    j0 = k2j(K_id, dof%m_K)

    ! K: defined at corners
    Kxx = dof%K11%data(i0,j0)
    Kyy = dof%K22%data(i0,j0)
    Kxy = dof%K12%data(i0,j0)

    call KCarte_to_KCanon(sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy)

    ! Perturbation
    if (cmp .eq. 0) then
      sigma1_sq = dabs(sigma1_sq + dof_SSD%K11%data(i0,j0)*randn() )
    elseif (cmp .eq. 1) then
      sigma2_sq = dabs(sigma2_sq + dof_SSD%K22%data(i0,j0)*randn() )
    else
      phi_K = phi_K + dof_SSD%K12%data(i0,j0)*randn()
    end if

    call KCanon_to_KCarte(Kxx, Kyy, Kxy, sigma1_sq, sigma2_sq, phi_K)

    dof%K11%data(i0,j0) = Kxx
    dof%K22%data(i0,j0) = Kyy
    dof%K12%data(i0,j0) = Kxy

  end subroutine propose_dof_Kcmp
  
  subroutine propose_dof_K(dof, dof_id, dof_SSD)
    type(dofdat), intent(inout) :: dof
    integer, intent(in) :: dof_id
    type(dofdat), intent(in) :: dof_SSD

    integer :: i0, j0, i, j
    real(kind = dp) :: sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy

    ! Define at corners
    if (dof%K11%type_id .ne. -1) &
      call abort_handle("Invalid type_id for K11", __FILE__, __LINE__)
      
    i0 = k2i(dof_id, dof%m_K)
    j0 = k2j(dof_id, dof%m_K)

    ! K: defined at corners
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

  subroutine allocate_theta(theta, dof)
    real(kind=dp), dimension(:), allocatable, intent(out) :: theta
    type(dofdat), intent(in) :: dof
    
    integer :: theta_len
    
    theta_len = 0
    theta_len = theta_len + size(dof%psi%data,1)*size(dof%psi%data,2)
    theta_len = theta_len + size(dof%K11%data,1)*size(dof%K11%data,2)
    theta_len = theta_len + size(dof%K22%data,1)*size(dof%K22%data,2)
    theta_len = theta_len + size(dof%K12%data,1)*size(dof%K12%data,2)
    
    allocate(theta(theta_len))

  end subroutine allocate_theta
  
  subroutine convert_dof_to_theta(theta, dof, rsc)
    real(kind=dp), dimension(:), intent(inout) :: theta
    type(dofdat), intent(in) :: dof
    real(kind = dp), intent(in) :: rsc

    integer :: i, j, k
    
    k = 0
    ! Psi
    do j = 1, size(dof%psi%data, 2)
      do i = 1, size(dof%psi%data, 1)
        k = k + 1
        theta(k) = dof%psi%data(i, j) *rsc
      end do
    end do
    ! K11
    do j = 1, size(dof%K11%data, 2)
      do i = 1, size(dof%K11%data, 1)
        k = k + 1
        theta(k) = dof%K11%data(i, j) *rsc
      end do
    end do
    ! K22
    do j = 1, size(dof%K22%data, 2)
      do i = 1, size(dof%K22%data, 1)
        k = k + 1
        theta(k) = dof%K22%data(i, j) *rsc
      end do
    end do
    ! K12
    do j = 1, size(dof%K12%data, 2)
      do i = 1, size(dof%K12%data, 1)
        k = k + 1
        theta(k) = dof%K12%data(i, j) *rsc
      end do
    end do
    
  end subroutine convert_dof_to_theta
  
  subroutine convert_theta_to_dof(dof, theta, sc)
    type(dofdat), intent(in) :: dof
    real(kind=dp), dimension(:), intent(in) :: theta
    real(kind = dp), intent(in) :: sc

    integer :: i, j, k
    
    k = 0
    ! Psi
    do j = 1, size(dof%psi%data, 2)
      do i = 1, size(dof%psi%data, 1)
        k = k + 1
        dof%psi%data(i, j) = theta(k) *sc
      end do
    end do
    ! K11
    do j = 1, size(dof%K11%data, 2)
      do i = 1, size(dof%K11%data, 1)
        k = k + 1
        dof%K11%data(i, j) = theta(k) *sc
      end do
    end do
    ! K22
    do j = 1, size(dof%K22%data, 2)
      do i = 1, size(dof%K22%data, 1)
        k = k + 1
        dof%K22%data(i, j) = theta(k) *sc
      end do
    end do
    ! K12
    do j = 1, size(dof%K12%data, 2)
      do i = 1, size(dof%K12%data, 1)
        k = k + 1
        dof%K12%data(i, j) = theta(k) *sc
      end do
    end do
    
  end subroutine convert_theta_to_dof
  
  subroutine compute_dJdm_Cid(dJdm, J, J_old, m, m_old, canon_id)
    real(kind=dp), dimension(:, :), intent(out) :: dJdm
    real(kind=dp), dimension(:), intent(in) :: J, J_old, m, m_old
    integer, intent(in) :: canon_id
    
    integer :: i
    real(kind=dp) :: dm
    
    dm = m(canon_id)-m_old(canon_id)
    do i = 1, size(J, 1)
      dJdm(i,canon_id) = (J(i) - J_old(i))/dm
    end do
    
  end subroutine compute_dJdm_Cid
  
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
