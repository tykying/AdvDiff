module advdiff_inference
  use advdiff_precision
  use advdiff_debug
  use advdiff_timing
  use advdiff_field
  use advdiff_timestep
  use advdiff_trajdata
  use advdiff_complib

  implicit none

!   public :: initialise_loglik, eval_cell_loglik, update_logPost
  public :: dofdat, IndFndat, priordat
  public :: allocate, deallocate, set
  public :: validate_IndFn
  public :: eval_INDk_loglik
  public :: reverse_dof_to_dof_inv
  public :: convert_dof_to_theta, convert_theta_to_dof
  public :: convert_dof_to_canon, convert_canon_to_dof
  public :: propose_dof_canon, propose_dof_EM
  public :: evaluate_logPost_OMP
  public :: tune_SSD

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
    integer, dimension(:), allocatable :: INDk_list
  end type IndFndat
  
  type priordat
    real(kind=dp) :: sigma_min, sigma_max
  end type priordat

  interface allocate
    module procedure allocate_dof, allocate_IndFn, allocate_dof_solver, &
                      allocate_theta, allocate_prior_param
  end interface allocate

  interface deallocate
    module procedure deallocate_dof, deallocate_IndFn, deallocate_prior_param
  end interface deallocate

  interface set
    module procedure set_dof
  end interface set
  
  interface print_info
    module procedure print_info_inference, print_info_rescaled, &
                      print_info_IndFn, print_info_RV
  end interface print_info

contains
  ! MPI
  pure integer function count_NInd(my_id, num_procs, NInd)
    ! Input: my_id
    ! Output: number of IndFunc attributed to my_id
    integer, intent(in) :: my_id, num_procs, NInd

    integer :: IndFunc

    count_NInd = 0
    do IndFunc = 1, NInd
      if (get_partition(IndFunc, num_procs) .eq. my_id) count_NInd = count_NInd + 1
    end do

  end function count_NInd

  pure integer function get_partition(k, num_procs)
    ! Input: index of a indicator function k
    !        total number of nodes = num_procs
    ! Output: if (my_id == output) then do the job
    integer, intent(in) :: k, num_procs

    ! Mod partition
    get_partition = mod(k, num_procs)

  end function get_partition
    
  ! Non-MPI below
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
    type_id = -1
    call allocate(dof%psi, m_psi, n_psi, 'psi', type_id) ! NOT include the boundaries
    call allocate(dof%K11, m_K, n_K, 'K11', type_id) ! Include the boundaries
    call allocate(dof%K22, m_K, n_K, 'K22', type_id)
    call allocate(dof%K12, m_K, n_K, 'K12', type_id)
    
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
  
  subroutine allocate_IndFn(InfFn, m_Ind, n_Ind, my_id, num_procs)
    type(IndFndat), intent(inout) :: InfFn
    integer, intent(in) :: m_Ind, n_Ind
    integer, optional, intent(in) :: my_id, num_procs
    
    integer :: NInd_total
    integer :: INDk, i

    ! Assume no overlapping
    InfFn%m_Ind = m_Ind
    InfFn%n_Ind = n_Ind
    NInd_total = m_Ind*n_Ind
    
    if (present(my_id) .and. present(num_procs)) then
      ! Parallelise IndFn
      InfFn%NInd = count_NInd(my_id, num_procs, NInd_total)
      write(6, "(a,i0,a,i0)") "my_id = ", my_id, "; InfFn%NInd = ", InfFn%NInd

      if (InfFn%NInd .gt. 0) then
        allocate(InfFn%INDk_list(InfFn%NInd))
        
        i = 0
        do INDk = 1, NInd_total
          if (get_partition(INDk, num_procs) .eq. my_id) then
            i = i + 1
            InfFn%INDk_list(i) = INDk
          end if
        end do
        
        if (i .ne. InfFn%NInd) then
          write(6, "(a,i0,i0)") "Inconsistent INDk_list: i, NInd = ", i, InfFn%NInd
          stop
        end if
      else
        allocate(InfFn%INDk_list(1))
        InfFn%INDk_list(1) = 0
      end if
    else
      InfFn%NInd = NInd_total
      allocate(InfFn%INDk_list(InfFn%NInd))
      
      do INDk = 1, InfFn%NInd
        InfFn%INDk_list(INDk) = INDk
      end do
    end if
    
  end subroutine allocate_IndFn
  
  subroutine deallocate_IndFn(InfFn)
    type(IndFndat), intent(inout) :: InfFn
    
    deallocate(InfFn%INDk_list)

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
!    real(kind=dp), parameter :: kappa_scale = 10000.0_dp            ! Layer 1
    real(kind=dp), parameter :: kappa_scale = 0.2_dp*10000.0_dp   ! Layer 2
    real(kind=dp), parameter :: psi_scale = 10.0_dp*1000.0_dp  ! Irrelevant
    
    if (present(sc)) then
      ssc = sc
    else
      ssc = 1.0_dp
    end if
    
    ! Initialise dof
    call zeros(dof%psi)
    call set(dof%K11, 0.5_dp*kappa_scale*ssc)
    call set(dof%K22, 0.5_dp*kappa_scale*ssc)
    call set(dof%K12, 0.0_dp*kappa_scale*ssc)

    call imprint_canon(dof)
    
  end subroutine init_dof

  subroutine set_dof(dof, dof_in)
    type(dofdat), intent(out) :: dof
    type(dofdat), intent(in) :: dof_in
    
    integer :: i

    call set(dof%psi, dof_in%psi)
    call set(dof%K11, dof_in%K11)
    call set(dof%K22, dof_in%K22)
    call set(dof%K12, dof_in%K12)

    dof%SlogPost = dof_in%SlogPost
    dof%occurrence = dof_in%occurrence
    dof%ndof = dof_in%ndof

    do i = 1, size(dof_in%canon)
      dof%canon(i) = dof_in%canon(i)
    end do
  end subroutine set_dof
  
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
  
  subroutine initialise_q_Gauss(q, Gauss_param, mesh)
    type(field), intent(inout) :: q
    real(kind=dp), dimension(3), intent(in) :: Gauss_param    !=(mean_x, mean_y, sigma)
    type(meshdat), intent(in) :: mesh

    integer :: i, j
    real(kind=dp) :: x, y, x0, y0, Sigma0
    real(kind=dp), parameter :: Pi = 4.0_dp * atan (1.0_dp)
    real(kind=dp) :: tmp
    
    x0 = Gauss_param(1)
    y0 = Gauss_param(2)
    Sigma0 = Gauss_param(3)
    
    do j = 1, mesh%n
      do i = 1, mesh%m
        x = fld_x(i, mesh%m, 2)
        y = fld_x(j, mesh%n, 2)
        
        tmp = dexp(-0.5_dp*((x-x0)**2+(y-y0)**2)/Sigma0**2) &
                            /(2.0_dp*Pi*Sigma0**2)
        
        if (dabs(tmp) .gt. 1D-16) then
          q%data(i,j) = tmp
        else
          q%data(i,j) = 0.0_dp
        end if
      end do
    end do
    
  end subroutine initialise_q_Gauss
  
  real(kind=dp) function eval_INDk_loglik(INDk, jumps, IndFn, mesh, dof_solver, h, nts)
    integer, intent(in) :: INDk
    type(jumpsdat), dimension(:), allocatable, intent(in) :: jumps
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

    eval_INDk_loglik = 0.0_dp

    ! Find out the cells corresponding to the indicator function
    call INDk_to_klist(klist, INDk, IndFn, mesh)
    call allocate(q, mesh%m, mesh%n, 'q', type_id=2)
    
    ! Solve FK equation
    call initialise_q(q, klist, mesh)
    call advdiff_q(q, dof_solver%psi, dof_solver%K11, dof_solver%K22, dof_solver%K12, h, nts)
    
    ! Evaluate likelihood
    do i = 1, size(klist, 1)
      k_i = klist(i)
      do jump = 1, jumps(k_i)%njumps
        ! Bilinear interpolations
        lik = eval_fldpt(q, mesh, jumps(k_i)%k_f(jump), jumps(k_i)%alpha_f(:,jump))
        ! Set the negative probability to be 1D-16
        lik = max(1D-16, lik)
        eval_INDk_loglik = eval_INDk_loglik + dlog(lik)
      end do
    end do

    call deallocate(q)
    deallocate(klist)
    
  end function eval_INDk_loglik
  
  subroutine allocate_dof_solver(dof_solver, mesh)
    type(dofdat), intent(out) :: dof_solver
    type(meshdat), intent(in) :: mesh
    
   ! Defined at corners (! 1= corner; 2= centres; *-1 = DOF)
    call allocate(dof_solver%psi, mesh%m, mesh%n, 'psi', type_id=1)

    call allocate(dof_solver%K11, mesh%m, mesh%n, 'K11', type_id=1)
    call allocate(dof_solver%K22, mesh%m, mesh%n, 'K22', type_id=1)
    call allocate(dof_solver%K12, mesh%m, mesh%n, 'K12', type_id=1)

    ! Other data members: to be assigned when interpolated
    
  end subroutine allocate_dof_solver
  
  subroutine intpl_dof_solver(dof_solver, dof, mesh)
    ! Interpolate from DOF to fields for the solver
    type(dofdat), intent(inout) :: dof_solver  ! Field to use for solver
    type(dofdat), intent(in) :: dof
    type(meshdat), intent(in) :: mesh
    
    integer :: m, n, Mx, My
    real(kind=dp), dimension(:, :), allocatable :: fldij ! Bilinear intepolation
    

    ! Quantities other than the field: same as dof
    dof_solver%m_psi = dof%m_psi
    dof_solver%n_psi = dof%n_psi
    dof_solver%m_K = dof%m_K
    dof_solver%n_K = dof%n_K

    dof_solver%ndof = dof%ndof
    dof_solver%occurrence = dof%occurrence
    dof_solver%SlogPost = dof%SlogPost
    
    allocate(dof_solver%canon(dof%ndof))
    dof_solver%canon = dof%canon
    
    ! fldij(i, j) = field(x_i, y_j) on computational domain
    ! DOF: can be coarse field(x_i, y_j), Fourier coefficients or anything
    
    ! Streamfunction: psi
#if EM_MEAN == 1
     dof_solver%psi%data = 0.0_dp
     dof_solver%psi%data(2:size(dof_solver%psi%data,1)-1, &
                         2:size(dof_solver%psi%data,2)-1) = dof%psi%data
#else
    ! N.B. DOF%psi = non-bondary values, assumed boundary value = 0
    m = dof%psi%m  ! = 15
    n = dof%psi%n

!   ! Bilinear intepolation
    allocate(fldij(m+2, n+2))
    fldij = 0.0_dp  ! Boundary condition for streamfunction
    fldij(2:(m+1), 2:(n+1)) = dof%psi%data
   
    Mx = mesh%m/(m+1)
    My = mesh%n/(n+1)
    ! Assumed periodic(= 0) at first and final row/column
    call bilinear_intpl(dof_solver%psi%data, fldij, Mx, My)
    deallocate(fldij)
#endif

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
        ! Not stictly negative but with a threshold -1D-14
        is_nneg = (is_nneg .and. (data(i,j) .ge. -1D-10))
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
!     real(kind=dp), parameter :: psi_scale = 100.0_dp*1000.0_dp  ! Layer 1
    real(kind=dp), parameter :: psi_scale = 10.0_dp*1000.0_dp  ! Layer 2
    
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
      canon_SSD(k) = 0.25_dp
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
  
  real(kind=dp) function evaluate_logPrior(dof_solver, prior_param)
    type(dofdat), intent(in) :: dof_solver
    type(priordat), intent(in) :: prior_param
    
    integer :: cid, list_i, list_f
    logical :: inrange
    
    inrange = .TRUE.

    ! Sigma 1 and Sigma 2
    list_i = dof_solver%m_psi*dof_solver%n_psi+1
    list_f = dof_solver%m_psi*dof_solver%n_psi+2*dof_solver%m_K*dof_solver%n_K
    
    do cid = list_i, list_f
      if ((dof_solver%canon(cid) .le. prior_param%sigma_min) .or. &
          (dof_solver%canon(cid) .gt. prior_param%sigma_max)) then
        inrange = .FALSE.
      end if
    end do
    
    if (inrange) then
      evaluate_logPrior = 0.0_dp
    else
      evaluate_logPrior = -1D15
    end if
    
    
  end function evaluate_logPrior
  
  subroutine allocate_prior_param(prior_param, sc)
    type(priordat), intent(inout) :: prior_param
    real(kind=dp), intent(in) :: sc
    
    prior_param%sigma_max = 1D5 * sc
    prior_param%sigma_min = 1.0_dp * sc
    
  end subroutine allocate_prior_param
    
  subroutine deallocate_prior_param(prior_param)
    type(priordat), intent(inout) :: prior_param
    
    prior_param%sigma_max = 0.0_dp
    prior_param%sigma_min = 0.0_dp
    
  end subroutine deallocate_prior_param
  
  real(kind=dp) function evaluate_logPost_OMP(prior_param, jumps, IndFn, mesh, dof, h, nts)
    type(priordat), intent(in) :: prior_param
    type(jumpsdat), dimension(:), allocatable, intent(in) :: jumps
    type(IndFndat), intent(in) :: IndFn
    type(meshdat), intent(in) :: mesh
    type(dofdat), intent(in) :: dof
    real(kind=dp), intent(in) :: h
    integer, intent(in) :: nts

    integer :: i, nts_adapted

    type(dofdat) :: dof_solver  ! Refined dof for solver
    real(kind=dp), dimension(:), allocatable :: loglik
    real(kind=dp) :: logPrior, dt_min, dt_arg
    
    call allocate(dof_solver, mesh)
    call intpl_dof_solver(dof_solver, dof, mesh)
    
    dt_arg = h/nts
    dt_min = dt_CFL(dof_solver%psi, 0.20_dp)
    
    nts_adapted = max(int(h/dt_min + 0.5_dp), nts)
    
    ! log-Prior
    logPrior =  evaluate_logPrior(dof_solver, prior_param)
    
    if (logPrior .gt. -1D-15) then
      ! log-likelihood
      allocate(loglik(IndFn%nIND))
      
      !$OMP PARALLEL DO
      do i = 1, IndFn%nIND
        loglik(i) = eval_INDk_loglik(IndFn%INDk_list(i), jumps, IndFn, mesh, dof_solver, h, nts_adapted)
      end do
      !$OMP END PARALLEL DO
    
      evaluate_logPost_OMP = sum(loglik) + logPrior

      deallocate(loglik)
    else
      evaluate_logPost_OMP = logPrior
    end if
    
    call deallocate(dof_solver)
    
  end function evaluate_logPost_OMP
  
  subroutine propose_dof_canon(dof, canon_SSD, canon_id, Z_rv)
    type(dofdat), intent(inout) :: dof
    real(kind = dp), dimension(:), intent(in) :: canon_SSD
    integer, intent(in) :: canon_id
    real(kind = dp), intent(in) :: Z_rv
    
    dof%canon(canon_id) = dof%canon(canon_id) + canon_SSD(canon_id)*Z_rv
    
    call convert_canon_to_dof(dof, dof%canon, canon_id)
    
  end subroutine propose_dof_canon
  
  subroutine propose_dof_EM(dof, canon_SSD, canon_id, Z_rv)
    type(dofdat), intent(inout) :: dof
    real(kind = dp), dimension(:), intent(in) :: canon_SSD
    integer, intent(in) :: canon_id
    real(kind = dp), intent(in) :: Z_rv

    integer :: canon_K
    
    ! No need to propose psi
    canon_K = canon_id + dof%m_psi*dof%n_psi
    
    dof%canon(canon_K) = dof%canon(canon_K) + canon_SSD(canon_K)*Z_rv
    
    call convert_canon_to_dof(dof, dof%canon, canon_K)
    
  end subroutine propose_dof_EM
  
  subroutine reverse_dof_to_dof_inv (dof_inv, dof)
    type(dofdat), intent(inout) :: dof_inv
    type(dofdat), intent(in) :: dof
    
    call set(dof_inv, dof)
    call scale(dof_inv%psi, -1.0_dp)
    
  end subroutine reverse_dof_to_dof_inv
  
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
  
  subroutine tune_SSD(canon_SSD, accept_ratio)
    real(kind=dp), dimension(:), intent(inout) :: canon_SSD
    real(kind=dp), dimension(:), intent(in) :: accept_ratio
    
    integer :: canon_id
    
    do canon_id = 1, size(canon_SSD, 1)
      if (accept_ratio(canon_id) .ge. 0.50_dp) then
        canon_SSD(canon_id) = canon_SSD(canon_id) * 1.50_dp
      elseif  (accept_ratio(canon_id) .le. 0.15_dp) then
        canon_SSD(canon_id) = canon_SSD(canon_id) * 0.50_dp
      end if
    end do
    
  end subroutine tune_SSD
    
  subroutine print_info_IndFn(IndFn)
    type(IndFndat), intent(in) :: IndFn
    
    write(6, "(a)") ""
    write(6, "(a)") " Indicator function info: "
    write(6, "(a,"//int_chr//",a,"//int_chr//")") &
      "|  IndFn%m_ind, n_ind", IndFn%m_ind, " , ", IndFn%n_ind
    write(6, "(a,"//int_chr//")") "|  IndFn%nIND", IndFn%nIND
    
    FLUSH(6)
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
  
  subroutine print_info_RV(RV)
    real(kind = dp), dimension(:), intent(in) :: RV
    
    write(6, "(a,"//sp_chr//")") "|  E(X) = ", sum(RV)/size(RV, 1)
    write(6, "(a,"//sp_chr//")") "|  Var(X) = ", sum(RV**2)/size(RV, 1)-(sum(RV)/size(RV, 1))**2
    
    FLUSH(6)
  end subroutine print_info_RV

end module advdiff_inference
