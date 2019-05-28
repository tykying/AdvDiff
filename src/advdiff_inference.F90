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
  public :: evaluate_loglik_IND, eval_INDk_loglik

  type dofdat
    integer :: m, n
    integer :: ndof
    type(field) :: psi, K11, K22, K12
  end type dofdat
  
  type IndFndat
    integer :: m_Ind, n_Ind
    integer :: NInd
  end type IndFndat

  type(timer), save :: loglik_timer, intpl_timer, reffld_timer

  interface allocate
    module procedure allocate_dof, allocate_rdof, allocate_IndFn
  end interface allocate

  interface deallocate
    module procedure deallocate_dof, deallocate_IndFn
  end interface deallocate

  interface set
    module procedure set_dof
  end interface set

  interface ij2k
    module procedure ij2k_dof
  end interface ij2k

  interface k2i
    module procedure k2i_dof
  end interface k2i

  interface k2j
    module procedure k2j_dof
  end interface k2j
  
  interface print_info
    module procedure print_info_inference, print_info_rescaled, print_info_IndFn
  end interface print_info

contains
  subroutine allocate_dof(dof, m, n)
    type(dofdat), intent(inout) :: dof
    integer, intent(in) :: m, n

    integer :: type_id ! 1-> corner; 0-> centres
    
    dof%m = m
    dof%n = n
    dof%ndof = m*n

    ! Defined at corners
    type_id = 1
    call allocate(dof%psi, m, n, 'psi', glayer=0, type_id=type_id)

    ! Defined at cell centre
    type_id = 0
!     write(6, *) "- Now switch to Kcorner -"
!     type_id = 1
    call allocate(dof%K11, m, n, 'K11', glayer=0, type_id=type_id)
    call allocate(dof%K22, m, n, 'K22', glayer=0, type_id=type_id)
    call allocate(dof%K12, m, n, 'K12', glayer=0, type_id=type_id)

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

  pure integer function ij2k_dof(i,j,dof)
    integer, intent(in) :: i, j
    type(dofdat), intent(in) :: dof

    ij2k_dof = (j-1)*dof%m + i
  end function ij2k_dof

  pure integer function k2i_dof(k,dof)
    integer, intent(in) :: k
    type(dofdat), intent(in) :: dof

    k2i_dof = mod(k-1, dof%m) + 1
  end function k2i_dof

  pure integer function k2j_dof(k,dof)
    integer, intent(in) :: k
    type(dofdat), intent(in) :: dof

    k2j_dof = (k-1)/ dof%m + 1
  end function k2j_dof

  ! For parallel computing
  pure integer function count_cell(my_id, num_procs, mesh)
    ! Input: my_id
    ! Output: number of cells attributed to my_id
    type(meshdat), intent(in) :: mesh
    integer, intent(in) :: my_id, num_procs

    integer :: cell

    count_cell = 0
    do cell = 1, mesh%ncell
      if (get_partition(cell, num_procs) .eq. my_id) count_cell = count_cell + 1
    end do

  end function count_cell

  pure integer function get_partition(k, num_procs)
    ! Input: index of a cell k
    !        total number of nodes = num_procs
    ! Output: if (my_id == output) then do the job
    integer, intent(in) :: k, num_procs

    ! Mod partition
    get_partition = mod(k, num_procs)

  end function get_partition

  subroutine evaluate_loglik_IND(loglik, jumps, IndFn, mesh, dof, h, nts, my_id, num_procs)
    real(kind=dp), dimension(:), intent(out) :: loglik
    type(jumpsdat), dimension(:), pointer, intent(in) :: jumps
    type(IndFndat), intent(in) :: IndFn
    type(meshdat), intent(in) :: mesh
    type(dofdat), intent(in) :: dof
    real(kind=dp), intent(in) :: h
    integer, intent(in) :: nts
    integer, optional, intent(in) :: my_id, num_procs

    integer :: INDk
    integer :: node_id, nproc

    ! Compatibility with single node
    if (present(my_id) .and. present(num_procs)) then
      node_id = my_id
      nproc = num_procs
    else
      node_id = 0
      nproc = 1
    end if

    do INDk = 1, IndFn%nIND
      if (get_partition(INDk, nproc) .eq. node_id) then
        loglik(INDk) = eval_INDk_loglik(INDk, jumps, IndFn, mesh, dof, h, nts)
      else
        loglik(INDk) = 0.0_dp
      end if
    end do

  end subroutine evaluate_loglik_IND
  
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
        klist(klist_len) = ij2k(i,j,mesh)
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
      call indicator_field(q, k2i(k_i, mesh), k2j(k_i, mesh))
    end do
    call scale(q, 1.0_dp/real(nzc, kind=dp))
  
  end subroutine initialise_q
  
  real(kind=dp) function eval_INDk_loglik(INDk, jumps, IndFn, mesh, dof, h, nts)
    integer, intent(in) :: INDk
    type(jumpsdat), dimension(:), pointer, intent(in) :: jumps
    type(IndFndat), intent(in) :: IndFn
    type(meshdat), intent(in) :: mesh
    type(dofdat), intent(in) :: dof
    real(kind=dp), intent(in) :: h
    integer, intent(in) :: nts
    integer, dimension(:), allocatable :: klist

    type(field) :: q
    integer :: jump
    real(kind=dp) :: lik
    integer :: k_i, i

    type(dofdat) :: rdof  ! Refined dof for computation

    call start(loglik_timer)

    eval_INDk_loglik = 0.0_dp

    ! Find out the cells corresponding to the indicator function
    call INDk_to_klist(klist, INDk, IndFn, mesh)
   
    call start(reffld_timer)
    call allocate(rdof, dof, mesh%m_reflv, mesh%n_reflv)
    call allocate(q, rdof%m, rdof%n, 'q', glayer=1,type_id=0)
    call stop(reffld_timer)
    
    ! Solve FK equation
    call initialise_q(q, klist, mesh)
    call advdiff_q(q, rdof%psi, rdof%K11, rdof%K22, rdof%K12, h, nts)
!       if (.not. is_nneg(q%data)) then
!         write(6, "(a,"//dp_chr//")") "Non-negative q!: ", minval(q%data)
!       end if

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
    call deallocate(rdof)

    deallocate(klist)
    
    call stop(loglik_timer)

  end function eval_INDk_loglik

  subroutine allocate_rdof(rdof, dof, m_reflv, n_reflv)
    type(dofdat), intent(out) :: rdof
    type(dofdat), intent(in) :: dof
    integer, intent(in) :: m_reflv, n_reflv
    
    rdof%m = dof%m * m_reflv
    rdof%n = dof%n * n_reflv
    
    call allocate(rdof%psi, dof%psi, m_reflv, n_reflv)
    call allocate(rdof%K11, dof%K11, m_reflv, n_reflv)
    call allocate(rdof%K22, dof%K22, m_reflv, n_reflv)
    call allocate(rdof%K12, dof%K12, m_reflv, n_reflv)

  end subroutine allocate_rdof

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

    !write(6, *) "propose_dof: Proposing two changes at a time"
    if (mod(dof_id, 2) .eq. 0) then
      call propose_dof_K(dof, dof_id, dof_SSD)
      call propose_dof_K(dof, dof%ndof-dof_id+1, dof_SSD)
    else
      call propose_dof_psi(dof, dof_id, dof_SSD)
      call propose_dof_psi(dof, dof%ndof-dof_id+1, dof_SSD)
      call imposeBC_cornerpsi(dof)
    end if
    
  end subroutine propose_dof

  subroutine propose_dof_K(dof, dof_id, dof_SSD)
    type(dofdat), intent(inout) :: dof
    integer, intent(in) :: dof_id
    type(dofdat), intent(in) :: dof_SSD

    integer :: i0, j0, i, j
    real(kind = dp) :: sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy

    i0 = k2i(dof_id, dof)
    j0 = k2j(dof_id, dof)

    if ( (dof%K11%type_id .eq. 0) .and. &
      (dof%K22%type_id .eq. 0) .and. (dof%K12%type_id .eq. 0) ) then
      ! K: defined at cell-centres
      Kxx = dof%K11%data(i0,j0)
      Kyy = dof%K22%data(i0,j0)
      Kxy = dof%K12%data(i0,j0)

      call KCarte_to_KCanon(sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy)

      ! Perturbation
      sigma1_sq = dabs(sigma1_sq + dof_SSD%K11%data(i0,j0)*randn() )
      sigma2_sq = dabs(sigma2_sq + dof_SSD%K22%data(i0,j0)*randn() )
      phi_K = phi_K + dof_SSD%K12%data(i0,j0)*randn()

      ! TODO: delete Adhoc
      !sigma2_sq = sigma1_sq
      !phi_K = 0.0_dp
      ! TODO: delete Adhoc

      call KCanon_to_KCarte(Kxx, Kyy, Kxy, sigma1_sq, sigma2_sq, phi_K)

      dof%K11%data(i0,j0) = Kxx
      dof%K22%data(i0,j0) = Kyy
      dof%K12%data(i0,j0) = Kxy
      
    elseif ( (dof%K11%type_id .eq. 1) .and. &
      (dof%K22%type_id .eq. 1) .and. (dof%K12%type_id .eq. 1) ) then
      ! K: defined at corners
      do j = j0, j0+1
      do i = i0, i0+1
        Kxx = dof%K11%data(i,j)
        Kyy = dof%K22%data(i,j)
        Kxy = dof%K12%data(i,j)

        call KCarte_to_KCanon(sigma1_sq, sigma2_sq, phi_K, Kxx, Kyy, Kxy)

        ! Perturbation
        sigma1_sq = dabs(sigma1_sq + dof_SSD%K11%data(i,j)*randn() )
        sigma2_sq = dabs(sigma2_sq + dof_SSD%K22%data(i,j)*randn() )
        phi_K = phi_K + dof_SSD%K12%data(i,j)*randn()

        call KCanon_to_KCarte(Kxx, Kyy, Kxy, sigma1_sq, sigma2_sq, phi_K)

        dof%K11%data(i,j) = Kxx
        dof%K22%data(i,j) = Kyy
        dof%K12%data(i,j) = Kxy
      end do
      end do
      
      ! Impose boundary condition
      call imposeBC_cornerK(dof)
    else
      call abort_handle("Invalid type_id for K", __FILE__, __LINE__)
    end if


  end subroutine propose_dof_K

  subroutine propose_dof_psi(dof, dof_id, dof_SSD)
    type(dofdat), intent(inout) :: dof
    integer, intent(in) :: dof_id
    type(dofdat), intent(in) :: dof_SSD

    integer :: i0, j0

    i0 = k2i(dof_id, dof)
    j0 = k2j(dof_id, dof)

    if (dof%psi%type_id .ne. 1) &
      call abort_handle("Invalid type_id for psi", __FILE__, __LINE__)

    ! Vary the psi values at the corners
    dof%psi%data(i0  ,j0  ) = dof%psi%data(i0  ,j0  ) + dof_SSD%psi%data(i0  ,j0  )*randn()
    dof%psi%data(i0+1,j0  ) = dof%psi%data(i0+1,j0  ) + dof_SSD%psi%data(i0+1,j0  )*randn()
    dof%psi%data(i0  ,j0+1) = dof%psi%data(i0  ,j0+1) + dof_SSD%psi%data(i0  ,j0+1)*randn()
    dof%psi%data(i0+1,j0+1) = dof%psi%data(i0+1,j0+1) + dof_SSD%psi%data(i0+1,j0+1)*randn()

  end subroutine propose_dof_psi

  subroutine imposeBC_cornerpsi(dof)
    type(dofdat), intent(inout) :: dof

    integer :: i, j, mp1, np1

    mp1 = size(dof%psi%data, 1)
    np1 = size(dof%psi%data, 2)

    do j = 1, np1
      dof%psi%data(1, j) = 0.0_dp
      dof%psi%data(mp1, j) = 0.0_dp
    end do

    do i = 1, mp1
      dof%psi%data(i, 1) = 0.0_dp
      dof%psi%data(i, np1) = 0.0_dp
    end do

  end subroutine imposeBC_cornerpsi

  ! TODO: cornerK
  subroutine imposeBC_cornerK(dof)
    type(dofdat), intent(inout) :: dof

    integer :: i, j, mp1, np1

    mp1 = size(dof%K11%data, 1)
    np1 = size(dof%K11%data, 2)

    do j = 1, np1
      dof%K11%data(1, j) = dof%K11%data(2, j)
      dof%K12%data(1, j) = dof%K12%data(2, j)
      dof%K22%data(1, j) = dof%K22%data(2, j)
      dof%K11%data(mp1, j) = dof%K11%data(mp1-1, j)
      dof%K12%data(mp1, j) = dof%K12%data(mp1-1, j)
      dof%K22%data(mp1, j) = dof%K22%data(mp1-1, j)
    end do

    do i = 1, mp1
      dof%K11%data(i, 1) = dof%K11%data(i, 2)
      dof%K12%data(i, 1) = dof%K12%data(i, 2)
      dof%K22%data(i, 1) = dof%K22%data(i, 2)
      dof%K11%data(i, np1) = dof%K11%data(i, np1-1)
      dof%K12%data(i, np1) = dof%K12%data(i, np1-1)
      dof%K22%data(i, np1) = dof%K22%data(i, np1-1)
    end do

  end subroutine imposeBC_cornerK

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
      "|  dof%m, n", dof%m, " , ", dof%n
    write(6, "(a,"//int_chr//")") "|  dof%ndof", dof%ndof
    write(6, "(a)") "| [type_id =0 -> centre; =1 -> corner]"
    write(6, "(a,"//int_chr//")") "|  psi%type_id", dof%psi%type_id
    write(6, "(a,"//int_chr//")") "|  K11%type_id", dof%K11%type_id
    write(6, "(a,"//int_chr//")") "|  K22%type_id", dof%K22%type_id
    write(6, "(a,"//int_chr//")") "|  K12%type_id", dof%K12%type_id
    
    
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
