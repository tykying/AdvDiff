module advdiff_trajdata

  use advdiff_precision
  use advdiff_debug
  use advdiff_timing
  use advdiff_field
  use advdiff_complib

  implicit none

  private

  public :: allocate, deallocate
  public :: jumpsdat, trajdat, meshdat, eval_fldpt
  public :: normalise_traj
  public :: traj2jumps, traj2jumps_inv, print_info
  public :: check_normalised_pos
  public :: read_uniform_h

  interface allocate
    module procedure allocate_trajdat, allocate_mesh, &
                     allocate_jumpsdat_ptr, allocate_jumpsdat
  end interface allocate

  interface deallocate
    module procedure deallocate_trajdat, deallocate_jumpsdat, deallocate_jumpsdat_ptr
  end interface deallocate

  interface eval_fldpt
    module procedure eval_fldpt_rect, eval_fldpt_rect_intpl
  end interface eval_fldpt

  interface print_info
    module procedure print_info_mesh, print_info_jumps, print_info_traj
  end interface print_info

  type jumpsdat
    integer :: njumps
    real(kind = dp), dimension(:, :), allocatable :: alpha_i
    real(kind = dp), dimension(:, :), allocatable :: alpha_f
    real(kind = dp), dimension(:), allocatable :: h
    integer, dimension(:), allocatable :: k_f
  end type jumpsdat

  type trajdat
    integer :: nparticles, nts0  ! nts0 = nts + 1 to include initial positions
    real(kind = dp), dimension(:, :), allocatable :: x
    real(kind = dp), dimension(:, :), allocatable :: y
    real(kind = dp), dimension(:, :), allocatable :: t
  end type trajdat

  type meshdat
    integer :: m, n, Ncell
  end type meshdat

contains
  subroutine allocate_trajdat(traj, nparticles, nts0)
    type(trajdat), intent(inout) :: traj
    integer, intent(in) :: nparticles, nts0

    traj%nparticles = nparticles
    traj%nts0 = nts0

    allocate(traj%x(nparticles, nts0))
    allocate(traj%y(nparticles, nts0))
    allocate(traj%t(nparticles, nts0))

  end subroutine allocate_trajdat

  subroutine deallocate_trajdat(traj)
    type(trajdat), intent(inout) :: traj

    deallocate(traj%x)
    deallocate(traj%y)
    deallocate(traj%t)

  end subroutine deallocate_trajdat

  subroutine allocate_jumpsdat_ptr(jumps, mesh)
    type(jumpsdat), dimension(:), allocatable, intent(inout) :: jumps
    type(meshdat), intent(in) :: mesh

    allocate(jumps(mesh%Ncell))

  end subroutine allocate_jumpsdat_ptr

  subroutine allocate_jumpsdat(jumps_c, njumps)
    type(jumpsdat), intent(inout) :: jumps_c
    integer, intent(in) :: njumps

    jumps_c%njumps = njumps
    allocate(jumps_c%alpha_i(2, njumps))
    allocate(jumps_c%alpha_f(2, njumps))
    allocate(jumps_c%h(njumps))
    allocate(jumps_c%k_f(njumps))

  end subroutine allocate_jumpsdat

  subroutine deallocate_jumpsdat_ptr(jumps)
    type(jumpsdat), dimension(:), allocatable, intent(inout) :: jumps
    integer :: cell
    
    do cell = 1, size(jumps,1)
      call deallocate(jumps(cell))
    end do
    deallocate(jumps)
    
  end subroutine deallocate_jumpsdat_ptr
  
  subroutine deallocate_jumpsdat(jumps_c)
    type(jumpsdat), intent(inout) :: jumps_c

    deallocate(jumps_c%alpha_i)
    deallocate(jumps_c%alpha_f)
    deallocate(jumps_c%h)
    deallocate(jumps_c%k_f)

  end subroutine deallocate_jumpsdat


  subroutine allocate_mesh(mesh, m, n)
    type(meshdat), intent(inout) :: mesh
    integer, intent(in) :: m, n
    
    mesh%m = m
    mesh%n = n
    mesh%Ncell = m*n
    
  end subroutine allocate_mesh

  pure real(kind = dp) function read_uniform_h(jumps)
    type(jumpsdat), dimension(:), allocatable, intent(in) :: jumps
    real(kind = dp) :: h
    integer :: cell
    
    ! mesh(1) may not contain any jump
    h = -1.0_dp
    cell = 1
    do while ( h < 0.0_dp )
      if (jumps(cell)%njumps .ge. 1) then
        h = jumps(cell)%h(1)
      else
        cell = cell + 1
      end if
    end do
    
    read_uniform_h = h
    
  end function read_uniform_h
  
  pure integer function CellInd(x,m)
    ! x: range from [0, 1]; bin_size = 1/m
    real(kind = dp), intent(in) :: x
    integer, intent(in) :: m
    real(kind=dp), parameter :: eps = 1D-10

    CellInd = floor(x*m) + 1
!     CellInd = ceiling(x*m+eps/m)  
    ! Note: Cannot use ceiling()-wrong in marginal case
  end function CellInd

  pure real(kind = dp) function alpha(x,m)
    ! x: range from [0, 1]; bin_size = 1/m
    real(kind = dp), intent(in) :: x
    integer, intent(in) :: m

    real(kind = dp) :: xm

    xm = x*m
    alpha = xm - floor(xm) - 0.5_dp
  end function alpha

  pure real(kind = dp) function eval_fldpt_rect(fld, mesh, k)
    type(field), intent(in) :: fld
    type(meshdat), intent(in) :: mesh
    integer, intent(in) :: k

    integer :: i, j

    i = k2i(k,mesh%m)
    j = k2j(k,mesh%m)

    ! Piecewise constant within the cell
    !eval_fldpt_rect = fld%data(i,j)
    eval_fldpt_rect = eval_field(fld, i, j)  ! in advdiff_field.F90

  end function eval_fldpt_rect

  pure real(kind = dp) function eval_fldpt_rect_intpl(fld, mesh, k, alpha_f)
    type(field), intent(in) :: fld
    type(meshdat), intent(in) :: mesh
    integer, intent(in) :: k
    real(kind = dp), dimension(2), intent(in) :: alpha_f

    integer :: i, j, i0, j0, i1, j1
    real(kind=dp) :: alpha_x, alpha_y

    i = k2i(k,mesh%m)
    j = k2j(k,mesh%m)

    if (alpha_f(1) .lt. 0.0_dp) then
      i0 = max(1, i-1)
      i1 = i
      alpha_x = alpha_f(1) + 1.0_dp
    else
      i0 = i
      i1 = min(mesh%m, i+1)
      alpha_x = alpha_f(1)
    end if

    if (alpha_f(2) .lt. 0.0_dp) then
      j0 = max(1, j-1)
      j1 = j
      alpha_y = alpha_f(2) + 1.0_dp
    else
      j0 = j
      j1 = min(mesh%n, j+1)
      alpha_y = alpha_f(2)
    end if

    ! Bilinear interpolation
    eval_fldpt_rect_intpl = &
        eval_field(fld,i0,j0)*(1.0_dp-alpha_x)*(1.0_dp-alpha_y) &
      + eval_field(fld,i1,j0)*alpha_x*(1.0_dp-alpha_y) &
      + eval_field(fld,i0,j1)*(1.0_dp-alpha_x)*alpha_y &
      + eval_field(fld,i1,j1)*alpha_x*alpha_y  ! in advdiff_field.F90

  end function eval_fldpt_rect_intpl

  subroutine normalise_traj(traj, m, n)
    type(trajdat), intent(inout) :: traj
    integer, intent(in) :: m, n
    integer :: t_ind, part

    ! Normalise the positions to [0,1]
    do t_ind = 1, size(traj%x, 2)
      do part = 1, size(traj%x, 1)
        traj%x(part, t_ind) = (traj%x(part, t_ind) - 1.0_dp)/m
        traj%y(part, t_ind) = (traj%y(part, t_ind) - 1.0_dp)/n
      end do
    end do

  end subroutine normalise_traj

  ! Convert trajectory (normalised to [0,1]) into jumps wrt mesh
  subroutine traj2jumps(jumps, traj, mesh)
    type(jumpsdat), dimension(:), intent(inout) :: jumps
    type(trajdat), intent(in) :: traj
    type(meshdat), intent(in) :: mesh

    integer, dimension(:), allocatable :: njumps_ptr
    integer, dimension(:,:), allocatable :: Jcell_ptr

    integer :: part, t_ind, cell
    integer :: i, j

    ! Determine number of jumps in each cell
    allocate(njumps_ptr(mesh%Ncell))
    allocate(Jcell_ptr(size(traj%x, 1), size(traj%x, 2)))

    do cell = 1, mesh%Ncell
      njumps_ptr(cell) = 0
    end do

    do part = 1, size(traj%x, 1)
      do t_ind = 1, size(traj%x, 2)
        i = CellInd(traj%x(part, t_ind),mesh%m)
        j = CellInd(traj%y(part, t_ind),mesh%n)
        cell = ij2k(i,j,mesh%m)

        Jcell_ptr(part, t_ind) = cell
        njumps_ptr(cell) = njumps_ptr(cell) + 1
      end do

      ! Remove the final count
      njumps_ptr(cell) = njumps_ptr(cell) - 1
    end do

    ! Allocate each jumps_c
    do cell = 1, mesh%Ncell
      call allocate(jumps(cell), njumps_ptr(cell))
    end do

    ! Reset njumps_ptr to be used as a counter
    do cell = 1, mesh%Ncell
      njumps_ptr(cell) = 0
    end do

    do part = 1, size(traj%x, 1)
      do t_ind = 1, (size(traj%x, 2)-1)
        cell = Jcell_ptr(part, t_ind)
        njumps_ptr(cell) = njumps_ptr(cell) + 1

        ! Assign jumpdats
        jumps(cell)%alpha_i(1, njumps_ptr(cell)) = alpha(traj%x(part, t_ind),mesh%m)
        jumps(cell)%alpha_i(2, njumps_ptr(cell)) = alpha(traj%y(part, t_ind),mesh%n)
        jumps(cell)%alpha_f(1, njumps_ptr(cell)) = alpha(traj%x(part, t_ind+1),mesh%m)
        jumps(cell)%alpha_f(2, njumps_ptr(cell)) = alpha(traj%y(part, t_ind+1),mesh%n)

        jumps(cell)%h(njumps_ptr(cell)) = traj%t(part, t_ind+1) - traj%t(part, t_ind)
        jumps(cell)%k_f(njumps_ptr(cell)) = Jcell_ptr(part, t_ind+1)
      end do
    end do

    deallocate(njumps_ptr)
    deallocate(Jcell_ptr)

  end subroutine traj2jumps

    ! Convert trajectory (normalised to [0,1]) into jumps wrt mesh
  subroutine traj2jumps_inv(jumps, traj, mesh)
    type(jumpsdat), dimension(:), intent(inout) :: jumps
    type(trajdat), intent(in) :: traj
    type(meshdat), intent(in) :: mesh

    integer, dimension(:), allocatable :: njumps_ptr
    integer, dimension(:,:), allocatable :: Jcell_ptr

    integer :: part, t_ind, cell
    integer :: i, j

    ! Determine number of jumps in each cell
    allocate(njumps_ptr(mesh%Ncell))
    allocate(Jcell_ptr(size(traj%x, 1), size(traj%x, 2)))

    do cell = 1, mesh%Ncell
      njumps_ptr(cell) = 0
    end do

    do part = 1, size(traj%x, 1)
      do t_ind = size(traj%x, 2), 1, -1  ! Backward in time
        i = CellInd(traj%x(part, t_ind),mesh%m)
        j = CellInd(traj%y(part, t_ind),mesh%n)
        cell = ij2k(i,j,mesh%m)

        Jcell_ptr(part, t_ind) = cell
        njumps_ptr(cell) = njumps_ptr(cell) + 1
      end do

      ! Remove the final count
      njumps_ptr(cell) = njumps_ptr(cell) - 1
    end do

    ! Allocate each jumps_c
    do cell = 1, mesh%Ncell
      call allocate(jumps(cell), njumps_ptr(cell))
    end do

    ! Reset njumps_ptr to be used as a counter
    do cell = 1, mesh%Ncell
      njumps_ptr(cell) = 0
    end do

    do part = 1, size(traj%x, 1)
      do t_ind = size(traj%x, 2), 2, -1 ! Backward in time
        cell = Jcell_ptr(part, t_ind)
        njumps_ptr(cell) = njumps_ptr(cell) + 1

        ! Assign jumpdats
        jumps(cell)%alpha_i(1, njumps_ptr(cell)) = alpha(traj%x(part, t_ind),mesh%m) ! Backward in time
        jumps(cell)%alpha_i(2, njumps_ptr(cell)) = alpha(traj%y(part, t_ind),mesh%n) ! Backward in time
        jumps(cell)%alpha_f(1, njumps_ptr(cell)) = alpha(traj%x(part, t_ind-1),mesh%m) ! Backward in time
        jumps(cell)%alpha_f(2, njumps_ptr(cell)) = alpha(traj%y(part, t_ind-1),mesh%n) ! Backward in time

        jumps(cell)%h(njumps_ptr(cell)) = traj%t(part, t_ind) - traj%t(part, t_ind-1)
        jumps(cell)%k_f(njumps_ptr(cell)) = Jcell_ptr(part, t_ind-1)
      end do
    end do

    deallocate(njumps_ptr)
    deallocate(Jcell_ptr)

  end subroutine traj2jumps_inv
  
  pure logical function check_normalised_pos(traj)
    type(trajdat), intent(in) :: traj
    integer :: t_ind, part
    real(kind=dp) :: x, y

    check_normalised_pos = .true.

    do t_ind = 1, size(traj%x, 2)
      do part = 1, size(traj%x, 1)
        x = traj%x(part, t_ind)
        y = traj%y(part, t_ind)

        if ((x < 0.0_dp) .or. (x > 1.0_dp)) then
          check_normalised_pos = .false.
        end if

        if ((y < 0.0_dp) .or. (y > 1.0_dp)) then
          check_normalised_pos = .false.
        end if
      end do
    end do
  end function check_normalised_pos

  subroutine print_info_mesh(mesh)
    type(meshdat), intent(in) :: mesh

    write(6, "(a)") ""
    write(6, "(a)") " Mesh info: "
    write(6, "(a,"//int_chr//",a,"//int_chr//")") &
      "|  mesh%m, n", mesh%m, " , ", mesh%n
    write(6, "(a,"//int_chr//")") "|  mesh%Ncell", mesh%Ncell
  end subroutine print_info_mesh

  subroutine print_info_traj(traj)
    type(trajdat), intent(in) :: traj
    integer :: part

    write(6, "(a)") ""
    write(6, "(a)") " Trajectory info: "
    do part = 1, size(traj%x, 1)
      write(6, "(a, i0, a,"//dp_chr//",a,"//dp_chr//",a,"//dp_chr//",a,"//dp_chr//")") &
        "|  p", part, ": x0, y0 =", traj%x(part, 1), ", ", traj%y(part, 1), &
         "; xf, yf =", traj%x(part, size(traj%x, 2)), ", ", traj%y(part, size(traj%x, 2))
    end do
  end subroutine print_info_traj

  subroutine print_info_jumps(jumps, mesh)
    type(jumpsdat), dimension(:), allocatable, intent(in) :: jumps
    type(meshdat), intent(in) :: mesh

    integer :: cell, jump

    write(6, "(a)") ""
    write(6, "(a)") " Cells with jumps: "
    do cell = 1, mesh%Ncell
      if (jumps(cell)%njumps > 0) then
        write(6, "(a, i0, a, i0)") "|  c", cell, ": njumps = ", jumps(cell)%njumps
      end if
      do jump = 1, jumps(cell)%njumps
        !write(6, "(a, i0, a, "//dp_chr//", a,"//dp_chr//", a, i0)") "j", jump, ": alpha_i = ", jumps(cell)%alpha_i(1, jump), ", ", jumps(cell)%alpha_i(2, jump), ", next in c", jumps(cell)%k_f(jump)
      end do
    end do
    !write(6, "(a)") "---------------------------------- "
    write(6, "(a)") ""

  end subroutine print_info_jumps

end module advdiff_trajdata
