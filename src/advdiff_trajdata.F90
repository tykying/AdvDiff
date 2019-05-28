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
  public :: initialise_traj, simulate_traj, timestep_traj
  public :: normalise_traj
  public :: traj2jumps, print_info
  public :: check_normalised_pos, reflective_BC
  public :: ij2k, k2i, k2j
  public :: read_uniform_h

  interface allocate
    module procedure allocate_trajdat, allocate_mesh, &
                     allocate_jumpsdat_ptr, allocate_jumpsdat
  end interface allocate

  interface deallocate
    module procedure deallocate_trajdat, deallocate_jumpsdat
  end interface deallocate

  interface eval_fldpt
    module procedure eval_fldpt_rect, eval_fldpt_rect_intpl
  end interface eval_fldpt

  interface print_info
    module procedure print_info_mesh, print_info_jumps, print_info_traj
  end interface print_info

  interface ij2k
    module procedure ij2k_mesh
  end interface ij2k

  interface k2i
    module procedure k2i_mesh
  end interface k2i

  interface k2j
    module procedure k2j_mesh
  end interface k2j

  type jumpsdat
    integer :: njumps
    real(kind = dp), dimension(:, :), pointer :: alpha_i
    real(kind = dp), dimension(:, :), pointer :: alpha_f
    real(kind = dp), dimension(:), pointer :: h
    integer, dimension(:), pointer :: k_f
  end type jumpsdat

  type trajdat
    integer :: nparticles, nts0  ! nts0 = nts + 1 to include initial positions
    real(kind = dp), dimension(:, :), pointer :: x
    real(kind = dp), dimension(:, :), pointer :: y
    real(kind = dp), dimension(:, :), pointer :: t
  end type trajdat

  type meshdat
    integer :: m, n, Ncell
    integer :: m_reflv, n_reflv   ! Related to DOF in inference.F90
  end type meshdat

contains
  subroutine allocate_trajdat(traj, nparticles, nts0)
    type(trajdat), intent(inout) :: traj
    integer, intent(in) :: nparticles, nts0

    traj%nparticles = nparticles
    traj%nts0 = nts0

    allocate(traj%x(traj%nparticles, traj%nts0))
    allocate(traj%y(traj%nparticles, traj%nts0))
    allocate(traj%t(traj%nparticles, traj%nts0))

  end subroutine allocate_trajdat

  subroutine deallocate_trajdat(traj)
    type(trajdat), intent(inout) :: traj

    deallocate(traj%x)
    deallocate(traj%y)
    deallocate(traj%t)

  end subroutine deallocate_trajdat

  subroutine allocate_jumpsdat_ptr(jumps, mesh)
    type(jumpsdat), dimension(:), pointer, intent(inout) :: jumps
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

  subroutine deallocate_jumpsdat(jumps_c)
    type(jumpsdat), intent(inout) :: jumps_c

    deallocate(jumps_c%alpha_i)
    deallocate(jumps_c%alpha_f)
    deallocate(jumps_c%h)
    deallocate(jumps_c%k_f)

  end subroutine deallocate_jumpsdat


  subroutine allocate_mesh(mesh, m, n, m_reflv, n_reflv)
    type(meshdat), intent(inout) :: mesh
    integer, intent(in) :: m, n
    integer, intent(in) :: m_reflv, n_reflv

    mesh%m = m*m_reflv
    mesh%n = n*n_reflv
    mesh%m_reflv = m_reflv
    mesh%n_reflv = n_reflv
    mesh%Ncell = mesh%m * mesh%n
    
  end subroutine allocate_mesh

  pure real(kind = dp) function read_uniform_h(jumps)
    type(jumpsdat), dimension(:), pointer, intent(in) :: jumps
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

    CellInd = ceiling(x*m)
  end function CellInd

  pure real(kind = dp) function alpha(x,m)
    ! x: range from [0, 1]; bin_size = 1/m
    real(kind = dp), intent(in) :: x
    integer, intent(in) :: m

    real(kind = dp) :: xm

    xm = x*m
    alpha = xm - floor(xm) - 0.5_dp
  end function alpha

  pure integer function ij2k_mesh(i,j,mesh)
    integer, intent(in) :: i, j
    type(meshdat), intent(in) :: mesh

    ij2k_mesh = (j-1)*mesh%m + i
  end function ij2k_mesh

  pure integer function k2i_mesh(k,mesh)
    integer, intent(in) :: k
    type(meshdat), intent(in) :: mesh

    k2i_mesh = mod(k-1, mesh%m) + 1
  end function k2i_mesh

  pure integer function k2j_mesh(k,mesh)
    integer, intent(in) :: k
    type(meshdat), intent(in) :: mesh

    k2j_mesh = (k-1)/ mesh%m + 1
  end function k2j_mesh

  pure real(kind = dp) function eval_fldpt_rect(fld, mesh, k)
    type(field), intent(in) :: fld
    type(meshdat), intent(in) :: mesh
    integer, intent(in) :: k

    integer :: i, j

    i = k2i(k,mesh)
    j = k2j(k,mesh)

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

    i = k2i(k,mesh)
    j = k2j(k,mesh)

    if (alpha_f(1) < 0.0_dp) then
      i0 = max(1, i-1)
      i1 = i
      alpha_x = alpha_f(1) + 1.0_dp
    else
      i0 = i
      i1 = min(mesh%m, i+1)
      alpha_x = alpha_f(1)
    end if

    if (alpha_f(2) < 0.0_dp) then
      j0 = max(1, j-1)
      j1 = j
      alpha_y = alpha_f(2) + 1.0_dp
    else
      j0 = j
      j1 = min(mesh%n, j+1)
      alpha_y = alpha_f(2)
    end if

    ! Piecewise constant within the cell
    !eval_fldpt_rect = fld%data(i,j)

    eval_fldpt_rect_intpl = &
        eval_field(fld,i0,j0)*(1.0_dp-alpha_x)*(1.0_dp-alpha_y) &
      + eval_field(fld,i1,j0)*alpha_x*(1.0_dp-alpha_y) &
      + eval_field(fld,i0,j1)*(1.0_dp-alpha_x)*alpha_y &
      + eval_field(fld,i1,j1)*alpha_x*alpha_y  ! in advdiff_field.F90

  end function eval_fldpt_rect_intpl

  subroutine initialise_traj(traj, npartc, nts, q, uniform)
    type(trajdat), intent(out) :: traj
    integer, intent(in) :: nts, npartc   ! nts + 1 = nts
    type(field), intent(in) :: q
    integer, optional, intent(in) :: uniform

    integer :: part, gl, i, j, nzc
    integer :: init_type = 0
    integer :: ii, jj
    integer :: npartcd
    gl = q%glayer

    ! count number of cells
    nzc = 0 ! non-zero cells
    do j= 1, q%n
      do i = 1, q%m
        if (q%data(i+gl, j+gl) > 0.0_dp) then
          nzc = nzc + 1
        end if
      end do
    end do

    write(6, "(a)") ""
    write(6, "(a)") " Initialise particle positions:"
    write(6, "(a, i0)") "|  number of non-zero cells = ", nzc
    call allocate(traj, nzc*npartc, nts+1)  !

    if (present(uniform)) init_type = uniform

    if (init_type .eq. 0) then
      ! Default: place at centre
      ! Place all particles at the cell centres
      nzc = 0
      do j= 1, q%n
        do i = 1, q%m
          if (q%data(i+gl, j+gl) > 0.0_dp) then
            do part = 1, npartc
              traj%x(part+nzc*npartc, 1) = real(i)+0.5_dp
              traj%y(part+nzc*npartc, 1) = real(j)+0.5_dp
              traj%t(part+nzc*npartc, 1) = 0.0_dp
            end do
            nzc = nzc + 1
          end if
        end do
      end do

    elseif (init_type .eq. 1) then
      ! Uniformly initialise particles inside a cell
      write(6, *) "Initialise particles uniformly in cells"
      nzc = 0
      npartcd = int(sqrt(real(npartc))+0.5_dp)
      if (abs(npartcd*npartcd - npartc) > 1e-14) then
        write(6, *) "npartc is not a sq number"
        stop 1
      end if
      do j= 1, q%n
        do i = 1, q%m
          if (q%data(i+gl, j+gl) > 0.0_dp) then
            ! Inside a cell with positive q
            part = 1
            do ii = 1, npartcd
              do jj = 1, npartcd
              ! Idea: partition [0,1] into 2*sqrt(npartc) edges
              ! Place particle only on odd number edges -> ensure uniformity across cells
                traj%x(part+nzc*npartc, 1) = real(i)+real(2*ii-1)/(2*npartcd)
                traj%y(part+nzc*npartc, 1) = real(j)+real(2*jj-1)/(2*npartcd)
                traj%t(part+nzc*npartc, 1) = 0.0_dp
                part = part + 1
              end do
            end do
            nzc = nzc + 1

          end if
        end do
      end do

    end if

    ! Initialise random number
    call init_random_seed()
  end subroutine initialise_traj

  subroutine timestep_traj(traj, dt, u_fld, v_fld, K11_fld, K22_fld, K12_fld)
    type(trajdat), intent(inout) :: traj
    real(kind = dp), intent(in) :: dt
    type(field), intent(in) :: u_fld, v_fld, K11_fld, K22_fld, K12_fld

    integer :: nparticles
    integer :: part, sub_ts
    real(kind = dp) ::  dt_sub
    integer :: i, j, nsub_ts

    real(kind = dp), dimension(:, :), pointer :: lu_fld, lv_fld
    real(kind = dp), dimension(:, :), pointer :: lK11_fld, lK22_fld, lK12_fld

    real(kind = dp) :: res_x, res_y, x, y
    real(kind = dp), dimension(2) :: dW
    real(kind = dp), dimension(2,2) :: K, Krt
    integer :: m, n

    lu_fld => u_fld%data
    lv_fld => v_fld%data
    lK11_fld => K11_fld%data
    lK12_fld => K12_fld%data
    lK22_fld => K22_fld%data

    m = K11_fld%m
    n = K11_fld%n

    ! Initialise particles
    nparticles = size(traj%x, 1)

    nsub_ts = 2  ! Number of sub-timestep
    dt_sub = dt/nsub_ts

    do part = 1, nparticles
      ! Copy previous position
      x = traj%x(part, 1)
      y = traj%y(part, 1)

      do sub_ts = 1, nsub_ts
       i = floor(x)
       j = floor(y)
       res_x = x - real(i)
       res_y = y - real(j)

       !! Advection, velocity = simple interpolation
       x = x + ((1.0_dp-res_x)*lu_fld(i, j) +res_x*lu_fld(i+1, j))*dt_sub
       y = y + ((1.0_dp-res_y)*lv_fld(i, j) +res_y*lv_fld(i, j+1))*dt_sub

       ! Diffuion
       K(1,1) = lK11_fld(i,j)
       K(1,2) = lK12_fld(i,j)
       K(2,1) = lK12_fld(i,j)
       K(2,2) = lK22_fld(i,j)

       call sqrt_matrix(Krt, K)

       ! TODO: use lapack for randn()
       dW = (/randn()*dsqrt(2.0*dt_sub), randn()*dsqrt(2.0*dt_sub)/)
       x = x + Krt(1,1)*dW(1) + Krt(1,2)*dW(2)
       y = y + Krt(2,1)*dW(1) + Krt(2,2)*dW(2)

       !! Impose reflective BC
       x = reflective_BC(x, m)
       y = reflective_BC(y, n)
      end do

      traj%x(part, 1) = x
      traj%y(part, 1) = y
      traj%t(part, 1) = traj%t(part, 1) + dt
    end do

  end subroutine timestep_traj

  subroutine simulate_traj(traj, t, u_fld, v_fld, K11_fld, K22_fld, K12_fld)
    type(trajdat), intent(inout) :: traj
    real(kind = dp), intent(in) :: t
    type(field), intent(in) :: u_fld, v_fld, K11_fld, K22_fld, K12_fld

    integer :: nparticles
    integer :: part, t_ind, sub_ts
    real(kind = dp) :: dt, dt_sub
    integer :: i, j, nsub_ts

    real(kind = dp), dimension(:, :), pointer :: lu_fld, lv_fld
    real(kind = dp), dimension(:, :), pointer :: lK11_fld, lK22_fld, lK12_fld

    real(kind = dp) :: res_x, res_y, x, y
    real(kind = dp), dimension(2) :: dW
    real(kind = dp), dimension(2,2) :: K, Krt
    integer :: m, n

    lu_fld => u_fld%data
    lv_fld => v_fld%data
    lK11_fld => K11_fld%data
    lK12_fld => K12_fld%data
    lK22_fld => K22_fld%data

    m = K11_fld%m
    n = K11_fld%n

    ! Initialise particles
    nparticles = size(traj%x, 1)

    dt = t/(size(traj%x, 2)-1)  ! nts = (size(traj%x, 2)-1)
    nsub_ts = 2  ! Number of sub-timestep
    dt_sub = dt/nsub_ts

    do t_ind = 2, size(traj%x, 2)
      do part = 1, nparticles
        ! Copy previous position
        x = traj%x(part, t_ind-1)
        y = traj%y(part, t_ind-1)

        do sub_ts = 1, nsub_ts
         i = floor(x)
         j = floor(y)
         res_x = x - real(i)
         res_y = y - real(j)

         !! Advection, velocity = simple interpolation
         x = x + ((1.0_dp-res_x)*lu_fld(i, j) +res_x*lu_fld(i+1, j))*dt_sub
         y = y + ((1.0_dp-res_y)*lv_fld(i, j) +res_y*lv_fld(i, j+1))*dt_sub

         !! Diffusion
         K(1,1) = lK11_fld(i,j)
         K(1,2) = lK12_fld(i,j)
         K(2,1) = lK12_fld(i,j)
         K(2,2) = lK22_fld(i,j)

         call sqrt_matrix(Krt, K)

         ! TODO: use lapack for randn()
         dW = (/randn()*dsqrt(2.0*dt_sub), randn()*dsqrt(2.0*dt_sub)/)
         x = x + Krt(1,1)*dW(1) + Krt(1,2)*dW(2)
         y = y + Krt(2,1)*dW(1) + Krt(2,2)*dW(2)

         !! Impose reflective BC
         x = reflective_BC(x, m)
         y = reflective_BC(y, n)
        end do

        traj%x(part, t_ind) = x
        traj%y(part, t_ind) = y
        traj%t(part, t_ind) = traj%t(part, t_ind-1) + dt

      end do
    end do

  end subroutine simulate_traj

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

  pure real(kind=dp) function reflective_BC(x, m)
    real(kind=dp), intent(in) :: x
    integer, intent(in) :: m
    real(kind=dp) :: x_tmp

    ! Reflect at the left: 1+(1-x)
    x_tmp = max(x, 2.0_dp-x)

    ! Reflect at the right: (m+1)-(x-(m+1))
    reflective_BC = min(x_tmp, 2.0_dp*(m+1)-x_tmp)
  end function

  ! Convert trajectory (normalised to [0,1]) into jumps wrt mesh
  subroutine traj2jumps(jumps, traj, mesh)
    type(jumpsdat), dimension(:), pointer, intent(inout) :: jumps
    type(trajdat), intent(in) :: traj
    type(meshdat), intent(in) :: mesh

    integer, dimension(:), pointer :: njumps_ptr
    integer, dimension(:,:), pointer :: Jcell_ptr

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
        cell = ij2k(i,j,mesh)

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

  pure real(kind=dp) function eval_loglik(jumps_c, likfn, mesh)
    type(jumpsdat), intent(in) :: jumps_c
    type(field), intent(in) :: likfn
    type(meshdat), intent(in) :: mesh

    integer :: jump

    eval_loglik = 0.0_dp
    do jump = 1, jumps_c%njumps
      eval_loglik = eval_loglik + dlog(eval_fldpt(likfn, mesh, jumps_c%k_f(jump)))
    end do

  end function eval_loglik

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
    write(6, "(a,"//int_chr//",a,"//int_chr//")") &
      "|  mesh%m_reflv, n_reflv", mesh%m_reflv, " , ", mesh%n_reflv
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
    type(jumpsdat), dimension(:), pointer, intent(in) :: jumps
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
