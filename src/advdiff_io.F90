module advdiff_io

  use advdiff_precision
  use advdiff_debug
  use advdiff_timing
  use advdiff_field
  use advdiff_trajdata
  use advdiff_inference

  implicit none

  private

  public :: read_header, read, write, write_txt, read_theta, read_QGfield

  interface read_header
    module procedure read_field_header, read_traj_header
  end interface read_header

  interface read
    module procedure read_field, read_traj, read_dof, read_ssd
  end interface read

  ! Interface: allowing overloading the 'write' function by the parameters input
  interface write
    module procedure write_field, write_dof, write_ssd
  end interface write
  
  interface write_txt
    module procedure write_array_dp1d, write_array_dp
  end interface write_txt


  type(timer), save :: write_field_timer
  public :: reset_io_timers, print_io_timers, write_theta

  integer, parameter :: input_unit = 10, output_unit = 11
  public :: input_unit, output_unit

contains

  subroutine read_field_header(fld, filename, t, index)
    type(field), intent(out) :: fld
    character(len = *), intent(in) :: filename
    real(kind = dp), intent(out) :: t
    integer, optional, intent(in) :: index

    character(len = 1) :: chr
    character(len = field_name_len) :: name
    character(len = max_line_len) :: format_chr, input_chr
    character(len = len_trim(filename) + 1 + max_int_len) :: lfilename
    integer :: i, m, n, type_id

    if(present(index)) then
      write(lfilename, "(a,a,i0)") trim(filename), "_", index
    else
      lfilename = trim(filename)
    end if

    open(unit = input_unit, file = trim(lfilename) // ".hdr", &
      & status = "old", action = "read")
    read(input_unit, "(a)") input_chr
    if(trim(input_chr) /= "serial") then
      call abort_handle("Invalid header", __FILE__, __LINE__)
    end if
    read(input_unit, "(a)") name
    read(input_unit, "(a)") input_chr
    i = scan(input_chr, " ")
    write(format_chr, "(a,i0,a)") "(i", i - 1, ",a1,"//int_chr//")"
    read(input_chr, trim(format_chr)) m, chr, n
    read(input_unit, "("//int_chr//")") type_id
    read(input_unit, "("//dp_chr//")") t
    close(input_unit)

    call allocate(fld, m, n, name, type_id=type_id)

  end subroutine read_field_header

  subroutine read_field(fld, filename, index)
    type(field), intent(inout) :: fld
    character(len = *), intent(in) :: filename
    integer, optional, intent(in) :: index

    character(len = len_trim(filename) + 1 + max_int_len) :: lfilename

    if(present(index)) then
      write(lfilename, "(a,a,i0)") trim(filename), "_", index
    else
      lfilename = trim(filename)
    end if

    open(unit = input_unit, file = trim(lfilename) // ".dat", &
      & status = "old", access = "stream", action = "read")
    read(input_unit) fld%data
    close(input_unit)

  end subroutine read_field

  subroutine write_field(fld, filename, t, index)
    type(field), intent(in) :: fld
    character(len = *), intent(in) :: filename
    real(kind = dp), intent(in) :: t
    integer, optional, intent(in) :: index

    character(len = len_trim(filename) + 1 + max_int_len) :: lfilename

    call start(write_field_timer)

    if(present(index)) then
      write(lfilename, "(a,a,i0)") trim(filename), "_", index
    else
      lfilename = trim(filename)
    end if

    open(unit = output_unit, file = trim(lfilename) // ".hdr", &
      & status = "replace", action = "write")
    write(output_unit, "(a)") "serial"
    write(output_unit, "(a)") trim(fld%name)
    write(output_unit, "(i0,a,i0)") fld%m, " ", fld%n
    write(output_unit, "(i0)") fld%type_id
    write(output_unit, "("//dp_chr//")") t
    close(output_unit)

    open(unit = output_unit, file = trim(lfilename) // ".dat", &
      & status = "replace", access = "stream", action = "write")
    write(output_unit) fld%data
    close(output_unit)

    call stop(write_field_timer)

  end subroutine write_field

  subroutine write_array_dp1d(array, filename)
    real(kind=dp), dimension(:), intent(in) :: array
    character(len = *), intent(in) :: filename
    character(len = 16) :: collen
    
    write(collen, "(i0)") size(array, 1)

    open(unit = output_unit, file = trim(filename) // ".txt", &
      & status = "replace", action = "write")
    write(output_unit, '('//collen//'('//dp_chr//', ","))') array
    close(output_unit)

  end subroutine write_array_dp1d
  
  subroutine write_array_dp(array, filename)
    real(kind=dp), dimension(:,:), intent(in) :: array
    character(len = *), intent(in) :: filename
    character(len = 16) :: collen
    integer :: i

    write(collen, "(i0)") size(array, 2)

    open(unit = output_unit, file = trim(filename) // ".txt", &
      & status = "replace", action = "write")
    do i = 1, size(array, 1)
      write(output_unit, '('//collen//'('//dp_chr//', ","))') array(i, :)
    end do
    close(output_unit)

  end subroutine write_array_dp

  subroutine read_traj_header(traj, filename, t, index)
    type(trajdat), intent(out) :: traj
    character(len = *), intent(in) :: filename
    real(kind = dp), intent(out) :: t
    integer, optional, intent(in) :: index

    character(len = 1) :: chr
    character(len = field_name_len) :: name
    character(len = max_line_len) :: format_chr, input_chr
    character(len = len_trim(filename) + 1 + max_int_len) :: lfilename
    integer :: i, nparticles, nts0, trajdim

    if(present(index)) then
      write(lfilename, "(a,a,i0)") trim(filename), "_", index
    else
      lfilename = trim(filename)
    end if

    open(unit = input_unit, file = trim(lfilename) // ".hdr", &
      & status = "old", action = "read")
    read(input_unit, "(a)") input_chr
    if(trim(input_chr) /= "serial") then
      call abort_handle("Invalid header", __FILE__, __LINE__)
    end if
    read(input_unit, "(a)") name
    read(input_unit, "(a)") input_chr
    i = scan(input_chr, " ")
    write(format_chr, "(a,i0,a)") "(i", i - 1, ",a1,"//int_chr//")"
    read(input_chr, trim(format_chr)) nparticles, chr, nts0
    read(input_unit, "("//int_chr//")") trajdim
    if (trajdim .ne. 2) then
      call abort_handle("Invalid header", __FILE__, __LINE__)
    end if
    read(input_unit, "("//dp_chr//")") t
    close(input_unit)

    call allocate(traj, nparticles, nts0)

  end subroutine read_traj_header

  subroutine read_traj(traj, filename, index)
    type(trajdat), intent(out) :: traj
    character(len = *), intent(in) :: filename
    integer, optional, intent(in) :: index

    character(len = len_trim(filename) + 1 + max_int_len) :: lfilename
    real(kind = dp) :: t

    if(present(index)) then
      write(lfilename, "(a,a,i0)") trim(filename), "_", index
    else
      lfilename = trim(filename)
    end if

    call read_traj_header(traj, filename, t, index)
    
    open(unit = input_unit, file = trim(lfilename) // "_x.dat", &
      & status = "old", access = "stream", action = "read")
    read(input_unit) traj%x
    close(input_unit)

    open(unit = input_unit, file = trim(lfilename) // "_y.dat", &
      & status = "old", access = "stream", action = "read")
    read(input_unit) traj%y
    close(input_unit)

    open(unit = input_unit, file = trim(lfilename) // "_t.dat", &
      & status = "old", access = "stream", action = "read")
    read(input_unit) traj%t
    close(input_unit)

  end subroutine read_traj

  
  subroutine write_dof(dof, filename, SlogPost, sc, index)
    type(dofdat), intent(in) :: dof
    character(len = *), intent(in) :: filename
    real(kind = dp), intent(in) :: SlogPost, sc
    integer, optional, intent(in) :: index


    type(dofdat) :: dof_rsc

    call allocate(dof_rsc, dof%m_psi, dof%n_psi, dof%m_K, dof%n_K)

    call set(dof_rsc, dof)

    call scale(dof_rsc%psi, 1.0_dp/sc)
    call scale(dof_rsc%K11, 1.0_dp/sc)
    call scale(dof_rsc%K22, 1.0_dp/sc)
    call scale(dof_rsc%K12, 1.0_dp/sc)

    if (present(index)) then
      call write(dof_rsc%K11, trim(filename)//"K11", SlogPost, index)
      call write(dof_rsc%K22, trim(filename)//"K22", SlogPost, index)
      call write(dof_rsc%K12, trim(filename)//"K12", SlogPost, index)

      call write(dof_rsc%psi, trim(filename)//"psi", SlogPost, index)
    else
      call write(dof_rsc%K11, trim(filename)//"K11", SlogPost)
      call write(dof_rsc%K22, trim(filename)//"K22", SlogPost)
      call write(dof_rsc%K12, trim(filename)//"K12", SlogPost)

      call write(dof_rsc%psi, trim(filename)//"psi", SlogPost)
    end if

    call deallocate(dof_rsc)

  end subroutine write_dof
  
  
  subroutine read_dof(dof, filename, sc, index)
    type(dofdat), intent(inout) :: dof
    character(len = *), intent(in) :: filename
    real(kind = dp), intent(in) :: sc
    integer, optional, intent(in) :: index
    
    real(kind = dp) :: SlogPost

    if (present(index)) then
      call read_header(dof%K11, trim(filename)//"K11", SlogPost, index)
      call read_header(dof%K22, trim(filename)//"K22", SlogPost, index)
      call read_header(dof%K12, trim(filename)//"K12", SlogPost, index)
      call read_header(dof%psi, trim(filename)//"psi", SlogPost, index)
      
      call read(dof%K11, trim(filename)//"K11", index)
      call read(dof%K22, trim(filename)//"K22", index)
      call read(dof%K12, trim(filename)//"K12", index)
      call read(dof%psi, trim(filename)//"psi", index)
    else
      call read_header(dof%K11, trim(filename)//"K11", SlogPost)
      call read_header(dof%K22, trim(filename)//"K22", SlogPost)
      call read_header(dof%K12, trim(filename)//"K12", SlogPost)
      call read_header(dof%psi, trim(filename)//"psi", SlogPost)

      call read(dof%K11, trim(filename)//"K11")
      call read(dof%K22, trim(filename)//"K22")
      call read(dof%K12, trim(filename)//"K12")
      call read(dof%psi, trim(filename)//"psi")
    end if
    
    ! TODO: Ad-hoc 
    dof%m_psi = dof%psi%m
    dof%n_psi = dof%psi%n
    dof%m_K = dof%K11%m
    dof%n_K = dof%K11%n
    dof%ndof = dof%m_psi*dof%n_psi + 3*dof%m_K*dof%n_K
    dof%SlogPost = SlogPost
    dof%occurrence = 1

    ! Rescale to computation domain
    call scale(dof%psi, sc)
    call scale(dof%K11, sc)
    call scale(dof%K22, sc)
    call scale(dof%K12, sc)
    
  end subroutine read_dof
  
  subroutine write_SSD(canon_SSD, filename)
    real(kind=dp), dimension(:), intent(in) :: canon_SSD
    character(len = *), intent(in) :: filename
        
    open(unit = output_unit, file = trim(filename) // ".dat", &
      & status = "replace", access = "stream", action = "write")
    write(output_unit) canon_SSD
    close(output_unit)

  end subroutine write_SSD
  
  subroutine read_SSD(canon_SSD, filename)
    real(kind=dp), dimension(:), intent(inout) :: canon_SSD
    character(len = *), intent(in) :: filename
        
    open(unit = input_unit, file = trim(filename) // ".dat", &
      & status = "old", access = "stream", action = "read")
    read(input_unit) canon_SSD
    close(input_unit)

  end subroutine read_SSD
  
  subroutine write_theta(dof, filename, sc, index)
    type(dofdat), intent(in) :: dof
    character(len = *), intent(in) :: filename
    real(kind=dp), optional, intent(in) :: sc
    integer, optional, intent(in) :: index
    
    character(len = len_trim(filename) + 1 + max_int_len) :: lfilename
    real(kind=dp), dimension(:), allocatable :: theta
    real(kind=dp) :: rsc
        
    if (present(sc)) then
      rsc = 1.0_dp/sc
    else
      rsc = 1.0_dp
    end if

    if(present(index)) then
      write(lfilename, "(a,a,i0)") trim(filename), "_", index
    else
      lfilename = trim(filename)
    end if

    open(unit = output_unit, file = trim(lfilename) // ".hdr", &
      & status = "replace", action = "write")
    write(output_unit, "(a)") "theta -- psi, K11, K22, K12 (m, n, type_id: 0center, 1corner) logPost, occurrence"
    write(output_unit, "(i0)") 4
    write(output_unit, "(i0,a,i0,a,i0)") size(dof%psi%data, 1), " ", &
                      size(dof%psi%data, 2), " ", dof%psi%type_id
    write(output_unit, "(i0,a,i0,a,i0)") size(dof%K11%data, 1), " ", &
                      size(dof%K11%data, 2), " ", dof%K11%type_id
    write(output_unit, "(i0,a,i0,a,i0)") size(dof%K22%data, 1), " ", &
                      size(dof%K22%data, 2), " ", dof%K22%type_id
    write(output_unit, "(i0,a,i0,a,i0)") size(dof%K12%data, 1), " ", &
                      size(dof%K12%data, 2), " ", dof%K12%type_id
    write(output_unit, "("//dp_chr//")") dof%SlogPost
    write(output_unit, "(i0)") dof%occurrence
    close(output_unit)
    
    call allocate(theta, dof)
    call convert_dof_to_theta(theta, dof, rsc)

    open(unit = output_unit, file = trim(lfilename) // ".dat", &
      & status = "replace", access = "stream", action = "write")
    write(output_unit) theta
    close(output_unit)
    
    deallocate(theta)
    
  end subroutine write_theta
  
  subroutine read_theta_header(dof, filename, index)
    type(dofdat), intent(out) :: dof
    character(len = *), intent(in) :: filename
    integer, optional, intent(in) :: index

    character(len = max_line_len) :: input_chr
    character(len = len_trim(filename) + 1 + max_int_len) :: lfilename
    integer :: m_psi, n_psi, type_id
    integer :: m_K, n_K
    real(kind=dp) :: SlogPost
    integer :: occurrence

    if(present(index)) then
      write(lfilename, "(a,a,i0)") trim(filename), "_", index
    else
      lfilename = trim(filename)
    end if
    
    open(unit = input_unit, file = trim(lfilename) // ".hdr", &
      & status = "old", action = "read")
    read(input_unit, "(a)") input_chr
    read(input_unit, "(a)") input_chr
    if(input_chr .ne. "4") then
      call abort_handle("Invalid header", __FILE__, __LINE__)
    end if
    read(input_unit, "(a)") input_chr
    read(input_chr, *) m_psi, n_psi, type_id
    read(input_unit, "(a)") input_chr
    read(input_chr, *) m_K, n_K, type_id
    read(input_unit, "(a)") input_chr
    read(input_chr, *) m_K, n_K, type_id
    read(input_unit, "(a)") input_chr
    read(input_chr, *) m_K, n_K, type_id
    read(input_unit, "("//dp_chr//")") SlogPost
    read(input_unit, "("//int_chr//")") occurrence
    close(input_unit)
    
    call allocate(dof, m_psi, n_psi, m_K, n_K)
    dof%ndof = m_psi*n_psi + 3*m_K*n_K
    dof%SlogPost = SlogPost
    dof%occurrence = occurrence
    
  end subroutine read_theta_header
  
  subroutine read_theta(dof, filename, sc, index)
    type(dofdat), intent(out) :: dof
    character(len = *), intent(in) :: filename
    real(kind=dp), optional, intent(in) :: sc
    integer, optional, intent(in) :: index

    character(len = len_trim(filename) + 1 + max_int_len) :: lfilename
    real(kind=dp), dimension(:), allocatable :: theta
    real(kind=dp) :: rsc

    if(present(index)) then
      write(lfilename, "(a,a,i0)") trim(filename), "_", index
    else
      lfilename = trim(filename)
    end if

    call read_theta_header(dof, filename, index)
    
    call allocate(theta, dof)
    theta = 0.0_dp
    
    write(6, *) "Reading: "// trim(lfilename) // ".dat"
    
    open(unit = input_unit, file = trim(lfilename) // ".dat", &
      & status = "old", access = "stream", action = "read")
    read(input_unit) theta
    close(input_unit)
    
    
    if(present(sc)) then
      rsc = sc
    else
      rsc = 1.0_dp
    end if
    
    call convert_theta_to_dof(dof, theta, rsc)
    call imprint_canon(dof)
    
    deallocate(theta)
    
  end subroutine read_theta
  
  subroutine read_QGfield_header(fld, t, nlayer, filename)
    type(field), intent(out) :: fld
    real(kind = dp), intent(out) :: t
    character(len = *), intent(in) :: filename
    integer, intent(out) :: nlayer

    character(len = 1) :: chr
    character(len = field_name_len) :: name
    character(len = max_line_len) :: format_chr, input_chr
    character(len = len_trim(filename) + 1 + max_int_len) :: lfilename
    integer :: i, m, n

    lfilename = trim(filename)

    open(unit = input_unit, file = trim(lfilename) // ".hdr", &
      & status = "old", action = "read")
    read(input_unit, "(a)") input_chr
    if(trim(input_chr) /= "serial") then
      call abort_handle("Invalid header", __FILE__, __LINE__)
    end if
    read(input_unit, "(a)") name
    read(input_unit, "(a)") input_chr
    i = scan(input_chr, " ")
    write(format_chr, "(a,i0,a)") "(i", i - 1, ",a1,"//int_chr//")"
    read(input_chr, trim(format_chr)) m, chr, n
    read(input_unit, "("//int_chr//")") nlayer
    if (nlayer .ne. 3) then
      call abort_handle("Invalid header: nlayer =/= 3", __FILE__, __LINE__)
    end if
    read(input_unit, "("//dp_chr//")") t
    close(input_unit)

    call allocate(fld, m-1, n-1, name, type_id=1)

  end subroutine read_QGfield_header

  subroutine read_QGfield(fld, t, filename, layer)
    type(field), intent(inout) :: fld
    real(kind=dp), intent(out) :: t
    character(len = *), intent(in) :: filename
    integer, intent(in) :: layer
    
    character(len = len_trim(filename) + 1 + max_int_len) :: lfilename
    
    real(kind=dp), dimension(:,:,:), allocatable :: tmp_fld
    integer :: nlayer

    lfilename = trim(filename)
    
    call read_QGfield_header(fld, t, nlayer, filename)
    
    allocate(tmp_fld(size(fld%data, 1), size(fld%data, 2), nlayer))
    
    open(unit = input_unit, file = trim(lfilename) // ".dat", &
      & status = "old", access = "stream", action = "read")
    read(input_unit) tmp_fld
    close(input_unit)
    
    fld%data = tmp_fld(:, :, layer)
    
    deallocate(tmp_fld)

  end subroutine read_QGfield
 

  subroutine reset_io_timers()
    call reset(write_field_timer)
  end subroutine reset_io_timers

  subroutine print_io_timers()
    write(6, "(a)") "I/O timers:"
    call print(write_field_timer, "write_field", prefix = "  ")
  end subroutine print_io_timers

end module advdiff_io
