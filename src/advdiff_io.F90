module advdiff_io

  use advdiff_precision
  use advdiff_debug
  use advdiff_timing
  use advdiff_parameters, only : input_unit, output_unit
  use advdiff_field
  use advdiff_trajdata
  use advdiff_inference

  implicit none

  private

  public :: read_header, read, write

  interface read_header
    module procedure read_field_header, read_traj_header, read_dof
  end interface read_header

  interface read
    module procedure read_field, read_traj, read_dof
  end interface read

  ! Interface: allowing overloading the 'write' function by the parameters input
  interface write
    module procedure write_field, write_traj, write_dof, write_array
  end interface write

  type(timer), save :: write_field_timer
  public :: reset_io_timers, print_io_timers

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
    integer :: i, m, n, glayer, type_id

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
    read(input_unit, "("//int_chr//")") glayer
    read(input_unit, "("//int_chr//")") type_id
    read(input_unit, "("//dp_chr//")") t
    close(input_unit)

    call allocate(fld, m, n, name, glayer=glayer, type_id=type_id)

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
    write(output_unit, "(i0)") fld%glayer
    write(output_unit, "(i0)") fld%type_id
    write(output_unit, "("//dp_chr//")") t
    close(output_unit)

    open(unit = output_unit, file = trim(lfilename) // ".dat", &
      & status = "replace", access = "stream", action = "write")
    write(output_unit) fld%data
    close(output_unit)

    call stop(write_field_timer)

  end subroutine write_field

  subroutine write_array(array, filename)
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

  end subroutine write_array

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

    if(present(index)) then
      write(lfilename, "(a,a,i0)") trim(filename), "_", index
    else
      lfilename = trim(filename)
    end if

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

  subroutine write_traj(traj, filename, t, index)
    type(trajdat), intent(in) :: traj
    character(len = *), intent(in) :: filename
    real(kind = dp), intent(in) :: t
    integer, optional, intent(in) :: index

    character(len = len_trim(filename) + 1 + max_int_len) :: lfilename

    if(present(index)) then
      write(lfilename, "(a,a,i0)") trim(filename), "_", index
    else
      lfilename = trim(filename)
    end if

    open(unit = output_unit, file = trim(lfilename) // ".hdr", &
      & status = "replace", action = "write")
    write(output_unit, "(a)") "serial"
    write(output_unit, "(a)") "traj"
    write(output_unit, "(i0,a,i0)") traj%nparticles, " ", traj%nts0
    write(output_unit, "(i0)") 2
    write(output_unit, "("//dp_chr//")") t
    close(output_unit)

    open(unit = output_unit, file = trim(lfilename) // "_x.dat", &
      & status = "replace", access = "stream", action = "write")
    write(output_unit) traj%x
    close(output_unit)

    open(unit = output_unit, file = trim(lfilename) // "_y.dat", &
      & status = "replace", access = "stream", action = "write")
    write(output_unit) traj%y
    close(output_unit)

    open(unit = output_unit, file = trim(lfilename) // "_t.dat", &
      & status = "replace", access = "stream", action = "write")
    write(output_unit) traj%t
    close(output_unit)

  end subroutine write_traj

  subroutine write_dof(dof, filename, SlogPost, sc, index)
    type(dofdat), intent(in) :: dof
    character(len = *), intent(in) :: filename
    real(kind = dp), intent(in) :: SlogPost, sc
    integer, optional, intent(in) :: index


    type(dofdat) :: dof_rsc

    call allocate(dof_rsc, dof%m, dof%n)

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
    dof%m = dof%psi%m
    dof%n = dof%psi%n
    dof%ndof = dof%m * dof%n

    ! Rescale to computation domain
    call scale(dof%psi, sc)
    call scale(dof%K11, sc)
    call scale(dof%K22, sc)
    call scale(dof%K12, sc)
    
  end subroutine read_dof

  subroutine reset_io_timers()
    call reset(write_field_timer)
  end subroutine reset_io_timers

  subroutine print_io_timers()
    write(6, "(a)") "I/O timers:"
    call print(write_field_timer, "write_field", prefix = "  ")
  end subroutine print_io_timers

end module advdiff_io