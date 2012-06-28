module GOLD_file_parser
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of GOLD.                                        *
!*                                                                     *
!* GOLD is free software; you can redistribute it and/or modify it and *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* GOLD is distributed in the hope that it will be useful, but WITHOUT *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, June 2005.                                     *
!*                                                                     *
!*    The subroutines here parse a set of input files for the value    *
!*  a named parameter and sets that parameter at run time.  Currently  *
!*  these files use the same format as the header file GOLD_memory.h.  *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_coms, only : root_PE, broadcast
use GOLD_error_handler, only : GOLD_error, FATAL, WARNING, GOLD_mesg
use GOLD_error_handler, only : is_root_pe, stdlog, stdout
use GOLD_time_manager, only : set_time, get_time, time_type, get_ticks_per_second
use GOLD_document, only : doc_param, doc_module, doc_init, doc_end

implicit none ; private

integer, parameter, public :: MAX_PARAM_FILES = 5 ! Maximum number of parameter files.
integer, parameter :: INPUT_STR_LENGTH = 120 ! Maximum linelength in parameter file.

! The all_PEs_read option should be eliminated with post-riga shared code.
logical :: all_PEs_read = .true.

type, private :: file_data_type ; private
  integer :: num_lines = 0
  character(len=INPUT_STR_LENGTH), pointer, dimension(:) :: line => NULL()
  logical,                         pointer, dimension(:) :: line_used => NULL()
end type file_data_type

type, public :: param_file_type ; private
  integer  :: nfiles = 0            ! The number of open files.
  integer  :: iounit(MAX_PARAM_FILES)   ! The unit number of an open file.
  character(len=200) :: filename(MAX_PARAM_FILES) ! The names of the open files.
  logical  :: NetCDF_file(MAX_PARAM_FILES)! If true, the input file is in NetCDF.
                                    ! This is not yet implemented.
  type(file_data_type) :: param_data(MAX_PARAM_FILES) ! Structures that contain 
                                    ! the valid data lines from the parameter
                                    ! files, enabling all subsequent reads of
                                    ! parameter data to occur internally.
  logical  :: report_unused = .false. ! If true, report any parameter
                                    ! lines that are not used in the run.
  logical  :: unused_params_fatal = .false.  ! If true, kill the run if there
                                    ! are any unused parameters.
  logical  :: log_to_stdout = .false. ! If true, all log messages are also
                                    ! sent to stdout.
  logical  :: log_open = .false.    ! True if the log file has been opened.
  integer  :: stdout, stdlog        ! The units from stdout() and stdlog().
  character(len=240) :: doc_file    ! A file where all run-time parameters, their
                                    ! settings and defaults are documented.
  logical  :: minimal_doc           ! If true, document only those run-time
                                    ! parameters that differ from defaults.
end type param_file_type

public read_param, open_param_file, close_param_file, log_param, log_version
public doc_param
public lowercase, uppercase

interface read_param
  module procedure read_param_int, read_param_real, read_param_logical, &
                   read_param_char, read_param_char_array, read_param_time, &
                   read_param_int_array, read_param_real_array
end interface
interface log_param
  module procedure log_param_int, log_param_real, log_param_logical, &
                   log_param_char, log_param_time, &
                   log_param_int_array, log_param_real_array
end interface

contains

subroutine open_param_file(filename, param_file, checkable)
  character(len=*),      intent(in)    :: filename
  type(param_file_type), intent(inout) :: param_file
  logical,     optional, intent(in)    :: checkable
  
  logical :: file_exists, unit_in_use, Netcdf_file, may_check
  integer :: ios, iounit, strlen, i

  may_check = .true. ; if (present(checkable)) may_check = checkable

  ! Check for non-blank filename
  strlen = len_trim(filename)
  if (strlen == 0) then
    call GOLD_error(FATAL, "open_param_file: Input file has not been specified.")
  endif

  ! Check that this file has not already been opened
  if (param_file%nfiles > 0) then
    inquire(file=trim(filename), number=iounit)
    if (iounit /= -1) then
      do i = 1, param_file%nfiles
        if (param_file%iounit(i) == iounit) then
          if (trim(param_file%filename(1)) /= trim(filename)) then
            call GOLD_error(FATAL, &
              "open_param_file: internal inconsistency! "//trim(filename)// &
              " is registered as open but has the wrong unit number!")
          else
            call GOLD_error(WARNING, &
              "open_param_file: file "//trim(filename)// &
              " has already been opened. This should NOT happen!"// &
              " Did you specify the same file twice in a namelist?")
            return
          endif ! filenames
        endif ! unit numbers
      enddo ! i
    endif
  endif

  ! Check that the file exists to readstdlog
  inquire(file=trim(filename), exist=file_exists)
  if (.not.file_exists) call GOLD_error(FATAL, &
      "open_param_file: Input file "// trim(filename)//" does not exist.")

  Netcdf_file = .false.
  if (strlen > 3) then
    if (filename(strlen-2:strlen) == ".nc") Netcdf_file = .true.
  endif

  if (Netcdf_file) &
    call GOLD_error(FATAL,"open_param_file: NetCDF files are not yet supported.")

  if (all_PEs_read .or. is_root_pe()) then
    ! Find an unused unit number.
    do iounit=10,512
      INQUIRE(iounit,OPENED=unit_in_use) ; if (.not.unit_in_use) exit
    enddo
    if (iounit >= 512) call GOLD_error(FATAL, &
        "open_param_file: No unused file unit could be found.")

    ! Open the parameter file.
    open(iounit, file=trim(filename), access='SEQUENTIAL', &
         form='FORMATTED', action='READ', position='REWIND', iostat=ios)
    if (ios /= 0) call GOLD_error(FATAL, "open_param_file: Error opening "// &
                                       trim(filename))
  else
    iounit = 1
  endif

  ! Store/register the unit and details
  i = param_file%nfiles + 1
  param_file%nfiles = i
  param_file%iounit(i) = iounit
  param_file%filename(i) = filename
  param_file%NetCDF_file(i) = Netcdf_file

  call GOLD_mesg("open_param_file: "// trim(filename)// &
                 " has been opened successfully.", 5)

  call populate_param_data(iounit, filename, param_file%param_data(i))

  call read_param(param_file,"SEND_LOG_TO_STDOUT",param_file%log_to_stdout)
  call read_param(param_file,"REPORT_UNUSED_PARAMS",param_file%report_unused)
  call read_param(param_file,"FATAL_UNUSED_PARAMS",param_file%unused_params_fatal)
  param_file%doc_file = " "
  call read_param(param_file,"DOCUMENT_FILE", param_file%doc_file)
  if (.not.may_check) then
    param_file%report_unused = .false.
    param_file%unused_params_fatal = .false.
  endif

  ! Open the log file.
  param_file%stdlog = stdlog() ; param_file%stdout = stdout()
  param_file%log_open = (stdlog() > 0)

  if (len_trim(param_file%doc_file) > 0) then
    param_file%minimal_doc = .false.
    call read_param(param_file, "MINIMAL_DOCUMENTATION", param_file%minimal_doc)
    call doc_init(param_file%doc_file, param_file%minimal_doc)
  endif

end subroutine open_param_file

subroutine close_param_file(param_file)
  type(param_file_type), intent(inout) :: param_file
  character(len=128) :: version = '$Id: GOLD_file_parser.F90,v 13.0.2.3.2.21 2011/09/19 16:10:35 Robert.Hallberg Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod   ! This module's name.
  integer :: i, n, num_unused

  ! Log the parameters for the parser.
  mod = "GOLD_file_parser"
  call log_version(param_file, mod, version, tagname)
  call log_param(param_file, mod, "SEND_LOG_TO_STDOUT", &
                        param_file%log_to_stdout, &
                 "If true, all log messages are also sent to stdout.", &
                 default=.false.)
  call log_param(param_file, mod, "REPORT_UNUSED_PARAMS", &
                        param_file%report_unused, &
                 "If true, report any parameter lines that are not used \n"//&
                 "in the run.", default=.false.)
  call log_param(param_file, mod, "FATAL_UNUSED_PARAMS", &
                        param_file%unused_params_fatal, &
                 "If true, kill the run if there are any unused \n"//&
                 "parameters.", default=.false.)
  call log_param(param_file, mod, "DOCUMENT_FILE", param_file%doc_file, &
                 "A file where all run-time parameters, their settings \n"//&
                 "and defaults are documented.", default=" ")
  if (len_trim(param_file%doc_file) > 0) &
    call log_param(param_file, mod, "MINIMAL_DOCUMENTATION", &
                   param_file%minimal_doc, &
                  "If true, document only those run-time parameters that \n"//&
                  "differ from their defaults.", default=.false.)

  num_unused = 0
  do i = 1, param_file%nfiles
    ! only root pe has the file open
    if (all_PEs_read .or. is_root_pe()) close(param_file%iounit(i))
    call GOLD_mesg("close_param_file: "// trim(param_file%filename(i))// &
                 " has been closed successfully.", 5)
    
    ! Check for unused lines.
    if (is_root_pe() .and. (param_file%report_unused .or. &
                            param_file%unused_params_fatal)) then
      do n=1,param_file%param_data(i)%num_lines
        if (.not.param_file%param_data(i)%line_used(n)) then
          num_unused = num_unused + 1
!          call GOLD_error(WARNING, "Unused line in "//trim(param_file%filename(i))//&
!                          " : "//trim(param_file%param_data(i)%line(n)))
          if (param_file%report_unused) &
            call GOLD_mesg("Unused line in "//trim(param_file%filename(i))// &
                            " : "//trim(param_file%param_data(i)%line(n)), 0)
        endif
      enddo
    endif

    param_file%iounit(i) = -1
    param_file%filename(i) = ""
    param_file%NetCDF_file(i) = .false.
    deallocate (param_file%param_data(i)%line)
    deallocate (param_file%param_data(i)%line_used)
  enddo

  if (is_root_pe() .and. (num_unused>0) .and. param_file%unused_params_fatal) &
    call GOLD_error(FATAL, "Run stopped because of unused parameter lines.")

  param_file%log_open = .false.
  call doc_end()

end subroutine close_param_file

subroutine populate_param_data(iounit, filename, param_data)
  integer,              intent(in) :: iounit
  character(len=*),     intent(in) :: filename
  type(file_data_type), intent(inout) :: param_data

  character(len=INPUT_STR_LENGTH) :: line
  integer  :: is, isd, isu, icom, verbose
  integer  :: last, num_lines
  logical  :: found_override, found_define, found_undef
  verbose = 1

  ! find the number of keyword lines in a parameter file
  ! allocate the space to hold the lines in param_data%line
  ! populate param_data%line with the keyword lines from parameter file

  if (iounit <= 0) return

  if (all_PEs_read .or. is_root_pe()) then
    ! rewind the parameter file
    rewind(iounit)

    ! count the number of valid entries in the parameter file
    num_lines = 0
    do while(.true.)
      read(iounit, '(a)', end=8, err=9) line
      last = len_trim(line)
      ! Eliminate either F90 or C comments from the line.
      icom = index(line(:last), "!") ; if (icom > 0) last = icom-1
      icom = index(line(:last), "/*") ; if (icom > 0) last = icom-1

      if (last < 1) cycle

      ! Detect keywords
      found_override = .false.; found_define = .false.; found_undef = .false.
      is = index(line(:last), "override" ); if (is > 0) found_override = .true.
      isd = index(line(:last), "define" ); if (isd > 0) found_define = .true.
      isu = index(line(:last), "undef" ); if (isu > 0) found_undef = .true.

      if (found_override .or. found_define .or. found_undef) num_lines = num_lines + 1 
    enddo ! while (.true.)

 8  continue

    ! allocate space to hold contents of the parameter file
    param_data%num_lines = num_lines
  endif  ! (is_root_pe())

! broadcast the number of valid entries in parameter file
  if (.not. all_PEs_read) then
    call GOLD_error(FATAL, "Code to read the param_file from only the root "//&
             "PE will not be enabled until the post-riga shared code is used.")
!    call broadcast(param_data%num_lines, root_pe())
  endif

! Set up the space for storing the actual lines.
  num_lines = param_data%num_lines
  allocate (param_data%line(num_lines))
  allocate (param_data%line_used(num_lines))
  param_data%line(:) = ' '
  param_data%line_used(:) = .false.

  if (all_PEs_read .or. is_root_pe()) then
    ! Read the actual lines.
    ! rewind the parameter file
    rewind(iounit)

    ! populate param_data%line
    num_lines = 0
    do while(.true.)
      read(iounit, '(a)', end=18, err=9) line
      last = len_trim(line)
      ! Eliminate either F90 or C comments from the line.
      icom = index(line(:last), "!") ; if (icom > 0) last = icom-1
      icom = index(line(:last), "/*") ; if (icom > 0) last = icom-1

      if (last < 1) cycle

      ! Detect keywords
      found_override = .false.; found_define = .false.; found_undef = .false.
      is = index(line(:last), "override" ); if (is > 0) found_override = .true.
      isd = index(line(:last), "define" ); if (isd > 0) found_define = .true.
      isu = index(line(:last), "undef" ); if (isu > 0) found_undef = .true.

      if (found_override .or. found_define .or. found_undef) then
        num_lines = num_lines + 1 
        param_data%line(num_lines) = line
      endif
    enddo ! while (.true.)

18  continue
    if (num_lines /= param_data%num_lines) &
      call GOLD_error(FATAL, 'GOLD_file_parser : Found different number of '// &
                      'valid lines on second reading of '//trim(filename))
  endif  ! (is_root_pe())

  ! broadcast the populated array param_data%line
  if (.not. all_PEs_read) then
    call GOLD_error(FATAL, "Code to read the param_file from only the root "//&
             "PE will not be enabled until the post-riga shared code is used.")
!   call broadcast(param_data%line, INPUT_STR_LENGTH, root_pe())
  endif

  return

9 call GOLD_error(FATAL, "GOLD_file_parser : "//&
                  "Error while reading file "//trim(filename))

end subroutine populate_param_data

subroutine read_param_int(param_file, varname, value, fail_if_missing)
  type(param_file_type), intent(in) :: param_file
  character(len=*), intent(in)      :: varname
  integer, intent(inout)            :: value
  logical, optional, intent(in)     :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical            :: found, defined

  call get_variable_line(param_file, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    read(value_string(1),*) value
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call GOLD_error(FATAL,'Unable to find variable '//trim(varname)// &
                             ' in any input files.')
      else
        call GOLD_error(FATAL,'Variable '//trim(varname)// &
                             ' found but not set in input files.')
      endif
    endif ; endif
  endif
end subroutine read_param_int

subroutine read_param_int_array(param_file, varname, value, fail_if_missing)
  type(param_file_type), intent(in) :: param_file
  character(len=*), intent(in)      :: varname
  integer, intent(inout)            :: value(:)
  logical, optional, intent(in)     :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical            :: found, defined

  call get_variable_line(param_file, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    read(value_string(1),*,end=991) value
 991 return
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call GOLD_error(FATAL,'Unable to find variable '//trim(varname)// &
                             ' in any input files.')
      else
        call GOLD_error(FATAL,'Variable '//trim(varname)// &
                             ' found but not set in input files.')
      endif
    endif ; endif
  endif
end subroutine read_param_int_array

subroutine read_param_real(param_file, varname, value, fail_if_missing)
  type(param_file_type), intent(in) :: param_file
  character(len=*), intent(in)      :: varname
  real, intent(inout)               :: value
  logical, optional, intent(in)     :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical            :: found, defined

  call get_variable_line(param_file, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    read(value_string(1),*) value
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call GOLD_error(FATAL,'Unable to find variable '//trim(varname)// &
                             ' in any input files.')
      else
        call GOLD_error(FATAL,'Variable '//trim(varname)// &
                             ' found but not set in input files.')
      endif
    endif ; endif
  endif
end subroutine read_param_real

subroutine read_param_real_array(param_file, varname, value, fail_if_missing)
  type(param_file_type), intent(in) :: param_file
  character(len=*), intent(in)      :: varname
  real, intent(inout)               :: value(:)
  logical, optional, intent(in)     :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.
  character(len=120) :: value_string(1)
  logical            :: found, defined

  call get_variable_line(param_file, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    read(value_string(1),*,end=991) value
 991 return
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call GOLD_error(FATAL,'Unable to find variable '//trim(varname)// &
                             ' in any input files.')
      else
        call GOLD_error(FATAL,'Variable '//trim(varname)// &
                             ' found but not set in input files.')
      endif
    endif ; endif
  endif
end subroutine read_param_real_array

subroutine read_param_char(param_file, varname, value, fail_if_missing)
  type(param_file_type), intent(in) :: param_file
  character(len=*), intent(in)      :: varname
  character(len=*), intent(inout)   :: value
  logical, optional, intent(in)     :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical            :: found, defined

  call get_variable_line(param_file, varname, found, defined, value_string)
  if (found) then
    value = trim(value_string(1))
  elseif (present(fail_if_missing)) then ; if (fail_if_missing) then
    call GOLD_error(FATAL,'Unable to find variable '//trim(varname)// &
                         ' in any input files.')
  endif ; endif

end subroutine read_param_char

subroutine read_param_char_array(param_file, varname, value, fail_if_missing)
  type(param_file_type), intent(in) :: param_file
  character(len=*), intent(in)      :: varname
  character(len=*), intent(inout)   :: value(:)
  logical, optional, intent(in)     :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.
  character(len=INPUT_STR_LENGTH) :: value_string(SIZE(value))
  logical            :: found, defined
  integer            :: i, i_out

  call get_variable_line(param_file, varname, found, defined, value_string)
  if (found) then
    i_out = 1
    do i=1,SIZE(value) ; if (len_trim(value_string(i)) > 0) then
      value(i_out) = trim(value_string(i)) ; i_out = i_out + 1
    endif ; enddo
    do i=i_out,SIZE(value) ; value(i) = " " ; enddo
  elseif (present(fail_if_missing)) then ; if (fail_if_missing) then
    call GOLD_error(FATAL,'Unable to find variable '//trim(varname)// &
                         ' in any input files.')
  endif ; endif

end subroutine read_param_char_array

subroutine read_param_logical(param_file, varname, value, fail_if_missing)
  type(param_file_type), intent(in) :: param_file
  character(len=*), intent(in)      :: varname
  logical, intent(inout)            :: value
  logical, optional, intent(in)     :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical            :: found, defined

  call get_variable_line(param_file, varname, found, defined, value_string)
  if (found) then
    value = defined
  elseif (present(fail_if_missing)) then ; if (fail_if_missing) then
    call GOLD_error(FATAL,'Unable to find variable '//trim(varname)// &
                         ' in any input files.')
  endif ; endif
end subroutine read_param_logical


subroutine read_param_time(param_file, varname, value, timeunit, fail_if_missing)
  type(param_file_type), intent(in) :: param_file
  character(len=*), intent(in)      :: varname
  type(time_type), intent(inout)    :: value
  real, intent(in)                  :: timeunit
  logical, optional, intent(in)     :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.  The unique argument
! to read time is the number of seconds to use as the unit of time being read.
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical            :: found, defined
  real               :: real_time
  integer            :: days, secs

  call get_variable_line(param_file, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    read( value_string(1), *) real_time
    days = int(real_time*(timeunit/86400.0))
    secs = int(floor((real_time*(timeunit/86400.0)-days)*86400.0 + 0.5))
    value = set_time(secs, days)
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call GOLD_error(FATAL,'Unable to find variable '//trim(varname)// &
                             ' in any input files.')
      else
        call GOLD_error(FATAL,'Variable '//trim(varname)// &
                             ' found but not set in input files.')
      endif
    endif ; endif
  endif
end subroutine read_param_time


subroutine get_variable_line(param_file, varname, found, defined, value_string)
  type(param_file_type), intent(in)  :: param_file
  character(len=*),      intent(in)  :: varname
  logical,               intent(out) :: found, defined
  character(len=*),      intent(out) :: value_string(:)

  character(len=INPUT_STR_LENGTH) :: line, val_str, lname, filename
  integer            :: is, id, isd, isu, icom, verbose, ipf
  integer            :: last, ival, oval, max_vals, count
  character(len=52)  :: set
  logical            :: found_override, found_define, found_undef
  logical            :: force_cycle, defined_in_line
  set = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
  verbose = 1

  ! Find the first instance (if any) where the named variable is found, and
  ! return variables indicating whether this variable is defined and the string
  ! that contains the value of this variable.
  found = .false.
  oval = 0; ival = 0;
  max_vals = SIZE(value_string)
  do is=1,max_vals ; value_string(is) = " " ; enddo

  paramfile_loop: do ipf = 1, param_file%nfiles
    filename = param_file%filename(ipf)

    do count = 1, param_file%param_data(ipf)%num_lines
      line = param_file%param_data(ipf)%line(count)
      last = len_trim(line)
      ! Eliminate either F90 or C comments from the line.
      icom = index(line(:last), "!") ; if (icom > 0) last = icom-1
      icom = index(line(:last), "/*") ; if (icom > 0) last = icom-1

      if (last < 1) cycle
      if (index(line(:last), trim(varname)) == 0) cycle

      ! Detect keywords
      found_override = .false.; found_define = .false.; found_undef = .false.
      is = index(line(:last), "override" ); if (is > 0) found_override = .true.
      isd = index(line(:last), "define" ); if (isd > 0) found_define = .true.
      isu = index(line(:last), "undef" ); if (isu > 0) found_undef = .true.

      ! Check for missing, mutually exclusive or incomplete keywords
      if (is_root_pe()) then
        if (.not. (found_define .or. found_undef)) &
               call GOLD_error(FATAL, "GOLD_file_parser : the parameter name '"// &
                 trim(varname)//"' was found without define or undef."// &
                 " Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
        if (found_define .and. found_undef) call GOLD_error(FATAL, &
                 "GOLD_file_parser : Both 'undef' and 'define' occur."// &
                 " Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
        if (found_override .and. .not. (found_define .or. found_undef)) &
               call GOLD_error(FATAL, "GOLD_file_parser : override was found "// &
                 " without a define or undef."// &
                 " Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
        if (found_override .and. (found_define .or. found_undef) &
            .and. (isd+isu<is)) &
               call GOLD_error(FATAL, "GOLD_file_parser : override was found "// &
                 " but was not the first keyword."// &
                 " Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
      endif

      ! Interpret the line and collect values, if any
      if (found_define) then
        ! Move starting pointer to first letter of defined name.
        is = isd + 5 + scan(line(isd+6:last), set)

        id = scan(line(is:last), ' ')  ! Find space between name and value
        if ( id == 0 ) then
          ! There is no space so the name is simply being defined.
          lname = trim(line(is:last))
          if (trim(lname) /= trim(varname)) cycle
          val_str = " "
        else
          ! There is a string or number after the name.
          lname = trim(line(is:is+id-1))
          if (trim(lname) /= trim(varname)) cycle
          val_str = trim(adjustl(line(is+id:last)))
          ! Remove starting and trailing quotes.
          id = index(val_str,ACHAR(34)) ; if (id > 0) val_str = val_str(id+1:)
          id = index(val_str,ACHAR(34)) ; if (id > 0) val_str = val_str(:id-1)
        endif
        found = .true. ; defined_in_line = .true.
      elseif (found_undef) then
        ! Move starting pointer to first letter of undefined name.
        is = isu + 4 + scan(line(isu+5:last), set)

        id = scan(line(is:last), ' ')  ! Find the first space after the name.
        if (id > 0) last = is + id - 1
        lname = trim(line(is:last))
        if (trim(lname) /= trim(varname)) cycle
        val_str = " "
        found = .true. ; defined_in_line = .false.
      else
        call GOLD_error(FATAL, "GOLD_file_parser: we should never reach this point")
      endif

      ! This line has now been used.
!     param_file%param_data(ipf)%line_used(count) = .true.
      ! Pathscale barfs on the above line
      call flag_line_as_read(param_file%param_data(ipf)%line_used,count)

      ! Detect inconsistencies
      force_cycle = .false.
      if (found_override .and. (oval >= max_vals)) then
        if (is_root_pe()) then
          if ((defined_in_line /= defined) .or. &
              (trim(val_str) /= trim(value_string(max_vals)))) then
            call GOLD_error(FATAL,"GOLD_file_parser : "//trim(varname)// &
                     " found with multiple inconsistent overrides."// &
                     " Line: '"//trim(line(:last))//"'"//&
                     " in file "//trim(filename)//" caused the model failure.")
          else
            call GOLD_error(WARNING,"GOLD_file_parser : "//trim(varname)// &
                     " over-ridden more times than is permitted."// &
                     " Line: '"//trim(line(:last))//"'"//&
                     " in file "//trim(filename)//" is being ignored.")
          endif
        endif
        force_cycle = .true.
      endif
      if (.not.found_override .and. (oval > 0)) then
        if (is_root_pe()) &
          call GOLD_error(WARNING,"GOLD_file_parser : "//trim(varname)// &
                   " has already been over-ridden."// &
                   " Line: '"//trim(line(:last))//"'"//&
                   " in file "//trim(filename)//" is being ignored.")
        force_cycle = .true.
      endif
      if (.not.found_override .and. (ival >= max_vals)) then
        if (is_root_pe()) then
          if ((defined_in_line /= defined) .or. &
              (trim(val_str) /= trim(value_string(max_vals)))) then
            call GOLD_error(FATAL,"GOLD_file_parser : "//trim(varname)// &
                     " found with multiple inconsistent definitions."// &
                     " Line: '"//trim(line(:last))//"'"//&
                     " in file "//trim(filename)//" caused the model failure.")
          else
            call GOLD_error(WARNING,"GOLD_file_parser : "//trim(varname)// &
                     " occurs more times than is permitted."// &
                     " Line: '"//trim(line(:last))//"'"//&
                     " in file "//trim(filename)//" is being ignored.")
          endif
        endif
        force_cycle = .true.
      endif
      if (force_cycle) cycle

      ! Store new values
      if (found_override) then
        oval = oval + 1
        value_string(oval) = trim(val_str)
        defined = defined_in_line
        if (verbose > 0 .and. ival > 0 .and. is_root_pe()) &
          call GOLD_error(WARNING,"GOLD_file_parser : "//trim(varname)// &
                 " over-ridden.  Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
      else ! (.not. found_overide)
        ival = ival + 1
        value_string(ival) = trim(val_str)
        defined = defined_in_line
        if (verbose > 1 .and. is_root_pe()) &
          call GOLD_error(WARNING,"GOLD_file_parser : "//trim(varname)// &
                 " set.  Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
      endif

    enddo ! param_file%param_data(ipf)%num_lines

  enddo paramfile_loop

end subroutine get_variable_line

subroutine flag_line_as_read(line_used,count)
  logical, dimension(:), pointer    :: line_used
  integer,               intent(in) :: count
  line_used(count) = .true.
end subroutine flag_line_as_read

function lowercase(input_string)
  character(len=*),     intent(in) :: input_string
  character(len=len(input_string)) :: lowercase
  integer, parameter :: co=iachar('a')-iachar('A') ! case offset
  integer :: k
!   This function returns a string in which all uppercase letters have been
! replaced by their lowercase counterparts.  It is loosely based on the
! lowercase function in mpp_util.F90.
  lowercase = input_string
  do k=1, len_trim(input_string)
    if (lowercase(k:k) >= 'A' .and. lowercase(k:k) <= 'Z') &
        lowercase(k:k) = achar(ichar(lowercase(k:k))+co)
  end do
end function lowercase

function uppercase(input_string)
  character(len=*),     intent(in) :: input_string
  character(len=len(input_string)) :: uppercase
  integer, parameter :: co=iachar('A')-iachar('a') ! case offset
  integer :: k
!   This function returns a string in which all lowercase letters have been
! replaced by their uppercase counterparts.  It is loosely based on the
! uppercase function in mpp_util.F90.
  uppercase = input_string
  do k=1, len_trim(input_string)
    if (uppercase(k:k) >= 'a' .and. uppercase(k:k) <= 'z') &
        uppercase(k:k) = achar(ichar(uppercase(k:k))+co)
  end do
end function uppercase

! The following subroutines write out to a log file.

subroutine log_version(param_file, modulename, version, tagname, desc)
  type(param_file_type), intent(in) :: param_file
  character(len=*), intent(in)      :: modulename
  character(len=*), intent(in)      :: version
  character(len=*), optional, intent(in) :: tagname, desc
! This subroutine writes the version of a module to a log file.
  character(len=240) :: mesg, myunits

!  write(mesg, '(a,": ",a)') trim(modulename), trim(version)
  if (present(tagname)) then
    mesg = trim(modulename)//": "//trim(version)//" Tag "//trim(tagname)
  else
    mesg = trim(modulename)//": "//trim(version)
  endif
  if (is_root_pe()) then
    if (param_file%log_open) write(param_file%stdlog,'(a)') trim(mesg)
    if (param_file%log_to_stdout) write(param_file%stdout,'(a)') trim(mesg)
  endif

  if (present(desc)) call doc_module(modulename, desc)

end subroutine log_version

subroutine log_param_int(param_file, modulename, varname, value, desc, units, &
                         default)
  type(param_file_type),      intent(in) :: param_file
  character(len=*),           intent(in) :: modulename
  character(len=*),           intent(in) :: varname
  integer,                    intent(in) :: value
  character(len=*), optional, intent(in) :: desc, units
  integer,          optional, intent(in) :: default
! This subroutine writes the value of an integer parameter to a log file,
! along with its name and the module it came from.
  character(len=240) :: mesg, myunits

  write(mesg, '("  ",a," ",a,": ",I)') trim(modulename), trim(varname), value
  if (is_root_pe()) then
    if (param_file%log_open) write(param_file%stdlog,'(a)') trim(mesg)
    if (param_file%log_to_stdout) write(param_file%stdout,'(a)') trim(mesg)
  endif

  myunits=" "; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) call doc_param(varname, desc, myunits, value, default)

end subroutine log_param_int

subroutine log_param_int_array(param_file, modulename, varname, value, desc, &
                               units, default)
  type(param_file_type),      intent(in) :: param_file
  character(len=*),           intent(in) :: modulename
  character(len=*),           intent(in) :: varname
  integer,                    intent(in) :: value(:)
  character(len=*), optional, intent(in) :: desc, units
  integer,          optional, intent(in) :: default
! This subroutine writes the value of an integer parameter to a log file,
! along with its name and the module it came from.
  character(len=1320) :: mesg, myunits

  write(mesg, '("  ",a," ",a,": ",A)') trim(modulename), trim(varname), trim(left_ints(value))
  if (is_root_pe()) then
    if (param_file%log_open) write(param_file%stdlog,'(a)') trim(mesg)
    if (param_file%log_to_stdout) write(param_file%stdout,'(a)') trim(mesg)
  endif

  myunits=" "; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) call doc_param(varname, desc, myunits, value, default)

end subroutine log_param_int_array

function left_int(i)
  character(len=19) :: left_int
  integer, intent(in) :: i
! Returns a character string of a left-formatted integer
! e.g. "123       "  (assumes 19 digit maximum)
  character(len=19) :: tmp
  write(tmp(1:19),'(I19)') i
  write(left_int(1:19),'(A)') adjustl(tmp)
end function left_int

function left_ints(i)
  integer, intent(in) :: i(:)
! Returns a character string of a comma-separated, compact formatted, integers
! e.g. "1, 2, 3, 4"
  character(len=1320) :: left_ints,tmp
  integer :: j
  write(left_ints(1:1320),'(A)') trim(left_int(i(1)))
  if (size(i)>1) then
    do j=2,size(i)
      tmp=left_ints
      write(left_ints(1:1320),'(A,", ",A)') trim(tmp),trim(left_int(i(j)))
    enddo
  endif
end function left_ints

subroutine log_param_real(param_file, modulename, varname, value, desc, units, &
                          default)
  type(param_file_type),      intent(in) :: param_file
  character(len=*),           intent(in) :: modulename
  character(len=*),           intent(in) :: varname
  real,                       intent(in) :: value
  character(len=*), optional, intent(in) :: desc, units
  real,             optional, intent(in) :: default
! This subroutine writes the value of a real parameter to a log file,
! along with its name and the module it came from.
  character(len=240) :: mesg, myunits

  write(mesg, '("  ",a," ",a,": ",ES19.12)') &
    trim(modulename), trim(varname), value
  if (is_root_pe()) then
    if (param_file%log_open) write(param_file%stdlog,'(a)') trim(mesg)
    if (param_file%log_to_stdout) write(param_file%stdout,'(a)') trim(mesg)
  endif

  myunits="not defined"; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) call doc_param(varname, desc, myunits, value, default)

end subroutine log_param_real

subroutine log_param_real_array(param_file, modulename, varname, value, desc, &
                                units, default)
  type(param_file_type),      intent(in) :: param_file
  character(len=*),           intent(in) :: modulename
  character(len=*),           intent(in) :: varname
  real,                       intent(in) :: value(:)
  character(len=*), optional, intent(in) :: desc, units
  real,             optional, intent(in) :: default
! This subroutine writes the value of a real parameter to a log file,
! along with its name and the module it came from.
  character(len=1320) :: mesg, myunits

 !write(mesg, '("  ",a," ",a,": ",ES19.12,99(",",ES19.12))') &
 !write(mesg, '("  ",a," ",a,": ",G,99(",",G))') &
 !  trim(modulename), trim(varname), value
  write(mesg, '("  ",a," ",a,": ",a)') &
    trim(modulename), trim(varname), trim(left_reals(value))
  if (is_root_pe()) then
    if (param_file%log_open) write(param_file%stdlog,'(a)') trim(mesg)
    if (param_file%log_to_stdout) write(param_file%stdout,'(a)') trim(mesg)
  endif

  myunits="not defined"; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) call doc_param(varname, desc, myunits, value, default)

end subroutine log_param_real_array

function left_real(r)
  character(len=30) :: left_real
  real, intent(in) :: r
! Returns a character string of a left-formatted real
! using either F or E format with trailing 0s removed
! e.g. "12.345    "
  character(len=30) :: tmp,expon
  integer :: i
 !write(0,*) 'In left_real, r=',r
  if ((abs(r)<1.e4.and.abs(r)>=1.e-3).or.r.eq.0.) then
    write(tmp(1:30),'(F,X)') r
  else
    write(tmp(1:30),'(ES,X)') r
  endif
  i=index(tmp,'E')
  ! Save exponent
  if (i>0) then
    write(expon(1:30),'(A)') trim(tmp(i:len_trim(tmp)))
    tmp(i:len_trim(tmp))=' '
  else
    expon=' '
  endif
  ! Left justify
  write(left_real(1:30),'(A)') adjustl(tmp)
  ! Strip trailing zeros
  do i=len_trim(left_real),2,-1
    if (left_real(i:i).ne.'0') exit
    left_real(i:i)=' '
  enddo
  ! Re-attach exponent
  write(left_real(i+1:len(left_real)),'(A)') trim(expon)
end function left_real

function left_reals(r)
  real, intent(in) :: r(:)
! Returns a character string of a comma-separated, compact formatted, reals
! e.g. "1., 2., 3., 5.E2"
  character(len=1320) :: left_reals,tmp
  integer :: j
  write(left_reals(1:1320),'(A)') trim(left_real(r(1)))
  if (size(r)>1) then
    do j=2,size(r)
      tmp=left_reals
      write(left_reals(1:1320),'(A,", ",A)') trim(tmp),trim(left_real(r(j)))
    enddo
  endif
end function left_reals

subroutine log_param_logical(param_file, modulename, varname, value, desc, &
                             units, default)
  type(param_file_type),      intent(in) :: param_file
  character(len=*),           intent(in) :: modulename
  character(len=*),           intent(in) :: varname
  logical,                    intent(in) :: value
  character(len=*), optional, intent(in) :: desc, units
  logical,          optional, intent(in) :: default
! This subroutine writes the value of a logical parameter to a log file,
! along with its name and the module it came from.
  character(len=240) :: mesg, myunits

  if (value) then
    write(mesg, '("  ",a," ",a,": True")') trim(modulename), trim(varname)
  else
    write(mesg, '("  ",a," ",a,": False")') trim(modulename), trim(varname)
  endif
  if (is_root_pe()) then
    if (param_file%log_open) write(param_file%stdlog,'(a)') trim(mesg)
    if (param_file%log_to_stdout) write(param_file%stdout,'(a)') trim(mesg)
  endif

  myunits="Boolean"; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) call doc_param(varname, desc, myunits, value, default)

end subroutine log_param_logical

subroutine log_param_char(param_file, modulename, varname, value, desc, units, &
                          default)
  type(param_file_type),      intent(in) :: param_file
  character(len=*),           intent(in) :: modulename
  character(len=*),           intent(in) :: varname
  character(len=*),           intent(in) :: value
  character(len=*), optional, intent(in) :: desc, units
  character(len=*), optional, intent(in) :: default
! This subroutine writes the value of a character string parameter to a log
! file, along with its name and the module it came from.
  character(len=240) :: mesg, myunits

  write(mesg, '("  ",a," ",a,": ",a)') &
    trim(modulename), trim(varname), trim(value)
  if (is_root_pe()) then
    if (param_file%log_open) write(param_file%stdlog,'(a)') trim(mesg)
    if (param_file%log_to_stdout) write(param_file%stdout,'(a)') trim(mesg)
  endif

  myunits=" "; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) call doc_param(varname, desc, myunits, value, default)

end subroutine log_param_char

subroutine log_param_time(param_file, modulename, varname, value, desc, units, &
                          default, timeunit)
  type(param_file_type),      intent(in) :: param_file
  character(len=*),           intent(in) :: modulename
  character(len=*),           intent(in) :: varname
  type(time_type),            intent(in) :: value
  character(len=*), optional, intent(in) :: desc, units
  type(time_type),  optional, intent(in) :: default
  real,             optional, intent(in) :: timeunit
! This subroutine writes the value of a time-type parameter to a log file,
! along with its name and the module it came from.
  real :: real_time, real_default
  logical :: use_timeunit = .false.
  character(len=240) :: mesg, myunits
  integer :: days, secs, ticks

  call get_time(value, secs, days, ticks)

  if (ticks == 0) then
    write(mesg, '("  ",a," ",a," (Time): ",i,":",i)') trim(modulename), &
       trim(varname), days, secs
  else
    write(mesg, '("  ",a," ",a," (Time): ",i,":",i,":",i)') trim(modulename), &
       trim(varname), days, secs, ticks
  endif
  if (is_root_pe()) then
    if (param_file%log_open) write(param_file%stdlog,'(a)') trim(mesg)
    if (param_file%log_to_stdout) write(param_file%stdout,'(a)') trim(mesg)
  endif

  if (present(desc)) then
    if (present(timeunit)) use_timeunit = (timeunit > 0.0)
    if (use_timeunit) then
      if (present(units)) then
        write(myunits(1:240),'(A)') trim(units)
      else
        if (abs(timeunit-1.0) < 0.01) then ; myunits = "seconds"
        elseif (abs(timeunit-3600.0) < 1.0) then ; myunits = "hours"
        elseif (abs(timeunit-86400.0) < 1.0) then ; myunits = "days"
        elseif (abs(timeunit-3.1e7) < 1.0e6) then ; myunits = "years"
        else ; write(myunits,'(es8.2," sec")') timeunit ; endif
      endif
      real_time = (86400.0/timeunit)*days + secs/timeunit
      if (ticks > 0) real_time = real_time + &
                           real(ticks) / (timeunit*get_ticks_per_second())
      if (present(default)) then
        call get_time(default, secs, days, ticks)
        real_default = (86400.0/timeunit)*days + secs/timeunit
        if (ticks > 0) real_default = real_default + &
                           real(ticks) / (timeunit*get_ticks_per_second())
        call doc_param(varname, desc, myunits, real_time, real_default)
      else
        call doc_param(varname, desc, myunits, real_time)
      endif
    else
      myunits='not defined'; if (present(units)) write(myunits(1:240),'(A)') trim(units)
      call doc_param(varname, desc, myunits, value, default)
    endif
  endif

end subroutine log_param_time

end module GOLD_file_parser
