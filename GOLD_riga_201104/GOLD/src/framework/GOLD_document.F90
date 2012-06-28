module GOLD_document
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
!*    The subroutines here provide hooks for document generation       *
!*  functions at various levels of granularity.                        *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_time_manager, only : time_type
use GOLD_error_handler, only : GOLD_error, FATAL, is_root_pe

implicit none ; private

public doc_param, doc_subroutine, doc_function, doc_module, doc_init, doc_end

interface doc_param
  module procedure doc_param_none, &
                   doc_param_logical, doc_param_logical_array, &
                   doc_param_int,     doc_param_int_array, &
                   doc_param_real,    doc_param_real_array, &
                   doc_param_char, &
                   doc_param_time
end interface

integer :: doc_unit = -1
logical :: minimal_doc = .false.

contains

! ----------------------------------------------------------------------

subroutine doc_param_none(varname, desc, units)
  character(len=*),           intent(in) :: varname, desc, units
! This subroutine handles parameter documentation with no value.
  integer :: numspc
  character(len=240) :: mesg

  if (is_root_pe() .and. (doc_unit > 0)) then
    numspc = max(1,32-8-len_trim(varname))
    mesg = "#define "//trim(varname)//repeat(" ",numspc)//"!"
    if (len_trim(units) > 0) mesg = trim(mesg)//"   ["//trim(units)//"]"

    write(doc_unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(desc)
  endif 
end subroutine doc_param_none

subroutine doc_param_logical(varname, desc, units, val, default)
  character(len=*),           intent(in) :: varname, desc, units
  logical,                    intent(in) :: val
  logical,          optional, intent(in) :: default
! This subroutine handles parameter documentation for logicals.
  integer :: numspc
  character(len=240) :: mesg

  if (is_root_pe() .and. (doc_unit > 0)) then
    if (present(default) .and. minimal_doc) then
      if (val == default) return
    endif

    numspc = max(1,32-8-len_trim(varname))
    if (val) then
      mesg = "#define "//trim(varname)//repeat(" ",numspc)//&
             "!   ["//trim(units)//"]"
    else
      mesg = "#undef  "//trim(varname)//repeat(" ",numspc)//&
             "!   ["//trim(units)//"]"
    endif

    if (present(default)) then
      if (default) then
        mesg = trim(mesg)//" default = True"
      else
        mesg = trim(mesg)//" default = False"
      endif
    endif

    write(doc_unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(desc)
  endif
end subroutine doc_param_logical

subroutine doc_param_logical_array(varname, desc, units, vals, default)
  character(len=*),           intent(in) :: varname, desc, units
  logical,                    intent(in) :: vals(:)
  logical,          optional, intent(in) :: default
! This subroutine handles parameter documentation for arrays of logicals.
  logical :: go_back
  integer :: i, numspc
  character(len=1280) :: mesg
  character(len=1240)  :: valstring

  if (is_root_pe() .and. (doc_unit > 0)) then
    if (present(default) .and. minimal_doc) then
      go_back = .true.
      do i=1,size(vals) ; if (vals(i) /= default) go_back = .false. ; enddo
      if (go_back) return
    endif

    if (vals(1)) then ; valstring = "True" ; else ; valstring = "False" ; endif
    do i=2,min(size(vals),128)
      if (vals(i)) then
        valstring = trim(valstring)//", True"
      else
        valstring = trim(valstring)//", False"
      endif
    enddo
    numspc = max(1,32-9-len_trim(varname)-len_trim(valstring))

    mesg = "#define "//trim(varname)//" "//trim(valstring)//&
           repeat(" ",numspc)//"!   ["//trim(units)//"]"

    if (present(default)) then
      if (default) then
        mesg = trim(mesg)//" default = True"
      else
        mesg = trim(mesg)//" default = False"
      endif
    endif

    write(doc_unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(desc)
  endif
end subroutine doc_param_logical_array

subroutine doc_param_int(varname, desc, units, val, default)
  character(len=*),           intent(in) :: varname, desc, units
  integer,                    intent(in) :: val
  integer,          optional, intent(in) :: default
! This subroutine handles parameter documentation for integers.
  integer :: numspc
  character(len=240) :: mesg
  character(len=32)  :: valstring

  if (is_root_pe() .and. (doc_unit > 0)) then
    if (present(default) .and. minimal_doc) then
      if (val == default) return
    endif

    valstring = int_string(val)
    numspc = max(1,32-9-len_trim(varname)-len_trim(valstring))

    mesg = "#define "//trim(varname)//" "//trim(valstring)//&
           repeat(" ",numspc)//"!"
    if (len_trim(units) > 0) mesg = trim(mesg)//"   ["//trim(units)//"]"

    if (present(default)) &
      mesg = trim(mesg)//" default = "//(trim(int_string(default)))

    write(doc_unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(desc)
  endif
end subroutine doc_param_int

subroutine doc_param_int_array(varname, desc, units, vals, default)
  character(len=*),           intent(in) :: varname, desc, units
  integer,                    intent(in) :: vals(:)
  integer,          optional, intent(in) :: default
! This subroutine handles parameter documentation for arrays of integers.
  logical :: go_back
  integer :: i, numspc
  character(len=1280) :: mesg
  character(len=1240)  :: valstring

  if (is_root_pe() .and. (doc_unit > 0)) then
    if (present(default) .and. minimal_doc) then
      go_back = .true.
      do i=1,size(vals) ; if (vals(i) /= default) go_back = .false. ; enddo
      if (go_back) return
    endif

    valstring = int_string(vals(1))
    do i=2,min(size(vals),128)
      valstring = trim(valstring)//", "//trim(int_string(vals(i)))
    enddo
    numspc = max(1,32-9-len_trim(varname)-len_trim(valstring))

    mesg = "#define "//trim(varname)//" "//trim(valstring)//&
           repeat(" ",numspc)//"!"
    if (len_trim(units) > 0) mesg = trim(mesg)//"   ["//trim(units)//"]"

    if (present(default)) &
      mesg = trim(mesg)//" default = "//(trim(int_string(default)))

    write(doc_unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(desc)
  endif

end subroutine doc_param_int_array

subroutine doc_param_real(varname, desc, units, val, default)
  character(len=*),           intent(in) :: varname, desc, units
  real,                       intent(in) :: val
  real,             optional, intent(in) :: default
! This subroutine handles parameter documentation for reals.
  integer :: numspc
  character(len=240) :: mesg
  character(len=32)  :: valstring

  if (is_root_pe() .and. (doc_unit > 0)) then
    if (present(default) .and. minimal_doc) then
      if (val == default) return
    endif

    valstring = real_string(val)
    numspc = max(1,32-9-len_trim(varname)-len_trim(valstring))

    mesg = "#define "//trim(varname)//" "//trim(valstring)//&
           repeat(" ",numspc)//"!   ["//trim(units)//"]"

    if (present(default)) &
      mesg = trim(mesg)//" default = "//trim(real_string(default))

    write(doc_unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(desc)
  endif
end subroutine doc_param_real

subroutine doc_param_real_array(varname, desc, units, vals, default)
  character(len=*),           intent(in) :: varname, desc, units
  real,                       intent(in) :: vals(:)
  real,             optional, intent(in) :: default
! This subroutine handles parameter documentation for arrays of reals.
  logical :: go_back
  integer :: i, numspc
  character(len=1280) :: mesg
  character(len=1240)  :: valstring

  if (is_root_pe() .and. (doc_unit > 0)) then
    if (present(default) .and. minimal_doc) then
      go_back = .true.
      do i=1,size(vals) ; if (vals(i) /= default) go_back = .false. ; enddo
      if (go_back) return
    endif

    valstring = real_string(vals(1))
    do i=2,min(size(vals),128)
      valstring = trim(valstring)//", "//trim(real_string(vals(i)))
    enddo

    numspc = max(1,32-9-len_trim(varname)-len_trim(valstring))

    mesg = "#define "//trim(varname)//" "//trim(valstring)//&
           repeat(" ",numspc)//"!   ["//trim(units)//"]"

    if (present(default)) &
      mesg = trim(mesg)//" default = "//trim(real_string(default))

    write(doc_unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(desc)
  endif

end subroutine doc_param_real_array

subroutine doc_param_char(varname, desc, units, val, default)
  character(len=*),           intent(in) :: varname, desc, units
  character(len=*),           intent(in) :: val
  character(len=*), optional, intent(in) :: default
! This subroutine handles parameter documentation for character strings.
  integer :: numspc
  character(len=240) :: mesg

  if (is_root_pe() .and. (doc_unit > 0)) then
    if (present(default) .and. minimal_doc) then
      if (trim(val) == trim(default)) return
    endif

    numspc = max(1,32-9-len_trim(varname)-len_trim(val))

    mesg = "#define "//trim(varname)//" "//trim(val)//&
           repeat(" ",numspc)//"!"
    if (len_trim(units) > 0) mesg = trim(mesg)//"   ["//trim(units)//"]"

    if (present(default)) &
      mesg = trim(mesg)//" default = "//trim(adjustl(default))

    write(doc_unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(desc)
  endif

end subroutine doc_param_char

subroutine doc_param_time(varname, desc, units, val, default)
  character(len=*),           intent(in) :: varname, desc, units
  type(time_type),            intent(in) :: val
  type(time_type),  optional, intent(in) :: default
! This subroutine handles parameter documentation for time-type variables.
!  ### This needs to be written properly!
  integer :: numspc
  character(len=240) :: mesg

  if (is_root_pe() .and. (doc_unit > 0)) then
    numspc = max(1,32-18-len_trim(varname))
    mesg = "#define "//trim(varname)//" Time-type"//repeat(" ",numspc)//"!"
    if (len_trim(units) > 0) mesg = trim(mesg)//"   ["//trim(units)//"]"

    write(doc_unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(desc)
  endif 

end subroutine doc_param_time

subroutine write_desc(desc, indent)
  character(len=*),           intent(in) :: desc
  integer,          optional, intent(in) :: indent
  character(len=240) :: mesg
  integer :: start_ind = 1, end_ind, indnt, tab

  indnt = 32 ; if (present(indent)) indnt = indent
  start_ind = 1
  do
    if (len_trim(desc(start_ind:)) < 1) exit

    end_ind = index(desc(start_ind:), "\n")

    if (end_ind > 0) then
      mesg = repeat(" ",indnt)//"! "//trim(desc(start_ind:start_ind+end_ind-2))
      do ; tab = index(mesg, "\t")
        if (tab == 0) exit
        mesg(tab:tab+1) = "  "
      enddo
      write(doc_unit, '(a)') trim(mesg)
      start_ind = start_ind+end_ind+1
    else
      mesg = repeat(" ",indnt)//"! "//trim(desc(start_ind:))
      do ; tab = index(mesg, "\t")
        if (tab == 0) exit
        mesg(tab:tab+1) = "  "
      enddo
      write(doc_unit, '(a)') trim(mesg)
      exit
    endif

  enddo
end subroutine write_desc

! ----------------------------------------------------------------------

function real_string(val)
  real, intent(in)  :: val
  character(len=32) :: real_string
  
  integer :: len, ind

  if ((abs(val) < 1.0e4) .and. (abs(val) >= 1.0e-3)) then
    write(real_string, '(F30.12)') val
    do
      len = len_trim(real_string)
      if ((len<2) .or. (real_string(len-1:len) == ".0") .or. &
          (real_string(len:len) /= "0")) exit
      real_string(len:len) = " "
    enddo
  elseif (val == 0) then
    real_string = "0.0"
  else
    write(real_string, '(ES21.15)') val
    do
      ind = index(real_string,"0E")
      if (ind == 0) exit
      if (real_string(ind-1:ind-1) == ".") exit
      real_string = real_string(1:ind-1)//real_string(ind+1:)
    enddo
  endif
  real_string = adjustl(real_string)
  
end function real_string

function int_string(val)
  integer, intent(in)  :: val
  character(len=24)    :: int_string
  write(int_string, '(i24)') val
  int_string = adjustl(int_string)
end function int_string

! ----------------------------------------------------------------------

subroutine doc_module(modname, desc)
  character(len=*),           intent(in) :: modname, desc
! This subroutine handles the module documentation
  character(len=240) :: mesg
  if (is_root_pe() .and. (doc_unit > 0)) then
    mesg = "    !  Parameters of module "//trim(modname)
    write(doc_unit, '(a)') trim(mesg)
    if (len_trim(desc) > 0) call write_desc(desc,8)
  endif 
end subroutine doc_module

subroutine doc_subroutine(modname, subname, desc)
  character(len=*),           intent(in) :: modname, subname, desc
! This subroutine handles the subroutine documentation
end subroutine doc_subroutine

subroutine doc_function(modname, fnname, desc)
  character(len=*),           intent(in) :: modname, fnname, desc
! This subroutine handles the function documentation
end subroutine doc_function

! ----------------------------------------------------------------------
  
subroutine doc_init(docfile, minimal)
  character(len=*),   intent(in) :: docfile
  logical,  optional, intent(in) :: minimal
! Argument: The name of the doc file.
  integer :: ios, unit
  logical :: opened, new_file

  if (present(minimal)) minimal_doc = minimal

  if (is_root_pe() .and. (len_trim(docfile) > 0) .and. (doc_unit<0)) then
    new_file = .true. ; if (doc_unit /= -1) new_file = .false.
    ! Find an unused unit number.
    do doc_unit=512,42,-1
      inquire( doc_unit, opened=opened)
      if (.not.opened) exit
    enddo

    if (opened) call GOLD_error(FATAL, &
        "doc_init failed to find an unused unit number.")

    if (new_file) then
      open(doc_unit, file=trim(docfile), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='REPLACE', iostat=ios)
    else ! This file is being reopened, and should be appended.
      open(doc_unit, file=trim(docfile), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='OLD', position='APPEND', iostat=ios)
    endif
    inquire(doc_unit, opened=opened)
    if ((.not.opened) .or. (ios /= 0)) then
      call GOLD_error(FATAL, "Failed to open doc file "//trim(docfile)//".")
    endif
  endif

end subroutine doc_init

subroutine doc_end

  if (doc_unit > 0) then
    close(doc_unit)
    doc_unit = -2
  endif

end subroutine doc_end

end module GOLD_document
