module GOLD_obsolete_params
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
!
!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, July 2010                                      *
!*                                                                     *
!*    The subroutine in this module checks whether obsolete parameters *
!*  are being used, and stops the model or issues a warning as         *
!*  appropriate.                                                       *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_error_handler, only : GOLD_error, FATAL, WARNING, is_root_pe
use GOLD_file_parser, only : read_param, log_version, param_file_type

implicit none ; private

#include <GOLD_memory.h>

public find_obsolete_params

contains

subroutine find_obsolete_params(param_file)
  type(param_file_type),       intent(in)    :: param_file
! Argument: param_file - A structure indicating the open file to parse for
!                        model parameter values.
  character(len=128) :: version = '$Id: GOLD_obsolete_params.F90,v 1.1.2.6 2011/08/10 15:46:02 Robert.Hallberg Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "find_obsolete_params"  ! This module's name.
  real    :: test_val
  integer :: test_int
  logical :: test_logic, test_logic2, test_logic3

  if (.not.is_root_pe()) return

  test_int = -1 ; call read_param(param_file,"NTSTEP",test_int)
  if (test_int /= -1) call GOLD_ERROR(FATAL, "find_obsolete_params: "// &
      "NTSTEP is an obsolete option.  Instead #define DT_THERM " // &
      "to set the thermodynamic time step.")

  test_logic = .false. ; call read_param(param_file,"JACOBIAN_PGF",test_logic)
  if (test_logic) call GOLD_ERROR(FATAL, "find_obsolete_params: "// &
         "JACOBIAN_PGF is an obsolete run-time flag.  Instead use "// &
         "#define ANALYTIC_FV_PGF.")
  test_logic = .true. ; call read_param(param_file,"JACOBIAN_PGF",test_logic)
  if (.not.test_logic) call GOLD_ERROR(WARNING, "find_obsolete_params: "// &
         "JACOBIAN_PGF is an obsolete run-time flag.")

  test_logic = .true. ; call read_param(param_file,"SADOURNY",test_logic)
  test_logic2 = .false. ; call read_param(param_file,"SADOURNY",test_logic2)
  if (test_logic  == test_logic2) call GOLD_ERROR(FATAL, "find_obsolete_params: "// &
         "SADOURNY is an obsolete run-time flag.  Use CORIOLIS_SCHEME instead.")

  ! Candidate for obsolescence?: BT_COR_SLOW_RATE
!  test_val = 0.0 ; call read_param(param_file,"BT_COR_SLOW_RATE",test_val)
!  if (test_val /= 0.0) call GOLD_ERROR(FATAL, "find_obsolete_params: "// &
!      "BT_COR_SLOW_RATE is an obsolete option.  Instead try using "// &
!      "#define BT_MASS_SOURCE_LIMIT 0.1.")
! test_val = 1.0 ; call read_param(param_file,"BT_COR_SLOW_RATE",test_val)
!  if (test_val /= 1.0) call GOLD_ERROR(WARNING, "find_obsolete_params: "// &
!      "BT_COR_SLOW_RATE is an obsolete option.")

  test_logic = .true. ; call read_param(param_file,"USE_LOCAL_PREF_CORRECT",test_logic)
  test_logic2 = .false. ; call read_param(param_file,"USE_LOCAL_PREF_CORRECT",test_logic2)
  if (test_logic  == test_logic2) then
    test_logic3 = .false. ; call read_param(param_file,"USE_LOCAL_PREF",test_logic3)
    if (test_logic3 == test_logic) then
      call GOLD_ERROR(WARNING,"find_obsolete_params: USE_LOCAL_PREF_CORRECT "// &
          "is an obsolete run-time flag, but is set consistently with USE_LOCAL_PREF.")
    else
      call GOLD_ERROR(FATAL,"find_obsolete_params: USE_LOCAL_PREF_CORRECT "// &
          "is an obsolete run-time flag.  Use [#define|#undef] USE_LOCAL_PREF instead.")
    endif
  endif

  test_logic = .false. ; call read_param(param_file,"RIGA_ITIDE_BUGS",test_logic)
  if (test_logic) call GOLD_ERROR(FATAL, "find_obsolete_params: "// &
         "RIGA_ITIDE_BUGS is an obsolete run-time flag, and should not be used.")
  test_logic = .true. ; call read_param(param_file,"RIGA_ITIDE_BUGS",test_logic)
  if (.not.test_logic) call GOLD_ERROR(WARNING, "find_obsolete_params: "// &
         "RIGA_ITIDE_BUGS is an obsolete run-time flag.")

  test_int = -1 ; call read_param(param_file,"ML_RADIATION_CODING",test_int)
  if (test_int == 1) call GOLD_ERROR(FATAL, "find_obsolete_params: "// &
    "ML_RADIATION_CODING is an obsolete option and the code previously "//&
    "used by setting it to 1 has been eliminated.")
  if (test_int /= -1) call GOLD_ERROR(WARNING, "find_obsolete_params: "// &
    "ML_RADIATION_CODING is an obsolete option.")


  ! Test for inconsistent parameter settings.
  test_logic = .false. ; test_logic2 = .false.
  call read_param(param_file,"SPLIT",test_logic)
  call read_param(param_file,"DYNAMIC_SURFACE_PRESSURE",test_logic2)
  if (test_logic2 .and. .not.test_logic) call GOLD_ERROR(FATAL, &
    "find_obsolete_params: #define DYNAMIC_SURFACE_PRESSURE is not yet "//&
    "implemented without #define SPLIT.")

  ! Write the file version number to the model log.
  call log_version(param_file, mod, version, tagname)

end subroutine find_obsolete_params

end module GOLD_obsolete_params
