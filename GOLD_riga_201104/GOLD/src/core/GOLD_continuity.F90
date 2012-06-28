module GOLD_continuity
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
!*  By Robert Hallberg and Alistair Adcroft, September 2006.           *
!*                                                                     *
!*    This file contains the driver routine which selects which        *
!*  continuity solver will be called, based on run-time input.         *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, vh                                    *
!*    j    x ^ x ^ x   At >:  u, uh                                    *
!*    j    > o > o >   At o:  h, hin                                   *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use GOLD_continuity_PPM, only : continuity_PPM, continuity_PPM_init
use GOLD_continuity_PPM, only : continuity_PPM_end, continuity_PPM_CS
use GOLD_continuity_Hallberg95, only : continuity_orig, continuity_orig_init
use GOLD_continuity_Hallberg95, only : continuity_orig_end, continuity_orig_CS
use GOLD_diag_mediator, only : time_type, diag_ptrs
use GOLD_domains, only : pass_var, pass_vector
use GOLD_error_handler, only : GOLD_error, FATAL, WARNING, is_root_pe
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_grid, only : ocean_grid_type
use GOLD_variables, only : ocean_OBC_type, BT_cont_type
use GOLD_vert_friction, only : vertvisc_remean

implicit none ; private

#include <GOLD_memory.h>

public continuity, continuity_init, continuity_end

integer :: id_clock_pass, id_clock_vertvisc

type, public :: continuity_CS ; private
  logical :: use_continuity_PPM  ! If true, use the PPM continuity solver.
  real    :: maxvel         ! Velocity components greater than maxvel,
                            ! in m s-1, are truncated.
  type(continuity_PPM_CS), pointer :: PPM_CSp => NULL()
  type(continuity_orig_CS), pointer :: orig_CSp => NULL()
end type continuity_CS

contains

subroutine continuity(u, v, hin, h, uh, vh, dt, G, CS, uhbt, vhbt, OBC, &
                      visc_rem_u, visc_rem_v, u_cor, v_cor, &
                      uhbt_aux, vhbt_aux, u_cor_aux, v_cor_aux, BT_cont)
  real, intent(in),  dimension(NXMEMQ_,NYMEM_,NZ_) :: u
  real, intent(in),  dimension(NXMEM_,NYMEMQ_,NZ_) :: v
  real, intent(in),  dimension(NXMEM_,NYMEM_,NZ_)  :: hin
  real, intent(out), dimension(NXMEM_,NYMEM_,NZ_)  :: h
  real, intent(out), dimension(NXMEMQ_,NYMEM_,NZ_) :: uh
  real, intent(out), dimension(NXMEM_,NYMEMQ_,NZ_) :: vh
  real, intent(in)                                 :: dt
  type(ocean_grid_type), intent(inout)             :: G
  type(continuity_CS), pointer                     :: CS
  real, intent(in), optional, dimension(NXMEMQ_,NYMEM_) :: uhbt
  real, intent(in), optional, dimension(NXMEM_,NYMEMQ_) :: vhbt
  type(ocean_OBC_type), pointer, optional          :: OBC
  real, intent(in), optional, dimension(NXMEMQ_,NYMEM_,NZ_) :: visc_rem_u
  real, intent(in), optional, dimension(NXMEM_,NYMEMQ_,NZ_) :: visc_rem_v
  real, intent(out), optional, dimension(NXMEMQ_,NYMEM_,NZ_) :: u_cor
  real, intent(out), optional, dimension(NXMEM_,NYMEMQ_,NZ_) :: v_cor
  real, intent(in), optional, dimension(NXMEMQ_,NYMEM_) :: uhbt_aux
  real, intent(in), optional, dimension(NXMEM_,NYMEMQ_) :: vhbt_aux
  real, intent(out), optional, dimension(NXMEMQ_,NYMEM_,NZ_) :: u_cor_aux
  real, intent(out), optional, dimension(NXMEM_,NYMEMQ_,NZ_) :: v_cor_aux
  type(BT_cont_type),                  pointer,     optional :: BT_cont
!    This subroutine time steps the layer thicknesses, using a monotonically
!  limit, directionally split PPM scheme, based on Lin (1994).

! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      hin - Initial layer thickness, in m.
!  (out)     h - Final layer thickness, in m.
!  (out)     uh - Volume flux through zonal faces = u*h*dy, m3 s-1.
!  (out)     vh - Volume flux through meridional faces = v*h*dx,
!                  in m3 s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 continuity_init.
!  (in, opt) uhbt - The summed volume flux through zonal faces, m3 s-1.
!  (in, opt) vhbt - The summed volume flux through meridional faces, m3 s-1.
!  (in, opt) OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in, opt) visc_rem_u - Both the fraction of the momentum originally in a
!  (in, opt) visc_rem_v - layer that remains after a time-step of viscosity,
!                         and the fraction of a time-step's worth of a
!                         barotropic acceleration that a layer experiences
!                         after viscosity is applied, in the zonal (_u) and
!                         meridional (_v) directions.  Nondimensional between
!                         0 (at the bottom) and 1 (far above the bottom).
!  (out, opt) u_cor - The zonal velocities that give uhbt as the depth-
!                     integrated transport, in m s-1.
!  (out, opt) v_cor - The meridional velocities that give vhbt as the
!                     depth-integrated transport, in m s-1.
!  (in, opt) uhbt_aux - A second set of summed volume fluxes through zonal
!  (in, opt) vhbt_aux - and meridional faces, both in m3 s-1.
!  (out, opt) u_cor_aux - The zonal and meridional velocities that give uhbt_aux
!  (out, opt) v_cor_aux - and vhbt_aux as the depth-integrated transports,
!                         both in m s-1.
!  (out, opt) BT_cont - A structure with elements that describe the effective
!                       open face areas as a function of barotropic flow.
  if (present(visc_rem_u) /= present(visc_rem_v)) call GOLD_error(FATAL, &
      "GOLD_continuity: Either both visc_rem_u and visc_rem_v or neither"// &
       " one must be present in call to continuity.")
  if (present(u_cor) /= present(v_cor)) call GOLD_error(FATAL, &
      "GOLD_continuity: Either both u_cor and v_cor or neither"// &
       " one must be present in call to continuity.")
  if (present(uhbt_aux) /= present(vhbt_aux)) call GOLD_error(FATAL, &
      "GOLD_continuity: Either both uhbt_aux and uhbt_aux or neither"// &
       " one must be present in call to continuity.")
  if (present(u_cor_aux) /= present(v_cor_aux)) call GOLD_error(FATAL, &
      "GOLD_continuity: Either both u_cor_aux and v_cor_aux or neither"// &
       " one must be present in call to continuity.")
  if (present(u_cor_aux) /= present(uhbt_aux)) call GOLD_error(FATAL, &
      "GOLD_continuity: u_cor_aux can only be calculated if uhbt_aux is"// &
      " provided, and uhbt_aux has no other purpose.  Include both arguments"//&
      " or neither.")

  if (CS%use_continuity_PPM) then
    call continuity_PPM(u, v, hin, h, uh, vh, dt, G, CS%PPM_CSp, uhbt, vhbt, OBC, &
                        visc_rem_u, visc_rem_v, u_cor, v_cor, &
                        uhbt_aux, vhbt_aux, u_cor_aux, v_cor_aux, BT_cont)
  else
    if (present(BT_cont)) then ; if (associated(BT_cont)) &
      call GOLD_error(FATAL, "GOLD_continuity: "//&
        "The use of the BT_cont_type is only available with CONTINUITY_PPM.")
    endif
    if (present(visc_rem_u) .and. present(u_cor)) then
      ! This call is necessary because continuity_orig lacks the ability to
      ! return the "corrected" velocities.
      call cpu_clock_begin(id_clock_vertvisc)
      call vertvisc_remean(u, v, hin, uhbt, vhbt, visc_rem_u, visc_rem_v, &
                           u_cor, v_cor, OBC, G, CS%maxvel)
      if (present(u_cor_aux)) &
        call vertvisc_remean(u, v, hin, uhbt_aux, vhbt_aux, visc_rem_u, visc_rem_v, &
                             u_cor_aux, v_cor_aux, OBC, G, CS%maxvel)
      call cpu_clock_end(id_clock_vertvisc)

      call cpu_clock_begin(id_clock_pass)
      call pass_vector(u_cor, v_cor, G%Domain)
      call cpu_clock_end(id_clock_pass)

      call continuity_orig(u_cor, v_cor, hin, h, uh, vh, dt, G, &
                           CS%orig_CSp, uhbt, vhbt, OBC)
    else
      call continuity_orig(u, v, hin, h, uh, vh, dt, G, CS%orig_CSp, uhbt, vhbt, OBC)
    endif
  endif

end subroutine continuity

subroutine continuity_init(Time, G, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ptrs), target, intent(inout) :: diag
  type(continuity_CS),     pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
  character(len=128) :: version = '$Id: GOLD_continuity.F90,v 13.0.2.2.2.6 2010/08/17 18:32:10 rwh Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "GOLD_continuity" ! This module's name.

  if (associated(CS)) then
    call GOLD_error(WARNING, "continuity_init called with associated control structure.")
    return
  endif
  allocate(CS)

  CS%use_continuity_PPM = .false.
  call read_param(param_file,"CONTINUITY_PPM",CS%use_continuity_PPM)
  if (.not.CS%use_continuity_PPM) then
    CS%maxvel = 1.0e6 ; call read_param(param_file,"MAXVEL",CS%maxvel)
  endif

  ! Write all relevant parameters to the model log.
  call log_version(param_file, mod, version, tagname)
  call log_param(param_file, mod, "CONTINUITY_PPM", CS%use_continuity_PPM, &
                 "If true, a positive-definite (or monotonic) piecewise \n"//&
                 "parabolic reconstruction is used for the continuity \n"//&
                 "solver.  The default should be changed to true.", &
                 default=.false.)
  if (.not.CS%use_continuity_PPM) &
    call log_param(param_file, mod, "MAXVEL", CS%maxvel, &
                 "The maximum velocity allowed before the velocity \n"//&
                 "components are truncated.", units="m s-1")

  if (CS%use_continuity_PPM) then
    call continuity_PPM_init(Time, G, param_file, diag, CS%PPM_CSp)
  else
    call continuity_orig_init(Time, G, param_file, diag, CS%orig_CSp)
  endif

  if (.not.CS%use_continuity_PPM) then
    id_clock_pass = cpu_clock_id('(Ocean continuity halo updates)', grain=CLOCK_ROUTINE)
    id_clock_vertvisc = cpu_clock_id('(Ocean continuity vertvisc remean)', grain=CLOCK_ROUTINE)
  endif

end subroutine continuity_init

subroutine continuity_end(CS)
  type(continuity_CS),     pointer       :: CS

  if (CS%use_continuity_PPM) then
    call continuity_PPM_end(CS%PPM_CSp)
  else
    call continuity_orig_end(CS%orig_CSp)
  endif

  deallocate(CS)

end subroutine continuity_end

end module GOLD_continuity
