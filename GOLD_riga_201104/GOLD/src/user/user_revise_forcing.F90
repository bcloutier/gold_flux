module user_revise_forcing
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
!*  This module provides a method for updating the forcing fluxes      *
!*  using user-written code without the need to duplicate the          *
!*  extensive code used to create or obtain the fluxes.                *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**
use GOLD_diag_mediator, only : post_data, query_averaging_enabled
use GOLD_diag_mediator, only : register_diag_field, diag_ptrs
use GOLD_domains, only : pass_var, pass_vector, AGRID
use GOLD_error_handler, only : GOLD_error, FATAL, WARNING, is_root_pe
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_grid, only : ocean_grid_type
use GOLD_io, only : file_exists, read_data, vardesc
use GOLD_restart, only : register_restart_field, GOLD_restart_CS
use GOLD_time_manager, only : time_type, operator(+), operator(/), get_time
use GOLD_tracer_flow_control, only : call_tracer_set_forcing
use GOLD_tracer_flow_control, only : tracer_flow_control_CS
use GOLD_variables, only : forcing, surface

implicit none ; private

public user_alter_forcing, user_revise_forcing_init

type, public :: user_revise_forcing_CS ; private
  real    :: cdrag               ! The quadratic bottom drag coefficient.
end type user_revise_forcing_CS

contains

subroutine user_alter_forcing(state, fluxes, day, G, CS)
  type(surface),            intent(in)    :: state
  type(forcing),            intent(inout) :: fluxes
  type(time_type),          intent(in)    :: day
  type(ocean_grid_type),    intent(in)    :: G
  type(user_revise_forcing_CS), pointer   :: CS
! This subroutine sets the surface wind stresses.
!
! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (out)     fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      day - Time of the fluxes.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure returned by a previous
!                 call to surface_forcing_init

end subroutine user_alter_forcing

subroutine user_revise_forcing_init(param_file,CS)
  type(param_file_type),     intent(in) :: param_file
  type(user_revise_forcing_CS), pointer :: CS

  character(len=128) :: version = '$Id: user_revise_forcing.F90,v 1.1.2.5 2008/08/20 06:13:55 rwh Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "user_revise_forcing" ! This module's name.

  call log_version(param_file, mod, version, tagname)

end subroutine user_revise_forcing_init

end module user_revise_forcing
