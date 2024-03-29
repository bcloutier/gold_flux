module user_initialization
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
!*  By Robert Hallberg, April 1994 - June 2002                         *
!*                                                                     *
!*    This subroutine initializes the fields for the simulations.      *
!*  The one argument passed to initialize, Time, is set to the         *
!*  current time of the simulation.  The fields which are initialized  *
!*  here are:                                                          *
!*    u - Zonal velocity in m s-1.                                     *
!*    v - Meridional velocity in m s-1.                                *
!*    h - Layer thickness in m.  (Must be positive.)                   *
!*    D - Basin depth in m.  (Must be positive.)                       *
!*    f - The Coriolis parameter, in s-1.                              *
!*    g - The reduced gravity at each interface, in m s-2.             *
!*    Rlay - Layer potential density (coordinate variable) in kg m-3.  *
!*  If TEMPERATURE is defined:                                         *
!*    T - Temperature in C.                                            *
!*    S - Salinity in psu.                                             *
!*  If BULKMIXEDLAYER is defined:                                      *
!*    Rml - Mixed layer and buffer layer potential densities in        *
!*          units of kg m-3.                                           *
!*  If SPONGE is defined:                                              *
!*    A series of subroutine calls are made to set up the damping      *
!*    rates and reference profiles for all variables that are damped   *
!*    in the sponge.                                                   *
!*  Any user provided tracer code is also first linked through this    *
!*  subroutine.                                                        *
!*                                                                     *
!*    Forcing-related fields (taux, tauy, buoy, ustar, etc.) are set   *
!*  in GOLD_surface_forcing.F90.                                       *
!*                                                                     *
!*    These variables are all set in the set of subroutines (in this   *
!*  file) USER_initialize_bottom_depth, USER_initialize_thickness,     *
!*  USER_initialize_velocity,  USER_initialize_temperature_salinity,   *
!*  USER_initialize_mixed_layer_density, USER_initialize_sponges,      *
!*  USER_set_coord, and USER_set_ref_profile.                          *
!*                                                                     *
!*    The names of these subroutines should be self-explanatory. They  *
!*  start with "USER_" to indicate that they will likely have to be    *
!*  modified for each simulation to set the initial conditions and     *
!*  boundary conditions.  Most of these take two arguments: an integer *
!*  argument specifying whether the fields are to be calculated        *
!*  internally or read from a NetCDF file; and a string giving the     *
!*  path to that file.  If the field is initialized internally, the    *
!*  path is ignored.                                                   *
!*                                                                     *
!*  Macros written all in capital letters are defined in GOLD_memory.h.*
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, f                                     *
!*    j+1  > o > o >   At ^:  v, tauy                                  *
!*    j    x ^ x ^ x   At >:  u, taux                                  *
!*    j    > o > o >   At o:  h, D, buoy, tr, T, S, Rml, ustar         *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_CompressComp, only : read_compress, register_compress
use GOLD_CompressComp, only : compress_init, Compress_CS
use GOLD_error_handler, only : GOLD_mesg, GOLD_error, FATAL, is_root_pe
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_grid, only : ocean_grid_type
use GOLD_io, only : close_file, create_file, fieldtype, file_exists
use GOLD_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
use GOLD_io, only : write_field, slasher, vardesc
use GOLD_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use GOLD_tracer, only : add_tracer_OBC_values, advect_tracer_CS
use GOLD_variables, only : thermo_var_ptrs, directories
use GOLD_variables, only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use GOLD_variables, only : OBC_FLATHER_E, OBC_FLATHER_W, OBC_FLATHER_N, OBC_FLATHER_S
use GOLD_EOS, only : calculate_density, calculate_density_derivs, EOS_type
implicit none ; private

#include <GOLD_memory.h>

public USER_set_coord, USER_initialize_topography, USER_initialize_thickness
public USER_initialize_velocity, USER_init_temperature_salinity
public USER_init_mixed_layer_density, USER_initialize_sponges
public USER_set_Open_Bdry_Conds, USER_set_ref_profile, USER_set_rotation

logical :: first_call = .true.

contains

subroutine USER_set_coord(Rlay, g_prime, G, param_file, eqn_of_state)
  real, dimension(:), intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  type(EOS_type),        pointer    :: eqn_of_state

  call GOLD_error(FATAL, &
   "USER_initialization.F90, USER_set_coord: " // &
   "Unmodified user routine called - you must edit the routine to use it")
  Rlay(:) = 0.0
  g_prime(:) = 0.0
  
  if (first_call) call write_user_log(param_file)

end subroutine USER_set_coord

subroutine USER_initialize_topography(D, G, param_file)
  real, dimension(NXMEM_,NYMEM_), intent(out) :: D
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  call GOLD_error(FATAL, &
   "USER_initialization.F90, USER_initialize_topography: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  D(:,:) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_topography

subroutine USER_initialize_thickness(h, G, param_file, T)
  real, intent(out), dimension(NXMEM_,NYMEM_, NZ_) :: h
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  real, intent(in), dimension(NXMEM_,NYMEM_, NZ_)  :: T
  call GOLD_error(FATAL, &
   "USER_initialization.F90, USER_initialize_thickness: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  h(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_thickness

subroutine USER_initialize_velocity(u, v, G, param_file)
  real, dimension(NXMEMQ_,NYMEM_, NZ_), intent(out) :: u
  real, dimension(NXMEM_,NYMEMQ_, NZ_), intent(out) :: v
  type(ocean_grid_type),                intent(in)  :: G
  type(param_file_type),                intent(in)  :: param_file
  call GOLD_error(FATAL, &
   "USER_initialization.F90, USER_initialize_velocity: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  u(:,:,1) = 0.0
  v(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_velocity

subroutine USER_init_temperature_salinity(T, S, G, param_file, eqn_of_state)
  real, dimension(NXMEM_,NYMEM_, NZ_), intent(out) :: T, S
  type(ocean_grid_type),               intent(in)  :: G
  type(param_file_type),               intent(in)  :: param_file
  type(EOS_type),                      pointer     :: eqn_of_state
  call GOLD_error(FATAL, &
   "USER_initialization.F90, USER_init_temperature_salinity: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  T(:,:,1) = 0.0
  S(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_init_temperature_salinity

subroutine USER_init_mixed_layer_density(Rml, G, param_file, use_temperature, &
                                         eqn_of_state, T, S, P_Ref)
  real, dimension(NXMEM_,NYMEM_, NZ_),          intent(out) :: Rml
  type(ocean_grid_type),                        intent(in)  :: G
  type(param_file_type),                        intent(in)  :: param_file
  logical,                                      intent(in)  :: use_temperature
  type(EOS_type),                      optional, pointer    :: eqn_of_state
  real, dimension(NXMEM_,NYMEM_, NZ_), optional, intent(in) :: T, S
  real,                                optional, intent(in) :: P_Ref
  call GOLD_error(FATAL, &
   "USER_initialization.F90, USER_init_mixed_layer_density: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  Rml(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_init_mixed_layer_density

subroutine USER_initialize_sponges(G, use_temperature, bulkmixedlayer, tv, &
                                   param_file, CSp, h)
  type(ocean_grid_type), intent(in) :: G
  logical, intent(in) :: use_temperature, bulkmixedlayer
  type(thermo_var_ptrs), intent(in) :: tv
  type(param_file_type), intent(in) :: param_file
  type(sponge_CS),       pointer    :: CSp
  real, intent(in), dimension(NXMEM_,NYMEM_, NZ_) :: h
  call GOLD_error(FATAL, &
   "USER_initialization.F90, USER_initialize_sponges: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_sponges

subroutine USER_set_Open_Bdry_Conds(OBC, tv, G, param_file, advect_tracer_CSp)
  type(ocean_OBC_type),  pointer    :: OBC
  type(thermo_var_ptrs), intent(in) :: tv
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  type(advect_tracer_CS), pointer   :: advect_tracer_CSp
  call GOLD_error(FATAL, &
   "USER_initialization.F90, USER_set_Open_Bdry_Conds: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_set_Open_Bdry_Conds

subroutine USER_set_ref_profile(new_sim, G, param_file, Compress_CSp, compress_file)
  logical, intent(in):: new_sim
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in) :: param_file
  type(Compress_CS), pointer   :: Compress_CSp
  character(len=*), intent(in) :: compress_file
  call GOLD_error(FATAL, &
   "USER_initialization.F90, USER_set_ref_profile: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_set_ref_profile

subroutine USER_set_rotation(G, param_file)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in) :: param_file
  call GOLD_error(FATAL, &
   "USER_initialization.F90, USER_set_rotation: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_set_rotation

subroutine write_user_log(param_file)
  type(param_file_type), intent(in) :: param_file

  character(len=128) :: version = '$Id: user_initialization.F90,v 1.1.2.11.2.10 2010/05/28 13:20:30 rwh Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "user_initialization" ! This module's name.

  call log_version(param_file, mod, version, tagname)
  first_call = .false.

end subroutine write_user_log

end module user_initialization
