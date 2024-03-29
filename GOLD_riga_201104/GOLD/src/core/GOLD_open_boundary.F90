module GOLD_open_boundary
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
!*  By Mehmet Ilicak and Robert Hallberg, 2010                         *
!*                                                                     *
!*    This module implements some aspects of internal open boundary    *
!*  conditions in GOLD.                                                *
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

use GOLD_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use GOLD_diag_mediator, only : diag_ptrs, time_type
use GOLD_domains, only : pass_var, pass_vector
use GOLD_error_handler, only : GOLD_mesg, GOLD_error, FATAL, WARNING
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_grid, only : ocean_grid_type
use GOLD_variables, only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use GOLD_variables, only : OBC_FLATHER_E, OBC_FLATHER_W, OBC_FLATHER_N, OBC_FLATHER_S

implicit none ; private

#include <GOLD_memory.h>

public Radiation_Open_Bdry_Conds, open_boundary_init, open_boundary_end

type, public :: open_boundary_CS ; private
  real :: gamma_uv ! The relative weighting for the baroclinic radiation
                   ! velocities (or speed of characteristics) at the
                   ! new time level (1) or the running mean (0) for velocities.
                   ! Valid values range from 0 to 1, with a default of 0.3.
  real :: gamma_h  ! The relative weighting for the baroclinic radiation
                   ! velocities (or speed of characteristics) at the
                   ! new time level (1) or the running mean (0) for thicknesses.
                   ! Valid values range from 0 to 1, with a default of 0.2. 
  real :: rx_max   ! The maximum magnitude of the baroclinic radiation
                   ! velocity (or speed of characteristics), in m s-1.  The
                   ! default value is 10 m s-1.
end type open_boundary_CS

integer :: id_clock_pass

contains

subroutine Radiation_Open_Bdry_Conds(OBC, u_new, u_old, v_new, v_old, &
                                     h_new, h_old, G, CS)
  type(ocean_OBC_type),                pointer       :: OBC
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(inout) :: u_new
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(in)    :: u_old
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(inout) :: v_new
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(in)    :: v_old    
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(inout) :: h_new
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(in)    :: h_old
  type(ocean_grid_type),               intent(inout) :: G
  type(open_boundary_CS),              pointer       :: CS

  real :: dhdt, dhdx, gamma_u, gamma_h, gamma_v
  real :: rx_max, ry_max ! coefficients for radiation
  real :: rx_new, rx_avg ! coefficients for radiation

  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(OBC)) return
  if (.not.(OBC%apply_OBC_u_flather_east .or. OBC%apply_OBC_u_flather_west .or. &
            OBC%apply_OBC_v_flather_north .or. OBC%apply_OBC_v_flather_south)) &
    return
  if (.not.associated(CS)) call GOLD_error(FATAL, &
         "GOLD_open_boundary: Module must be initialized before it is used.")

  gamma_u = CS%gamma_uv ; gamma_v = CS%gamma_uv ; gamma_h = CS%gamma_h
  rx_max = CS%rx_max ; ry_max = CS%rx_max

  if (OBC%apply_OBC_u_flather_east .or. OBC%apply_OBC_u_flather_west) then   
    do k=1,nz ; do j=js,je ; do I=is-1,ie ; if (OBC%OBC_mask_u(I,j)) then
      if (OBC%OBC_kind_u(I,j) == OBC_FLATHER_E) then
        dhdt = u_old(I-1,j,k)-u_new(I-1,j,k) !old-new
        dhdx = u_new(I-1,j,k)-u_new(I-2,j,k) !in new time backward sasha for I-1
        rx_new = 0.0
        if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
        rx_avg = (1.0-gamma_u)*OBC%rx_old_u(I,j,k) + gamma_u*rx_new   
        OBC%rx_old_u(I,j,k) = rx_avg
        u_new(I,j,k) = (u_old(I,j,k) + rx_avg*u_new(I-1,j,k)) / (1.0+rx_avg) 

    !   dhdt = h_old(I,j,k)-h_new(I,j,k) !old-new
    !   dhdx = h_new(I,j,k)-h_new(I-1,j,k) !in new time
    !   rx_new = 0.0
    !   if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
    !   rx_avg = (1.0-gamma_h)*OBC%rx_old_h(I,j,k) + gamma_h*rx_new  
    !   OBC%rx_old_h(I,j,k) = rx_avg	  
    !	  h_new(I+1,j,k) = (h_old(I+1,j,k) + rx_avg*h_new(I,j,k)) / (1.0+rx_avg) !original	  
      endif     
      if (OBC%OBC_kind_u(I,j) == OBC_FLATHER_W) then
        dhdt = u_old(I+1,j,k)-u_new(I+1,j,k) !old-new
        dhdx = u_new(I+1,j,k)-u_new(I+2,j,k) !in new time backward sasha for I+1
        rx_new = 0.0
        if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
        rx_avg = (1.0-gamma_u)*OBC%rx_old_u(I,j,k) + gamma_u*rx_new   
        OBC%rx_old_u(I,j,k) = rx_avg
        u_new(I,j,k) = (u_old(I,j,k) + rx_avg*u_new(I+1,j,k)) / (1.0+rx_avg)

    !   dhdt = h_old(I+1,j,k)-h_new(I+1,j,k) !old-new
    !   dhdx = h_new(I+1,j,k)-h_new(I+2,j,k) !in new time 
    !   rx_new = 0.0
    !   if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
    !   rx_avg = (1.0-gamma_h)*OBC%rx_old_h(I,j,k) + gamma_h*rx_new
    !   OBC%rx_old_h(I,j,k) = rx_avg
    !   h_new(I,j,k) = (h_old(I,j,k) + rx_avg*h_new(I+1,j,k)) / (1.0+rx_avg) !original      
      endif     
    endif ; enddo ; enddo ; enddo      
  endif

  if (OBC%apply_OBC_v_flather_north .or. OBC%apply_OBC_v_flather_south) then   
    do k=1,nz ; do J=js-1,je ; do i=is,ie ; if (OBC%OBC_mask_v(i,J)) then
      if (OBC%OBC_kind_v(i,J) == OBC_FLATHER_N) then
        dhdt = v_old(i,J-1,k)-v_new(i,J-1,k) !old-new
        dhdx = v_new(i,J-1,k)-v_new(i,J-2,k) !in new time backward sasha for J-1
        rx_new = 0.0
        if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
        rx_avg = (1.0-gamma_v)*OBC%ry_old_v(i,J,k) + gamma_v*rx_new   
        OBC%ry_old_v(i,J,k) = rx_avg
        v_new(i,J,k) = (v_old(I,j,k) + rx_avg*v_new(i,J-1,k)) / (1.0+rx_avg) 

    !   dhdt = h_old(i,J,k)-h_new(i,J,k) !old-new
    !   dhdx = h_new(i,J,k)-h_new(i,J-1,k) !in new time
    !   rx_new = 0.0
    !   if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
    !   rx_avg = (1.0-gamma_h)*OBC%ry_old_h(i,J,k) + gamma_h*rx_new 	  
    !   OBC%ry_old_h(i,J,k) = rx_avg	  
    !	  h_new(i,J+1,k) = (h_old(i,J+1,k) + rx_avg*h_new(i,J,k)) / (1.0+rx_avg) !original	  
      endif

      if (OBC%OBC_kind_v(i,J) == OBC_FLATHER_S) then
        dhdt = v_old(i,J+1,k)-v_new(i,J+1,k) !old-new
        dhdx = v_new(i,J+1,k)-v_new(i,J+2,k) !in new time backward sasha for J+1
        rx_new = 0.0
        if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
        rx_avg = (1.0-gamma_v)*OBC%ry_old_v(i,J,k) + gamma_v*rx_new   
        OBC%ry_old_v(i,J,k) = rx_avg
        v_new(i,J,k) = (v_old(I,j,k) + rx_avg*v_new(i,J+1,k)) / (1.0+rx_avg) 

    !   dhdt = h_old(i,J+1,k)-h_new(i,J+1,k) !old-new
    !   dhdx = h_new(i,J+1,k)-h_new(i,J+2,k) !in new time
    !   rx_new = 0.0
    !   if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
    !   rx_avg = (1.0-gamma_h)*OBC%ry_old_h(i,J,k) + gamma_h*rx_new 	  
    !   OBC%ry_old_h(i,J,k) = rx_avg	  
    !	  h_new(i,J,k) = (h_old(i,J,k) + rx_avg*h_new(i,J+1,k)) / (1.0+rx_avg) !original	  
      endif

    endif ; enddo ; enddo ; enddo      
  endif

  call cpu_clock_begin(id_clock_pass)
  call pass_vector(u_new, v_new, G%Domain)
  call pass_var(h_new, G%Domain)
  call cpu_clock_end(id_clock_pass)

end subroutine Radiation_Open_Bdry_Conds

subroutine open_boundary_init(Time, G, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ptrs), target, intent(inout) :: diag
  type(open_boundary_CS),  pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
  character(len=128) :: version = '$Id: GOLD_open_boundary.F90,v 1.1.2.5 2010/08/13 21:54:13 rwh Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "GOLD_open_boundary" ! This module's name.
  logical :: flather_east, flather_west, flather_north, flather_south

  if (associated(CS)) then
    call GOLD_error(WARNING, "continuity_init called with associated control structure.")
    return
  endif

  flather_east = .false.  ; flather_west = .false. 
  flather_north = .false. ; flather_south = .false.
  call read_param(param_file, "APPLY_OBC_U_FLATHER_EAST", flather_east)
  call read_param(param_file, "APPLY_OBC_U_FLATHER_WEST", flather_west)
  call read_param(param_file, "APPLY_OBC_V_FLATHER_NORTH", flather_north)
  call read_param(param_file, "APPLY_OBC_V_FLATHER_SOUTH", flather_south)

  if (.not.(flather_east .or. flather_west .or. flather_north .or. &
            flather_south)) return

  allocate(CS)
  CS%gamma_uv = 0.3 ; CS%gamma_h = 0.2 ; CS%rx_max = 10.0

  call read_param(param_file, "OBC_RADIATION_MAX", CS%rx_max)
  call read_param(param_file, "OBC_RAD_VEL_WT", CS%gamma_uv)
  call read_param(param_file, "OBC_RAD_THICK_WT", CS%gamma_h)

  call log_version(param_file, mod, version, tagname)
  call log_param(param_file, mod, "APPLY_OBC_U_FLATHER_EAST", flather_east, &
                 "If true, some zonal velocity points use Flather open \n"//&
                 "boundary conditions on the east side of the ocean.", &
                 default=.false.)
  call log_param(param_file, mod, "APPLY_OBC_U_FLATHER_WEST", flather_west, &
                 "If true, some zonal velocity points use Flather open \n"//&
                 "boundary conditions on the west side of the ocean.", &
                 default=.false.)
  call log_param(param_file, mod, "APPLY_OBC_V_FLATHER_NORTH", flather_north, &
                 "If true, some meridional velocity points use Flather \n"//&
                 "open boundary conditions on the north side of the ocean.", &
                 default=.false.)
  call log_param(param_file, mod, "APPLY_OBC_V_FLATHER_SOUTH", flather_south, &
                 "If true, some meridional velocity points use Flather \n"//&
                 "open boundary conditions on the north side of the ocean.", &
                 default=.false.)
  call log_param(param_file, mod, "OBC_RADIATION_MAX", CS%rx_max, &
                 "The maximum magnitude of the baroclinic radiation \n"//&
                 "velocity (or speed of characteristics).  This is only \n"//&
                 "used if one of the APPLY_OBC_[UV]_FLATHER_... is true.", &
                 units="m s-1", default=10.0)
  call log_param(param_file, mod, "OBC_RAD_VEL_WT", CS%gamma_uv, &
                 "The relative weighting for the baroclinic radiation \n"//&
                 "velocities (or speed of characteristics) at the new \n"//&
                 "time level (1) or the running mean (0) for velocities. \n"//&
                 "Valid values range from 0 to 1. This is only used if \n"//&
                 "one of the APPLY_OBC_[UV]_FLATHER_...  is true.", &
                 units="nondim",  default=0.3)
  call log_param(param_file, mod, "OBC_RAD_THICK_WT", CS%gamma_h, &
                 "The relative weighting for the baroclinic radiation \n"//&
                 "velocities (or speed of characteristics) at the new \n"//&
                 "time level (1) or the running mean (0) for thicknesses. \n"//&
                 "Valid values range from 0 to 1. This is only used if \n"//&
                 "one of the APPLY_OBC_[UV]_FLATHER_...  is true.", &
                 units="nondim",  default=0.3)

  id_clock_pass = cpu_clock_id('(Ocean OBC halo updates)', grain=CLOCK_ROUTINE)

end subroutine open_boundary_init

subroutine open_boundary_end(CS)
  type(open_boundary_CS), pointer :: CS
  deallocate(CS)
end subroutine open_boundary_end

end module GOLD_open_boundary
