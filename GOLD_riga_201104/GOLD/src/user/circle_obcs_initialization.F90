module circle_obcs_initialization
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

!***********************************************************************
!*                                                                     *
!*  The module configures the model for the "circle_obcs" experiment.  *
!*  circle_obcs = Test of Open Boundary Conditions for an anomaly      *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use GOLD_error_handler, only : GOLD_mesg, GOLD_error, FATAL, is_root_pe
use GOLD_file_parser, only : read_param, open_param_file, param_file_type
use GOLD_grid, only : ocean_grid_type
use GOLD_tracer, only : add_tracer_OBC_values, advect_tracer_CS
use GOLD_variables, only : thermo_var_ptrs, directories, ocean_OBC_type
use GOLD_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use GOLD_file_parser, only : log_param, log_version
implicit none ; private

#include <GOLD_memory.h>

public circle_obcs_initialize_thickness

contains

subroutine circle_obcs_initialize_thickness(h, G, param_file)
  real, intent(out), dimension(NXMEM_,NYMEM_, NZ_) :: h
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

!  This subroutine initializes layer thicknesses for the circle_obcs experiment
  real :: e0(SZK_(G))     ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: e_pert(SZK_(G)) ! Interface height perturbations, positive     !
                          ! upward, in m.                                !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  real :: max_depth  ! The minimum depth in m.
  character(len=128) :: version = '$Id: circle_obcs_initialization.F90,v 1.1.2.3 2011/09/19 16:10:36 Robert.Hallberg Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod   ! This module's name.
  character(len=60) :: axis_units
  integer :: i, j, k, is, ie, js, je, nz
  real :: diskrad, rad, xCenter, xRadius, south_lat, west_lon, len_lat, len_lon, lonC, latC

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call GOLD_mesg("  circle_obcs_initialization.F90, circle_obcs_initialize_thickness: setting thickness", 5)

  ! Parameters read by cartesian grid initialization
  south_lat = 0.0 ; call read_param(param_file, "SOUTHLAT", south_lat, .true.)
  call read_param(param_file, "LENLAT", len_lat, .true.)
  west_lon = 0.0 ; call read_param(param_file, "WESTLON", west_lon)
  call read_param(param_file, "LENLON", len_lon, .true.)
  call read_param(param_file, "MAXIMUM_DEPTH", max_depth, .true.)
  diskrad = 0.0 ; call read_param(param_file, "DISK_RADIUS", diskrad, .true.)

  do k=1,nz
    e0(k) = -max_depth * real(k-1) / real(nz)
  enddo

  ! Uniform thicknesses for base state
  do j=js,je ; do i=is,ie                        !
    eta1D(nz+1) = -1.0*G%D(i,j)
    do k=nz,1,-1
      eta1D(k) = e0(k)
      if (eta1D(k) < (eta1D(k+1) + G%Angstrom_z)) then
        eta1D(k) = eta1D(k+1) + G%Angstrom_z
        h(i,j,k) = G%Angstrom_z
      else
        h(i,j,k) = eta1D(k) - eta1D(k+1)
      endif
    enddo
  enddo ; enddo

  ! Perturb base state by circular anomaly in center
  k=Nz
  latC = south_lat + 0.5*len_lat
  lonC = west_lon + 0.5*len_lon
  do j=js,je ; do i=is,ie
    rad = sqrt((G%geolonh(i,j)-lonC)**2+(G%geolath(i,j)-latC)**2)/(diskrad)
    ! if (rad <= 6.*diskrad) h(i,j,k) = h(i,j,k)+10.0*exp( -0.5*( rad**2 ) )
    rad = min( rad, 1. ) ! Flatten outside radius of diskrad
    rad = rad*(2.*asin(1.)) ! Map 0-1 to 0-pi
    h(i,j,k) = h(i,j,k) + 10.0*0.5*(1.+cos(rad)) ! cosine bell
    h(i,j,k-1) = h(i,j,k-1) - 10.0*0.5*(1.+cos(rad)) ! cosine bell
  enddo ; enddo

  axis_units = 'degrees' ; call read_param(param_file,"AXIS_UNITS",axis_units)

  mod = "circle_obcs_initialize_thickness"
  call log_version(param_file, mod, version, tagname, "")
  call log_param(param_file, mod, "MAXIMUM_DEPTH", max_depth, &
                 "The maximum depth of the ocean.", units="m")
  call log_param(param_file, mod, "SOUTHLAT", south_lat, &
                 "The southern latitude of the domain or the equivalent \n"//&
                 "starting value for the y-axis.", units=axis_units, default=0.)
  call log_param(param_file, mod, "LENLAT", len_lat, &
                 "The latitudinal or y-direction length of the domain.", &
                 units=axis_units)
  call log_param(param_file, mod, "WESTLON", west_lon, &
                 "The western longitude of the domain or the equivalent \n"//&
                 "starting value for the x-axis.", units=axis_units, &
                 default=0.0)
  call log_param(param_file, mod, "LENLON", len_lon, &
                 "The longitudinal or x-direction length of the domain.", &
                 units=axis_units)
  call log_param(param_file, mod, "DISK_RADIUS", diskrad, &
                 "The radius of the initially elevated disk in the \n"//&
                 "circle_obcs test case.", units=axis_units)
end subroutine circle_obcs_initialize_thickness

end module circle_obcs_initialization
