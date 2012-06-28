module benchmark_initialization
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
!*  The module configures the model for the "DOME" experiment.         *
!*  DOME = Dense Overflow and Mixing Experiment?                       *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use GOLD_error_handler, only : GOLD_mesg, GOLD_error, FATAL, is_root_pe
use GOLD_file_parser, only : read_param, open_param_file, param_file_type
use GOLD_file_parser, only : log_param, log_version
use GOLD_grid, only : ocean_grid_type
use GOLD_tracer, only : add_tracer_OBC_values, advect_tracer_CS
use GOLD_variables, only : thermo_var_ptrs, directories
use GOLD_variables, only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use GOLD_EOS, only : calculate_density, calculate_density_derivs, EOS_type
implicit none ; private

#include <GOLD_memory.h>

public benchmark_initialize_topography
public benchmark_initialize_thickness
public benchmark_init_temperature_salinity

contains

! -----------------------------------------------------------------------------
subroutine benchmark_initialize_topography(D, G, param_file)
  real, intent(out), dimension(NXMEM_,NYMEM_) :: D
  type(ocean_grid_type), intent(in)           :: G
  type(param_file_type), intent(in)           :: param_file
! Arguments: D          - the bottom depth in m. Intent out.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine sets up the DOME topography
  real :: min_depth, max_depth ! The minimum and maximum depths in m.
  real :: PI                   ! 3.1415926... calculated as 4*atan(1)
  real :: D0                   ! A constant to make the maximum     !
                               ! basin depth MAXIMUM_DEPTH.         !
  real :: expdecay = 400000.0  ! A decay scale of associated with   !
                               ! the sloping boundaries, in m.      !
  real :: Dedge = 100.0        ! The depth in m at the basin edge.  !
  real :: south_lat, west_lon, len_lon, len_lat, x, y
  character(len=128) :: version = '$Id: benchmark_initialization.F90,v 1.1.2.4 2011/09/19 16:10:36 Robert.Hallberg Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "benchmark_initialize_topography" ! This subroutine's name.
  character(len=60)  :: axis_units
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call GOLD_mesg("  benchmark_initialization.F90, benchmark_initialize_topography: setting topography", 5)

  min_depth = 0.0 ; call read_param(param_file, "MINIMUM_DEPTH", min_depth)
  call read_param(param_file, "MAXIMUM_DEPTH", max_depth, .true.)
!  Calculate the depth of the bottom.                                !

  ! These expressions force rounding of approximate values in
  ! consistent way.
  south_lat = 0.0 ; call read_param(param_file, "SOUTHLAT", south_lat)
  call read_param(param_file, "LENLAT", len_lat, .true.)
  west_lon = 0.0 ; call read_param(param_file, "WESTLON", west_lon)
  call read_param(param_file, "LENLON", len_lon, .true.)

  PI = 4.0*atan(1.0)
  D0 = max_depth / 0.5;

  do i=is,ie ; do j=js,je
    x=(G%geolonh(i,j)-west_lon)/len_lon
    y=(G%geolath(i,j)-south_lat)/len_lat
!  This sets topography that has a reentrant channel to the south.
    D(i,j) = -D0 * ( y*(1.0 + 0.6*cos(4.0*PI*x)) &
                   + 0.75*exp(-6.0*y) &
                   + 0.05*cos(10.0*PI*x) - 0.7 )
    if (D(i,j) > max_depth) D(i,j) = max_depth
    if (D(i,j) < min_depth) D(i,j) = 0.
  enddo ; enddo

  axis_units = 'degrees' ; call read_param(param_file,"AXIS_UNITS",axis_units)
  call log_version(param_file, mod, version, tagname, "")
  call log_param(param_file, mod, "MAXIMUM_DEPTH", max_depth, &
                 "The maximum depth of the ocean.", units="m")
  call log_param(param_file, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)
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

end subroutine benchmark_initialize_topography
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine benchmark_initialize_thickness(h, G, param_file, eqn_of_state, P_ref)
  real, intent(out), dimension(NXMEM_,NYMEM_, NZ_) :: h
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  type(EOS_type),        pointer    :: eqn_of_state
  real,                                intent(in)  :: P_Ref
! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      eqn_of_state - integer that selects the equatio of state
!  (in)      P_Ref - The coordinate-density reference pressure in Pa.

!  This subroutine initializes layer thicknesses for the DOME experiment
  real :: e0(SZK_(G)+1)   ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: e_pert(SZK_(G)) ! Interface height perturbations, positive     !
                          ! upward, in m.                                !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  real, dimension(SZK_(G)) :: T0, pres, S0, rho_guess, drho, drho_dT, drho_dS
  real :: max_depth, efold, minarg, Tabyss, Tint, arg
  character(len=40)  :: mod = "benchmark_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, nz, itt, nkml

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  nkml = 2

  call GOLD_mesg("  benchmark_initialization.F90, benchmark_initialize_thickness: setting thickness", 5)

  call read_param(param_file, "MAXIMUM_DEPTH", max_depth, .true.)

! <!-- This block calculates T0(k) for purposes of 
    do k=1,nz
      pres(k) = P_Ref ; S0(k) = 35.0
    enddo
    T0(1) = 29.0

    call calculate_density(T0(1),S0(1),pres(1),rho_guess(1),1,1,eqn_of_state)
    call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,1,eqn_of_state)

! A first guess of the layers' temperatures.
    do k=1,nz
      T0(k) = T0(1) + (G%Rlay(k) - rho_guess(1)) / drho_dT(1)
    enddo

! Refine the guesses for each layer.
    do itt = 1,6
      call calculate_density(T0,S0,pres,rho_guess,1,nz,eqn_of_state)
      call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,nz,eqn_of_state)
      do k=1,nz
        T0(k) = T0(k) + (G%Rlay(k) - rho_guess(k)) / drho_dT(k)
      enddo
    enddo
! -->

    do k=1,nz
      e_pert(k) = 0.0
      e0(k) = -1.0*max_depth * real(k-1) / real(nz)
    enddo

    efold = 800.0
    minarg = exp(-(10.0*max_depth-50.0)/efold)
    Tabyss = 0.1
    do j=js,je ; do i=is,ie

!    This sets e_pert, the upward interface height perturbations
!  relative to e0(k).
!      do k=1, NZ ; e_pert(k) = e_here(k) - e0(k) ; enddo
!  This puts in a broad thermocline, based on surface and abyssal
!  temperatures, and an e-folding scale of 800 m.
      do k=2,nkml+1
        e0(k) = -50.0 + EPSILON * (nkml + 1 - k)
      enddo
      do k=nkml+2,nz
        Tint = 0.5*(T0(k-1)+T0(k))
        if (Tint >= T0(1)) then
          e0(k) = e0(k-1)-EPSILON
        else
          if (T0(1)<=Tabyss) e0(k) = -(13.0*max_depth) ! This shouldn't happen!
          arg = (Tint-Tabyss) / (T0(1)-Tabyss)
          if (arg < minarg)then
            e0(k) = -10.0 * max_depth
          else
            e0(k) = -50.0 + efold * log(arg)
          endif
        endif
      enddo

!  The remainder of this subroutine should not be changed.           !

!    This sets the initial thickness (in m) of the layers.  The      !
!  thicknesses are set to insure that: 1.  each layer is at least    !
!  EPSILON thick, and 2.  the interfaces are where they should be    !
!  based on the resting depths and interface height perturbations,   !
!  as long at this doesn't interfere with 1.                         !
      eta1D(nz+1) = -1.0*G%D(i,j)
      do k=nz,1,-1
        eta1D(k) = e0(k) + e_pert(k)
        if (eta1D(k) < (eta1D(k+1) + EPSILON)) then
          eta1D(k) = eta1D(k+1) + EPSILON
          h(i,j,k) = EPSILON
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
    enddo ; enddo

  call log_param(param_file, mod, "MAXIMUM_DEPTH", max_depth, &
                 "The maximum depth of the ocean.", units="m")
end subroutine benchmark_initialize_thickness
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine benchmark_init_temperature_salinity(T, S, G, param_file, &
               eqn_of_state, P_Ref)
  real, dimension(NXMEM_,NYMEM_, NZ_), intent(out) :: T, S
  type(ocean_grid_type),               intent(in)  :: G
  type(param_file_type),               intent(in)  :: param_file
  type(EOS_type),                      pointer     :: eqn_of_state
  real,                                intent(in)  :: P_Ref
!  This function puts the initial layer temperatures and salinities  !
! into T(:,:,:) and S(:,:,:).                                        !

! Arguments: T - The potential temperature that is being initialized.
!  (out)     S - The salinity that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      eqn_of_state - integer that selects the equatio of state
!  (in)      P_Ref - The coordinate-density reference pressure in Pa.
  real :: T0(SZK_(G)), S0(SZK_(G))
  real :: pres(SZK_(G))      ! Reference pressure in kg m-3.             !
  real :: drho_dT(SZK_(G))   ! Derivative of density with temperature in !
                        ! kg m-3 K-1.                               !
  real :: drho_dS(SZK_(G))   ! Derivative of density with salinity in    !
                        ! kg m-3 PSU-1.                             !
  real :: rho_guess(SZK_(G)) ! Potential density at T0 & S0 in kg m-3.   !
  real    :: PI        ! 3.1415926... calculated as 4*atan(1)
  real :: lat, len_lat, south_lat
  character(len=40)  :: mod = "benchmark_init_temperature_salinity" ! This subroutine's name.
  character(len=60)  :: axis_units
  integer :: i, j, k, is, ie, js, je, nz, itt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  south_lat = 0.0 ; call read_param(param_file, "SOUTHLAT", south_lat)
  call read_param(param_file, "LENLAT", len_lat, .true.)

    do k=1,nz
      pres(k) = P_Ref ; S0(k) = 35.0
    enddo

    T0(1) = 29.0
    call calculate_density(T0(1),S0(1),pres(1),rho_guess(1),1,1,eqn_of_state)
    call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,1,eqn_of_state)

! A first guess of the layers' temperatures.                         !
    do k=1,nz
      T0(k) = T0(1) + (G%Rlay(k) - rho_guess(1)) / drho_dT(1)
    enddo

! Refine the guesses for each layer.                                 !
    do itt = 1,6
      call calculate_density(T0,S0,pres,rho_guess,1,nz,eqn_of_state)
      call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,nz,eqn_of_state)
      do k=1,nz
        T0(k) = T0(k) + (G%Rlay(k) - rho_guess(k)) / drho_dT(k)
      enddo
    enddo

    do k=1,nz ; do i=is,ie ; do j=js,je
      T(i,j,k) = T0(k)
      S(i,j,k) = S0(k)
    enddo ; enddo ; enddo
    PI = 4.0*atan(1.0)
    do i=is,ie ; do j=js,je
      lat = 0.5*(G%geolatq(i,j)+G%geolatq(i,j-1))
      T(i,j,1) = 15.0 - 12.0*cos(PI*(lat-south_lat)/(len_lat))
      T(i,j,2) = T(i,j,1)
    enddo ; enddo

  axis_units = 'degrees' ; call read_param(param_file,"AXIS_UNITS",axis_units)
  call log_param(param_file, mod, "SOUTHLAT", south_lat, &
                 "The southern latitude of the domain or the equivalent \n"//&
                 "starting value for the y-axis.", units=axis_units)
  call log_param(param_file, mod, "LENLAT", len_lat, &
                 "The latitudinal or y-direction length of the domain.", &
                 units=axis_units)
end subroutine benchmark_init_temperature_salinity
! -----------------------------------------------------------------------------

end module benchmark_initialization
