module user_surface_forcing
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
!*  Rewritten by Robert Hallberg, June 2009                            *
!*                                                                     *
!*    This file contains the subroutines that a user should modify to  *
!*  to set the surface wind stresses and fluxes of buoyancy or         *
!*  temperature and fresh water.  They are called when the run-time    *
!*  parameters WIND_CONFIG or BUOY_CONFIG are set to "USER".  The      *
!*  standard version has simple examples, along with run-time error    *
!*  messages that will cause the model to abort if this code has not   *
!*  been modified.  This code is intended for use with relatively      *
!*  simple specifications of the forcing.  For more complicated forms, *
!*  it is probably a good idea to read the forcing from input files    *
!*  using "file" for WIND_CONFIG and BUOY_CONFIG.                      *
!*                                                                     *
!*    USER_wind_forcing should set the surface wind stresses (taux and *
!*  tauy) perhaps along with the surface friction velocity (ustar).    *
!*                                                                     *
!*    USER_buoyancy forcing is used to set the surface buoyancy        *
!*  forcing, which may include a number of fresh water flux fields     *
!*  (evap, liq_precip, froz_precip, liq_runoff, froz_runoff, and       *
!*  virt_precip) and the surface heat fluxes (sw, lw, latent and sens) *
!*  if temperature and salinity are state variables, or it may simply  *
!*  be the buoyancy flux if it is not.  This routine also has coded a  *
!*  restoring to surface values of temperature and salinity.           *
!*                                                                     *
!*  Macros written all in capital letters are defined in GOLD_memory.h.*
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, tauy                                  *
!*    j    x ^ x ^ x   At >:  u, taux                                  *
!*    j    > o > o >   At o:  h, fluxes.                               *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**
use GOLD_diag_mediator, only : post_data, query_averaging_enabled
use GOLD_diag_mediator, only : register_diag_field, diag_ptrs
use GOLD_domains, only : pass_var, pass_vector, AGRID
use GOLD_error_handler, only : GOLD_error, FATAL, WARNING, is_root_pe
use GOLD_file_parser, only : read_param, param_file_type, log_version, log_param
use GOLD_grid, only : ocean_grid_type
use GOLD_io, only : file_exists, read_data
use GOLD_time_manager, only : time_type, operator(+), operator(/), get_time
use GOLD_tracer_flow_control, only : call_tracer_set_forcing
use GOLD_tracer_flow_control, only : tracer_flow_control_CS
use GOLD_variables, only : forcing, surface

implicit none ; private

public USER_wind_forcing, USER_buoyancy_forcing, USER_surface_forcing_init

type, public :: user_surface_forcing_CS ; private
  !   This control structure should be used to store any run-time variables
  ! associated with the user-specified forcing.  It can be readily modified
  ! for a specific case, and because it is private there will be no changes
  ! needed in other code (although they will have to be recompiled).
  !   The variables in the cannonical example are used for some common
  ! cases, but do not need to be used.

  logical :: use_temperature ! If true, temperature and salinity are used as
                             ! state variables.
  logical :: bulkmixedlayer  ! If true, a refined bulk mixed layer is used.
  logical :: restorebuoy     ! If true, use restoring surface buoyancy forcing.
  real :: Rho0               !   The density used in the Boussinesq
                             ! approximation, in kg m-3.
  real :: C_p                !   The heat capacity of sea water in J kg-1 K-1.
  real :: G_Earth            !   The gravitational acceleration in m s-2.
  real :: Flux_const         !   The restoring rate at the surface, in m s-1.
  real :: gust_const         !   A constant unresolved background gustiness
                             ! that contributes to ustar, in Pa.

  type(diag_ptrs), pointer :: diag ! A pointer to a structure that is used when
                             ! writing diagnostics.
end type user_surface_forcing_CS

contains

subroutine USER_wind_forcing(state, fluxes, day, G, CS)
  type(surface),                 intent(inout) :: state
  type(forcing),                 intent(inout) :: fluxes
  type(time_type),               intent(in)    :: day
  type(ocean_grid_type),         intent(inout) :: G
  type(user_surface_forcing_CS), pointer       :: CS

!   This subroutine sets the surface wind stresses, fluxes%taux and fluxes%tauy.
! These are the stresses in the direction of the model grid (i.e. the same
! direction as the u- and v- velocities.)  They are both in Pa.
!   In addition, this subroutine can be used to set the surface friction
! velocity, fluxes%ustar, in m s-1. This is needed with a bulk mixed layer.
!
! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (out)     fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      day - Time of the fluxes.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure returned by a previous
!                 call to user_surface_forcing_init

  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq

  !   When modifying the code, comment out this error message.  It is here
  ! so that the original (unmodified) version is not accidentally used.
  call GOLD_error(FATAL, "User_wind_surface_forcing: " // &
     "User forcing routine called without modification." )

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq

  ! Allocate the forcing arrays, if necessary.
  call alloc_if_needed(fluxes%taux, Isdq, Iedq, jsd, jed)
  call alloc_if_needed(fluxes%tauy, isd, ied, Jsdq, Jedq)
  if (CS%bulkmixedlayer) call alloc_if_needed(fluxes%ustar, isd, ied, jsd, jed)

  !  Set the surface wind stresses, in units of Pa.  A positive taux
  !  accelerates the ocean to the (pseudo-)east.

  !  The i-loop extends to is-1 so that taux can be used later in the
  ! calculation of ustar - otherwise the lower bound would be Isq.
  do j=js,je ; do I=is-1,Ieq
    fluxes%taux(I,j) = G%umask(I,j) * 0.0  ! Change this to the desired expression.
  enddo ; enddo
  do J=js-1,Jeq ; do i=is,ie
    fluxes%tauy(i,J) = G%vmask(i,J) * 0.0  ! Change this to the desired expression.
  enddo ; enddo

  !    Set the surface friction velocity, in units of m s-1.  ustar
  !  is always positive.
  if (associated(fluxes%ustar)) then ; do j=js,je ; do i=is,ie
    !  This expression can be changed if desired, but need not be.
    fluxes%ustar(i,j) = G%hmask(i,j) * sqrt(CS%gust_const/CS%Rho0 + &
       sqrt(0.5*(fluxes%taux(I-1,j)**2 + fluxes%taux(I,j)**2) + &
            0.5*(fluxes%tauy(i,J-1)**2 + fluxes%tauy(i,J)**2))/CS%Rho0)
  enddo ; enddo ; endif

end subroutine USER_wind_forcing

subroutine USER_buoyancy_forcing(state, fluxes, day, dt, G, CS)
  type(surface),                 intent(inout) :: state
  type(forcing),                 intent(inout) :: fluxes
  type(time_type),               intent(in)    :: day
  real,                          intent(in)    :: dt
  type(ocean_grid_type),         intent(in)    :: G
  type(user_surface_forcing_CS), pointer       :: CS

!    This subroutine specifies the current surface fluxes of buoyancy or
!  temperature and fresh water.  It may also be modified to add
!  surface fluxes of user provided tracers.

!    When temperature is used, there are long list of fluxes that need to be
!  set - essentially the same as for a full coupled model, but most of these
!  can be simply set to zero.  The net fresh water flux should probably be
!  set in fluxes%evap and fluxes%liq_precip, with any salinity restoring
!  appearing in fluxes%virt_precip, and the other water flux components
!  (froz_precip, liq_runoff and froz_runoff) left as arrays full of zeros. 
!  Evap is usually negative and precip is usually positive.  All heat fluxes
!  are in W m-2 and positive for heat going into the ocean.  All fresh water
!  fluxes are in kg m-2 s-1 and positive for water moving into the ocean.

! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (out)     fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      day_start - Start time of the fluxes.
!  (in)      day_interval - Length of time over which these fluxes
!                           will be applied.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure returned by a previous
!                 call to user_surface_forcing_init

  real :: Temp_restore   ! The temperature that is being restored toward, in C.
  real :: Salin_restore  ! The salinity that is being restored toward, in PSU.
  real :: density_restore  ! The potential density that is being restored
                         ! toward, in kg m-3.
  real :: rhoXcp ! The mean density times the heat capacity, in J m-3 K-1.
  real :: buoy_rest_const  ! A constant relating density anomalies to the
                           ! restoring buoyancy flux, in m5 s-3 kg-1.

  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  !   When modifying the code, comment out this error message.  It is here
  ! so that the original (unmodified) version is not accidentally used.
  call GOLD_error(FATAL, "User_buoyancy_surface_forcing: " // &
    "User forcing routine called without modification." )

  ! Allocate and zero out the forcing arrays, as necessary.  This portion is
  ! usually not changed.
  if (CS%use_temperature) then
    call alloc_if_needed(fluxes%evap, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%liq_precip, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%froz_precip, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%liq_runoff, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%froz_runoff, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%virt_precip, isd, ied, jsd, jed)

    call alloc_if_needed(fluxes%sw, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%lw, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%latent, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%sens, isd, ied, jsd, jed)
  else ! This is the buoyancy only mode.
    call alloc_if_needed(fluxes%buoy, isd, ied, jsd, jed)
  endif


  ! MODIFY THE CODE IN THE FOLLOWING LOOPS TO SET THE BUOYANCY FORCING TERMS.

  if ( CS%use_temperature ) then
    ! Set whichever fluxes are to be used here.  Any fluxes that
    ! are always zero do not need to be changed here.
    do j=js,je ; do i=is,ie
      ! Fluxes of fresh water through the surface are in units of kg m-2 s-1
      ! and are positive downward - i.e. evaporation should be negative.
      fluxes%evap(i,j) = -0.0 * G%hmask(i,j)
      fluxes%liq_precip(i,j) = 0.0 * G%hmask(i,j)

      ! virt_precip will be set later, if it is needed for salinity restoring.
      fluxes%virt_precip(i,j) = 0.0

      !   Heat fluxes are in units of W m-2 and are positive into the ocean.
      fluxes%lw(i,j) = 0.0 * G%hmask(i,j)
      fluxes%latent(i,j) = 0.0 * G%hmask(i,j)
      fluxes%sens(i,j) = 0.0 * G%hmask(i,j)
      fluxes%sw(i,j) = 0.0 * G%hmask(i,j)
    enddo ; enddo
  else ! This is the buoyancy only mode.
    do j=js,je ; do i=is,ie
      !   fluxes%buoy is the buoyancy flux into the ocean in m2 s-3.  A positive
      ! buoyancy flux is of the same sign as heating the ocean.
      fluxes%buoy(i,j) = 0.0 * G%hmask(i,j)
    enddo ; enddo
  endif

  if (CS%restorebuoy) then
    if (CS%use_temperature) then
      call alloc_if_needed(fluxes%heat_restore, isd, ied, jsd, jed)
      !   When modifying the code, comment out this error message.  It is here
      ! so that the original (unmodified) version is not accidentally used.
      call GOLD_error(FATAL, "User_buoyancy_surface_forcing: " // &
        "Temperature and salinity restoring used without modification." )

      rhoXcp = CS%Rho0 * CS%C_p
      do j=js,je ; do i=is,ie
        !   Set Temp_restore and Salin_restore to the temperature (in C) and
        ! salinity (in PSU) that are being restored toward.
        Temp_restore = 0.0
        Salin_restore = 0.0

        fluxes%heat_restore(i,j) = (G%hmask(i,j) * (rhoXcp * CS%Flux_const)) * &
            (Temp_restore - state%SST(i,j)) 
        fluxes%virt_precip(i,j) = - (G%hmask(i,j) * (CS%Rho0*CS%Flux_const)) * &
            ((Salin_restore - state%SSS(i,j)) / &
             (0.5 * (Salin_restore + state%SSS(i,j))))
      enddo ; enddo
    else
      !   When modifying the code, comment out this error message.  It is here
      ! so that the original (unmodified) version is not accidentally used.
      call GOLD_error(FATAL, "User_buoyancy_surface_forcing: " // &
        "Buoyancy restoring used without modification." )

      ! The -1 is because density has the opposite sign to buoyancy.
      buoy_rest_const = -1.0 * (CS%G_Earth * CS%Flux_const) / CS%Rho0
      do j=js,je ; do i=is,ie
       !   Set density_restore to an expression for the surface potential
       ! density in kg m-3 that is being restored toward.
        density_restore = 1030.0

        fluxes%buoy(i,j) = G%hmask(i,j) * buoy_rest_const * &
                          (density_restore - state%Rml(i,j))
      enddo ; enddo
    endif
  endif                                             ! end RESTOREBUOY

end subroutine USER_buoyancy_forcing

subroutine alloc_if_needed(ptr, isd, ied, jsd, jed)
  ! If ptr is not associated, this routine allocates it with the given size
  ! and zeros out its contents.  This is equivalent to safe_alloc_ptr in
  ! GOLD_diag_mediator, but is here so as to be completely transparent.
  real, pointer :: ptr(:,:)
  integer :: isd, ied, jsd, jed
  if (.not.ASSOCIATED(ptr)) then
    allocate(ptr(isd:ied,jsd:jed))
    ptr(:,:) = 0.0
  endif
end subroutine alloc_if_needed

subroutine USER_surface_forcing_init(Time, G, param_file, diag, CS)
  type(time_type),               intent(in) :: Time
  type(ocean_grid_type),         intent(in) :: G
  type(param_file_type),         intent(in) :: param_file
  type(diag_ptrs), target,       intent(in) :: diag
  type(user_surface_forcing_CS), pointer    :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure used for writing diagnostics.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module

  character(len=128) :: version = '$Id: user_surface_forcing.F90,v 1.1.2.2.2.8 2010/09/30 22:53:40 rwh Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "user_surface_forcing" ! This module's name.

  if (associated(CS)) then
    call GOLD_error(WARNING, "USER_surface_forcing_init called with an associated "// &
                             "control structure.")
    return
  endif
  allocate(CS)
  CS%diag => diag

  ! Read in any run-time parameters.
  call read_param(param_file,"TEMPERATURE",CS%use_temperature,.true.)
  CS%bulkmixedlayer = .false. ; CS%restorebuoy = .false.
  call read_param(param_file,"BULKMIXEDLAYER",CS%bulkmixedlayer)
  call read_param(param_file,"RESTOREBUOY",CS%restorebuoy)

  call read_param(param_file,"G_EARTH",CS%G_Earth,.true.)
  if (CS%use_temperature) call read_param(param_file,"C_P",CS%C_p,.true.)
  call read_param(param_file,"RHO_0",CS%Rho0,.true.)

  CS%gust_const = 0.02 ! default value is 0.02 Pa.
  call read_param(param_file,"GUST_CONST",CS%gust_const)
  if (CS%restorebuoy) then
    call read_param(param_file,"FLUXCONST",CS%Flux_const, .true.)
    ! Convert CS%Flux_const from m day-1 to m s-1.
    CS%Flux_const = CS%Flux_const / 86400.0
  endif

  call log_version(param_file, mod, version, tagname, "")
  call log_param(param_file, mod, "TEMPERATURE", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state \n"//&
                 "variables.")
  call log_param(param_file, mod, "BULKMIXEDLAYER", CS%bulkmixedlayer, &
                 "If defined, use a refined Kraus-Turner-like bulk mixed \n"//&
                 "layer with transitional buffer layers.")
  call log_param(param_file, mod, "G_EARTH", CS%G_Earth, &
                 "The gravitational acceleration of the Earth.", units="m s-2")
  call log_param(param_file, mod, "GUST_CONST", CS%gust_const, &
                 "The background gustiness in the winds.", units="Pa", &
                 default=0.02)
  call log_param(param_file, mod, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3")
  if (CS%use_temperature) call log_param(param_file, mod, "C_P", CS%C_p, &
                 "The heat capacity of sea water, approximated as a \n"//&
                 "constant. This is only used if TEMPERATURE is true.", &
                 units="J kg-1 K-1")
  call log_param(param_file, mod, "RESTOREBUOY", CS%restorebuoy, &
                 "If true, the buoyancy fluxes drive the model back \n"//&
                 "toward some specified surface state with a rate \n"//&
                 "given by FLUXCONST.", default= .false.)
  if (CS%restorebuoy) call log_param(param_file,mod, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes \n"//&
                 "to the relative surface anomalies (akin to a piston \n"//&
                 "velocity).  Note the non-MKS units.", units="m day-1")

end subroutine USER_surface_forcing_init

end module user_surface_forcing
