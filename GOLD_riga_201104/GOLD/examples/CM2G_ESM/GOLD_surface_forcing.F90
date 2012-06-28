module GOLD_surface_forcing
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of GOLD.                                         *
!*                                                                     *
!* GOLD is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* GOLD is distributed in the hope that it will be useful, but WITHOUT  *
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
!*  By Robert Hallberg, May 2004                                       *
!*                                                                     *
!*    This program contains the subroutines that transform the surface *
!*  forcing fields for the coupled model and manage the output of      *
!*  these fields.                                                      *
!*                                                                     *
!*  Variables written all in capital letters are defined in GOLD_memory.h.*
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
use GOLD_diag_mediator, only : post_data, query_averaging_enabled, diag_ptrs
use GOLD_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use GOLD_domains, only : pass_vector, global_field_sum, BITWISE_EXACT_SUM
use GOLD_error_handler, only : GOLD_error, WARNING, NOTE, is_root_pe
use GOLD_file_parser, only : read_param, param_file_type
use GOLD_io, only : slasher, write_version_number
use GOLD_variables, only : ocean_grid_type, forcing, surface
!   Forcing is a structure containing pointers to the forcing fields
! which may be used to drive GOLD.  All fluxes are positive downward.
!   Surface is a structure containing pointers to various fields that
! may be used describe the surface state of GOLD.
!   ice_ocean_boundary_type is a structure corresponding to forcing, but with
! the elements, units, and conventions that exactly conform to the use for
! MOM-based coupled models.

use coupler_types_mod, only : coupler_2d_bc_type
use time_interp_external_mod, only : init_external_field, time_interp_external, &
                                     time_interp_external_init
use fms_mod, only : read_data

implicit none ; private

#include <GOLD_memory.h>

public convert_IOB_to_fluxes, surface_forcing_init, average_forcing

type, public :: surface_forcing_CS ; private
  character(len=8) :: wind_stagger ! 'A', 'B', or 'C' to indicate the staggering
                                 ! of the winds that are being provided in calls
                                 ! to update_ocean_model.
  real    :: cdrag               ! The quadratic bottom drag coefficient.
  logical :: use_temperature ! If true, temperature and salinity are used as
                             ! state variables.
  real :: Rho0               !   The density used in the Boussinesq
                             ! approximation, in kg m-3.
  real :: max_p_surf         !   The maximum surface pressure that can be
                             ! exerted by the atmosphere and floating sea-ice,
                             ! in Pa.  This is needed because the FMS coupling
                             ! structure does not limit the water that can be
                             ! frozen out of the ocean and the ice-ocean heat
                             ! fluxes are treated explicitly.
  real :: Flux_const         !   The restoring rate at the surface, in m s-1.
  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
                             ! ocean diagnostic fields and control variables.
  character(len=200) :: inputdir ! The directory where NetCDF input files are.
  integer :: id_taux = -1, id_tauy = -1, id_ustar = -1
  integer :: id_PminusE = -1, id_evap = -1, id_precip = -1
  integer :: id_liq_precip = -1, id_froz_precip = -1
  integer :: id_Net_Heating = -1, id_sw = -1, id_sw_pen = -1, id_LwLatSens = -1, id_buoy = -1
  integer :: id_psurf = -1, id_saltflux = -1, id_TKE_tidal = -1
  integer :: id_srestore = -1  ! An id number for time_interp_external.
  logical :: read_TKE_tidal = .false.
  real, pointer, dimension(:,:) :: TKE_tidal
  real, pointer, dimension(:,:) :: ustar_tidal

end type surface_forcing_CS

type, public :: ice_ocean_boundary_type
  real, pointer, dimension(:,:) :: u_flux =>NULL()   ! i-direction wind stress (Pa)
  real, pointer, dimension(:,:) :: v_flux =>NULL()   ! j-direction wind stress (Pa)
  real, pointer, dimension(:,:) :: t_flux =>NULL()   ! sensible heat flux (W/m2)
  real, pointer, dimension(:,:) :: q_flux =>NULL()   ! specific humidity flux (kg/m2/s)
  real, pointer, dimension(:,:) :: salt_flux =>NULL()! salt flux (kg/m2/s)
  real, pointer, dimension(:,:) :: lw_flux =>NULL()  ! long wave radiation (w/m2)
  real, pointer, dimension(:,:) :: sw_flux_vis_dir => NULL() ! direct visible sw radiation (w/m2)
  real, pointer, dimension(:,:) :: sw_flux_vis_dif => NULL() ! diffuse visible sw radiation (w/m2)
  real, pointer, dimension(:,:) :: sw_flux_nir_dir => NULL() ! direct Near InfraRed sw radiation (w/m2)
  real, pointer, dimension(:,:) :: sw_flux_nir_dif => NULL() ! diffuse Near InfraRed sw radiation (w/m2)
  real, pointer, dimension(:,:) :: lprec =>NULL()    ! mass flux of liquid precip (kg/m2/s)
  real, pointer, dimension(:,:) :: fprec =>NULL()    ! mass flux of frozen precip (kg/m2/s)
  real, pointer, dimension(:,:) :: runoff =>NULL()   ! mass flux of liquid runoff (kg/m2/s)
  real, pointer, dimension(:,:) :: calving =>NULL()  ! mass flux of frozen runoff (kg/m2/s)
  real, pointer, dimension(:,:) :: p =>NULL()        ! pressure of overlying ice and atmosphere
                                                     ! on ocean surface (Pa)
  integer :: xtype                                   ! REGRID, REDIST or DIRECT
  type(coupler_2d_bc_type)      :: fluxes            ! A structure that may contain an
                                                     ! array of named fields used for
                                                     ! passive tracer fluxes.
end type ice_ocean_boundary_type

contains

subroutine convert_IOB_to_fluxes(IOB, fluxes, index_bounds, Time, G, CS, state, restore_salt)
  use GOLD_constants, only : hlv, hlf
  use GOLD_domains, only : BGRID_NE
  type(ice_ocean_boundary_type), intent(in), target :: IOB
  type(forcing),              intent(inout) :: fluxes
  integer, dimension(4),      intent(in)    :: index_bounds
  type(time_type),            intent(in)    :: Time
  type(ocean_grid_type),      intent(inout) :: G
  type(surface_forcing_CS),   pointer       :: CS
  type(surface),              intent(in)    :: state
  logical, optional,          intent(in)    :: restore_salt
! This subroutine translates the Ice_ocean_boundary_type into a
! GOLD forcing type, including changes of units, sign conventions,
! and puting the fields into arrays with GOLD-standard halos.
! Arguments: IOB - the ice-ocean boundary type, containing all the fluxes used
!                  to drive the ocean in a coupled model.
!  (out)     fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      index_bounds - the i- and j- size of the arrays in IOB.
!  (in)      Time - The time of the fluxes, used for interpolating the salinity
!                   to the right time, when it is being restored.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure returned by a previous
!                 call to surface_forcing_init.
!  (in)      state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (in)      restore_salt - if true, salinity is restored to a target value.
  real ALLOCABLE, save :: taux_at_q(NXMEMQ_,NYMEMQ_), tauy_at_q(NXMEMQ_,NYMEMQ_)
  real :: data_srestore(SZ1_(G%D),SZ2_(G%D)), PmE_adj(SZ1_(G%D),SZ2_(G%D))
  real :: work_sum(SZ1_(G%D),SZ2_(G%D))
  real, save :: area_surf
  real :: PmE_adj_total, salt_damp_factor
  real :: IRho0, taux2, tauy2
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, i0, j0
  integer :: isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq
  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd
  logical, save :: first_call = .true.
  logical :: restore_salinity

  isc_bnd = index_bounds(1) ; iec_bnd = index_bounds(2)
  jsc_bnd = index_bounds(3) ; jec_bnd = index_bounds(4)
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq

  IRho0 = 1.0/CS%Rho0

  restore_salinity = .false.
  if (present(restore_salt)) restore_salinity = restore_salt

  if (first_call) then
    call safe_alloc_ptr(fluxes%taux,Isdq,Iedq,jsd,jed)      ; fluxes%taux(:,:) = 0.0
    call safe_alloc_ptr(fluxes%tauy,isd,ied,Jsdq,Jedq)      ; fluxes%tauy(:,:) = 0.0
    call safe_alloc_ptr(fluxes%ustar,isd,ied,jsd,jed)       ; fluxes%ustar(:,:) = 0.0
    call safe_alloc_ptr(fluxes%PminusE,isd,ied,jsd,jed)     ; fluxes%PminusE(:,:) = 0.0
    call safe_alloc_ptr(fluxes%sw,isd,ied,jsd,jed)          ; fluxes%sw(:,:) = 0.0
    call safe_alloc_ptr(fluxes%sw_vis_dir,isd,ied,jsd,jed); fluxes%sw_vis_dir(:,:) = 0.0
    call safe_alloc_ptr(fluxes%sw_vis_dif,isd,ied,jsd,jed); fluxes%sw_vis_dif(:,:) = 0.0
    call safe_alloc_ptr(fluxes%sw_nir_dir,isd,ied,jsd,jed); fluxes%sw_nir_dir(:,:) = 0.0
    call safe_alloc_ptr(fluxes%sw_nir_dif,isd,ied,jsd,jed); fluxes%sw_nir_dif(:,:) = 0.0
    call safe_alloc_ptr(fluxes%sw_pen,isd,ied,jsd,jed)      ; fluxes%sw_pen(:,:) = 0.0    
    call safe_alloc_ptr(fluxes%lw_lat_sens,isd,ied,jsd,jed) ; fluxes%lw_lat_sens(:,:) = 0.0
    call safe_alloc_ptr(fluxes%p_surf,isd,ied,jsd,jed)      ; fluxes%p_surf(:,:) = 0.0
    call safe_alloc_ptr(fluxes%salt_flux,isd,ied,jsd,jed)   ; fluxes%salt_flux(:,:) = 0.0
    call safe_alloc_ptr(fluxes%TKE_tidal,isd,ied,jsd,jed)   ; fluxes%TKE_tidal(:,:) = 0.0
    call safe_alloc_ptr(fluxes%ustar_tidal,isd,ied,jsd,jed) ; fluxes%ustar_tidal(:,:) = 0.0
    call safe_alloc_ptr(fluxes%river,isd,ied,jsd,jed)      ; fluxes%river(:,:) = 0.0        

    ALLOC(taux_at_q(Isdq:Iedq,Jsdq:Jedq)) ; taux_at_q(:,:) = 0.0
    ALLOC(tauy_at_q(Isdq:Iedq,Jsdq:Jedq)) ; tauy_at_q(:,:) = 0.0
    first_call = .false.
    if ((CS%wind_stagger(1:1) == 'b') .or. (CS%wind_stagger(1:1) == 'B')) then
      if (is_root_pe()) call GOLD_error(NOTE,"B-grid input wind stagger")
    else
      if (is_root_pe()) call GOLD_error(NOTE,"C-grid input wind stagger")
    endif

    do j=js,je ; do i=is,ie
      work_sum(i,j) = G%DXDYh(i,j) * G%hmask(i,j)
    enddo ; enddo
    area_surf = global_field_sum(G%Domain%mpp_domain, work_sum, BITWISE_EXACT_SUM)

    do j=js-2,je+2 ; do i=is-2,ie+2

       if (CS%read_TKE_tidal) then
           fluxes%TKE_tidal(i,j) = CS%TKE_tidal(i,j)
           fluxes%ustar_tidal(i,j) = CS%ustar_tidal(i,j)
       endif

    ! This is a hard-coded work put in only at Gibraltar.
      if ((abs(G%geolath(i,j)-35.5) < 0.5) .and. (abs(G%geolonh(i,j)+5.5) < 0.5)) then
        fluxes%TKE_tidal(i,j) = 0.4 ! 2 W/m-2, but in a box 10x wider than Gibraltar.
        fluxes%ustar_tidal(i,j) = sqrt(CS%cdrag) * 1.1 ! Scales as TKE_tidal^1/3
      endif
      if ((abs(G%geolath(i,j)-35.5) < 1.5) .and. (abs(G%geolonh(i,j)+6.5) < 0.5)) then
        fluxes%TKE_tidal(i,j) = 0.02 ! .01 W/m-2 across region.
        fluxes%ustar_tidal(i,j) = sqrt(CS%cdrag) * 0.19
      endif
      if ((abs(G%geolath(i,j)-35.5) < 0.5) .and. (abs(G%geolonh(i,j)+6.5) < 0.5)) then
        fluxes%TKE_tidal(i,j) = 0.05 ! .05 W/m-2 across region.
        fluxes%ustar_tidal(i,j) = sqrt(CS%cdrag) * 0.25
      endif
    enddo ; enddo

  endif

  if (restore_salinity) then
    salt_damp_factor = CS%Flux_const

    call time_interp_external(CS%id_srestore,Time,data_srestore)
    do j=js,je ; do i=is,ie
      if (G%hmask(i,j) > 0.5) then
        pme_adj(i,j) = salt_damp_factor * &
          (state%SSS(i,j) - data_srestore(i,j))/(state%SSS(i,j))
      else
        pme_adj(i,j) = 0.0
      endif
      work_sum(i,j) = G%DXDYh(i,j) * pme_adj(i,j)
    enddo ; enddo
    PmE_adj_total = global_field_sum(G%Domain%mpp_domain,work_sum(:,:), &
                                     BITWISE_EXACT_SUM)/area_surf
  endif

  i0 = is - isc_bnd ; j0 = js - jsc_bnd
  do j=js,je ; do i=is,ie
    if ((CS%wind_stagger(1:1) == 'b') .or. (CS%wind_stagger(1:1) == 'B')) then
      if (ASSOCIATED(IOB%u_flux)) taux_at_q(i,j) = IOB%u_flux(i-i0,j-j0)
      if (ASSOCIATED(IOB%v_flux)) tauy_at_q(i,j) = IOB%v_flux(i-i0,j-j0)
    else ! C-grid wind stresses.
      if (ASSOCIATED(IOB%u_flux)) fluxes%taux(i,j) = IOB%u_flux(i-i0,j-j0)
      if (ASSOCIATED(IOB%v_flux)) fluxes%tauy(i,j) = IOB%v_flux(i-i0,j-j0)
    endif

    fluxes%PminusE(i,j) = 0.0
    if (ASSOCIATED(IOB%lprec)) &
      fluxes%PminusE(i,j) = fluxes%PminusE(i,j) + IOB%lprec(i-i0,j-j0)
    if (ASSOCIATED(IOB%fprec)) &
      fluxes%PminusE(i,j) = fluxes%PminusE(i,j) + IOB%fprec(i-i0,j-j0)
    if (ASSOCIATED(IOB%q_flux)) &
      fluxes%PminusE(i,j) = fluxes%PminusE(i,j) - IOB%q_flux(i-i0,j-j0)
    if (ASSOCIATED(IOB%runoff)) then
         fluxes%PminusE(i,j) = fluxes%PminusE(i,j) + IOB%runoff(i-i0,j-j0)
         fluxes%river(i,j) = IOB%runoff(i-i0,j-j0) * G%hmask(i,j) * Irho0 ! convert from kg m-2 s-1 to m s-1
! NOTE: runoff is contained in PmE. Stored as a separate array so that GOLD_mixed_layer
! can add an additional source of TKE to mix river outflow to simulate estuaries.
    else
      fluxes%river(i,j) = 0.0
    endif
    if (ASSOCIATED(IOB%calving)) then
      fluxes%PminusE(i,j) = fluxes%PminusE(i,j) + IOB%calving(i-i0,j-j0)
      fluxes%river(i,j) = fluxes%river(i,j) + IOB%calving(i-i0,j-j0) * G%hmask(i,j) * Irho0
    endif
    ! Convert from kg m-2 s-1 to m s-1.
    fluxes%PminusE(i,j) = G%hmask(i,j) * Irho0 * fluxes%PminusE(i,j)
    ! The restoring fluxes are already in m s-1.
    if (restore_salinity) fluxes%PminusE(i,j) = &
      fluxes%PminusE(i,j) + (pme_adj(i,j) - PmE_adj_total)*G%hmask(i,j)

    fluxes%LW_lat_sens(i,j) = 0.0
    if (ASSOCIATED(IOB%lw_flux)) &
      fluxes%LW_lat_sens(i,j) = fluxes%LW_lat_sens(i,j) + IOB%lw_flux(i-i0,j-j0)
    if (ASSOCIATED(IOB%t_flux)) &
      fluxes%LW_lat_sens(i,j) = fluxes%LW_lat_sens(i,j) - IOB%t_flux(i-i0,j-j0)
    if (ASSOCIATED(IOB%fprec)) &
      fluxes%LW_lat_sens(i,j) = fluxes%LW_lat_sens(i,j) - IOB%fprec(i-i0,j-j0)*hlf
    if (ASSOCIATED(IOB%calving)) &
      fluxes%LW_lat_sens(i,j) = fluxes%LW_lat_sens(i,j) - IOB%calving(i-i0,j-j0)*hlf
    if (ASSOCIATED(IOB%q_flux)) &
      fluxes%LW_lat_sens(i,j) = fluxes%LW_lat_sens(i,j) - IOB%q_flux(i-i0,j-j0)*hlv
    fluxes%LW_lat_sens(i,j) = G%hmask(i,j) * fluxes%LW_lat_sens(i,j)

    if (ASSOCIATED(IOB%sw_flux_vis_dir)) &
      fluxes%sw_vis_dir(i,j) = G%hmask(i,j) * IOB%sw_flux_vis_dir(i-i0,j-j0)
      fluxes%sw_vis_dif(i,j) = G%hmask(i,j) * IOB%sw_flux_vis_dif(i-i0,j-j0)
      fluxes%sw_nir_dir(i,j) = G%hmask(i,j) * IOB%sw_flux_nir_dir(i-i0,j-j0)
      fluxes%sw_nir_dif(i,j) = G%hmask(i,j) * IOB%sw_flux_nir_dif(i-i0,j-j0)
      fluxes%sw(i,j) = fluxes%sw_vis_dir(i,j) + fluxes%sw_vis_dif(i,j) + &
        fluxes%sw_nir_dir(i,j) + fluxes%sw_nir_dif(i,j)
    if (ASSOCIATED(IOB%salt_flux)) &
      fluxes%salt_flux(i,j) = -G%hmask(i,j) * Irho0 * IOB%salt_flux(i-i0,j-j0)
  enddo ; enddo

  if (ASSOCIATED(IOB%p) .and. (CS%max_p_surf >= 0.0)) then
    do j=js,je ; do i=is,ie
      fluxes%p_surf(i,j) = G%hmask(i,j) * MIN(IOB%p(i-i0,j-j0),CS%max_p_surf)
    enddo ; enddo
  elseif (ASSOCIATED(IOB%p)) then
    do j=js,je ; do i=is,ie
      fluxes%p_surf(i,j) = G%hmask(i,j) * IOB%p(i-i0,j-j0)
    enddo ; enddo
  endif

  if ((CS%wind_stagger(1:1) == 'b') .or. (CS%wind_stagger(1:1) == 'B')) then
    call pass_vector(taux_at_q,tauy_at_q,G%Domain,stagger=BGRID_NE)

    do j=js,je ; do I=Isq,Ieq
      fluxes%taux(I,j) = 0.0
      If ((G%qmask(I,J) + G%qmask(I,J-1)) > 0) &
        fluxes%taux(I,j) = (G%qmask(I,J)*taux_at_q(I,J) + G%qmask(I,J-1)*taux_at_q(I,J-1)) / &
            (G%qmask(I,J) + G%qmask(I,J-1))
    enddo ; enddo

    do J=Jsq,Jeq ; do i=is,ie
      fluxes%tauy(i,J) = 0.0
      if ((G%qmask(I,J) + G%qmask(I-1,J)) > 0) &
        fluxes%tauy(i,J) = (G%qmask(I,J)*tauy_at_q(I,J) + G%qmask(I-1,J)*tauy_at_q(I-1,J)) / &
            (G%qmask(I,J) + G%qmask(I-1,J))
    enddo ; enddo

      ! ustar is required for GOLD's mixed layer formulation.  The background value
      ! of 0.02 Pa is a relatively small value intended to give reasonable behavior
      ! in regions of very weak winds.
    do j=js,je ; do i=is,ie
      if (((G%qmask(I,J) + G%qmask(I-1,J-1)) + (G%qmask(I,J-1) + G%qmask(I-1,J))) > 0) then
        fluxes%ustar(i,j) = sqrt(0.02*Irho0 + Irho0*sqrt( &
          ((G%qmask(I,J)*(taux_at_q(I,J)**2 + tauy_at_q(I,J)**2) + &
            G%qmask(I-1,J-1)*(taux_at_q(I-1,J-1)**2 + tauy_at_q(I-1,J-1)**2)) + &
           (G%qmask(I,J-1)*(taux_at_q(I,J-1)**2 + tauy_at_q(I,J-1)**2) + &
            G%qmask(I-1,J)*(taux_at_q(I-1,J)**2 + tauy_at_q(I-1,J)**2)) ) / &
          ((G%qmask(I,J) + G%qmask(I-1,J-1)) + (G%qmask(I,J-1) + G%qmask(I-1,J))) ))
      else
        fluxes%ustar(i,j) = sqrt(0.02*Irho0)
      endif
    enddo ; enddo
  else ! C-grid wind stresses.
    call pass_vector(fluxes%taux,fluxes%tauy,G%Domain)
    do j=js,je ; do i=is,ie
      taux2 = 0.0
      if ((G%umask(I-1,j) + G%umask(I,j)) > 0) &
        taux2 = (G%umask(I-1,j)*fluxes%taux(I-1,j)**2 + &
                 G%umask(I,j)*fluxes%taux(I,j)**2) / (G%umask(I-1,j) + G%umask(I,j))

      tauy2 = 0.0
      if ((G%vmask(i,J-1) + G%vmask(i,J)) > 0) &
        tauy2 = (G%vmask(i,J-1)*fluxes%tauy(i,J-1)**2 + &
                 G%vmask(i,J)*fluxes%tauy(i,J)**2) / (G%vmask(i,J-1) + G%vmask(i,J))

      fluxes%ustar(i,j) = sqrt(0.02*Irho0 + Irho0*sqrt(taux2 + tauy2) )
    enddo ; enddo
  endif

  !   At a later time, it might prove valuable to translate this array into the
  ! index space of the ocean model, rather than leaving it in the (haloless)
  ! arrays that come in from the surface forcing.
  fluxes%tr_fluxes => IOB%fluxes

end subroutine convert_IOB_to_fluxes


subroutine average_forcing(fluxes, dt, G, CS)
  type(forcing),         intent(in) :: fluxes
  real,                  intent(in) :: dt
  type(ocean_grid_type), intent(in) :: G
  type(surface_forcing_CS), pointer :: CS
!   This subroutine offers forcing fields for time averaging.  These
! fields must first be registered in surface_forcing_init (below).
! This subroutine will typically not be modified, except when new
! forcing fields are added.
!
! Arguments: fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields are unallocated.
!  (in)      dt - The amount of time over which to average.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure returned by a previous
!                 call to surface_forcing_init.

  if (query_averaging_enabled(CS%diag)) then
    if ((CS%id_taux > 0) .and. ASSOCIATED(fluxes%taux)) &
      call post_data(CS%id_taux, fluxes%taux, CS%diag)
    if ((CS%id_tauy > 0) .and. ASSOCIATED(fluxes%tauy)) &
      call post_data(CS%id_tauy, fluxes%tauy, CS%diag)
    if ((CS%id_ustar > 0) .and. ASSOCIATED(fluxes%ustar)) &
      call post_data(CS%id_ustar, fluxes%ustar, CS%diag)

    if ((CS%id_PminusE > 0) .and. ASSOCIATED(fluxes%PminusE)) &
      call post_data(CS%id_PminusE, fluxes%PminusE, CS%diag)
    if ((CS%id_evap > 0) .and. ASSOCIATED(fluxes%evap)) &
      call post_data(CS%id_evap, fluxes%evap, CS%diag)
    if ((CS%id_precip > 0) .and. ASSOCIATED(fluxes%precip)) &
      call post_data(CS%id_precip, fluxes%precip, CS%diag)
    if ((CS%id_liq_precip > 0) .and. ASSOCIATED(fluxes%liq_precip)) &
      call post_data(CS%id_liq_precip, fluxes%liq_precip, CS%diag)
    if ((CS%id_froz_precip > 0) .and. ASSOCIATED(fluxes%froz_precip)) &
      call post_data(CS%id_froz_precip, fluxes%froz_precip, CS%diag)

    if ((CS%id_Net_Heating > 0) .and. ASSOCIATED(fluxes%Net_Heating)) &
      call post_data(CS%id_Net_Heating, fluxes%Net_Heating, CS%diag)
    if ((CS%id_sw > 0) .and. ASSOCIATED(fluxes%sw)) &
      call post_data(CS%id_sw, fluxes%sw, CS%diag)
    if ((CS%id_sw_pen > 0) .and. ASSOCIATED(fluxes%sw_pen)) &
      call post_data(CS%id_sw_pen, fluxes%sw_pen, CS%diag)
    if ((CS%id_LwLatSens > 0) .and. ASSOCIATED(fluxes%lw_lat_sens)) &
      call post_data(CS%id_LwLatSens, fluxes%lw_lat_sens, CS%diag)

    if ((CS%id_psurf > 0) .and. ASSOCIATED(fluxes%p_surf)) &
      call post_data(CS%id_psurf, fluxes%p_surf, CS%diag)
    if ((CS%id_saltflux > 0) .and. ASSOCIATED(fluxes%salt_flux)) &
      call post_data(CS%id_saltflux, fluxes%salt_flux, CS%diag)
    if ((CS%id_TKE_tidal > 0) .and. ASSOCIATED(fluxes%TKE_tidal)) &
      call post_data(CS%id_TKE_tidal, fluxes%TKE_tidal, CS%diag)

    if ((CS%id_buoy > 0) .and. ASSOCIATED(fluxes%buoy)) &
      call post_data(CS%id_buoy, fluxes%buoy, CS%diag)
  endif

end subroutine average_forcing

subroutine surface_forcing_init(Time, G, param_file, diag, CS, restore_salt)
  type(time_type),          intent(in) :: Time
  type(ocean_grid_type),    intent(in) :: G
  type(param_file_type),    intent(in) :: param_file
  type(diag_ptrs), target,  intent(in) :: diag
  type(surface_forcing_CS), pointer    :: CS
  logical, optional,       intent(in) :: restore_salt
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      restore_salt - If present and true, salinity restoring will be
!                           applied in this model.
  real :: C1_3
  character(len=200) :: filename              ! The name of an input file.
  character(len=128) :: version = '$Id: GOLD_surface_forcing.F90,v 1.1.2.1 2008/01/24 15:52:14 wga Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  integer :: isd,ied,jsd,jed
  
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  
  if (associated(CS)) then
    call GOLD_error(WARNING, "surface_forcing_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  call write_version_number (version, tagname)

  CS%diag => diag
  CS%inputdir = "." ; call read_param(param_file,"INPUTDIR",CS%inputdir)
  CS%inputdir = slasher(CS%inputdir)
  call read_param(param_file,"CDRAG",CS%cdrag,.true.)
  call read_param(param_file,"TEMPERATURE",CS%use_temperature,.true.)
  call read_param(param_file,"RHO_0",CS%Rho0,.true.)
  CS%max_p_surf = -1.0 ; call read_param(param_file,"MAX_P_SURF",CS%max_p_surf)
  CS%wind_stagger = 'C'
  call read_param(param_file,"WIND_STAGGER",CS%wind_stagger)
  if (restore_salt) then
    call read_param(param_file,"FLUXCONST",CS%Flux_const, .true.)
    ! Convert CS%Flux_const from m day-1 to m s-1.
    CS%Flux_const = CS%Flux_const / 86400.0
  endif

! ==MJH    
  CS%read_TKE_tidal = .false.
  call read_param(param_file,"READ_TKE_TIDAL",CS%read_TKE_tidal)
  call read_param(param_file,"TKE_TIDAL_FILE",filename)  
  
  if (CS%read_TKE_tidal) then
      
      call safe_alloc_ptr(CS%TKE_tidal,isd,ied,jsd,jed)   ; CS%TKE_tidal(:,:) = 0.0
      call safe_alloc_ptr(CS%ustar_tidal,isd,ied,jsd,jed) ; CS%ustar_tidal(:,:) = 0.0
      
      filename = trim(CS%inputdir) // trim(filename)
      call read_data(filename,'tke_tidal',CS%TKE_tidal,domain=G%domain%mpp_domain,timelevel=1)
      C1_3 = 1.0/3.0
      CS%ustar_tidal(:,:)=sqrt(CS%cdrag)*CS%TKE_tidal(:,:)**C1_3
      
  endif
! ==MJH  
  
  call time_interp_external_init

  CS%id_taux = register_diag_field('ocean_model', 'taux', G%axesu1, Time, &
        'Zonal Wind Stress', 'Pascal')
  CS%id_tauy = register_diag_field('ocean_model', 'tauy', G%axesv1, Time, &
        'Meridional Wind Stress', 'Pascal')
    CS%id_ustar = register_diag_field('ocean_model', 'ustar', G%axesh1, Time, &
        'Surface friction velocity', 'meter second-1')

  if (CS%use_temperature) then
    CS%id_PminusE = register_diag_field('ocean_model', 'PmE', G%axesh1, Time, &
          'Precipitation minus Evaporation', 'meter second-1')
    CS%id_evap = register_diag_field('ocean_model', 'evap', G%axesh1, Time, &
          'Evaporation at ocean surface (usually negative)', 'meter second-1')
    CS%id_precip = register_diag_field('ocean_model', 'precip', G%axesh1, Time, &
          'Precipitation into ocean', 'meter second-1')
    CS%id_froz_precip = register_diag_field('ocean_model', 'froz_precip', G%axesh1, Time, &
          'Frozen Precipitation into ocean', 'meter second-1')
    CS%id_liq_precip = register_diag_field('ocean_model', 'liq_precip', G%axesh1, Time, &
          'Liquid Precipitation into ocean', 'meter second-1')

    CS%id_Net_Heating = register_diag_field('ocean_model', 'Net_Heat', G%axesh1, Time, &
          'Net Surface Heating of Ocean', 'K meter second-1')
    CS%id_sw = register_diag_field('ocean_model', 'SW', G%axesh1, Time, &
        'Shortwave radiation flux into ocean', 'Watt meter-2')
    CS%id_sw_pen = register_diag_field('ocean_model', 'SW_pen', G%axesh1, Time, &
        'Penetrating Shortwave radiation flux into ocean', 'Watt meter-2')
    CS%id_LwLatSens = register_diag_field('ocean_model', 'LwLatSens', G%axesh1, Time, &
          'Combined longwave, latent, and sensible heating', 'Watt meter-2')

    CS%id_psurf = register_diag_field('ocean_model', 'p_surf', G%axesh1, Time, &
          'Pressure at ice-ocean or atmosphere-ocean interface', 'Pascal')
    CS%id_saltflux = register_diag_field('ocean_model', 'salt_flux', G%axesh1, Time, &
          'Salt flux into ocean at surface', 'PSU meter second-1')
    CS%id_TKE_tidal = register_diag_field('ocean_model', 'TKE_tidal', G%axesh1, Time, &
          'Tidal source of BBL mixing', 'Watt meter-2')
  else
    CS%id_buoy = register_diag_field('ocean_model', 'buoy', G%axesh1, Time, &
          'Buoyancy forcing', 'meter2 second-3')
  endif

  if (present(restore_salt)) then ; if (restore_salt) then
    filename = trim(CS%inputdir) // "salt_restore.nc"
    CS%id_srestore = init_external_field(filename,'salt',domain=G%Domain%mpp_domain)
  endif ; endif

end subroutine surface_forcing_init


subroutine surface_forcing_end(CS, fluxes)
  type(surface_forcing_CS), pointer       :: CS
  type(forcing), optional,  intent(inout) :: fluxes
! Arguments:  CS - A pointer to the control structure returned by a previous
!                  call to surface_forcing_init, it will be deallocated here.
!  (inout)    fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.

  if (present(fluxes)) then
    if (associated(fluxes%taux))        deallocate(fluxes%taux)
    if (associated(fluxes%tauy))        deallocate(fluxes%tauy)
    if (associated(fluxes%ustar))       deallocate(fluxes%ustar)
    if (associated(fluxes%buoy))        deallocate(fluxes%buoy)
    if (associated(fluxes%Net_heating)) deallocate(fluxes%Net_heating)
    if (associated(fluxes%sw))          deallocate(fluxes%sw)
    if (associated(fluxes%sw_vis_dir))  deallocate(fluxes%sw_vis_dir)
    if (associated(fluxes%sw_vis_dif))  deallocate(fluxes%sw_vis_dif)
    if (associated(fluxes%sw_nir_dir))  deallocate(fluxes%sw_nir_dir)
    if (associated(fluxes%sw_nir_dif))  deallocate(fluxes%sw_nir_dif)
    if (associated(fluxes%sw_pen))      deallocate(fluxes%sw_pen)
    if (associated(fluxes%lw))          deallocate(fluxes%lw)
    if (associated(fluxes%latent))      deallocate(fluxes%latent)
    if (associated(fluxes%sens))        deallocate(fluxes%sens)
    if (associated(fluxes%lw_lat_sens)) deallocate(fluxes%lw_lat_sens)
    if (associated(fluxes%PminusE))     deallocate(fluxes%PminusE)
    if (associated(fluxes%evap))        deallocate(fluxes%evap)
    if (associated(fluxes%precip))      deallocate(fluxes%precip)
    if (associated(fluxes%liq_precip))  deallocate(fluxes%liq_precip)
    if (associated(fluxes%froz_precip)) deallocate(fluxes%froz_precip)
    if (associated(fluxes%p_surf))      deallocate(fluxes%p_surf)
    if (associated(fluxes%salt_flux))   deallocate(fluxes%salt_flux)
    if (associated(fluxes%TKE_tidal))   deallocate(fluxes%TKE_tidal)
    if (associated(fluxes%ustar_tidal)) deallocate(fluxes%ustar_tidal)
    ! Deallocate any elements of fluxes%tr_fluxes.
    if (associated(fluxes%tr_fluxes))   deallocate(fluxes%tr_fluxes)
  endif

  if (associated(CS)) deallocate(CS)
  CS => NULL()

end subroutine surface_forcing_end

end module GOLD_surface_forcing
