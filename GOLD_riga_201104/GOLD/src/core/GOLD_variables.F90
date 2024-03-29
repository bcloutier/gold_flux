module GOLD_variables

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

use GOLD_domains, only : GOLD_domain_type, get_domain_extent
use GOLD_checksums, only : hchksum, qchksum, uchksum, vchksum
use GOLD_error_handler, only : GOLD_error, FATAL
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_grid, only : ocean_grid_type, GOLD_variables_init => GOLD_grid_init
use GOLD_io, only : vardesc
use GOLD_EOS, only : EOS_type

use coupler_types_mod, only : coupler_2d_bc_type

implicit none ; private

#include <GOLD_memory.h>

public GOLD_variables_init, GOLD_forcing_chksum, GOLD_thermovar_chksum
public ocean_grid_type, vardesc, alloc_BT_cont_type, dealloc_BT_cont_type

type, public :: p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d
type, public :: p2d
  real, dimension(:,:), pointer :: p => NULL()
end type p2d

type, public :: directories
  character(len=120) :: &
    restart_input_dir = ' ',& ! The directory to read restart and input files.
    restart_output_dir = ' ',&! The directory into which to write restart files.
    output_directory = ' ', & ! The directory to use to write the model output.
    input_filename  = ' '     ! A string that indicates the input files or how
                              ! the run segment should be started.
end type directories

!   The following structure contains pointers to the forcing fields
! which may be used to drive GOLD.  All fluxes are positive downward.
! Pointers to unused fluxes should be set to NULL.
type, public :: forcing
  real, pointer, dimension(:,:) :: &
    taux => NULL(), &       ! The zonal wind stress, in Pa.
    tauy => NULL(), &       ! The meridional wind stress, in Pa.
    ustar => NULL(), &      ! The surface friction velocity, in units of m s-1.
    buoy => NULL(), &       ! The buoyancy flux into the ocean in m2 s-3.

    sw => NULL(), &         ! The shortwave heat flux into the ocean, in W m-2.
    sw_vis_dir => NULL(), & ! The visible, direct shortwave heat flux into the
                            ! ocean, in W m-2.
    sw_vis_dif => NULL(), & ! The visible, diffuse shortwave heat flux into the
                            ! ocean, in W m-2.
    sw_nir_dir => NULL(), & ! The near-IR, direct shortwave heat flux into the
                            ! ocean, in W m-2.
    sw_nir_dif => NULL(), & ! The near-IR, diffuse shortwave heat
                            ! flux into the ocean, in W m-2.
    lw => NULL(), &         ! The longwave heat flux into the ocean, in W m-2.
                            ! This field is typically negative.
    latent => NULL(), &     ! The latent heat flux into the ocean, in W m-2.
                            ! This field is typically negative.
    sens => NULL(), &       ! The sensible heat flux into the ocean, in W m-2.
                            ! This field is typically negative.
    heat_restore => NULL(), & ! The heat flux into the ocean from temperature
                            ! restoring, in W m-2.
                            ! This field is typically negative.

    evap => NULL(), &       ! The negative of the fresh water flux out of the
                            ! ocean, in kg m-2 s-1.
    liq_precip => NULL(), & ! The liquid water flux into the ocean, in
                            ! kg m-2 s-1.
    froz_precip => NULL(), &  ! The frozen water flux into the ocean,
                            ! in kg m-2 s-1.
    virt_precip => NULL(), & ! The virtual water flux into the ocean associated
                            ! with salinity restoring, in kg m-2 s-1.
    runoff_hflx => NULL(), & ! Heat flux associated with liq_runoff in W m-2.
    calving_hflx => NULL(), & ! Heat flux associated with froz_runoff in W m-2.

    p_surf_full => NULL(), & ! The pressure at the top ocean interface, in Pa.
                             ! If there is sea-ice, this is at the interface
                             ! between the ice and ocean.
    p_surf => NULL(), &      ! The pressure at the top ocean interface, in Pa,
                             ! as used to drive the ocean model.  If p_surf is
                             ! limited, this may be smaller than p_surf_full,
                             ! otherwise they are the same.
    salt_flux => NULL(), &   ! The net salt flux into the ocean in kg Salt m-2 s-1.
    TKE_tidal => NULL(), &   ! The tidal source of energy driving mixing in the
                             ! bottom boundary layer, in W m-2.
    ustar_tidal => NULL(), & ! The tidal contribution to bottom ustar, in m s-1.
    liq_runoff => NULL(), &  ! Mass of river runoff in units of kg m-2 s-1.
    froz_runoff => NULL(), & ! Mass of calving in units of kg m-2 s-1.

    ustar_shelf => NULL(), & ! The friction velocity under ice-shelves in m s-1.
                             ! This was calculated by the ocean the previous
                             ! time step.
    frac_shelf_h => NULL(), &! Fractional ice shelf coverage of h-, u-, and v-
    frac_shelf_u => NULL(), &! cells, nondimensional from 0 to 1. These are only
    frac_shelf_v => NULL(), &! associated if ice shelves are enabled, and are
                             ! exactly 0 away from shelves or on land.
    rigidity_ice_u => NULL(),& ! The depth-integrated lateral viscosity of
    rigidity_ice_v => NULL()   ! ice shelves at u- or v-points, in m3 s-1.
  type(coupler_2d_bc_type), pointer :: tr_fluxes  => NULL()
                                            ! A structure that may contain an
                                            ! array of named fields used for
                                            ! passive tracer fluxes.
       !!! NOTE: ALL OF THE ARRAYS IN TR_FLUXES USE THE COUPLER'S INDEXING
       !!!       CONVENTION AND HAVE NO HALOS!  THIS IS DONE TO CONFORM TO
       !!!       THE TREATMENT IN MOM4, BUT I DON'T LIKE IT!
end type forcing

!   The following structure contains pointers to various fields
! which may be used describe the surface state of GOLD, and which
! will be returned to a the calling program
type, public :: surface
  real, pointer, dimension(:,:) :: &
    SST => NULL(), &     ! The sea surface temperature in C.
    SSS => NULL(), &     ! The sea surface salinity in psu.
    Rml => NULL(), &     ! The mixed layer density in kg m-3.
    u => NULL(), &       ! The mixed layer zonal velocity in m s-1.
    v => NULL(), &       ! The mixed layer meridional velocity in m s-1.
    Hml => NULL(), &     ! The mixed layer depth in m.
    ocean_mass => NULL(), &  ! The total mass of the ocean in kg m-2.
    ocean_heat => NULL(), &  ! The total heat content of the ocean in C kg m-2.
    ocean_salt => NULL(), &  ! The total salt content of the ocean in kgSalt m-2.
    taux_shelf => NULL(), &  ! The zonal and meridional stresses on the ocean
    tauy_shelf => NULL(), &  ! under shelves, in Pa.
    frazil => NULL(), &  ! The energy needed to heat the ocean column to the
                         ! freezing point over the call to step_GOLD, in J m-2.
    salt_deficit => NULL(), & ! The salt needed to maintain the ocean column
                         ! at a minimum salinity of 0.01 PSU over the call to
                         ! step_GOLD, in kgSalt m-2.
    TempxPmE => NULL(), &  ! The net inflow of water into the ocean times
                         ! the temperature at which this inflow occurs during
                         ! the call to step_GOLD, in deg C kg m-2.
                         !   This should be prescribed in the forcing fields,
                         ! but as it often is not, this is a useful heat budget
                         ! diagnostic.
    internal_heat => NULL() , & ! Any internal or geothermal heat sources that
                         ! are applied to the ocean integrated over the call
                         ! to step_GOLD, in deg C kg m-2.
    sea_lev => NULL()    ! The sea level in m.  If a reduced  surface gravity is
                         ! used, that is compensated for in sea_lev.
  type(coupler_2d_bc_type), pointer :: tr_fields  => NULL()
                                          ! A structure that may contain an
                                          ! array of named fields describing
                                          ! tracer-related quantities.
       !!! NOTE: ALL OF THE ARRAYS IN TR_FIELDS USE THE COUPLER'S INDEXING
       !!!       CONVENTION AND HAVE NO HALOS!  THIS IS DONE TO CONFORM TO
       !!!       THE TREATMENT IN MOM4, BUT I DON'T LIKE IT!
end type surface

!   The following structure contains pointers to an assortment of
! thermodynamic fields that may be available, including potential
! temperature, salinity and mixed layer density.
type, public :: thermo_var_ptrs
!   If allocated, the following variables have nz layers.
  real, pointer :: T(:,:,:) => NULL()   ! Potential temperature in C.
  real, pointer :: S(:,:,:) => NULL()   ! Salnity in psu.
  type(EOS_type), pointer :: eqn_of_state => NULL() ! Type that indicates the
                                        ! equation of state to use.
!   If allocated, the Rml has nz layers, but only nk_Rml are typically used.
  real, pointer :: Rml(:,:,:) => NULL() !   The mixed and buffer layer
                                        ! coordinate-reference potential
                                        ! densities in kg m-3.
  integer :: nk_Rml = 0                 !   The number of layers in Rml that
                                        ! are actually variable.
  real :: P_Ref                         !   The coordinate-density reference
                                        ! pressure in Pa.  This is the
                                        ! pressure used to calculate Rml
                                        ! from T and S when eqn_of_state
                                        ! is associated.
  real, pointer, dimension(:,:) :: &
    Hml => NULL(), &     !   The surface mixed layer depth in m.                        

!  These arrays are accumulated fluxes for communication with other components.
    frazil => NULL(), &  !   The energy needed to heat the ocean column to the
                         ! freezing point since calculate_surface_state was
                         ! last called, in units of J m-2.
    salt_deficit => NULL(), & !   The salt needed to maintain the ocean column
                         ! at a minumum salinity of 0.01 PSU since the last time
                         ! that calculate_surface_state was called, in units
                         ! of gSalt m-2.
    TempxPmE => NULL(), &!   The net inflow of water into the ocean times the
                         ! temperature at which this inflow occurs since the
                         ! last call to calculate_surface_state, in units of
                         ! deg C kg m-2. This should be prescribed in the
                         ! forcing fields, but as it often is not, this is a
                         ! useful heat budget diagnostic.
    internal_heat => NULL() ! Any internal or geothermal heat sources that
                         ! have been applied to the ocean since the last call to
                         ! calculate_surface_state, in units of deg C kg m-2.
end type thermo_var_ptrs

type, public :: ocean_internal_state
! This structure contains pointers to all of the prognostic variables allocated
! here and in GOLD.F90.  It is useful for sending them for diagnostics, and for
! preparation for ensembles later on.  All variables have the same names as the
! local (public) variables they refer to in GOLD.F90.
  real, pointer, dimension(:,:,:,:) :: &
    u => NULL(), v => NULL(), h => NULL()
  real, pointer, dimension(:,:,:) :: &
    uh => NULL(), vh => NULL(), CAu => NULL(), CAv => NULL(), &
    PFu  => NULL(), PFv => NULL(), diffu => NULL(), diffv => NULL(), &
    T => NULL(), S => NULL(), pbce => NULL(), eta => NULL(), &
    u_accel_bt => NULL(), v_accel_bt => NULL(), u_av => NULL(), v_av => NULL()
end type ocean_internal_state

type, public :: vertvisc_type
!   This structure contains vertical viscosities, drag coefficients, and
! related fields.
  real, pointer, dimension(:,:) :: &
    bbl_thick_u => NULL(), & ! The bottom boundary layer thickness at the zonal
    bbl_thick_v => NULL(), & ! and meridional velocity points, in m.
    kv_bbl_u => NULL(), &    ! The bottom boundary layer viscosity at the zonal
    kv_bbl_v => NULL(), &    ! and meridional velocity points, in m2 s-1.
    ustar_BBL => NULL(), &   ! The turbulence velocity in the bottom boundary
                             ! layer at h points, in m s-1.
    TKE_BBL => NULL(), &     ! A term related to the bottom boundary layer
                             ! source of turbulent kinetic energy, currently
                             ! in units of m3 s-3, but will later be changed
                             ! to W m-2.
    taux_shelf => NULL(), &  ! The zonal and meridional stresses on the ocean
    tauy_shelf => NULL(), &  ! under shelves, in Pa.
    tbl_thick_shelf_u => NULL(), & ! Thickness of the viscous top boundary
    tbl_thick_shelf_v => NULL(), & ! layer under ice shelves, in m.
    kv_tbl_shelf_u => NULL(), &  ! Viscosity in the viscous top boundary
    kv_tbl_shelf_v => NULL(), &  ! layer under ice shelves, in m2 s-1.
    nkml_visc_u => NULL(), & ! The number of layers in the viscous surface
    nkml_visc_v => NULL()    ! mixed layer.  These are not integers because
                             ! there may be fractional layers.  This is done
                             ! in terms of layers, not depth, to facilitate
                             ! the movement of the viscous boundary layer with
                             ! the flow.
  real, pointer, dimension(:,:,:) :: &
    Ray_u => NULL(), &  ! The Rayleigh drag velocity to be applied to each layer
    Ray_v => NULL(), &  ! at u- and v-points, in m s-1.
    Kd_extra_T => NULL(), & ! The extra diffusivities of temperature and
    Kd_extra_S => NULL(), & ! salinity due to double diffusion at the interfaces
                        ! relative to the diffusivity of density, in m2 s-1.  
                        ! One of these is always 0.  Kd_extra_S is positive for
                        ! salt fingering; Kd_extra_T is positive for double
                        ! diffusive convection.  These are only allocated if
                        ! DOUBLE_DIFFUSION is true.
    Kd_turb => NULL(), &! The turbulent diapycnal diffusivity at the interfaces
                        ! between each layer, in m2 s-1.
    TKE_turb => NULL()  ! The turbulent kinetic energy per unit mass defined
                        ! at the interfaces between each layer, in m2 s-2.
end type vertvisc_type

integer, parameter, public :: OBC_NONE = 0, OBC_SIMPLE = 1, OBC_WALL = 2
integer, parameter, public :: OBC_FLATHER_E = 4, OBC_FLATHER_W = 5
integer, parameter, public :: OBC_FLATHER_N = 6, OBC_FLATHER_S = 7

type, public :: ocean_OBC_type
! This structure is used to apply specified open boundary conditions.
!!!!!!!!!!!!!!!!!!!!!!! Mehmet !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical :: apply_OBC_u_flather_east = .false.  ! If true, some zonal velocity
  logical :: apply_OBC_u_flather_west = .false.  ! points in the local domain use flather open
                                                 ! boundary conditions.
  logical :: apply_OBC_v_flather_north = .false.  ! If true, some meridional velocity
  logical :: apply_OBC_v_flather_south = .false.  ! points in the local domain use flather open
                                                 ! boundary conditions.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical :: apply_OBC_u = .false.  ! If true, some zonal or meridional velocity
  logical :: apply_OBC_v = .false.  ! points in the local domain use open
                                    ! boundary conditions.
  logical, pointer, dimension(:,:) :: &
    OBC_mask_u => NULL(), & ! These arrays are true at zonal or meridional
    OBC_mask_v => NULL()    ! velocity points that have prescribed open boundary
                            ! conditions.
  integer, pointer, dimension(:,:) :: &
    OBC_kind_u => NULL(), & ! These arrays indicate the kind of open boundary
    OBC_kind_v => NULL()    ! conditions that are to be applied at the u and v
                            ! points, and can be OBC_NONE, OBC_SIMPLE, OBC_WALL,
                            ! or one of OBC_FLATHER_[EWNS].  Generally these
                            ! should be consistent with OBC_mask_[uv], with
                            ! OBC_mask_[uv] .false. for OBC_kind_[uv] = NONE
                            ! and true for all other values.
!!!!!!!!!!!!!!!!!!!!!!! Mehmet !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
  ! The following apply at points with OBC_kind_[uv] = OBC_FLATHER_x.
  real, pointer, dimension(:,:,:) :: &
    rx_old_u => NULL(), &  ! The rx_old_u value for radition coeff for u-velocity in x-direction
    ry_old_v => NULL(), &  ! The ry_old_v value for radition coeff for v-velocity in y-direction
    rx_old_h => NULL(), &  ! The rx_old_h value for radition coeff for layer thickness h in x-direction
    ry_old_h => NULL()     ! The ry_old_h value for radition coeff for layer thickness h in y-direction

  !   The following can be used to specify the outer-domain values of the
  ! surface height and barotropic velocity.  If these are not allocated, the
  ! default with Flather boundary conditions is the same as if they were
  ! filled with zeros.  With simple OBCs, these should not be allocated.
  real, pointer, dimension(:,:) :: &
    ubt_outer => NULL(), &    ! The u-velocity in the outer domain, in m s-1.
    vbt_outer => NULL(), &    ! The v-velocity in the outer domain, in m s-1.
    eta_outer_u => NULL(), &  ! The sea surface height anomaly or water column
    eta_outer_v => NULL()     ! mass anomaly in the outer domain in m or kg m-2.

  ! The following apply at points with OBC_kind_[uv] = OBC_SIMPLE.
  real, pointer, dimension(:,:,:) :: &
    u => NULL(), &  ! The prescribed values of the zonal (u) or meridional (v)
    v => NULL(), &  ! velocities at OBC points, in m s-1.
    uh => NULL(), & ! The prescribed values of the zonal (uh) or meridional (vh)
    vh => NULL()    ! volume transports at OBC points, in m3 s-1.
end type ocean_OBC_type

type, public :: optics_type
  integer :: nbands    ! The number of penetrating bands of SW radiation.
  real, pointer, dimension(:,:,:,:) :: &
    opacity_band => NULL()  ! The SW optical depth per unit thickness, m-1.
                       ! The number of bands of radiation is the most rapidly
                       ! varying (first) index.
  real, pointer, dimension(:,:,:) :: &
    SW_pen_band => NULL()  ! The shortwave radiation at the surface in each of
                       ! the nbands bands that penetrates beyond the surface,
                       ! in W m-2. The most rapidly varying dimension is the
                       ! band.
  real, pointer, dimension(:) :: &
    min_wavelength_band => NULL(), & ! The range of wavelengths in each band of
    max_wavelength_band => NULL()    ! penetrating shortwave radiation, in nm.
end type optics_type

type, public :: BT_cont_type
  real, pointer, dimension(:,:) :: &
    FA_u_EE => NULL(), &  ! The FA_u_XX variables are the effective open face
    FA_u_E0 => NULL(), &  ! areas for barotropic transport through the zonal
    FA_u_W0 => NULL(), &  ! faces, all in H m, with the XX indicating where
    FA_u_WW => NULL(), &  ! the transport is from, with _EE drawing from points
                          ! far to the east, _E0 from points nearby from the
                          ! east, _W0 nearby from the west, and _WW from far to
                          ! the west.
    uBT_WW => NULL(), &   ! uBT_WW is the barotropic velocity, in m s-1, beyond
                          ! which the marginal open face area is FA_u_WW.
                          ! uBT_EE must be non-negative.
    uBT_EE => NULL(), &   ! uBT_EE is the barotropic velocity, in m s-1, beyond
                          ! which the marginal open face area is FA_u_EE.
                          ! uBT_EE must be non-positive.
    FA_v_NN => NULL(), &  ! The FA_v_XX variables are the effective open face
    FA_v_N0 => NULL(), &  ! areas for barotropic transport through the meridional
    FA_v_S0 => NULL(), &  ! faces, all in H m, with the XX indicating where
    FA_v_SS => NULL(), &  ! the transport is from, with _NN drawing from points
                          ! far to the north, _N0 from points nearby from the
                          ! north, _S0 nearby from the south, and _SS from far
                          ! to the south.
    vBT_SS => NULL(), &   ! vBT_SS is the barotropic velocity, in m s-1, beyond
                          ! which the marginal open face area is FA_v_SS.
                          ! vBT_NN must be non-negative.
    vBT_NN => NULL()      ! vBT_NN is the barotropic velocity, in m s-1, beyond
                          ! which the marginal open face area is FA_v_NN.
                          ! vBT_NN must be non-positive.
end type BT_cont_type

contains


subroutine alloc_BT_cont_type(BT_cont, G)
  type(BT_cont_type),    pointer    :: BT_cont
  type(ocean_grid_type), intent(in) :: G

  integer :: isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq

  if (associated(BT_cont)) call GOLD_error(FATAL, &
    "alloc_BT_cont_type called with an associated BT_cont_type pointer.")

  allocate(BT_cont)
  allocate(BT_cont%FA_u_WW(Isdq:Iedq,jsd:jed)) ; BT_cont%FA_u_WW(:,:) = 0.0
  allocate(BT_cont%FA_u_W0(Isdq:Iedq,jsd:jed)) ; BT_cont%FA_u_W0(:,:) = 0.0
  allocate(BT_cont%FA_u_E0(Isdq:Iedq,jsd:jed)) ; BT_cont%FA_u_E0(:,:) = 0.0
  allocate(BT_cont%FA_u_EE(Isdq:Iedq,jsd:jed)) ; BT_cont%FA_u_EE(:,:) = 0.0
  allocate(BT_cont%uBT_WW(Isdq:Iedq,jsd:jed))  ; BT_cont%uBT_WW(:,:) = 0.0
  allocate(BT_cont%uBT_EE(Isdq:Iedq,jsd:jed))  ; BT_cont%uBT_EE(:,:) = 0.0

  allocate(BT_cont%FA_v_SS(isd:ied,Jsdq:Jedq)) ; BT_cont%FA_v_SS(:,:) = 0.0
  allocate(BT_cont%FA_v_S0(isd:ied,Jsdq:Jedq)) ; BT_cont%FA_v_S0(:,:) = 0.0
  allocate(BT_cont%FA_v_N0(isd:ied,Jsdq:Jedq)) ; BT_cont%FA_v_N0(:,:) = 0.0
  allocate(BT_cont%FA_v_NN(isd:ied,Jsdq:Jedq)) ; BT_cont%FA_v_NN(:,:) = 0.0
  allocate(BT_cont%vBT_SS(isd:ied,Jsdq:Jedq))  ; BT_cont%vBT_SS(:,:) = 0.0
  allocate(BT_cont%vBT_NN(isd:ied,Jsdq:Jedq))  ; BT_cont%vBT_NN(:,:) = 0.0

end subroutine alloc_BT_cont_type

subroutine dealloc_BT_cont_type(BT_cont)
  type(BT_cont_type), pointer :: BT_cont

  if (.not.associated(BT_cont)) return

  deallocate(BT_cont%FA_u_WW) ; deallocate(BT_cont%FA_u_W0)
  deallocate(BT_cont%FA_u_E0) ; deallocate(BT_cont%FA_u_EE)
  deallocate(BT_cont%uBT_WW)  ; deallocate(BT_cont%uBT_EE)

  deallocate(BT_cont%FA_v_SS) ; deallocate(BT_cont%FA_v_S0)
  deallocate(BT_cont%FA_v_N0) ; deallocate(BT_cont%FA_v_NN)
  deallocate(BT_cont%vBT_SS)  ; deallocate(BT_cont%vBT_NN)

  deallocate(BT_cont)

end subroutine dealloc_BT_cont_type


subroutine GOLD_forcing_chksum(mesg, fluxes, G, haloshift)
  character(len=*),                    intent(in) :: mesg
  type(forcing),                       intent(in) :: fluxes
  type(ocean_grid_type),               intent(in) :: G
  integer, optional,                   intent(in) :: haloshift
!   This subroutine writes out chksums for the model's basic state variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m.
!  (in)      uh - Volume flux through zonal faces = u*h*dy, m3 s-1.
!  (in)      vh - Volume flux through meridional faces = v*h*dx, in m3 s-1.
!  (in)      G - The ocean's grid structure.
  integer :: is, ie, js, je, nz, hshift
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  hshift=1; if (present(haloshift)) hshift=haloshift

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  if (associated(fluxes%taux)) &
    call uchksum(fluxes%taux, mesg//" fluxes%taux",G,haloshift=1)
  if (associated(fluxes%tauy)) &
    call vchksum(fluxes%tauy, mesg//" fluxes%tauy",G,haloshift=1)
  if (associated(fluxes%ustar)) &
    call hchksum(fluxes%ustar, mesg//" fluxes%ustar",G,haloshift=1)
  if (associated(fluxes%buoy)) &
    call hchksum(fluxes%buoy, mesg//" fluxes%buoy ",G,haloshift=hshift)
  if (associated(fluxes%sw)) &
    call hchksum(fluxes%sw, mesg//" fluxes%sw",G,haloshift=hshift)
  if (associated(fluxes%sw_vis_dir)) &
    call hchksum(fluxes%sw_vis_dir, mesg//" fluxes%sw_vis_dir",G,haloshift=hshift)
  if (associated(fluxes%sw_vis_dif)) &
    call hchksum(fluxes%sw_vis_dif, mesg//" fluxes%sw_vis_dif",G,haloshift=hshift)
  if (associated(fluxes%sw_nir_dir)) &
    call hchksum(fluxes%sw_nir_dir, mesg//" fluxes%sw_nir_dir",G,haloshift=hshift)
  if (associated(fluxes%sw_nir_dif)) &
    call hchksum(fluxes%sw_nir_dif, mesg//" fluxes%sw_nir_dif",G,haloshift=hshift)
  if (associated(fluxes%lw)) &
    call hchksum(fluxes%lw, mesg//" fluxes%lw",G,haloshift=hshift)
  if (associated(fluxes%latent)) &
    call hchksum(fluxes%latent, mesg//" fluxes%latent",G,haloshift=hshift)
  if (associated(fluxes%sens)) &
    call hchksum(fluxes%sens, mesg//" fluxes%sens",G,haloshift=hshift)
  if (associated(fluxes%evap)) &
    call hchksum(fluxes%evap, mesg//" fluxes%evap",G,haloshift=hshift)
  if (associated(fluxes%liq_precip)) &
    call hchksum(fluxes%liq_precip, mesg//" fluxes%liq_precip",G,haloshift=hshift)
  if (associated(fluxes%froz_precip)) &
    call hchksum(fluxes%froz_precip, mesg//" fluxes%froz_precip",G,haloshift=hshift)
  if (associated(fluxes%virt_precip)) &
    call hchksum(fluxes%virt_precip, mesg//" fluxes%virt_precip",G,haloshift=hshift)
  if (associated(fluxes%p_surf)) &
    call hchksum(fluxes%p_surf, mesg//" fluxes%p_surf",G,haloshift=hshift)
  if (associated(fluxes%salt_flux)) &
    call hchksum(fluxes%salt_flux, mesg//" fluxes%salt_flux",G,haloshift=hshift)
  if (associated(fluxes%TKE_tidal)) &
    call hchksum(fluxes%TKE_tidal, mesg//" fluxes%TKE_tidal",G,haloshift=hshift)
  if (associated(fluxes%ustar_tidal)) &
    call hchksum(fluxes%ustar_tidal, mesg//" fluxes%ustar_tidal",G,haloshift=hshift)
  if (associated(fluxes%liq_runoff)) &
    call hchksum(fluxes%liq_runoff, mesg//" fluxes%liq_runoff",G,haloshift=hshift)
  if (associated(fluxes%froz_runoff)) &
    call hchksum(fluxes%froz_runoff, mesg//" fluxes%froz_runoff",G,haloshift=hshift)
end subroutine GOLD_forcing_chksum

subroutine GOLD_thermovar_chksum(mesg, tv, G)
  character(len=*),                    intent(in) :: mesg
  type(thermo_var_ptrs),               intent(in) :: tv
  type(ocean_grid_type),               intent(in) :: G
!   This subroutine writes out chksums for the model's basic state variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m.
!  (in)      uh - Volume flux through zonal faces = u*h*dy, m3 s-1.
!  (in)      vh - Volume flux through meridional faces = v*h*dx, in m3 s-1.
!  (in)      G - The ocean's grid structure.
  integer :: is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  if (associated(tv%T)) &
    call hchksum(tv%T, mesg//" tv%T",G)
  if (associated(tv%S)) &
    call hchksum(tv%S, mesg//" tv%S",G)
  if (associated(tv%Rml)) &
    call hchksum(tv%Rml, mesg//" tv%Rml",G)
  if (associated(tv%frazil)) &
    call hchksum(tv%frazil, mesg//" tv%frazil",G)
  if (associated(tv%salt_deficit)) &
    call hchksum(tv%salt_deficit, mesg//" tv%salt_deficit",G)
  if (associated(tv%TempxPmE)) &
    call hchksum(tv%TempxPmE, mesg//" tv%TempxPmE",G)
end subroutine GOLD_thermovar_chksum

end module GOLD_variables
