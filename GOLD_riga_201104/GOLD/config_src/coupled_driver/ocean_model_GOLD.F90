module ocean_model_mod
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

!-----------------------------------------------------------------------
!
! This is the top level module for the ocean model.  It contains routines for
! initialization, termination and update of ocean model state.  This
! particular version wraps all of the calls for GOLD in the calls that had
! been used for MOM4.
!
! <CONTACT EMAIL="Robert.Hallberg@noaa.gov"> Robert Hallberg
! </CONTACT>
!
!<OVERVIEW>
! This code is a stop-gap wrapper of the GOLD code to enable it to be called
! in the same way as MOM4.
!</OVERVIEW>

use GOLD, only : initialize_GOLD, step_GOLD, GOLD_control_struct, GOLD_end
use GOLD, only : calculate_surface_state
use GOLD_constants, only : CELSIUS_KELVIN_OFFSET
use GOLD_diag_mediator, only : enable_averaging, disable_averaging
use GOLD_domains, only : pass_vector, AGRID, BGRID_NE, CGRID_NE
use GOLD_error_handler, only : GOLD_error, FATAL, WARNING, is_root_pe
use GOLD_file_parser, only : read_param, close_param_file, param_file_type
use GOLD_file_parser, only : log_param, log_version
use GOLD_grid, only : ocean_grid_type
use GOLD_initialization, only : get_GOLD_Input
use GOLD_io, only : close_file, file_exists, read_data, write_version_number
use GOLD_restart, only : save_restart
use GOLD_sum_output, only : write_energy, accumulate_net_input
use GOLD_sum_output, only : GOLD_sum_output_init, sum_output_CS
use GOLD_surface_forcing, only : surface_forcing_init, convert_IOB_to_fluxes
use GOLD_surface_forcing, only : average_forcing, ice_ocn_bnd_type_chksum
use GOLD_surface_forcing, only : ice_ocean_boundary_type, surface_forcing_CS
use GOLD_surface_forcing, only : forcing_save_restart
use GOLD_time_manager, only : time_type, get_time, set_time, operator(>)
use GOLD_time_manager, only : operator(+), operator(-), operator(*), operator(/)
use GOLD_time_manager, only : operator(/=)
use GOLD_tracer_flow_control, only : call_tracer_register, tracer_flow_control_init
use GOLD_variables, only : forcing, surface, directories

use coupler_types_mod, only : coupler_2d_bc_type
use mpp_domains_mod, only : domain2d, mpp_get_layout, mpp_get_global_domain
use mpp_domains_mod, only : mpp_define_domains, mpp_get_compute_domain, mpp_get_data_domain
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux
use fms_mod, only : stdout
use mpp_mod, only : mpp_chksum

#ifdef _USE_TOPAZ
use GOLD_OCEAN_TOPAZ_MOD,only : TOPAZ_coupler_flux_init
#endif
#ifdef _USE_GENERIC_TRACER
use GOLD_generic_tracer, only : GOLD_generic_flux_init
#endif

implicit none ; private

#include <GOLD_memory.h>

public ocean_model_init, ocean_model_end, update_ocean_model
public ocean_model_save_restart, Ocean_stock_pe
public ice_ocean_boundary_type
public ocean_model_init_sfc, ocean_model_flux_init
public ocean_model_restart
public ice_ocn_bnd_type_chksum
public ocean_public_type_chksum
public    ocean_model_data_get
interface ocean_model_data_get
   module procedure ocean_model_data1D_get 
   module procedure ocean_model_data2D_get 
end interface



! For communication with FMS coupler
type, public ::  ocean_public_type
  type(domain2d) :: Domain       ! The domain for the surface fields.
  logical :: is_ocean_pe         ! .true. on processors that run the ocean model.
  character(len=32) :: instance_name = "" ! A name that can be used to identify
                                 ! this instance of an ocean model, for example
                                 ! in ensembles when writing messages.
  integer, pointer, dimension(:) :: pelist => NULL()   ! The list of ocean PEs.
  logical, pointer, dimension(:,:) :: maskmap =>NULL() ! A pointer to an array
                    ! indicating which logical processors are actually used for
                    ! the ocean code. The other logical processors would be all
                    ! land points and are not assigned to actual processors.
                    ! This need not be assigned if all logical processors are used.

  integer :: stagger = BGRID_NE  ! The staggering relative to the tracer points
                                 ! points of the two velocity components. Valid
                                 ! entries include AGRID, BGRID_NE, CGRID_NE,
                                 ! BGRID_SW, and CGRID_SW, corresponding to the
                                 ! community-standard Arakawa notation.  (These
                                 ! are named integers taken from mpp_parameter_mod.)
                                 ! Following MOM, this is BGRID_NE by default.
  real, pointer, dimension(:,:)  :: t_surf =>NULL()  ! SST on t-cell (degrees Kelvin)
  real, pointer, dimension(:,:)  :: s_surf =>NULL()  ! SSS on t-cell (psu)
  real, pointer, dimension(:,:)  :: u_surf =>NULL()  ! i-velocity the points indicated
                                                     ! by velocity_stagger. (m/s)
  real, pointer, dimension(:,:)  :: v_surf =>NULL()  ! j-velocity the points indicated
                                                     ! by velocity_stagger. (m/s)
  real, pointer, dimension(:,:)  :: sea_lev =>NULL() ! Sea level in m after correction
                                                     ! for surface pressure, i.e.
                                                     ! dzt(1) + eta_t + patm/rho0/grav (m)
  real, pointer, dimension(:,:)  :: frazil  =>NULL() ! Accumulated heating (Joules/m^2)
                                                     ! from frazil formation in the ocean.
  type(coupler_2d_bc_type)       :: fields    ! A structure that may contain an
                                              ! array of named tracer-related fields.
  integer                        :: avg_kount ! Used for accumulating averages of this type.
  integer, dimension(3)          :: axes = 0  ! Axis numbers that are available
                                              ! for I/O using this surface data.

  real, pointer, dimension(:,:)  :: area =>NULL() ! cell area of the ocean surface.

end type ocean_public_type


type, public :: ocean_state_type ; private
  ! This type is private, and can therefore vary between different ocean models.
  ! All information entire ocean state may be contained here, although it is not
  ! necessary that this is implemented with all models.
  logical :: is_ocean_PE = .false.  ! True if this is an ocean PE.
  type(time_type) :: Time    ! The ocean model's time and master clock.
  type(time_type) :: restint ! The time between saves of the restart file.
  type(time_type) :: restart_time ! The next time to write restart files.
  integer :: Restart_control ! An integer that is bit-tested to determine whether
                             ! incremental restart files are saved and whether they
                             ! have a time stamped name.  +1 (bit 0) for generic
                             ! files and +2 (bit 1) for time-stamped files.  A
                             ! restart file is saved at the end of a run segment
                             ! unless Restart_control is negative.
  type(time_type) :: energysavedays  ! The interval between writing the energies
                                     ! and other integral quantities of the run.
  type(time_type) :: write_energy_time ! The next time to write to the energy file.

  integer :: nstep = 0        ! The number of calls to update_ocean.
  integer :: m_last           ! The last time level calculated by step_GOLD.

  logical :: restore_salinity ! If true, the coupled GOLD driver adds a term to
                              ! restore salinity to a specified value.
  real :: press_to_z          ! A conversion factor between pressure and ocean
                              ! depth in m, usually 1/(rho_0*g), in m Pa-1.
  real    :: C_p              !   The heat capacity of seawater, in J K-1 kg-1.

  type(directories) :: dirs   ! A structure containing several relevant directory paths.
  type(forcing)   :: fluxes   ! A structure containing pointers to
                              ! the ocean forcing fields.
  type(surface)   :: state    ! A structure containing pointers to
                              ! the ocean surface state fields.
  type(ocean_grid_type), pointer :: grid => NULL() ! A pointer to a structure
                              ! containing metrics and related information.
  type(GOLD_control_struct), pointer :: GOLD_CSp => NULL()
  type(surface_forcing_CS),  pointer :: forcing_CSp => NULL()
  type(sum_output_CS),       pointer :: sum_output_CSp => NULL()
end type ocean_state_type

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_model_init">
!
! <DESCRIPTION>
! Initialize the ocean model.
! </DESCRIPTION>
!
subroutine ocean_model_init(Ocean_sfc, OS, Time_init, Time_in)
  type(ocean_public_type), target, intent(inout) :: Ocean_sfc
  type(ocean_state_type),          pointer       :: OS
  type(time_type),                 intent(in)    :: Time_init
  type(time_type),                 intent(in)    :: Time_in
!   This subroutine initializes both the ocean state and the ocean surface type.
! Because of the way that indicies and domains are handled, Ocean_sfc must have
! been used in a previous call to initialize_ocean_type.

! Arguments: Ocean_sfc - A structure containing various publicly visible ocean
!                    surface properties after initialization, this is intent(out).
!  (out,private) OS - A structure whose internal contents are private
!                    to ocean_model_mod that may be used to contain all
!                    information about the ocean's interior state.
!  (in)      Time_init - The start time for the coupled model's calendar.
!  (in)      Time_in - The time at which to initialize the ocean model.
  real :: Time_unit   ! The time unit in seconds for RESTINT and ENERGYSAVEDAYS.
  real :: Rho0        ! The Boussinesq ocean density, in kg m-3.
  real :: G_Earth     ! The gravitational acceleration in m s-2.
  character(len=128) :: version = '$Id: ocean_model_GOLD.F90,v 13.0.2.9.2.16 2011/07/21 17:12:42 Robert.Hallberg Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "ocean_model_init"  ! This module's name.
  integer :: secs, days
  type(param_file_type) :: param_file

  if (associated(OS)) then
    call GOLD_error(WARNING, "ocean_model_init called with an associated "// &
                    "ocean_state_type structure. Model is already initialized.")
    return
  endif
  allocate(OS)

  OS%is_ocean_pe = Ocean_sfc%is_ocean_pe
  if (.not.OS%is_ocean_pe) return

  OS%state%tr_fields => Ocean_sfc%fields
  OS%Time = Time_in
  call initialize_GOLD(OS%Time, param_file, OS%dirs, OS%GOLD_CSp, Time_in)
  OS%grid => OS%GOLD_CSp%grid

  OS%Restart_control = 1 ; OS%restint = Time_in + set_time(0, 2000000000)
  OS%energysavedays = set_time(0,1) ; OS%restore_salinity = .false.
  call read_param(param_file,"RESTART_CONTROL",OS%Restart_control)
  Time_unit = 86400.0 ;  call read_param(param_file,"TIMEUNIT",Time_unit)
  call read_param(param_file,"RESTINT",OS%restint,Time_unit)
  call read_param(param_file,"ENERGYSAVEDAYS",OS%energysavedays,Time_unit)
  call read_param(param_file,"RESTORE_SALINITY",OS%restore_salinity)

  call read_param(param_file,"RHO_0",Rho0,.true.)
  call read_param(param_file,"G_EARTH",G_Earth,.true.)
  OS%press_to_z = 1.0/(Rho0*G_Earth)
  call read_param(param_file, "C_P", OS%C_p, .true.)

  call surface_forcing_init(Time_in, OS%grid, param_file, OS%GOLD_CSp%diag, &
                            OS%forcing_CSp, OS%restore_salinity)
  call GOLD_sum_output_init(OS%grid, param_file, OS%dirs%output_directory, &
                           OS%GOLD_CSp%ntrunc, Time_init, OS%sum_output_CSp)

  call write_energy(OS%GOLD_CSp%u(:,:,:,1), OS%GOLD_CSp%v(:,:,:,1), OS%GOLD_CSp%h(:,:,:,1), &
                    OS%GOLD_CSp%tv, OS%Time, 0, OS%grid, OS%sum_output_CSp, &
                    OS%GOLD_CSp%tracer_flow_CSp)

  ! restart_time is the next integral multiple of restint.
  OS%restart_time = Time_init + OS%restint * &
      (1 + (OS%Time - Time_init) / OS%restint)

  ! write_energy_time is the next integral multiple of energysavedays.
  OS%write_energy_time = Time_init + OS%energysavedays * &
      (1 + (OS%Time - Time_init) / OS%energysavedays)

  call initialize_ocean_public_type(OS%grid%Domain%mpp_domain,Ocean_sfc)

!  call convert_state_to_ocean_type(state, Ocean_sfc, OS%grid)

  call log_version(param_file, mod, version, tagname, "")
  call log_param(param_file, mod, "RESTART_CONTROL", OS%Restart_control, &
                 "An integer whose bits encode which restart files are \n"//&
                 "written. Add 2 (bit 1) for a time-stamped file, and odd \n"//&
                 "(bit 0) for a non-time-stamped file.  A restart file \n"//&
                 "will be saved at the end of the run segment for any \n"//&
                 "non-negative value.", default=1)
  call log_param(param_file, mod, "TIMEUNIT", Time_unit, &
                 "The time unit for DAYMAX, ENERGYSAVEDAYS, and RESTINT.", &
                 units="s", default=86400.0)
  call log_param(param_file, mod, "RESTINT", OS%restint, &
                 "The interval between saves of the restart file in units \n"//&
                 "of TIMEUNIT.  Use a value that is larger than DAYMAX to \n"//&
                 "not save incremental restart files  within a run.  Use \n"//&
                 "0 not to save restart files at all.  The default is to \n"//&
                 "use a larger value than DAYMAX.", timeunit=Time_unit)
  call log_param(param_file, mod, "ENERGYSAVEDAYS",OS%energysavedays, &
                 "The interval in units of TIMEUNIT between saves of the \n"//&
                 "energies of the run and other globally summed diagnostics.", &
                 default=set_time(0,1), timeunit=Time_unit)
  call log_param(param_file, mod, "RESTORE_SALINITY",OS%restore_salinity, &
                 "If true, the coupled driver will add a globally-balanced \n"//&
                 "fresh-water flux that drives sea-surface salinity \n"//&
                 "toward specified values.", default=.false.)
  call log_param(param_file, mod, "RHO_0",Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3")
  call log_param(param_file, mod, "G_EARTH",G_Earth, &
                 "The gravitational acceleration of the Earth.", units="m s-2")
  call log_param(param_file, mod, "C_P", OS%C_p, &
                 "The heat capacity of sea water, approximated as a \n"//&
                 "constant. This is only used if TEMPERATURE is true.", &
                 units="J kg-1 K-1")

  call close_param_file(param_file)
 
  if (is_root_pe()) &
    write(*,'(/12x,a/)') '======== COMPLETED GOLD INITIALIZATION ========'

  OS%m_last = 1

end subroutine ocean_model_init
! </SUBROUTINE> NAME="ocean_model_init"


!#######################################################################
! <SUBROUTINE NAME="update_ocean_model">
!
! <DESCRIPTION>
! Update in time the ocean model fields.  This code wraps the call to step_GOLD
! with MOM4's call.
! </DESCRIPTION>
!

subroutine update_ocean_model(Ice_ocean_boundary, OS, Ocean_sfc, &
                              time_start_update, Ocean_coupling_time_step)
  type(ice_ocean_boundary_type), intent(inout) :: Ice_ocean_boundary
  type(ocean_state_type),        pointer       :: OS
  type(ocean_public_type),       intent(inout) :: Ocean_sfc
  type(time_type), intent(in)                  :: time_start_update
  type(time_type), intent(in)                  :: Ocean_coupling_time_step
!   This subroutine uses the forcing in Ice_ocean_boundary to advance the
! ocean model's state from the input value of Ocean_state (which must be for
! time time_start_update) for a time interval of Ocean_coupling_time_step,
! returning the publicly visible ocean surface properties in Ocean_sfc and
! storing the new ocean properties in Ocean_state.

! Arguments: Ice_ocean_boundary - A structure containing the various forcing
!                                 fields coming from the ice. It is intent in.
!  (inout)   Ocean_state - A structure containing the internal ocean state.
!  (out)     Ocean_sfc - A structure containing all the publicly visible ocean
!                        surface fields after a coupling time step.
!  (in)      time_start_update - The time at the beginning of the update step.
!  (in)      Ocean_coupling_time_step - The amount of time over which to advance
!                                       the ocean.

! Note: although several types are declared intent(inout), this is to allow for
!   the possibility of halo updates and to keep previously allocated memory.
!   In practice, Ice_ocean_boundary is intent in, Ocean_state is private to
!   this module and intent inout, and Ocean_sfc is intent out.
  type(time_type) :: Master_time ! This allows step_GOLD to temporarily change
                                 ! the time that is seen by internal modules.
  type(time_type) :: Time1       ! The value of the ocean model's time at the
                                 ! start of a call to step_GOLD.
  integer :: index_bnds(4)       ! The computational domain index bounds in the
                                 ! ice-ocean boundary type.

  real :: time_step         ! The time step of a call to step_GOLD in seconds.
  integer :: secs, days

  call get_time(Ocean_coupling_time_step, secs, days)
  time_step = 86400.0*real(days) + real(secs)
  
  if (time_start_update /= OS%Time) then
    call GOLD_error(WARNING, "update_ocean_model: internal clock does not "//&
                             "with time_start_update argument.")
  endif
  if (.not.associated(OS)) then
    call GOLD_error(FATAL, "update_ocean_model called with an unassociated "// &
                    "ocean_state_type structure. ocean_model_init must be "//  &
                    "called first to allocate this structure.")
    return
  endif

!  Translate Ice_ocean_boundary into fluxes.
  call mpp_get_compute_domain(Ocean_sfc%Domain, index_bnds(1), index_bnds(2), &
                              index_bnds(3), index_bnds(4))
  call convert_IOB_to_fluxes(Ice_ocean_boundary, OS%fluxes, index_bnds, OS%Time, &
                             OS%grid, OS%forcing_CSp, OS%state, OS%restore_salinity)
  Master_time = OS%Time ; Time1 = OS%Time

  OS%m_last = step_GOLD(OS%fluxes, OS%state, Time1, time_step, OS%GOLD_CSp)

  OS%Time = Master_time + Ocean_coupling_time_step
  OS%nstep = OS%nstep + 1

  call enable_averaging(time_step, OS%Time, OS%GOLD_CSp%diag)
  call average_forcing(OS%fluxes, time_step, OS%grid, OS%forcing_CSp)
  call accumulate_net_input(OS%fluxes, OS%state, time_step, OS%grid, OS%sum_output_CSp)
  call disable_averaging(OS%GOLD_CSp%diag)

!  See if it is time to write out the energy.
  if (OS%Time + ((Ocean_coupling_time_step)/2) > OS%write_energy_time) then
    call write_energy(OS%GOLD_CSp%u(:,:,:,OS%m_last), OS%GOLD_CSp%v(:,:,:,OS%m_last), &
                      OS%GOLD_CSp%h(:,:,:,OS%m_last), OS%GOLD_CSp%tv, OS%Time, OS%nstep, &
                      OS%grid, OS%sum_output_CSp, OS%GOLD_CSp%tracer_flow_CSp)
    OS%write_energy_time = OS%write_energy_time + OS%energysavedays
  endif

! Translate state into Ocean.
!  call convert_state_to_ocean_type(OS%state, Ocean_sfc, OS%grid, &
!                                   Ice_ocean_boundary%p, OS%press_to_z)
  call convert_state_to_ocean_type(OS%state, Ocean_sfc, OS%grid)

end subroutine update_ocean_model
! </SUBROUTINE> NAME="update_ocean_model"

!#######################################################################
! <SUBROUTINE NAME="ocean_model_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine ocean_model_restart(OS, timestamp)
   type(ocean_state_type),        pointer :: OS
   character(len=*), intent(in), optional :: timestamp

   if (BTEST(OS%Restart_control,1)) then
     call save_restart(OS%dirs%restart_output_dir,OS%Time,OS%m_last, OS%grid, &
                       OS%GOLD_CSp%restart_CSp, .true.)
     call forcing_save_restart(OS%forcing_CSp, OS%grid, OS%Time, &
                               OS%dirs%restart_output_dir, .true.)
   endif
   if (BTEST(OS%Restart_control,0)) then
     call save_restart(OS%dirs%restart_output_dir,OS%Time,OS%m_last, OS%grid, &
                       OS%GOLD_CSp%restart_CSp)
     call forcing_save_restart(OS%forcing_CSp, OS%grid, OS%Time, &
                               OS%dirs%restart_output_dir)
   endif
  
end subroutine ocean_model_restart
! </SUBROUTINE> NAME="ocean_model_restart"

!#######################################################################
! <SUBROUTINE NAME="ocean_model_end">
!
! <DESCRIPTION>
! Close down the ocean model
! </DESCRIPTION>
!
subroutine ocean_model_end(Ocean_sfc, Ocean_state, Time)
  type(ocean_public_type),           intent(inout) :: Ocean_sfc
  type(ocean_state_type),            pointer       :: Ocean_state
  type(time_type),                   intent(in)    :: Time
!   This subroutine terminates the model run, saving the ocean state in a
! restart file and deallocating any data associated with the ocean.

! Arguments: Ocean_sfc - An ocean_public_type structure that is to be
!                        deallocated upon termination.
!  (inout)   Ocean_state - A pointer to the structure containing the internal
!                          ocean state to be deallocated upon termination.
!  (in)      Time - The model time, used for writing restarts.

  call ocean_model_save_restart(Ocean_state, Time)
  call GOLD_end(Ocean_state%GOLD_CSp)
end subroutine ocean_model_end
! </SUBROUTINE> NAME="ocean_model_end"


subroutine ocean_model_save_restart(OS, Time, directory, filename_suffix)
  type(ocean_state_type),     pointer    :: OS
  type(time_type),            intent(in) :: Time
  character(len=*), optional, intent(in) :: directory
  character(len=*), optional, intent(in) :: filename_suffix
! Arguments: Ocean_state - A structure containing the internal ocean state (in).
!  (in)      Time - The model time at this call.  This is needed for mpp_write calls.
!  (in, opt) directory - An optional directory into which to write these restart files.
!  (in, opt) filename_suffix - An optional suffix (e.g., a time-stamp) to append
!                              to the restart file names.

! Note: This is a new routine - it will need to exist for the new incremental
!   checkpointing.  It will also be called by ocean_model_end, giving the same
!   restart behavior as now in FMS.
  character(len=200) :: restart_dir

  if (present(directory)) then ; restart_dir = directory
  else ; restart_dir = OS%dirs%restart_output_dir ; endif

  call save_restart(restart_dir, Time, OS%m_last, OS%grid, &
                    OS%GOLD_CSp%restart_CSp)
  call forcing_save_restart(OS%forcing_CSp, OS%grid, Time, &
                            restart_dir)

end subroutine ocean_model_save_restart

!#######################################################################

subroutine initialize_ocean_public_type(input_domain, Ocean_sfc)
  type(domain2D), intent(in)             :: input_domain
  type(ocean_public_type), intent(inout) :: Ocean_sfc
  integer :: xsz, ysz, layout(2)
  ! ice-ocean-boundary fields are always allocated using absolute indicies
  ! and have no halos.
  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd

  call mpp_get_layout(input_domain,layout)
  call mpp_get_global_domain(input_domain, xsize=xsz, ysize=ysz)
  call mpp_define_domains((/1,xsz,1,ysz/),layout,Ocean_sfc%Domain)

  call mpp_get_compute_domain(Ocean_sfc%Domain, isc_bnd, iec_bnd, &
                              jsc_bnd, jec_bnd)

  allocate ( Ocean_sfc%t_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%s_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%u_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%v_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%sea_lev(isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%area   (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%frazil (isc_bnd:iec_bnd,jsc_bnd:jec_bnd))

  Ocean_sfc%t_surf  = 0.0  ! time averaged sst (Kelvin) passed to atmosphere/ice model
  Ocean_sfc%s_surf  = 0.0  ! time averaged sss (psu) passed to atmosphere/ice models
  Ocean_sfc%u_surf  = 0.0  ! time averaged u-current (m/sec) passed to atmosphere/ice models
  Ocean_sfc%v_surf  = 0.0  ! time averaged v-current (m/sec)  passed to atmosphere/ice models
  Ocean_sfc%sea_lev = 0.0  ! time averaged thickness of top model grid cell (m) plus patm/rho0/grav
  Ocean_sfc%frazil  = 0.0  ! time accumulated frazil (J/m^2) passed to ice model
  Ocean_sfc%area    = 0.0

end subroutine initialize_ocean_public_type

subroutine convert_state_to_ocean_type(state, Ocean_sfc, G, patm, press_to_z)
  type(surface),           intent(inout) :: state
  type(ocean_public_type), target, intent(inout) :: Ocean_sfc
  type(ocean_grid_type),   intent(inout) :: G
  real,          optional, intent(in)    :: patm(:,:)
  real,          optional, intent(in)    :: press_to_z
! This subroutine translates the coupler's ocean_data_type into GOLD's
! surface state variable.  This may eventually be folded into the GOLD
! code that calculates the surface state in the first place.
! Note the offset in the arrays because the ocean_data_type has no
! halo points in its arrays and always uses absolute indicies.
  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd
  integer :: i, j, i0, j0, is, ie, js, je
  real :: IgR0

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  call pass_vector(state%u,state%v,G%Domain)

  call mpp_get_compute_domain(Ocean_sfc%Domain, isc_bnd, iec_bnd, &
                              jsc_bnd, jec_bnd)
  if (present(patm)) then
    ! Check that the inidicies in patm are (isc_bnd:iec_bnd,jsc_bnd:jec_bnd).
    if (.not.present(press_to_z)) call GOLD_error(FATAL, &
        'convert_state_to_ocean_type: press_to_z must be present if patm is.')
  endif

  i0 = is - isc_bnd ; j0 = js - jsc_bnd
  do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
    Ocean_sfc%t_surf(i,j) = state%SST(i+i0,j+j0) + CELSIUS_KELVIN_OFFSET
    Ocean_sfc%s_surf(i,j) = state%SSS(i+i0,j+j0)
    Ocean_sfc%u_surf(i,j) = G%qmask(i+i0,j+j0)*0.5*(state%u(i+i0,j+j0)+state%u(i+i0,j+j0+1))
    Ocean_sfc%v_surf(i,j) = G%qmask(i+i0,j+j0)*0.5*(state%v(i+i0,j+j0)+state%v(i+i0+1,j+j0))
    Ocean_sfc%sea_lev(i,j) = state%sea_lev(i+i0,j+j0)
    if (present(patm)) &
      Ocean_sfc%sea_lev(i,j) = Ocean_sfc%sea_lev(i,j) + patm(i,j) * press_to_z
    Ocean_sfc%frazil(i,j) = state%frazil(i+i0,j+j0)
    Ocean_sfc%area(i,j)   =  G%DXDYh(i+i0,j+j0)  
  enddo ; enddo

  if (.not.associated(state%tr_fields,Ocean_sfc%fields)) &
    call GOLD_error(FATAL,'state%tr_fields is not pointing to Ocean_sfc%fields')
  
end subroutine convert_state_to_ocean_type


!#######################################################################
! <SUBROUTINE NAME="ocean_model_init_sfc">
!
! <DESCRIPTION>
!   This subroutine extracts the surface properties from the ocean's internal
! state and stores them in the ocean type returned to the calling ice model.
! It has to be separate from the ocean_initialization call because the coupler
! module allocates the space for some of these variables.
! </DESCRIPTION>

subroutine ocean_model_init_sfc(OS, Ocean_sfc)
  type(ocean_state_type),  pointer       :: OS
  type(ocean_public_type), intent(inout) :: Ocean_sfc

  call calculate_surface_state(OS%state, OS%GOLD_CSp%u(:,:,:,1), &
           OS%GOLD_CSp%v(:,:,:,1), OS%GOLD_CSp%h(:,:,:,1), OS%GOLD_CSp%ave_ssh,&
           OS%grid, OS%GOLD_CSp)

  call convert_state_to_ocean_type(OS%state, Ocean_sfc, OS%grid)

end subroutine ocean_model_init_sfc
! </SUBROUTINE NAME="ocean_model_init_sfc">

!   I have no idea what the following subroutine was intended to do, but it
! appears to be necessary for compiling the Memphis coupled model.  In the MOM
! version, it makes a call to something in the "OCEAN_TPM_MOD_FLUX_INIT".  I
! believe that the equivalent calls are already taken care of inside of
! GOLD_initialize, and given that it has no arguments, in the GOLD paradigm,
! it can do nothing. I think that it should be elminated.  -RWH

! UPDATE May, 8 2007
! added aof_set_couple_flux function calls to this routine. These calls
! are only made by non Ocean PEs and only when CFCs are being used (though
! this condition may expand. These calls are normally done during 
! ocean_model_init. The problem is that they are only done by Ocean PEs. All
! processors need to make this call in order to get the number of ocean/atmos
! fluxes correct in the coupler framework. 
! The trick now is to let the atmos pes know when there is no CFCs and not 
! to call aof_set_coupler_flux
!WGA

subroutine ocean_model_flux_init(OS)
  type(ocean_state_type),  pointer       :: OS
  integer :: dummy
  character(len=128) :: default_ice_restart_file = 'ice_ocmip2_cfc.res.nc'
  character(len=128) :: default_ocean_restart_file = 'ocmip2_cfc.res.nc'
  character(len=40)  :: mod = "ocean_model_flux_init"  ! This module's name.

  type(param_file_type) :: param_file
  type(directories) :: dirs_tmp  ! A structure containing several relevant directory paths.
  logical :: use_OCMIP_CFCs = .false.
  logical :: use_GOLD_TOPAZ = .false.
  logical :: use_GOLD_generic_tracer = .false.

  call get_GOLD_Input(param_file, dirs_tmp)

  if(.not.associated(OS)) then
    call read_param(param_file,"USE_OCMIP2_CFC",use_OCMIP_CFCs)
    if (use_OCMIP_CFCs)then
      dummy = aof_set_coupler_flux('cfc_11_flux', &
        flux_type = 'air_sea_gas_flux', implementation = 'ocmip2', &
        param = (/ 9.36e-07, 9.7561e-06 /), &
        ice_restart_file = default_ice_restart_file, &
        ocean_restart_file = default_ocean_restart_file,  &
        caller = "register_OCMIP2_CFC")
      dummy = aof_set_coupler_flux('cfc_12_flux', &
        flux_type = 'air_sea_gas_flux', implementation = 'ocmip2', &
        param = (/ 9.36e-07, 9.7561e-06 /), &
        ice_restart_file = default_ice_restart_file, &
        ocean_restart_file = default_ocean_restart_file, &
        caller = "register_OCMIP2_CFC")
    endif 
  endif

  call read_param(param_file,"USE_TOPAZ",use_GOLD_TOPAZ)
  call read_param(param_file,"USE_generic_tracer",use_GOLD_generic_tracer)

  call log_param(param_file, mod, "USE_OCMIP2_CFC", use_OCMIP_CFCs)
  call log_param(param_file, mod, "USE_TOPAZ", use_GOLD_TOPAZ)
  call log_param(param_file, mod, "USE_generic_tracer", use_GOLD_generic_tracer)

  call close_param_file(param_file)
   

  if (use_GOLD_TOPAZ) then
#ifdef _USE_TOPAZ
    call TOPAZ_coupler_flux_init
#else
    call GOLD_error(FATAL, &
       "call_tracer_register: use_GOLD_TOPAZ=.true. BUT not compiled with _USE_TOPAZ")
#endif
  endif

  if (use_GOLD_generic_tracer) then
#ifdef _USE_GENERIC_TRACER
    call GOLD_generic_flux_init
#else
    call GOLD_error(FATAL, &
       "call_tracer_register: use_GOLD_generic_tracer=.true. BUT not compiled with _USE_GENERIC_TRACER")
#endif
  endif
    
end subroutine ocean_model_flux_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Ocean_stock_pe - returns stocks of heat, water, etc. for conservation checks.!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine Ocean_stock_pe(OS, index, value, time_index)
  use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT,ISTOCK_SALT
  type(ocean_state_type), pointer     :: OS
  integer,                intent(in)  :: index
  real,                   intent(out) :: value
  integer,      optional, intent(in)  :: time_index
! Arguments: OS - A structure containing the internal ocean state.
!  (in)      index - Index of conservation quantity of interest.
!  (in)      value -  Sum returned for the conservation quantity of interest.
!  (in,opt)  time_index - Index for time level to use if this is necessary.

  real :: to_heat, to_mass, to_salt, PSU_to_kg ! Conversion constants.
  integer :: i, j, k, m, is, ie, js, je, nz, ind

  value = 0.0
  if (.not.associated(OS)) return
  if (.not.OS%is_ocean_pe) return
  
  is = OS%grid%isc ; ie = OS%grid%iec
  js = OS%grid%jsc ; je = OS%grid%jec ; nz = OS%grid%ke
  m = OS%m_last

  select case (index)
    case (ISTOCK_WATER)
      ! Return the mass of fresh water in the ocean on this PE in kg.
      to_mass = OS%grid%H_to_kg_m2
      if (OS%grid%Boussinesq) then
        do k=1,nz ; do j=js,je ; do i=is,ie ; if (OS%grid%hmask(i,j) > 0.5) then
          value = value + to_mass*(OS%GOLD_CSp%h(i,j,k,m) * OS%grid%DXDYh(i,j))
        endif ; enddo ; enddo ; enddo
      else
        ! In non-Boussinesq mode, the mass of salt needs to be subtracted.
        PSU_to_kg = 1.0e-3
        do k=1,nz ; do j=js,je ; do i=is,ie ; if (OS%grid%hmask(i,j) > 0.5) then
          value = value + to_mass * ((1.0 - PSU_to_kg*OS%GOLD_CSp%tv%S(i,j,k))*&
                                  (OS%GOLD_CSp%h(i,j,k,m) * OS%grid%DXDYh(i,j)))
        endif ; enddo ; enddo ; enddo
      endif
    case (ISTOCK_HEAT)
      ! Return the heat content of the ocean on this PE in J.
      to_heat = OS%grid%H_to_kg_m2 * OS%C_p
      do k=1,nz ; do j=js,je ; do i=is,ie ; if (OS%grid%hmask(i,j) > 0.5) then
        value = value + (to_heat * OS%GOLD_CSp%tv%T(i,j,k)) * &
                        (OS%GOLD_CSp%h(i,j,k,m)*OS%grid%DXDYh(i,j))
      endif ; enddo ; enddo ; enddo
    case (ISTOCK_SALT)
      ! Return the mass of the salt in the ocean on this PE in kg.
      ! The 1000 converts salinity in PSU to salt in kg kg-1.
      to_salt = OS%grid%H_to_kg_m2 / 1000.0
      do k=1,nz ; do j=js,je ; do i=is,ie ; if (OS%grid%hmask(i,j) > 0.5) then
        value = value + (to_salt * OS%GOLD_CSp%tv%S(i,j,k)) * &
                        (OS%GOLD_CSp%h(i,j,k,m)*OS%grid%DXDYh(i,j))
      endif ; enddo ; enddo ; enddo
    case default ; value = 0.0
  end select

end subroutine Ocean_stock_pe

subroutine ocean_model_data2D_get(OS,Ocean, name, array2D,isc,jsc)
  use GOLD_constants, only : CELSIUS_KELVIN_OFFSET
  type(ocean_state_type),     pointer    :: OS
  type(ocean_public_type),    intent(in) :: Ocean
  character(len=*)          , intent(in) :: name
  real, dimension(isc:,jsc:), intent(out):: array2D
  integer                   , intent(in) :: isc,jsc

  integer :: g_isc, g_iec, g_jsc, g_jec,g_isd, g_ied, g_jsd, g_jed, i, j

  if (.not.associated(OS)) return
  if (.not.OS%is_ocean_pe) return
  
! The problem is %DXDYh is on GOLD domain but Ice_Ocean_Boundary%... is on mpp domain.
! We want to return the GOLD data on the mpp (compute) domain
! Get Gold domain extents 
  call mpp_get_compute_domain(OS%grid%Domain%mpp_domain, g_isc, g_iec, g_jsc, g_jec)
  call mpp_get_data_domain   (OS%grid%Domain%mpp_domain, g_isd, g_ied, g_jsd, g_jed)

  g_isc = g_isc-g_isd+1 ; g_iec = g_iec-g_isd+1 ; g_jsc = g_jsc-g_jsd+1 ; g_jec = g_jec-g_jsd+1
 

  select case(name)
  case('area')
     array2D(isc:,jsc:) = OS%grid%DXDYh(g_isc:g_iec,g_jsc:g_jec)
  case('mask')     
     array2D(isc:,jsc:) = OS%grid%hmask(g_isc:g_iec,g_jsc:g_jec)
!OR same result
!     do j=g_jsc,g_jec; do i=g_isc,g_iec
!        array2D(isc+i-g_isc,jsc+j-g_jsc) = OS%grid%hmask(i,j)
!     enddo; enddo
  case('t_surf')
     array2D(isc:,jsc:) = Ocean%t_surf(isc:,jsc:)-CELSIUS_KELVIN_OFFSET
  case('t_pme')
     array2D(isc:,jsc:) = Ocean%t_surf(isc:,jsc:)-CELSIUS_KELVIN_OFFSET
  case('t_runoff')
     array2D(isc:,jsc:) = Ocean%t_surf(isc:,jsc:)-CELSIUS_KELVIN_OFFSET
  case('t_calving')
     array2D(isc:,jsc:) = Ocean%t_surf(isc:,jsc:)-CELSIUS_KELVIN_OFFSET
  case('btfHeat')
     array2D(isc:,jsc:) = 0
  case default
     call GOLD_error(FATAL,'get_ocean_grid_data2D: unknown argument name='//name)
  end select
  

end subroutine ocean_model_data2D_get

subroutine ocean_model_data1D_get(OS,Ocean, name, value)
  type(ocean_state_type),     pointer    :: OS
  type(ocean_public_type),    intent(in) :: Ocean
  character(len=*)          , intent(in) :: name
  real                      , intent(out):: value
  
  if (.not.associated(OS)) return
  if (.not.OS%is_ocean_pe) return

  select case(name)
  case('c_p')
     value = OS%C_p
  case default
     call GOLD_error(FATAL,'get_ocean_grid_data1D: unknown argument name='//name)
  end select
  

end subroutine ocean_model_data1D_get

subroutine ocean_public_type_chksum(id, timestep, ocn)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ocean_public_type), intent(in) :: ocn
    integer ::   n,m, outunit 
 
    outunit = stdout()

    write(outunit,*) "BEGIN CHECKSUM(ocean_type):: ", id, timestep
    write(outunit,100) 'ocean%t_surf   ',mpp_chksum(ocn%t_surf )
    write(outunit,100) 'ocean%s_surf   ',mpp_chksum(ocn%s_surf )
    write(outunit,100) 'ocean%u_surf   ',mpp_chksum(ocn%u_surf )
    write(outunit,100) 'ocean%v_surf   ',mpp_chksum(ocn%v_surf )
    write(outunit,100) 'ocean%sea_lev  ',mpp_chksum(ocn%sea_lev)
    write(outunit,100) 'ocean%frazil   ',mpp_chksum(ocn%frazil )

    do n = 1, ocn%fields%num_bcs  !{
       do m = 1, ocn%fields%bc(n)%num_fields  !{
          write(outunit,101) 'ocean%',trim(ocn%fields%bc(n)%name), &
               trim(ocn%fields%bc(n)%field(m)%name), &
               mpp_chksum(ocn%fields%bc(n)%field(m)%values)
       enddo  !} m
    enddo  !} n
101 FORMAT("   CHECKSUM::",A6,a,'%',a," = ",Z20)


100 FORMAT("   CHECKSUM::",A20," = ",Z20)
end subroutine ocean_public_type_chksum

end module ocean_model_mod
