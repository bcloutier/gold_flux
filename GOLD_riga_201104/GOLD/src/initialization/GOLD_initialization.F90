module GOLD_initialization
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
!*  file): initialize_topography, initialize_thickness,                *
!*  initialize_velocity, initialize_temp_sal,                          *
!*  initialize_mixed_layer_density, initialize_sponges,                *
!*  set_coordinate, and set_ref_profile.                               *
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
use GOLD_checksums, only : hchksum, qchksum, uchksum, vchksum, chksum
use GOLD_domains, only : pass_var, pass_vector, sum_across_PEs
use GOLD_error_handler, only : GOLD_mesg, GOLD_error, FATAL, WARNING, is_root_pe
use GOLD_file_parser, only : read_param, open_param_file, param_file_type
use GOLD_file_parser, only : log_param, log_version
use GOLD_grid, only : ocean_grid_type
use GOLD_interface_heights, only : find_eta
use GOLD_io, only : close_file, create_file, fieldtype, file_exists
use GOLD_io, only : open_file, read_data, read_axis_data, SINGLE_FILE, MULTIPLE
use GOLD_io, only : slasher, vardesc, write_field
use GOLD_grid_initialize, only : initialize_masks, set_grid_metrics
use GOLD_restart, only : restore_state, GOLD_restart_CS
use GOLD_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use GOLD_time_manager, only : time_type, set_time
use GOLD_tracer, only : add_tracer_OBC_values, advect_tracer_CS
use GOLD_variables, only : thermo_var_ptrs, directories, ocean_OBC_type
use GOLD_variables, only : OBC_NONE, OBC_SIMPLE, OBC_FLATHER_E, OBC_FLATHER_W
use GOLD_variables, only : OBC_FLATHER_N, OBC_FLATHER_S
use GOLD_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use GOLD_EOS, only : int_specific_vol_dp
use user_initialization, only : user_set_coord, user_initialize_topography
use user_initialization, only : user_initialize_thickness, user_initialize_velocity
use user_initialization, only : user_init_temperature_salinity
use user_initialization, only : user_set_Open_Bdry_Conds, user_initialize_sponges
use DOME_initialization, only : DOME_initialize_thickness
use DOME_initialization, only : DOME_initialize_topography
use DOME_initialization, only : DOME_set_Open_Bdry_Conds
use DOME_initialization, only : DOME_initialize_sponges
use benchmark_initialization, only : benchmark_initialize_thickness
use benchmark_initialization, only : benchmark_initialize_topography
use benchmark_initialization, only : benchmark_init_temperature_salinity
use circle_obcs_initialization, only : circle_obcs_initialize_thickness
implicit none ; private

#include <GOLD_memory.h>

public GOLD_initialize, Get_GOLD_Input, GOLD_initialize_topography

! This structure is to simplify communication with the calling code.
type, public :: GOLD_initialization_struct
  type(advect_tracer_CS), pointer :: advect_tracer_CSp => NULL()
  type(sponge_CS), pointer :: sponge_CSp => NULL()
  type(Compress_CS), pointer :: compress_CSp => NULL()
  type(ocean_OBC_type), pointer :: OBC => NULL()
end type GOLD_initialization_struct

contains

! -----------------------------------------------------------------------------
subroutine GOLD_initialize(u, v, h, tv, Time, G, PF, dirs, &
                              restart_CS, CS, Time_in)
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(out)   :: u
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(out)   :: v
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(out)   :: h
  type(thermo_var_ptrs),               intent(inout) :: tv
  type(time_type),                     intent(inout) :: Time
  type(ocean_grid_type),               intent(inout) :: G
  type(param_file_type),               intent(in)    :: PF
  type(directories),                   intent(in)    :: dirs
  type(GOLD_restart_CS),                pointer       :: restart_CS
  type(GOLD_initialization_struct),     intent(inout) :: CS
  type(time_type), optional,           intent(in) :: Time_in
! Arguments: u  - Zonal velocity, in m s-1.
!  (out)     v  - Meridional velocity, in m s-1.
!  (out)     h  - Layer thickness, in m.
!  (out)     tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (out)     Time    - Time at the start of the run segment.
!  (inout)   G       - The ocean's grid structure.
!  (in)      PF      - A structure indicating the open file to parse for
!                      model parameter values.
!  (in)      dirs    - A structure containing several relevant directory paths.
!  (inout)   restart_CS - A pointer to the restart control structure.
!  (inout)   CS      - A structure of pointers to be exchanged with GOLD.F90.
!  (in)      Time_in - Time at the start of the run segment. Time_in overrides
!                      any value set for Time.

  character(len=200) :: filename   ! The name of an input file.
  character(len=200) :: filename2  ! The name of an input files.
  character(len = 200) :: inputdir ! The directory where NetCDF input files are.
  character(len=200) :: config
  logical :: FROM_FILE = .true., NOT_FROM_FILE = .false.
  logical :: new_sim, write_geom
  logical :: use_temperature, bulkmixedlayer, use_sponge
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  logical :: depress_sfc ! If true, remove the mass that would be displaced
                         ! by a large surface pressure, such as with an ice
                         ! sheet.
  logical :: Analytic_FV_PGF, obsol_test
  logical :: apply_OBC_u, apply_OBC_v
!!!!!!!!!!!! Mehmet !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical :: apply_OBC_u_flather_east, apply_OBC_u_flather_west
  logical :: apply_OBC_v_flather_north, apply_OBC_v_flather_south
  integer :: i,j,k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  logical :: convert
  type(EOS_type), pointer :: eos => NULL()
  logical :: debug = .FALSE. ! indicates whether to write debugging output
  character(len=128) :: version = '$Id: GOLD_initialization.F90,v 13.0.2.33.2.53 2011/10/18 20:48:30 Alistair.Adcroft Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "GOLD_initialization" ! This module's name.
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq

  call GOLD_mesg(" GOLD_initialization.F90, GOLD_initialize: subroutine entered", 3)
  call log_version(PF, mod, version, tagname)

  new_sim = .false.
  if ((dirs%input_filename(1:1) == 'n') .and. &
      (LEN_TRIM(dirs%input_filename) == 1)) new_sim = .true.

  inputdir = "." ;  call read_param(PF, "INPUTDIR", inputdir)
  inputdir = slasher(inputdir)
  call read_param(PF, "TEMPERATURE", use_temperature, .true.)
  bulkmixedlayer = .false.; call read_param(PF,"BULKMIXEDLAYER",bulkmixedlayer)
  call read_param(PF, "DEBUG", debug)

  use_EOS = associated(tv%eqn_of_state)
  if (use_EOS) eos => tv%eqn_of_state

  call log_param(PF, mod, "TEMPERATURE", use_temperature, &
                 "If true, Temperature and salinity are used as state \n"//&
                 "variables.")
  call log_param(PF, mod, "BULKMIXEDLAYER", bulkmixedlayer, &
                 "If true, use a Kraus-Turner-like bulk mixed layer \n"//&
                 "with transitional buffer layers.", default=.false.)
  call log_param(PF, mod, "DEBUG", debug)

! ====================================================================
!    Initialize fields that are time invariant - metrics, topography,
!  masks, vertical coordinate, Coriolis parameter, a reference T & S
!  profile for compensating for compressibility.
! ====================================================================

! Set-up the layer densities, G%Rlay, and reduced gravities, G%g_prime.
  call read_param(PF, "COORD_CONFIG", config, .true.)
  call log_param(PF, mod, "COORD_CONFIG", config, &
                 "This specifies how layers are to be defined: \n"//&
                 " \t file - read coordinate information from the file \n"//&
                 " \t\t specified by (COORD_FILE).\n"//&
                 " \t ts_ref - use reference temperature and salinity \n"//&
                 " \t\t (T_REF and S_REF) to determine surface density \n"//&
                 " \t\t and GINT calculate internal densities. \n"//&
                 " \t gprime - use reference density (RHO_0) for surface \n"//&
                 " \t\t density and GINT calculate internal densities. \n"//&
                 " \t ts_profile - use temperature and salinity profiles \n"//&
                 " \t\t (read from COORD_FILE) to set layer densities. \n"//&
                 " \t USER - call a user modified routine.")
  select case ( trim(config) )
    case ("gprime")
      call set_coord_from_gprime(G%Rlay, G%g_prime, G, PF)
    case ("layer_ref")
      call set_coord_from_layer_density(G%Rlay, G%g_prime, G, PF)
    case ("ts_ref")
      call set_coord_from_ts_ref(G%Rlay, G%g_prime, G, PF, eos, tv%P_Ref)
    case ("ts_profile")
      call set_coord_from_TS_profile(G%Rlay, G%g_prime, G, PF, eos, tv%P_Ref)
    case ("file")
      call set_coord_from_file(G%Rlay, G%g_prime, G, PF)
    case ("USER")
      call user_set_coord(G%Rlay, G%g_prime, G, PF, eos)
    case default ;call GOLD_error(FATAL,"GOLD_initialize: "// &
      "Unrecognized coordinate setup"//trim(config))
  end select
  if (debug) call chksum(G%Rlay, "GOLD_initialize: Rlay ", 1, nz)
  if (debug) call chksum(G%g_prime, "GOLD_initialize: g_prime ", 1, nz)

! Set up the parameters of the physical domain (i.e. the grid), G
  call set_grid_metrics(G, PF)

! Set up the bottom depth, G%D either analytically or from file
  call GOLD_initialize_topography(G%D, G, PF)

!    This call sets seamasks that prohibit flow over any point with  !
!  a bottom that is shallower than min_depth from PF.                !
  call initialize_masks(G, PF)
  if (debug) then
    call hchksum(G%D, 'GOLD_initialize: depth ', G, haloshift=1)
    call hchksum(G%hmask, 'GOLD_initialize: hmask ', G)
    call uchksum(G%umask, 'GOLD_initialize: umask ', G)
    call vchksum(G%vmask, 'GOLD_initialize: vmask ', G)
    call qchksum(G%qmask, 'GOLD_initialize: qmask ', G)
  endif

! Modulate geometric scales according to geography.
  config = "none" ; call read_param(PF, "CHANNEL_CONFIG", config)
  call log_param(PF, mod, "CHANNEL_CONFIG", config, &
                 "A string that determines which set of channels are \n"//&
                 "restricted to specific  widths. The default (none) \n"//&
                 "has no partially restricted channels.  The other \n"//&
                 "option that is currently coded is (global_1deg).", &
                 default="none")
  select case ( trim(config) )
    case ("none")
    case ("global_1deg") ; call reset_face_lengths_global_1deg(G, PF)
    case default ; call GOLD_error(FATAL, "GOLD_initialize: "// &
      "Unrecognized channel configuration "//trim(config))
  end select

!   This call sets the topography at velocity points.
  if (G%bathymetry_at_vel) then
    config = "max" ; call read_param(PF, "VELOCITY_DEPTH_CONFIG", config)
    call log_param(PF, mod, "VELOCITY_DEPTH_CONFIG", config, &
                   "A string that determines how the topography is set at \n"//&
                   "velocity points. This may be 'min' or 'max'.", &
                   default="max")
    select case ( trim(config) )
      case ("max") ; call set_velocity_depth_max(G)
      case ("min") ; call set_velocity_depth_min(G)
      case default ; call GOLD_error(FATAL, "GOLD_initialize: "// &
        "Unrecognized velocity depth configuration "//trim(config))
    end select
  endif

!    Calculate the value of the Coriolis parameter at the latitude   !
!  of the q grid points, in s-1.
  call read_param(PF, "ROTATION", config, .true.)
  call log_param(PF, mod, "ROTATION", config, &
                 "This specifies how the Coriolis parameter is specified: \n"//&
                 " \t 2omegasinlat - Use twice the planetary rotation rate \n"//&
                 " \t\t times the sine of latitude.\n"//&
                 " \t betaplane - Use a beta-plane or f-plane. \n"//&
                 " \t USER - call a user modified routine.")
  select case (trim(config))
    case ("2omegasinlat"); call set_rotation_planetary(G%f, G, PF)
    case ("beta"); call set_rotation_beta_plane(G%f, G, PF)
    case ("betaplane"); call set_rotation_beta_plane(G%f, G, PF)
   !case ("nonrotating") ! Note from AJA: Missing case?
    case default ; call GOLD_error(FATAL,"GOLD_initialize: "// &
      "Unrecognized rotation setup "//trim(config))
  end select
  if (debug) then
    call qchksum(G%f, "GOLD_initialize: f ", G)
  endif

! Write out all of the grid data used by this run.
  write_geom=.true. ; call read_param(PF,"ALWAYS_WRITE_GEOM",write_geom)
  if (write_geom .or. new_sim) then
    call write_ocean_geometry_file(G, PF, dirs%output_directory)
    call write_vertgrid_file(G, PF, dirs%output_directory)
  endif

!    Set up the reference compressibility for the simulation, either
!  by specifying a profile or by reading a climatological state and
!  doing a 3-D best smoothed fit.
  Analytic_FV_PGF = .false.
  call read_param(PF, "ANALYTIC_FV_PGF", Analytic_FV_PGF)
  call log_param(PF, mod, "ANALYTIC_FV_PGF", Analytic_FV_PGF)
  if (use_EOS .and. (.not.Analytic_FV_PGF)) then
    call compress_init(G, PF, CS%compress_CSp)
    call set_ref_profile(G, PF, CS%compress_CSp)
  endif

!====================================================================
!    Initialize temporally evolving fields, either as initial
!  conditions or by reading them from a restart (or saves) file.
!====================================================================

  if (new_sim) then
!  This block initializes all of the fields internally.              !
    call GOLD_mesg("Run initialized internally.", 3)

    if (present(Time_in)) Time = Time_in
    ! Otherwise leave Time at its input value.

!   Initialize thickness, h
    call read_param(PF, "THICKNESS_CONFIG", config, .true.)
    call log_param(PF, mod, "THICKNESS_CONFIG", config, &
                 "A string that determines how the initial layer \n"//&
                 "thicknesses are specified for a new run: \n"//&
                 " \t file - read interface heights from the file specified \n"//&
                 " \t thickness_file - read thicknesses from the file specified \n"//&
                 " \t\t by (THICKNESS_FILE).\n"//&
                 " \t uniform - uniform thickness layers evenly distributed \n"//&
                 " \t\t between the surface and MAXIMUM_DEPTH. \n"//&
                 " \t DOME - use a slope and channel configuration for the \n"//&
                 " \t\t DOME sill-overflow test case. \n"//&
                 " \t benchmark - use the benchmark test case thicknesses. \n"//&
                 " \t search - search a density profile for the interface \n"//&
                 " \t\t densities. This is not yet implemented. \n"//&
                 " \t circle_obcs - the circle_obcs test case is used. \n"//&
                 " \t USER - call a user modified routine.")
    select case (trim(config))
      case ("file");    call initialize_thickness_from_file(h, G, PF, .false.)
      case ("thickness_file"); call initialize_thickness_from_file(h, G, PF, .true.)
      case ("uniform"); call initialize_thickness_uniform(h, G, PF)
      case ("DOME");    call DOME_initialize_thickness(h, G, PF)
      case ("benchmark"); call benchmark_initialize_thickness(h, G, PF, &
                                   tv%eqn_of_state, tv%P_Ref)
      case ("search");  call initialize_thickness_search
      case ("circle_obcs"); call circle_obcs_initialize_thickness(h, G, PF)
      case ("USER");    call user_initialize_thickness(h, G, PF,tv%T)
      case default ; call GOLD_error(FATAL,  "GOLD_initialize: "//&
             "Unrecognized layer thickness configuration "//trim(config))
    end select
    call pass_var(h, G%Domain)
    if (debug) call hchksum(h, "GOLD_initialize: h ", G)

!   Initialize velocity components, u and v
    config = "zero" ; call read_param(PF, "VELOCITY_CONFIG", config)
    call log_param(PF, mod, "VELOCITY_CONFIG", config, &
                 "A string that determines how the initial velocities \n"//&
                 "are specified for a new run: \n"//&
                 " \t file - read velocities from the file specified \n"//&
                 " \t\t by (VELOCITY_FILE). \n"//&
                 " \t zero - the fluid is initially at rest. \n"//&
                 " \t uniform - the flow is uniform (determined by\n"//&
                 " \t\t paremters TORUS_U and TORUS_V).\n"//&
                 " \t USER - call a user modified routine.")
    select case (trim(config))
      case ("file"); call initialize_velocity_from_file(u, v, G, PF)
      case ("zero"); call initialize_velocity_zero(u, v, G, PF)
      case ("uniform"); call initialize_velocity_uniform(u, v, G, PF)
      case ("circular"); call initialize_velocity_circular(u, v, G, PF)
      case ("USER"); call user_initialize_velocity(u, v, G, PF)
      case default ; call GOLD_error(FATAL,  "GOLD_initialize: "//&
             "Unrecognized velocity configuration "//trim(config))
    end select

    call pass_vector(u, v, G%Domain)
    if (debug) call uchksum(u, "GOLD_initialize: u ", G)
    if (debug) call vchksum(v, "GOLD_initialize: v ", G)

!   Initialize temperature and salinity (T and S)
    if ( use_temperature ) then
      call read_param(PF, "TS_CONFIG", config, .true.)
      call log_param(PF, mod, "TS_CONFIG", config, &
                 "A string that determines how the initial tempertures \n"//&
                 "and salinities are specified for a new run: \n"//&
                 " \t file - read velocities from the file specified \n"//&
                 " \t\t by (TS_FILE). \n"//&
                 " \t fit - find the temperatures that are consistent with \n"//&
                 " \t\t the layer densities and salinity S_REF. \n"//&
                 " \t TS_profile - use temperature and salinity profiles \n"//&
                 " \t\t (read from TS_FILE) to set layer densities. \n"//&
                 " \t benchmark - use the benchmark test case T & S. \n"//&
                 " \t USER - call a user modified routine.")
      select case (trim(config))
        case ("fit"); call initialize_temp_salt_fit(tv%T, tv%S, G, PF, eos, tv%P_Ref)
        case ("file"); call initialize_temp_salt_from_file(tv%T, tv%S, G, PF)
        case ("benchmark"); call benchmark_init_temperature_salinity(tv%T, tv%S, &
                                     G, PF, eos, tv%P_Ref)
        case ("TS_profile") ; call initialize_temp_salt_from_profile(tv%T, tv%S, G, PF)
        case ("USER"); call user_init_temperature_salinity(tv%T, tv%S, G, PF, eos)
        case default ; call GOLD_error(FATAL,  "GOLD_initialize: "//&
               "Unrecognized Temp & salt configuration "//trim(config))
      end select
      if (debug) then
        call hchksum(tv%T, "GOLD_initialize: T ", G)
        call hchksum(tv%S, "GOLD_initialize: S ", G)
      endif
    endif

!   The mixed layer density is typically determined from T & S.       !
    if (bulkmixedlayer .and. .not.use_eos) then
      if (.not.associated(tv%Rml)) call GOLD_error(FATAL, &
        "GOLD_initialize: tv%Rml is not assoc'd in a run with a bulk ML.")

      if (use_temperature) then
        call initialize_mixed_layer_density(tv%Rml, G, PF, .true., eos, &
                                            tv%T, tv%S, tv%P_Ref)
      else
        call initialize_mixed_layer_density(tv%Rml, G, PF, .false.)
      endif
      if (debug) call hchksum(tv%Rml, "GOLD_initialize: Rml ", G)
    endif

!   Optionally convert the thicknesses from m to kg m-2.  This is particularly
! useful in a non-Boussinesq model.
    convert = .false. ; call read_param(PF, "CONVERT_THICKNESS_UNITS", convert)
    call log_param(PF, mod, "CONVERT_THICKNESS_UNITS", convert, &
                 "If true,  convert the thickness initial conditions from \n"//&
                 "units of m to kg m-2 or vice versa, depending on whether \n"//&
                 "BOUSSINESQ is defined. This does not apply if a restart \n"//&
                 "file is read.", default=.false.)
    if (convert) call convert_thickness(h, G, PF, tv)

!  Remove the mass that would be displaced by an ice shelf or inverse barometer.
    depress_sfc = .false.
    call read_param(PF, "DEPRESS_INITIAL_SURFACE", depress_sfc)
    call log_param(PF, mod, "DEPRESS_INITIAL_SURFACE", depress_sfc, &
                 "If true,  depress the initial surface to avoid huge \n"//&
                 "tsunamis when a large surface pressure is applied.", &
                 default=.false.)
    if (depress_sfc) call depress_surface(h, G, PF, tv)

  else ! Previous block for new_sim=.T., this block restores state
!    This line calls a subroutine that reads the initial conditions  !
!  from a previously generated file.                                 !
    call restore_state(dirs%input_filename, dirs%restart_input_dir, Time, &
                       G, restart_CS)
    if (present(Time_in)) Time = Time_in
  endif

  use_sponge = .false. ; call read_param(PF,"SPONGE",use_sponge)
  call log_param(PF, mod, "SPONGE", use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. \n"//&
                 "The exact location and properties of those sponges are \n"//&
                 "specified via SPONGE_CONFIG.", default=.false.)
  if ( use_sponge ) then
! The 3 arguments here are (1) a flag indicating whether the sponge  !
! values are to be read from files, (2) the name of a file containing!
! the state toward which the model is damped, and (3) the file in    !
! which the 2-D damping rate field can be found.                     !
    call read_param(PF, "SPONGE_CONFIG", config, .true.)
    call log_param(PF, mod, "SPONGE_CONFIG", config, &
                 "A string that sets how the sponges are configured: \n"//&
                 " \t file - read sponge properties from the file \n"//&
                 " \t\t specified by (SPONGE_FILE).\n"//&
                 " \t DOME - use a slope and channel configuration for the \n"//&
                 " \t\t DOME sill-overflow test case. \n"//&
                 " \t USER - call a user modified routine.", default="file")
    select case (trim(config))
      case ("DOME"); call DOME_initialize_sponges(G, tv, PF, CS%sponge_CSp)
      case ("USER"); call user_initialize_sponges(G, use_temperature, &
                              bulkmixedlayer, tv, PF, CS%sponge_CSp, h)
      case ("file"); call initialize_sponges_file(G, use_temperature, &
                              bulkmixedlayer, tv, PF, CS%sponge_CSp)
      case default ; call GOLD_error(FATAL,  "GOLD_initialize: "//&
             "Unrecognized sponge configuration "//trim(config))
    end select
  endif

! This subroutine call sets optional open boundary conditions.
  apply_OBC_u = .false. ; call read_param(PF, "APPLY_OBC_U", apply_OBC_u)
  apply_OBC_v = .false. ; call read_param(PF, "APPLY_OBC_V", apply_OBC_v)
  call log_param(PF, mod, "APPLY_OBC_U", apply_OBC_u, &
                 "If true, open boundary conditions may be set at some \n"//&
                 "u-points, with the configuration controlled by OBC_CONFIG", &
                 default=.false.)
  call log_param(PF, mod, "APPLY_OBC_V", apply_OBC_v, &
                 "If true, open boundary conditions may be set at some \n"//&
                 "v-points, with the configuration controlled by OBC_CONFIG", &
                 default=.false.)
  if (apply_OBC_u .or. apply_OBC_v) then 
    call read_param(PF, "OBC_CONFIG", config, .true.)
    call log_param(PF, mod, "OBC_CONFIG", config, &
                 "A string that sets how the open boundary conditions are \n"//&
                 " configured: \n"//&
                 " \t DOME - use a slope and channel configuration for the \n"//&
                 " \t\t DOME sill-overflow test case. \n"//&
                 " \t USER - call a user modified routine.", default="file")
    if (trim(config) == "DOME") then
      call DOME_set_Open_Bdry_Conds(CS%OBC, tv, G, PF, CS%advect_tracer_CSp)
    elseif (trim(config) == "USER") then
      call user_set_Open_Bdry_Conds(CS%OBC, tv, G, PF, CS%advect_tracer_CSp)
    else
      call GOLD_error(FATAL, "The open boundary conditions specified by "//&
              "OBC_CONFIG = "//trim(config)//" have not been fully implemented.")
      call set_Open_Bdry_Conds(CS%OBC, tv, G, PF, CS%advect_tracer_CSp)
    endif
  endif

  apply_OBC_u_flather_east = .false. ; apply_OBC_u_flather_west = .false.
  apply_OBC_v_flather_north = .false. ; apply_OBC_v_flather_south = .false.
  call read_param(PF,"APPLY_OBC_U_FLATHER_EAST",apply_OBC_u_flather_east)
  call read_param(PF,"APPLY_OBC_U_FLATHER_WEST",apply_OBC_u_flather_west)
  call read_param(PF,"APPLY_OBC_V_FLATHER_NORTH",apply_OBC_v_flather_north)
  call read_param(PF,"APPLY_OBC_V_FLATHER_SOUTH",apply_OBC_v_flather_south)

  call log_param(PF, mod, "APPLY_OBC_U_FLATHER_EAST", apply_OBC_u_flather_east,&
                 "Apply a Flather open boundary condition on the eastern \n"//&
                 "side of the global domain", default=.false.)
  call log_param(PF, mod, "APPLY_OBC_U_FLATHER_WEST", apply_OBC_u_flather_west,&
                 "Apply a Flather open boundary condition on the western \n"//&
                 "side of the global domain", default=.false.)
  call log_param(PF, mod, "APPLY_OBC_V_FLATHER_NORTH", apply_OBC_v_flather_north,&
                 "Apply a Flather open boundary condition on the northern \n"//&
                 "side of the global domain", default=.false.)
  call log_param(PF, mod, "APPLY_OBC_V_FLATHER_SOUTH", apply_OBC_v_flather_south,&
                 "Apply a Flather open boundary condition on the southern \n"//&
                 "side of the global domain", default=.false.)
  if (apply_OBC_u_flather_east .or. apply_OBC_u_flather_west .or. apply_OBC_v_flather_north .or. apply_OBC_v_flather_south) then
    call set_Flather_Bdry_Conds(CS%OBC, tv, G, PF, CS%advect_tracer_CSp)
  endif

  call GOLD_mesg(" GOLD_initialization.F90, GOLD_initialize: complete", 3)

end subroutine GOLD_initialize
! -----------------------------------------------------------------------------


! -----------------------------------------------------------------------------
subroutine set_coord_from_gprime(Rlay, g_prime, G, param_file)
  real, dimension(:), intent(out)   :: Rlay, g_prime
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real :: g_int   ! Reduced gravities across the internal interfaces, in m s-2.
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  character(len=40)  :: mod = "set_coord_from_gprime" ! This subroutine's name.
  integer :: k, nz
  nz = G%ke

  call GOLD_mesg("  GOLD_initialization.F90, set_coord_from_gprime: setting coordinate", 5)

  g_fs = G%g_Earth ; call read_param(param_file, "GFS", g_fs)
  call read_param(param_file, "GINT", g_int, .true.)
  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = g_int ; enddo
  Rlay(1) = G%Rho0
  do k=2,nz ; Rlay(k) = Rlay(k-1) + g_prime(k)*(G%Rho0/G%g_Earth) ; enddo

  call log_param(param_file, mod, "GFS" , g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call log_param(param_file, mod, "GINT", g_int, &
                 "The reduced gravity across internal interfaces.", &
                 units="m s-2")

end subroutine set_coord_from_gprime
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_from_layer_density(Rlay, g_prime, G, param_file)
  real, dimension(:), intent(out)   :: Rlay, g_prime
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  real :: Rlay_Ref! The surface layer's target density, in kg m-3.
  real :: RLay_range ! The range of densities, in kg m-3.
  character(len=40)  :: mod = "set_coord_from_layer_density" ! This subroutine's name.
  integer :: k, nz
  nz = G%ke

  call GOLD_mesg("  GOLD_initialization.F90, set_coord_from_layer_density: setting coordinate", 5)

  call read_param(param_file, "RLAY_REF", Rlay_Ref, .true.)
  RLay_range = 2.0 ; call read_param(param_file, "RLAY_RANGE", Rlay_range)
  g_fs = G%g_Earth ; call read_param(param_file, "GFS", g_fs)
  g_prime(1) = g_fs
  Rlay(1) = Rlay_Ref
  do k=2,nz
     Rlay(k) = Rlay(k-1) + RLay_range/(real(nz-1))
  enddo
!    These statements set the interface reduced gravities.           !
  do k=2,nz
     g_prime(k) = (G%g_Earth/G%Rho0) * (Rlay(k) - Rlay(k-1))
  enddo

  call log_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call log_param(param_file, mod, "RLAY_REF", Rlay_Ref, &
                 "The reference potential density used for layer 1.", &
                 units="kg m-3")
  call log_param(param_file, mod, "RLAY_RANGE", Rlay_range, &
                 "The range of reference potential densities in the layers.", &
                 units="kg m-3", default=2.0)
end subroutine set_coord_from_layer_density
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_from_TS_ref(Rlay, g_prime, G, param_file, eqn_of_state, &
                                 P_Ref)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(param_file_type), intent(in)  :: param_file
  type(EOS_type),        pointer     :: eqn_of_state
  real,                  intent(in)  :: P_Ref
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      eqn_of_state - integer selecting the equation of state
!  (in)      P_Ref - The coordinate-density reference pressure in Pa.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real :: T_ref   ! Reference temperature
  real :: S_ref   ! Reference salinity
  real :: g_int   ! Reduced gravities across the internal interfaces, in m s-2.
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  character(len=40)  :: mod = "set_coord_from_TS_ref" ! This subroutine's name.
  integer :: k, nz
  nz = G%ke

  call GOLD_mesg("  GOLD_initialization.F90, set_coord_from_TS_ref: setting coordinate", 5)

  call read_param(param_file,"T_REF",T_Ref,.true.)
  call read_param(param_file,"S_REF",S_Ref,.true.)

  g_fs = G%g_Earth ; call read_param(param_file, "GFS", g_fs)
                                      !
!    These statements set the interface reduced gravities.           !
  g_prime(1) = g_fs
  call read_param(param_file,"GINT",g_int,.true.)
  do k=2,nz ; g_prime(k) = g_int ; enddo

!    The uppermost layer's density is set here.  Subsequent layers'  !
!  densities are determined from this value and the g values.        !
!        T0 = 28.228 ; S0 = 34.5848 ; Pref = P_Ref
  call calculate_density(T_ref, S_ref, P_ref, Rlay(1), 1,1, eqn_of_state)

!    These statements set the layer densities.                       !
  do k=2,nz ; Rlay(k) = Rlay(k-1) + g_prime(k)*(G%Rho0/G%g_Earth) ; enddo

  call log_param(param_file, mod, "T_REF", T_Ref, &
                 "The initial temperature of the lightest layer.", units="degC")
  call log_param(param_file, mod, "S_REF", S_Ref, &
                 "The initial salinities.", units="PSU")
  call log_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call log_param(param_file, mod, "GINT", g_int, &
                 "The reduced gravity across internal interfaces.", &
                 units="m s-2")

end subroutine set_coord_from_TS_ref
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_from_TS_profile(Rlay, g_prime, G, param_file, &
                                     eqn_of_state, P_Ref)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(param_file_type), intent(in)  :: param_file
  type(EOS_type),        pointer     :: eqn_of_state
  real,                  intent(in)  :: P_Ref
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      eqn_of_state - integer that selects equation of state
!  (in)      P_Ref - The coordinate-density reference pressure in Pa.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real, dimension(SZK_(G)) :: T0, S0,  Pref
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  integer :: k, nz
  character(len=40)  :: mod = "set_coord_from_TS_profile" ! This subroutine's name.
  character(len=200) :: filename,coord_file,inputdir ! Strings for file/path
  nz = G%ke

  call GOLD_mesg("  GOLD_initialization.F90, set_coord_from_TS_profile: setting coordinate", 5)

  g_fs = G%g_Earth ; call read_param(param_file, "GFS", g_fs)

  inputdir = "." ; call read_param(param_file,"INPUTDIR",inputdir)
  inputdir = slasher(inputdir)
  call read_param(param_file, "COORD_FILE", coord_file, .true.)
  filename = trim(inputdir)//trim(coord_file)
  call read_data(filename,"PTEMP",T0(:),domain=G%Domain%mpp_domain)
  call read_data(filename,"SALT",S0(:),domain=G%Domain%mpp_domain)

  if (.not.file_exists(filename)) call GOLD_error(FATAL, &
      " set_coord_from_TS_profile: Unable to open " //trim(filename))
!    These statements set the interface reduced gravities.           !
  g_prime(1) = g_fs
  do k=1,nz ; Pref(k) = P_ref ; enddo
  call calculate_density(T0, S0, Pref, Rlay, 1,nz,eqn_of_state)
  do k=2,nz; g_prime(k) = (G%g_Earth/G%Rho0) * (Rlay(k) - Rlay(k-1)); enddo

  call log_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call log_param(param_file, mod, "COORD_FILE", coord_file, &
                 "The file from which the coordinate temperatures and \n"//&
                 "salnities are read.")
  call log_param(param_file, mod, "INPUTDIR/COORD_FILE", filename)

end subroutine set_coord_from_TS_profile
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_from_file(Rlay, g_prime, G, param_file)
  real, dimension(:), intent(out)   :: Rlay, g_prime
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  integer :: k, nz
  character(len=40)  :: mod = "set_coord_from_file" ! This subroutine's name.
  character(len=40)  :: coord_var
  character(len=200) :: filename,coord_file,inputdir ! Strings for file/path
  nz = G%ke

  call GOLD_mesg("  GOLD_initialization.F90, set_coord_from_file: setting coordinate", 5)

  inputdir = "." ; call read_param(param_file, "INPUTDIR", inputdir)
  inputdir = slasher(inputdir)
  call read_param(param_file,"COORD_FILE", coord_file, .true.)
  g_fs = G%g_Earth ; call read_param(param_file, "GFS", g_fs)
  coord_var = "Layer" ; call read_param(param_file, "COORD_VAR", coord_var)

  filename = trim(inputdir)//trim(coord_file)
  if (.not.file_exists(filename)) call GOLD_error(FATAL, &
      " set_coord_from_file: Unable to open "//trim(filename))

  call read_axis_data(filename, coord_var, Rlay)
  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = (G%g_Earth/G%Rho0) * (Rlay(k) - Rlay(k-1)) ; enddo
  do k=1,nz ; if (g_prime(k) <= 0.0) then
    call GOLD_error(FATAL, "GOLD_initialization set_coord_from_file: "//&
       "Zero or negative g_primes read from variable "//"Layer"//" in file "//&
       trim(filename))
  endif ; enddo

  call log_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call log_param(param_file, mod, "COORD_FILE", coord_file, &
                 "The file from which the coordinate densities are read.")
  call log_param(param_file, mod, "COORD_VAR", coord_var, &
                 "The variable in COORD_FILE that is to be used for the \n"//&
                 "coordinate densities.", default="Layer")
  call log_param(param_file, mod, "INPUTDIR/COORD_FILE", filename)
end subroutine set_coord_from_file
! -----------------------------------------------------------------------------

subroutine GOLD_initialize_topography(D, G, PF)
  real, intent(out), dimension(NXMEM_,NYMEM_) :: D
  type(ocean_grid_type), intent(in)           :: G
  type(param_file_type), intent(in)           :: PF
! Arguments: D  - the bottom depth in m. Intent out.
!  (in)      G  - The ocean's grid structure.
!  (in)      PF - A structure indicating the open file to parse for
!                         model parameter values.

!  This subroutine makes the appropriate call to set up the bottom depth.
!  This is a separate subroutine so that it can be made public and shared with
!  the ice-sheet code or other components.
! Set up the bottom depth, G%D either analytically or from file
  character(len=40)  :: mod = "GOLD_initialize_topography" ! This subroutine's name.
  character(len=200) :: config

  call read_param(PF, "TOPO_CONFIG", config, .true.)
  call log_param(PF, mod, "TOPO_CONFIG", config, &
                 "This specifies how bathymetry is specified: \n"//&
                 " \t file - read bathymetric information from the file \n"//&
                 " \t\t specified by (TOPO_FILE).\n"//&
                 " \t flat - flat bottom set to MAXIMUM_DEPTH. \n"//&
                 " \t bowl - an analytically specified bown-shaped basin \n"//&
                 " \t\t ranging between MAXIMUM_DEPTH and MINIMUM_DEPTH. \n"//&
                 " \t spoon - a similar shape to 'bowl', but with an vertical \n"//&
                 " \t\t wall at the southern face. \n"//&
                 " \t halfpipe - a zonally uniform channel with a half-sine \n"//&
                 " \t\t profile in the meridional direction. \n"//&
                 " \t DOME - use a slope and channel configuration for the \n"//&
                 " \t\t DOME sill-overflow test case. \n"//&
                 " \t benchmark - use the benchmark test case topography. \n"//&
                 " \t USER - call a user modified routine.")
  select case ( trim(config) )
    case ("file");  call initialize_topography_from_file(D, G, PF)
    case ("flat");  call initialize_topography_flat(D, G, PF)
 !   case ("flat");     call initialize_topography_named(D, G, PF, config)
    case ("spoon");    call initialize_topography_named(D, G, PF, config)
    case ("bowl");     call initialize_topography_named(D, G, PF, config)
    case ("halfpipe"); call initialize_topography_named(D, G, PF, config)
    case ("DOME");  call DOME_initialize_topography(D, G, PF)
    case ("benchmark");  call benchmark_initialize_topography(D, G, PF)
    case ("USER");  call user_initialize_topography(D, G, PF)
    case default ;  call GOLD_error(FATAL,"GOLD_initialize: "// &
      "Unrecognized topography setup "//trim(config))
  end select
  if (trim(config) .ne. "DOME") then
    call limit_topography(D, G, PF)
  endif
  
end subroutine GOLD_initialize_topography

! -----------------------------------------------------------------------------
subroutine initialize_topography_from_file(D, G, param_file )
  real, intent(out), dimension(NXMEM_,NYMEM_) :: D
  type(ocean_grid_type), intent(in)           :: G
  type(param_file_type), intent(in)           :: param_file
! Arguments: D          - the bottom depth in m. Intent out.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

!  This subroutine reads depths from a file and puts it into D(:,:) in m.
  character(len=200) :: filename,topo_file,inputdir ! Strings for file/path
  character(len=200) :: grid_config                 ! Grid format
  character(len=200) :: topo_varname, topo_var_dflt ! Variable name in file
  character(len=40)  :: mod = "initialize_topography_from_file" ! This subroutine's name.

  call GOLD_mesg("  GOLD_initialization.F90, initialize_topography_from_file: reading topography", 5)

  inputdir = "." ; call read_param(param_file,"INPUTDIR",inputdir)
  inputdir = slasher(inputdir)
  call read_param(param_file,"GRID_CONFIG",grid_config)
  topo_file = "topog.nc" ; call read_param(param_file,"TOPO_FILE",topo_file)
  select case (trim(grid_config))
    case ("mosaic");
      ! The Mosaic file provides the link to the topography
      ! which we ought to read from the mosaic at this point. TBD. (AJA)
      topo_var_dflt = "depth"
    case default;
      topo_var_dflt = 'depth_t'
  end select
  topo_varname = topo_var_dflt
  call read_param(param_file,"TOPO_VARNAME",topo_varname)

  filename = trim(inputdir)//trim(topo_file)
  if (.not.file_exists(filename, G%Domain)) call GOLD_error(FATAL, &
       " initialize_topography_from_file: Unable to open "//trim(filename))

  call read_data(filename,trim(topo_varname),D,domain=G%Domain%mpp_domain)

  call log_param(param_file, mod, "TOPO_FILE", topo_file, &
                 "The file from which the bathymetry is read.", &
                 default="topog.nc")
  call log_param(param_file, mod, "INPUTDIR/TOPO_FILE", filename)
  call log_param(param_file, mod, "TOPO_VARNAME", topo_varname, &
                 "The name of the bathymetry variable in TOPO_FILE.", &
                 default=topo_var_dflt)
end subroutine initialize_topography_from_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_topography_flat(D, G, param_file)
  real, intent(out), dimension(NXMEM_,NYMEM_) :: D
  type(ocean_grid_type), intent(in)           :: G
  type(param_file_type), intent(in)           :: param_file
! Arguments: D          - the bottom depth in m. Intent out.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine places the bottom depth in m into D(:,:), shaped flat.
  real :: max_depth ! The maximum depths in m.
  character(len=40)  :: mod = "initialize_topography_flat" ! This subroutine's name.
  integer :: i, j

  call GOLD_mesg("  GOLD_initialization.F90, initialize_topography_flat: setting topography", 5)

  call read_param(param_file, "MAXIMUM_DEPTH", max_depth, .true.)

  do i=G%isc,G%iec ; do j=G%jsc,G%jec; D(i,j) = max_depth; enddo ; enddo

  call log_param(param_file, mod, "MAXIMUM_DEPTH", max_depth, &
                 "The maximum depth of the ocean.", units="m")

end subroutine initialize_topography_flat
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_topography_named(D, G, param_file, topog_config)
  real, intent(out), dimension(NXMEM_,NYMEM_) :: D
  type(ocean_grid_type), intent(in)           :: G
  type(param_file_type), intent(in)           :: param_file
  character(len=*),      intent(in)           :: topog_config
! Arguments: D          - the bottom depth in m. Intent out.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      topog_config - The name of an idealized topographic configuration.

! This subroutine places the bottom depth in m into D(:,:), shaped in a spoon
  real :: min_depth, max_depth ! The minimum and maximum depths in m.
  real :: PI                   ! 3.1415926... calculated as 4*atan(1)
  real :: D0                   ! A constant to make the maximum     !
                               ! basin depth MAXIMUM_DEPTH.         !
  real :: expdecay             ! A decay scale of associated with   !
                               ! the sloping boundaries, in m.      !
  real :: Dedge                ! The depth in m at the basin edge.  !
  real :: south_lat, west_lon, len_lon, len_lat, Rad_earth
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed, xhalo, yhalo
  character(len=40)  :: mod = "initialize_topography_named" ! This subroutine's name.
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call GOLD_mesg("  GOLD_initialization.F90, initialize_topography_named: "//&
                 "setting topography "//trim(topog_config), 5)

  min_depth = 0.0 ; call read_param(param_file, "MINIMUM_DEPTH", min_depth)
  call read_param(param_file, "MAXIMUM_DEPTH", max_depth, .true.)

  Dedge = 100.0 ; call read_param(param_file, "EDGE_DEPTH", Dedge)
  expdecay = 400000.0 ; call read_param(param_file, "TOPOG_SLOPE_SCALE", expdecay)

  ! These expressions force rounding of approximate values in
  ! consistent way.
  xhalo = G%isc-G%isd ; yhalo = G%jsc-G%jsd
  south_lat = ANInt(G%gridlatq(yhalo)*1024.0)/1024.0
  west_lon = ANInt(G%gridlonq(xhalo)*1024.0)/1024.0
  len_lon = ANInt((G%gridlonq(G%Domain%nxtot+xhalo)-G%gridlonq(xhalo))*1024.0)/1024.0
  len_lat = ANInt((G%gridlatq(G%Domain%nytot+yhalo)-G%gridlatq(yhalo))*1024.0)/1024.0

  call read_param(param_file,"SOUTHLAT",south_lat)
  call read_param(param_file,"LENLAT",len_lat)
  call read_param(param_file,"WESTLON",west_lon)
  call read_param(param_file,"LENLON",len_lon)
  call read_param(param_file,"RAD_EARTH",Rad_Earth,.true.)

  PI = 4.0*atan(1.0)

  if (trim(topog_config) == "flat") then
    do i=is,ie ; do j=js,je ; D(i,j) = max_depth ; enddo ; enddo
  elseif (trim(topog_config) == "spoon") then
    D0 = (max_depth - Dedge) / &
             ((1.0 - exp(-0.5*len_lat*Rad_earth*PI/(180.0 *expdecay))) * &
              (1.0 - exp(-0.5*len_lat*Rad_earth*PI/(180.0 *expdecay))))
    do i=is,ie ; do j=js,je
  !  This sets a bowl shaped (sort of) bottom topography, with a       !
  !  maximum depth of max_depth.                                   !
      D(i,j) =  Dedge + D0 * &
             (sin(PI * (G%geolonh(i,j) - (west_lon)) / len_lon) * &
           (1.0 - exp((G%geolath(i,j) - (south_lat+len_lat))*Rad_earth*PI / &
                      (180.0*expdecay)) ))
    enddo ; enddo
  elseif (trim(topog_config) == "spoon") then
    D0 = (max_depth - Dedge) / &
             ((1.0 - exp(-0.5*len_lat*Rad_earth*PI/(180.0 *expdecay))) * &
              (1.0 - exp(-0.5*len_lat*Rad_earth*PI/(180.0 *expdecay))))

  !  This sets a bowl shaped (sort of) bottom topography, with a
  !  maximum depth of max_depth.
    do i=is,ie ; do j=js,je
      D(i,j) =  Dedge + D0 * &
             (sin(PI * (G%geolonh(i,j) - west_lon) / len_lon) * &
             ((1.0 - exp(-(G%geolath(i,j) - south_lat)*Rad_Earth*PI/ &
                          (180.0*expdecay))) * &
             (1.0 - exp((G%geolath(i,j) - (south_lat+len_lat))* &
                         Rad_Earth*PI/(180.0*expdecay)))))
    enddo ; enddo
  elseif (trim(topog_config) == "halfpipe") then
    D0 = max_depth - Dedge
    do i=is,ie ; do j=js,je
      D(i,j) =  Dedge + D0 * ABS(sin(PI*(G%geolath(i,j) - south_lat)/len_lat))
    enddo ; enddo
  else
    call GOLD_error(FATAL,"initialize_topography_named: "// &
      "Unrecognized topography name "//trim(topog_config))
  endif

  ! This is here just for safety.  Hopefully it doesn't do anything.
  do i=is,ie ; do j=js,je
    if (D(i,j) > max_depth) D(i,j) = max_depth
    if (D(i,j) < min_depth) D(i,j) = 0.5*min_depth
  enddo ; enddo

  call log_param(param_file, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)
  call log_param(param_file, mod, "MAXIMUM_DEPTH", max_depth, &
                 "The maximum depth of the ocean.", units="m")
  if (trim(topog_config) /= "flat") then
    call log_param(param_file, mod, "EDGE_DEPTH", Dedge, &
                   "The depth at the edge of one of the named topographies.", &
                   units="m", default=100.0)
    call log_param(param_file, mod, "SOUTHLAT", south_lat, &
                   "The southern latitude of the domain.", units="degrees")
    call log_param(param_file, mod, "LENLAT", len_lat, &
                   "The latitudinal length of the domain.", units="degrees")
    call log_param(param_file, mod, "WESTLON", west_lon, &
                   "The western longitude of the domain.", units="degrees")
    call log_param(param_file, mod, "LENLON", len_lon, &
                   "The longitudinal length of the domain.", units="degrees")
    call log_param(param_file, mod, "RAD_EARTH", Rad_Earth, &
                   "The radius of the Earth.", units="m")
    call log_param(param_file, mod, "TOPOG_SLOPE_SCALE", Dedge, &
                   "The exponential decay scale used in defining some of \n"//&
                   "the named topographies.", units="m", default=400000.0)
  endif

end subroutine initialize_topography_named
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine limit_topography(D, G, param_file)
  real, intent(inout), dimension(NXMEM_,NYMEM_) :: D
  type(ocean_grid_type), intent(in)             :: G
  type(param_file_type), intent(in)             :: param_file
! Arguments: D          - the bottom depth in m. Intent in/out.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine ensures that    min_depth < D(x,y) < max_depth
  integer :: i, j
  character(len=40)  :: mod = "limit_topography" ! This subroutine's name.
  real :: min_depth, max_depth

  call GOLD_mesg("  GOLD_initialization.F90, limit_topography: limiting topography", 5)

  min_depth = 0.0 ; call read_param(param_file,"MINIMUM_DEPTH",min_depth)
  max_depth = 1.0e6 ; call read_param(param_file,"MAXIMUM_DEPTH",max_depth)

! Make sure that min_depth < D(x,y) < max_depth
  do j=G%jsd,G%jed ; do i=G%isd,G%ied
    D(i,j) = min( max( D(i,j), 0.5*min_depth ), max_depth )
  enddo ; enddo

  call log_param(param_file, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)
  call log_param(param_file, mod, "MAXIMUM_DEPTH", max_depth, &
                 "The maximum depth of the ocean.", units="m", default=1.0e6)
end subroutine limit_topography
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_rotation_planetary(f, G, param_file)
  real, dimension(NXMEMQ_,NYMEMQ_), intent(out) :: f
  type(ocean_grid_type),            intent(in)  :: G
  type(param_file_type),            intent(in)  :: param_file
! Arguments: f          - Coriolis parameter (vertical component) in s^-1
!     (in)   G          - grid type
!     (in)   param_file - parameter file type

! This subroutine sets up the Coriolis parameter for a sphere
  integer :: I, J
  real    :: PI, omega

  call GOLD_mesg("  GOLD_initialization.F90, set_rotation_planetary: setting f (Coriolis)", 5)

  call read_param(param_file,"OMEGA",omega,.true.)
  PI = 4.0*atan(1.0)

  do I=G%Isdq,G%Iedq ; do J=G%Jsdq,G%Jedq
    f(I,J) = ( 2.0 * omega ) * sin( ( PI * G%geolatq(I,J) ) / 180.)
  enddo ; enddo

  call log_param(param_file, "set_rotation_planetary", "OMEGA", omega, &
                 "The rotation rate of the earth.", units="s-1")
end subroutine set_rotation_planetary
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_rotation_beta_plane(f, G, param_file)
  real, dimension(NXMEMQ_,NYMEMQ_), intent(out) :: f
  type(ocean_grid_type),            intent(in)  :: G
  type(param_file_type),            intent(in)  :: param_file
! Arguments: f          - Coriolis parameter (vertical component) in s^-1
!     (in)   G          - grid type
!     (in)   param_file - parameter file type

! This subroutine sets up the Coriolis parameter for a beta-plane
  integer :: I, J
  real    :: f_0, beta, y_scl, Rad_Earth, PI
  character(len=40)  :: mod = "set_rotation_beta_plane" ! This subroutine's name.
  character(len=200) :: axis_units

  call GOLD_mesg("  GOLD_initialization.F90, set_rotation_beta_plane: setting f (Coriolis)", 5)

  f_0 = 0.0 ; call read_param(param_file,"F_0",f_0)
  beta = 0.0 ; call read_param(param_file,"BETA",beta)
  axis_units = "degrees"; call read_param(param_file,"AXIS_UNITS",axis_units)

  PI = 4.0*atan(1.0)
  select case (axis_units(1:1))
    case ("d")
      call read_param(param_file,"RAD_EARTH",Rad_Earth,.true.)
      y_scl = Rad_Earth/PI
    case ("k"); y_scl = 1.E3
    case ("m"); y_scl = 1.
    case ("c"); y_scl = 1.E-2
    case default ; call GOLD_error(FATAL, &
      " set_rotation_beta_plane: unknown AXIS_UNITS = "//trim(axis_units))
  end select

  do I=G%Isdq,G%Iedq ; do J=G%Jsdq,G%Jedq
    f(I,J) = f_0 + beta * ( G%geolatq(I,J) * y_scl )
  enddo ; enddo

  call log_param(param_file, mod, "F_0", f_0, &
                 "The reference value of the Coriolis parameter with the \n"//&
                 "betaplane option.", units="s-1", default=0.0)
  call log_param(param_file, mod, "BETA", beta, &
                 "The northward gradient of the Coriolis parameter with \n"//&
                 "the betaplane option.", units="m-1 s-1", default=0.0)
  call log_param(param_file, mod, "AXIS_UNITS", axis_units)
end subroutine set_rotation_beta_plane
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_thickness_from_file(h, G, param_file, file_has_thickness)
  real, intent(out), dimension(NXMEM_,NYMEM_, NZ_) :: h
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  logical,               intent(in) :: file_has_thickness
! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      file_has_thickness - If true, this file contains thicknesses;
!                                 otherwise it contains interface heights.

!  This subroutine reads the layer thicknesses from file.
  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1)
  integer :: inconsistent = 0
  real :: dilate     ! The amount by which each layer is dilated to agree
                     ! with the bottom depth and free surface height, nondim.
  logical :: correct_thickness = .false.
  character(len=40)  :: mod = "initialize_thickness_from_file" ! This subroutine's name.
  character(len=200) :: filename, thickness_file, inputdir, mesg ! Strings for file/path
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call GOLD_mesg("  GOLD_initialization.F90, initialize_thickness_from_file: reading thickness", 5)

  inputdir = "." ; call read_param(param_file,"INPUTDIR",inputdir)
  inputdir = slasher(inputdir)
  call read_param(param_file,"THICKNESS_FILE",thickness_file,.true.)
  call read_param(param_file,"ADJUST_THICKNESS",correct_thickness)

  filename = trim(inputdir)//trim(thickness_file)
  if (.not.file_exists(filename, G%Domain)) call GOLD_error(FATAL, &
         " initialize_thickness_from_file: Unable to open "//trim(filename))

  if (file_has_thickness) then
    call read_data(filename,"h",h(:,:,:),domain=G%Domain%mpp_domain)
  else
    call read_data(filename,"eta",eta(:,:,:),domain=G%Domain%mpp_domain)

    if (correct_thickness) then
      ! All mass below the bottom removed if the topography is shallower than
      ! the input file would indicate.  G%D is positive downward,
      ! eta is negative downward.
      do j=js,je ; do i=is,ie
        if (-eta(i,j,nz+1) > G%D(i,j) + 0.1) eta(i,j,nz+1) = -G%D(i,j)
      enddo ; enddo
    endif

    do k=nz,1,-1 ; do j=js,je ; do i=is,ie
      if (eta(i,j,k) < (eta(i,j,k+1) + G%Angstrom_z)) then
        eta(i,j,k) = eta(i,j,k+1) + G%Angstrom_z
        h(i,j,k) = G%Angstrom_z
      else
        h(i,j,k) = eta(i,j,k) - eta(i,j,k+1)
      endif
    enddo ; enddo ; enddo

  !  Check for consistency between the interface heights and topography.!
    if (correct_thickness) then
      do j=js,je ; do i=is,ie
        !   The whole column is dilated to accomodate deeper topography than
        ! the input file would indicate.
        if (-eta(i,j,nz+1) < G%D(i,j) - 0.1) then
          dilate = (eta(i,j,1)+G%D(i,j)) / (eta(i,j,1)-eta(i,j,nz+1))
          do k=1,nz ; h(i,j,k) = h(i,j,k) * dilate ; enddo
        endif
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        if (abs(eta(i,j,nz+1) + G%D(i,j)) > 1.0) inconsistent = inconsistent + 1
      enddo ; enddo
      call sum_across_PEs(inconsistent)

      if ((inconsistent > 0) .and. (is_root_pe())) then
        write(mesg,'("Thickness initial conditions are inconsistent ",'// &
                 '"with topography in ",I5," places.")') inconsistent
        call GOLD_error(WARNING, mesg)
      endif
    endif

    call log_param(param_file, mod, "ADJUST_THICKNESS", correct_thickness, &
                 "If true, all mass below the bottom removed if the \n"//&
                 "topography is shallower than the thickness input file \n"//&
                 "would indicate.", default=.false.)
  endif
  call log_param(param_file, mod, "THICKNESS_FILE", thickness_file, &
                 "The name of the thickness file.")
  call log_param(param_file, mod, "INPUTDIR/THICKNESS_FILE", filename)
end subroutine initialize_thickness_from_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_thickness_uniform(h, G, param_file)
  real, intent(out), dimension(NXMEM_,NYMEM_, NZ_) :: h
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file

! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

!  This subroutine initializes the layer thicknesses to be uniform.
  real :: e0(SZK_(G))     ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  real :: max_depth ! The maximum depths in m.
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call GOLD_mesg("  GOLD_initialization.F90, initialize_thickness_uniform: setting thickness", 5)

  call read_param(param_file,"MAXIMUM_DEPTH",max_depth,.true.)
  do k=1,nz
    e0(k) = -max_depth * real(k-1) / real(nz)
  enddo

  do j=js,je ; do i=is,ie
!    This sets the initial thickness (in m) of the layers.  The      !
!  thicknesses are set to insure that: 1.  each layer is at least an !
!  Angstrom thick, and 2.  the interfaces are where they should be   !
!  based on the resting depths and interface height perturbations,   !
!  as long at this doesn't interfere with 1.                         !
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

  call log_param(param_file, "initialize_thickness_uniform", "MAXIMUM_DEPTH", &
                 max_depth, "The maximum depth of the ocean.", units="m")
end subroutine initialize_thickness_uniform
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_thickness_search
! search density space for location of layers
  call GOLD_error(FATAL,"  GOLD_initialization.F90, initialize_thickness_search: NOT IMPLEMENTED")
end subroutine initialize_thickness_search
! -----------------------------------------------------------------------------

subroutine convert_thickness(h, G, param_file, tv)
  real, intent(inout), dimension(NXMEM_,NYMEM_, NZ_) :: h
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  type(thermo_var_ptrs), intent(in) :: tv
! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    p_top, p_bot
  real :: dz_geo(SZ1_(h),SZ2_(h))      ! The change in geopotential height
                                       ! across a layer, in m2 s-2.
  real :: rho(SZI_(G))
  real :: I_gEarth
  logical :: Boussinesq, Bulkmixedlayer
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkml, nkbl
  integer :: itt, max_itt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  max_itt = 10
  Boussinesq = .true. ; call read_param(param_file,"BOUSSINESQ",Boussinesq)
  I_gEarth = 1.0 / G%g_Earth

  if (Boussinesq) then
    call GOLD_error(FATAL,"Not yet converting thickness with Boussinesq approx.")
  else
    if (associated(tv%eqn_of_state)) then
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        p_bot(i,j) = 0.0 ; p_top(i,j) = 0.0
      enddo ; enddo
      do k=1,nz
        do j=js,je
          do i=is,ie ; p_top(i,j) = p_bot(i,j) ; enddo
          call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_top(:,j), rho, &
                                 is, ie-is+1, tv%eqn_of_state)
          do i=is,ie
            p_bot(i,j) = p_top(i,j) + G%g_Earth * h(i,j,k) * rho(i)
          enddo
        enddo

        do itt=1,max_itt
          call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p_top, p_bot, &
                                   0.0, G, tv%eqn_of_state, dz_geo)
          if (itt < max_itt) then ; do j=js,je
            call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_bot(:,j), rho, &
                                   is, ie-is+1, tv%eqn_of_state)
            ! Use Newton's method to correct the bottom value.
            !   The hydrostatic equation is linear to such a
            ! high degree that no bounds-checking is needed.
            do i=is,ie
              p_bot(i,j) = p_bot(i,j) + rho(i) * (G%g_Earth*h(i,j,k) - dz_geo(i,j))
            enddo
          enddo ; endif
        enddo

        do j=js,je ; do i=is,ie
          h(i,j,k) = (p_bot(i,j) - p_top(i,j)) * G%kg_m2_to_H * I_gEarth
        enddo ; enddo
      enddo
    else
      bulkmixedlayer = .false.
      call read_param(param_file,"BULKMIXEDLAYER",bulkmixedlayer)
      if (bulkmixedlayer) then
        nkml = 1 ; call read_param(param_file,"NKML",nkml)
        nkbl = 1 ; call read_param(param_file,"NKBL",nkbl)
        do k=1,nkml+nkbl ; do j=js,je ; do i=is,ie
          h(i,j,k) = h(i,j,k) * tv%Rml(i,j,k) * G%kg_m2_to_H
        enddo ; enddo ; enddo
        do k=nkml+nkbl+1,nz ; do j=js,je ; do i=is,ie
          h(i,j,k) = h(i,j,k) * G%Rlay(k) * G%kg_m2_to_H
        enddo ; enddo ; enddo
      else
        do k=1,nz ; do j=js,je ; do i=is,ie
          h(i,j,k) = h(i,j,k) * G%Rlay(k) * G%kg_m2_to_H
        enddo ; enddo ; enddo
      endif
    endif
  endif

end subroutine convert_thickness

subroutine depress_surface(h, G, param_file, tv)
  real, dimension(NXMEM_,NYMEM_, NZ_), intent(inout) :: h
  type(ocean_grid_type),               intent(in)    :: G
  type(param_file_type),               intent(in)    :: param_file
  type(thermo_var_ptrs),               intent(in)    :: tv
! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

  real, dimension(SZI_(G),SZJ_(G)) :: &
    eta_sfc  ! The free surface height that the model should use, in m.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: &
    eta  ! The free surface height that the model should use, in m.
  real :: dilate  ! A ratio by which layers are dilated, nondim.
  real :: scale_factor ! A scaling factor for the eta_sfc values that are read
                       ! in, which can be used to change units, for example.
  character(len=200) :: inputdir, eta_srf_file ! Strings for file/path
  character(len=200) :: filename, eta_srf_var  ! Strings for file/path
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! Read the surface height (or pressure) from a file.

  inputdir = "." ; call read_param(param_file,"INPUTDIR",inputdir)
  inputdir = slasher(inputdir)
  call read_param(param_file,"SURFACE_HEIGHT_IC_FILE",eta_srf_file,.true.)
  call read_param(param_file,"SURFACE_HEIGHT_IC_VAR",eta_srf_var,.true.)
  filename = trim(inputdir)//trim(eta_srf_file)
  call read_data(filename,eta_srf_var,eta_sfc,domain=G%Domain%mpp_domain)
  scale_factor = 1.0
  call read_param(param_file,"SURFACE_HEIGHT_IC_SCALE", scale_factor)

  if (scale_factor /= 1.0) then ; do j=js,je ; do i=is,ie
    eta_sfc(i,j) = eta_sfc(i,j) * scale_factor
  enddo ; enddo ; endif

  ! Convert thicknesses to interface heights.
  call find_eta(h, tv, G%g_Earth, G, eta)
  
  do j=js,je ; do i=is,ie ; if (G%hmask(i,j) > 0.0) then
!    if (eta_sfc(i,j) < eta(i,j,nz+1)) then
      ! Issue a warning?
!    endif
    if (eta_sfc(i,j) > eta(i,j,1)) then
      ! Dilate the water column to agree, but only up to 10-fold.
      if (eta_sfc(i,j) - eta(i,j,nz+1) > 10.0*(eta(i,j,1) - eta(i,j,nz+1))) then
        dilate = 10.0
        call GOLD_error(WARNING, "Free surface height dilation attempted "//&
               "to exceed 10-fold.")
      else
        dilate = (eta_sfc(i,j) - eta(i,j,nz+1)) / (eta(i,j,1) - eta(i,j,nz+1))
      endif
      do k=1,nz ; h(i,j,k) = h(i,j,k) * dilate ; enddo
    elseif (eta(i,j,1) > eta_sfc(i,j)) then
      ! Remove any mass that is above the target free surface.
      do k=1,nz
        if (eta(i,j,k) <= eta_sfc(i,j)) exit
        if (eta(i,j,k+1) >= eta_sfc(i,j)) then
          h(i,j,k) = G%Angstrom
        else
          h(i,j,k) = max(G%Angstrom, h(i,j,k) * &
              (eta_sfc(i,j) - eta(i,j,k+1)) / (eta(i,j,k) - eta(i,j,k+1)) )
        endif
      enddo
    endif
  endif ; enddo ; enddo

  call log_param(param_file, "depress_surface", "SURFACE_HEIGHT_IC_FILE", &
         eta_srf_file,"The initial condition file for the surface height.")
  call log_param(param_file,"depress_surface", "SURFACE_HEIGHT_IC_VAR", &
         eta_srf_var, "The initial condition variable for the surface height.")
  call log_param(param_file, "depress_surface", "SURFACE_HEIGHT_IC_SCALE", &
    scale_factor, "A scaling factor to convert SURFACE_HEIGHT_IC_VAR into \n"//&
                 "units of m", units="variable", default=1.0)
end subroutine depress_surface

! -----------------------------------------------------------------------------
subroutine initialize_velocity_from_file(u, v, G, param_file)
  real, dimension(NXMEMQ_,NYMEM_, NZ_), intent(out) :: u
  real, dimension(NXMEM_,NYMEMQ_, NZ_), intent(out) :: v
  type(ocean_grid_type),                intent(in)  :: G
  type(param_file_type),                intent(in)  :: param_file
! Arguments: u - The zonal velocity that is being initialized.
!  (out)     v - The meridional velocity that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file -  parameter file type

!   This subroutine reads the initial velocity components from file
  character(len=200) :: filename,velocity_file,inputdir ! Strings for file/path

  call GOLD_mesg("  GOLD_initialization.F90, initialize_velocity_from_file: reading u and v", 5)

  inputdir = "." ; call read_param(param_file, "INPUTDIR", inputdir)
  inputdir = slasher(inputdir)
  call read_param(param_file,"VELOCITY_FILE", velocity_file, .true.)

!  Read the velocities from a netcdf file.                           !
  filename = trim(inputdir)//trim(velocity_file)
  if (.not.file_exists(filename, G%Domain)) call GOLD_error(FATAL, &
         " initialize_velocity_from_file: Unable to open "//trim(filename))

  call read_data(filename,"u",u(:,:,:),domain=G%Domain%mpp_domain)
  call read_data(filename,"v",v(:,:,:),domain=G%Domain%mpp_domain)

  call log_param(param_file, "initialize_velocity_from_file", &
                 "VELOCITY_FILE", velocity_file, &
                 "The name of the velocity initial condition file.")
  call log_param(param_file, "initialize_velocity_from_file", &
                 "INPUTDIR/VELOCITY_FILE", filename)
end subroutine initialize_velocity_from_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_velocity_zero(u, v, G, param_file)
  real, dimension(NXMEMQ_,NYMEM_, NZ_), intent(out) :: u
  real, dimension(NXMEM_,NYMEMQ_, NZ_), intent(out) :: v
  type(ocean_grid_type),                intent(in)  :: G
  type(param_file_type),                intent(in)  :: param_file
! Arguments: u - The zonal velocity that is being initialized.
!  (out)     v - The meridional velocity that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file -  parameter file type

!   This subroutine sets the initial velocity components to zero
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq


  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u(I,j,k) = 0.0
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v(i,J,k) = 0.0
  enddo ; enddo ; enddo

end subroutine initialize_velocity_zero
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_velocity_uniform(u, v, G, param_file)
  real, dimension(NXMEMQ_,NYMEM_, NZ_), intent(out) :: u
  real, dimension(NXMEM_,NYMEMQ_, NZ_), intent(out) :: v
  type(ocean_grid_type),                intent(in)  :: G
  type(param_file_type),                intent(in)  :: param_file
! Arguments: u - The zonal velocity that is being initialized.
!  (out)     v - The meridional velocity that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file -  parameter file type

!   This subroutine sets the initial velocity components to uniform
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  real    :: initial_u_const, initial_v_const
  character(len=200) :: mod ! Strings
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq

  call read_param(param_file,"INITIAL_U_CONST",initial_u_const,.true.)
  call read_param(param_file,"INITIAL_V_CONST",initial_v_const,.true.)

  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u(I,j,k) = initial_u_const
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v(i,J,k) = initial_v_const
  enddo ; enddo ; enddo

  mod = "initialize_velocity_uniform"
  call log_param(param_file, mod, "INITIAL_U_CONST", initial_u_const, &
                 "A initial uniform value for the zonal flow.", units="m/s")
  call log_param(param_file, mod, "INITIAL_V_CONST", initial_v_const, &
                 "A initial uniform value for the meridional flow.", units="m/s")
end subroutine initialize_velocity_uniform
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_velocity_circular(u, v, G, param_file)
  real, dimension(NXMEMQ_,NYMEM_, NZ_), intent(out) :: u
  real, dimension(NXMEM_,NYMEMQ_, NZ_), intent(out) :: v
  type(ocean_grid_type),                intent(in)  :: G
  type(param_file_type),                intent(in)  :: param_file
! Arguments: u - The zonal velocity that is being initialized.
!  (out)     v - The meridional velocity that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file -  parameter file type

!   This subroutine sets the initial velocity components to be circular with
! no flow at edges of domain and center.
  character(len=200) :: mod ! Strings
  real :: south_lat, len_lat, west_lon, len_lon, circular_max_u
  real :: dpi, psi1, psi2
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq

  call read_param(param_file,"SOUTHLAT",south_lat)
  call read_param(param_file,"LENLAT",len_lat)
  call read_param(param_file,"WESTLON",west_lon)
  call read_param(param_file,"LENLON",len_lon)
  circular_max_u=0.; call read_param(param_file,"CIRCULAR_MAX_U",circular_max_u)
  dpi=acos(0.0)*2.0 ! pi

  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    psi1 = my_psi(I,j)
    psi2 = my_psi(I,j-1)
    u(I,j,k) = (psi1-psi2)/G%dy_u(I,j)! *(circular_max_u*len_lon/(2.0*dpi))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    psi1 = my_psi(i,J)
    psi2 = my_psi(i-1,J)
    v(i,J,k) = (psi2-psi1)/G%dx_v(i,J)! *(circular_max_u*len_lon/(2.0*dpi))
  enddo ; enddo ; enddo

  mod = "initialize_velocity_circular"
  call log_param(param_file, mod, "CIRCULAR_MAX_U", circular_max_u, &
                 "The amplitude of zonal flow from which to scale the\n"// &
                 "circular stream function (m/s).", &
                 units="m/s", default=0.)

  contains

  real function my_psi(ig,jg) ! in-line function
    integer :: ig, jg
    real :: x, y, r
    x = 2.0*(G%geolonq(ig,jg)-west_lon)/len_lon-1.0  ! -1<x<1
    y = 2.0*(G%geolatq(ig,jg)-south_lat)/len_lat-1.0 ! -1<y<1
    r = sqrt( x**2 + y**2 ) ! Circulat stream fn nis fn of radius only
    r = min(1.0,r) ! Flatten stream function in corners of box
    my_psi = 0.5*(1.0 - cos(dpi*r))
    my_psi = my_psi * (circular_max_u*len_lon*1e3/dpi) ! len_lon is in km
  end function my_psi

end subroutine initialize_velocity_circular
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_temp_salt_from_file(T, S, G, param_file)
  real, dimension(NXMEM_,NYMEM_, NZ_), intent(out) :: T, S
  type(ocean_grid_type),               intent(in)  :: G
  type(param_file_type),               intent(in)  :: param_file
!  This function puts the initial layer temperatures and salinities  !
! into T(:,:,:) and S(:,:,:).                                        !

! Arguments: T - The potential temperature that is being initialized.
!  (out)     S - The salinity that is being initialized.
!  (in)      from_file - .true. if the variables that are set here are to
!                        be read from a file; .false. to be set internally.
!  (in)      filename - The name of the file to read.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
  character(len=200) :: filename, ts_file, salt_file, inputdir, mod ! Strings for file/path
  character(len=64)  :: temp_var, salt_var

  call GOLD_mesg("  GOLD_initialization.F90, initialize_temp_salt_from_file: reading T and S", 5)

  inputdir = "." ; call read_param(param_file,"INPUTDIR",inputdir)
  inputdir = slasher(inputdir)
  call read_param(param_file,"TS_FILE",ts_file,.true.)

! Read the temperatures and salinities from a netcdf file.           !
  filename = trim(inputdir)//trim(ts_file)
  if (.not.file_exists(filename, G%Domain)) call GOLD_error(FATAL, &
     " initialize_temp_salt_from_file: Unable to open "//trim(filename))
  temp_var = "PTEMP" ; salt_var = "SALT"
  call read_param(param_file,"TEMP_IC_VAR",temp_var)
  call read_param(param_file,"SALT_IC_VAR",salt_var)

  call read_data(filename,temp_var,T(:,:,:),domain=G%Domain%mpp_domain)

  salt_file = ts_file
  call read_param(param_file,"SALT_FILE",salt_file)
  filename = trim(inputdir)//trim(ts_file)
  if (.not.file_exists(filename, G%Domain)) call GOLD_error(FATAL, &
     " initialize_temp_salt_from_file: Unable to open "//trim(filename))
  call read_data(filename,salt_var,S(:,:,:),domain=G%Domain%mpp_domain)

  mod = "initialize_temp_salt_from_file"
  call log_param(param_file, mod, "TS_FILE", ts_file, &
                 "The initial condition file for temperature.")
  call log_param(param_file, mod, "SALT_FILE", salt_file, &
                 "The initial condition file for salinity.", default=trim(ts_file))
  call log_param(param_file, mod, "TEMP_IC_VAR", temp_var, &
                 "The initial condition variable for potential temperature.", &
                 default="PTEMP")
  call log_param(param_file, mod, "SALT_IC_VAR", temp_var, &
                 "The initial condition variable for salinity.", default="SALT")
  call log_param(param_file, mod, "INPUTDIR/TS_FILE", filename)
end subroutine initialize_temp_salt_from_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_temp_salt_from_profile(T, S, G, param_file)
  real, dimension(NXMEM_,NYMEM_, NZ_), intent(out) :: T, S
  type(ocean_grid_type),               intent(in)  :: G
  type(param_file_type),               intent(in)  :: param_file
!  This function puts the initial layer temperatures and salinities  !
! into T(:,:,:) and S(:,:,:).                                        !

! Arguments: T - The potential temperature that is being initialized.
!  (out)     S - The salinity that is being initialized.
!  (in)      from_file - .true. if the variables that are set here are to
!                        be read from a file; .false. to be set internally.
!  (in)      filename - The name of the file to read.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
  real, dimension(SZK_(G)) :: T0, S0
  integer :: i, j, k
  character(len=200) :: filename, ts_file, inputdir ! Strings for file/path

  call GOLD_mesg("  GOLD_initialization.F90, initialize_temp_salt_from_file: reading T and S", 5)

  inputdir = "." ; call read_param(param_file,"INPUTDIR",inputdir)
  inputdir = slasher(inputdir)
  call read_param(param_file,"TS_FILE",ts_file,.true.)

! Read the temperatures and salinities from a netcdf file.           !
  filename = trim(inputdir)//trim(ts_file)
  if (.not.file_exists(filename)) call GOLD_error(FATAL, &
     " initialize_temp_salt_from_profile: Unable to open "//trim(filename))

  call read_data(filename,"PTEMP",T0(:),domain=G%Domain%mpp_domain)
  call read_data(filename,"SALT", S0(:),domain=G%Domain%mpp_domain)

  do k=1,G%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    T(i,j,k) = T0(k) ; S(i,j,k) = S0(k)
  enddo ; enddo ; enddo

  call log_param(param_file, "initialize_temp_salt_from_profile", &
                 "TS_FILE", ts_file, &
                 "The initial condition file for temperature.")
  call log_param(param_file, "initialize_temp_salt_from_profile", &
                 "INPUTDIR/TS_FILE", filename)
end subroutine initialize_temp_salt_from_profile
! -----------------------------------------------------------------------------


! -----------------------------------------------------------------------------
subroutine initialize_temp_salt_fit(T, S, G, param_file, eqn_of_state, P_Ref)
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
  real :: T_Ref         ! Reference Temperature
  real :: S_Ref         ! Reference Salinity
  real :: pres(SZK_(G))      ! An array of the reference pressure in Pa.
  real :: drho_dT(SZK_(G))   ! Derivative of density with temperature in kg m-3 K-1.                              !
  real :: drho_dS(SZK_(G))   ! Derivative of density with salinity in kg m-3 PSU-1.                             !
  real :: rho_guess(SZK_(G)) ! Potential density at T0 & S0 in kg m-3.
  character(len=40)  :: mod = "initialize_temp_salt_fit" ! This subroutine's name.
  integer :: i, j, k, itt, nz
  nz = G%ke

  call GOLD_mesg("  GOLD_initialization.F90, initialize_temp_salt_fit: setting T and S", 5)

  call read_param(param_file,"T_REF",T_Ref,.true.)
  call read_param(param_file,"S_REF",S_Ref,.true.)
  do k=1,nz
    pres(k) = P_Ref ; S0(k) = S_Ref
  enddo
  T0(1) = T_Ref

  call calculate_density(T0(1),S0(1),pres(1),rho_guess(1),1,1,eqn_of_state)
  call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,1,eqn_of_state)

! A first guess of the layers' temperatures.                         !
  do k=nz,1,-1
    T0(k) = T0(1) + (G%Rlay(k) - rho_guess(1)) / drho_dT(1)
  enddo

! Refine the guesses for each layer.                                 !
  do itt=1,6
    call calculate_density(T0,S0,pres,rho_guess,1,nz,eqn_of_state)
    call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,nz,eqn_of_state)
    do k=1,nz
      T0(k) = T0(k) + (G%Rlay(k) - rho_guess(k)) / drho_dT(k)
    enddo
  enddo

  do k=1,nz ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    T(i,j,k) = T0(k) ; S(i,j,k) = S0(k)
  enddo ; enddo ; enddo

  call log_param(param_file, mod, "T_REF", T_Ref, &
                 "A reference temperature used in initialization.", units="degC")
  call log_param(param_file, mod, "S_REF", S_Ref, &
                 "A reference salinity used in initialization.", units="PSU")
end subroutine initialize_temp_salt_fit
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_mixed_layer_density(Rml, G, param_file, use_temperature, &
                                          eqn_of_state, T, S, P_Ref)
  real, dimension(NXMEM_,NYMEM_, NZ_),          intent(out) :: Rml
  type(ocean_grid_type),                        intent(in)  :: G
  type(param_file_type),                        intent(in)  :: param_file
  logical,                                      intent(in)  :: use_temperature
  type(EOS_type),                      optional, pointer    :: eqn_of_state
  real, dimension(NXMEM_,NYMEM_, NZ_), optional, intent(in) :: T, S
  real,                                optional, intent(in) :: P_Ref
! Set the initial mixed layer and buffer layer densities.            !

  real :: pres(SZ1_(Rml)) ! An array of the reference pressure in Pa.
  character(len=40)  :: mod = "initialize_mixed_layer_density" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, ied, nkml, nkbl
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; ied = G%ied

  nkml = 1 ; call read_param(param_file,"NKML",nkml)
  nkbl = 1 ; call read_param(param_file,"NKBL",nkbl)
  if (use_temperature) then
    pres(:) = P_Ref
    if (.not.present(T) .or. .not.present(S) .or. .not.present(eqn_of_state)) &
      call GOLD_error(FATAL, " initialize_mixed_layer_density: "// &
          "T, S, and eqn_of_state must be present with use_temperature true.")

    do k=1,nkml+nkbl ; do j=js-1,je
      call calculate_density(T(:,j,k),S(:,j,k),pres,Rml(:,j,k),is-1,ie-is+2,eqn_of_state)
    enddo ; enddo
  else
    ! Specify some analytic expression for mixed layer density.
    ! The following is just an example.
    do k=1,nkml+nkbl ; do j=js-1,je ; do i=is-1,ie
      Rml(i,j,k) = G%Rlay(nkml+nkbl)
    enddo ; enddo ; enddo
  endif
  do k=nkml+nkbl+1,G%ke ; do j=js-1,je ; do i=is-1,ie
    Rml(i,j,k) = G%Rlay(k)
  enddo ; enddo ; enddo

  call log_param(param_file, mod, "NKML", nkml, &
                 "The number of sublayers within the mixed layer if \n"//&
                 "BULKMIXEDLAYER is true.", units="nondim", default=1)
  call log_param(param_file, mod, "NKBL", nkbl, &
                 "The number of layers that are used as variable density \n"//&
                 "buffer layers if BULKMIXEDLAYER is true.", units="nondim", &
                 default=1)
end subroutine initialize_mixed_layer_density
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_sponges_file(G, use_temperature, bulkmixedlayer, tv, &
                                   param_file, CSp)
  type(ocean_grid_type), intent(in) :: G
  logical, intent(in) :: use_temperature, bulkmixedlayer
  type(thermo_var_ptrs), intent(in) :: tv
  type(param_file_type), intent(in) :: param_file
  type(sponge_CS),       pointer    :: CSp
!   This subroutine sets the inverse restoration time (Idamp), and   !
! the values towards which the interface heights and an arbitrary    !
! number of tracers should be restored within each sponge. The       !
! interface height is always subject to damping, and must always be  !
! the first registered field.                                        !

! Arguments: from_file - .true. if the variables that are used here are to
!                        be read from a file; .false. to be set internally.
!  (in)      filename - The name of the file to read for all fields
!                       except the inverse damping rate.
!  (in)      damp_file - The name of the file from which to read the
!                        inverse damping rate.
!  (in)      G - The ocean's grid structure.
!  (in)      use_temperature - If true, T & S are state variables.
!  (in)      bulkmixedlayer - If true, a bulk mixed layer is used.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CSp - A pointer that is set to point to the control structure
!                  for this module

  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! The target interface heights, in m.
  real, dimension (SZI_(G),SZJ_(G),SZK_(G)) :: &
    tmp, tmp2 ! A temporary array for tracers.
  real, dimension (SZI_(G),SZJ_(G)) :: &
    tmp_2d ! A temporary array for tracers.
  real :: Idamp(SZI_(G),SZJ_(G))    ! The inverse damping rate, in s-1.
  real :: pres(SZI_(G))     ! An array of the reference pressure, in Pa.

  integer :: i, j, k, is, ie, js, je, nz
  character(len=40) :: potemp_var, salin_var, Idamp_var, eta_var
  character(len=200) :: damping_file, state_file  ! Strings for filenames
  character(len=200) :: filename, inputdir ! Strings for file/path and path.
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  pres(:) = 0.0 ; eta(:,:,:) = 0.0 ; tmp(:,:,:) = 0.0 ; Idamp(:,:) = 0.0

  inputdir = "." ; call read_param(param_file, "INPUTDIR", inputdir)
  inputdir = slasher(inputdir)
  call read_param(param_file,"SPONGE_DAMPING_FILE", damping_file, .true.)
  state_file = damping_file
  call read_param(param_file,"SPONGE_STATE_FILE", state_file)

  potemp_var = "PTEMP" ; salin_var="SALT" ; eta_var = "ETA" ; Idamp_var = "IDAMP"
  call read_param(param_file, "SPONGE_PTEMP_VAR", potemp_var)
  call read_param(param_file, "SPONGE_SALT_VAR", salin_var)
  call read_param(param_file, "SPONGE_ETA_VAR", eta_var)
  call read_param(param_file, "SPONGE_IDAMP_VAR", Idamp_var)

  filename = trim(inputdir)//trim(damping_file)
  if (.not.file_exists(filename, G%Domain)) &
    call GOLD_error(FATAL, " initialize_sponges: Unable to open "//trim(filename))

  call read_data(filename,"Idamp",Idamp(:,:), domain=G%Domain%mpp_domain)

! Now register all of the fields which are damped in the sponge.     !
! By default, momentum is advected vertically within the sponge, but !
! momentum is typically not damped within the sponge.                !

  filename = trim(inputdir)//trim(state_file)
  if (.not.file_exists(filename, G%Domain)) &
    call GOLD_error(FATAL, " initialize_sponges: Unable to open "//trim(filename))

!  The first call to set_up_sponge_field is for the interface height.!
  call read_data(filename, eta_var, eta(:,:,:), domain=G%Domain%mpp_domain)

  do j=js,je ; do i=is,ie
    eta(i,j,nz+1) = -G%D(i,j)
  enddo ; enddo
  do k=nz,1,-1 ; do j=js,je ; do i=is,ie
    if (eta(i,j,k) < (eta(i,j,k+1) + G%Angstrom_z)) &
      eta(i,j,k) = eta(i,j,k+1) + G%Angstrom_z
  enddo ; enddo ; enddo
! Set the inverse damping rates so that the model will know where to !
! apply the sponges, along with the interface heights.               !
  call initialize_sponge(Idamp, eta, G, param_file, CSp)

!   Now register all of the tracer fields which are damped in the    !
! sponge. By default, momentum is advected vertically within the     !
! sponge, but momentum is typically not damped within the sponge.    !

  if ( bulkmixedlayer ) then
!   The second call to set_up_sponge_field must be for Rml if        !
! BULKMIXEDLAYER is defined. The remaining calls can be in any order.!
! Only a single layer's worth of the Rml reference is used, even if  !
! there are multiple parts of the mixed layer (i.e. nkml>1).         !
    if ( use_temperature ) then
      do i=is-1,ie ; pres(i) = tv%P_Ref ; enddo

      call read_data(filename, potemp_var, tmp(:,:,:), domain=G%Domain%mpp_domain)
      call read_data(filename, salin_var, tmp2(:,:,:), domain=G%Domain%mpp_domain)

      do j=js,je
        call calculate_density(tmp(:,j,1),tmp2(:,j,1),pres,tmp_2d(:,j),is,ie-is+1,tv%eqn_of_state)
        do i=is,ie ; tmp(i,j,1) = tmp_2d(i,j) ; enddo
      enddo

    else
      call read_data(filename,"Rml",tmp(:,:,1), domain=G%Domain%mpp_domain)
    endif

    call set_up_sponge_field(tmp,tv%Rml,1,CSp)
  endif

!  The remaining calls to set_up_sponge_field can be in any order.   !
  if ( use_temperature ) then
    call read_data(filename, potemp_var, tmp(:,:,:), domain=G%Domain%mpp_domain)
    call set_up_sponge_field(tmp, tv%T, nz, CSp)
    call read_data(filename, salin_var, tmp(:,:,:), domain=G%Domain%mpp_domain)
    call set_up_sponge_field(tmp, tv%S, nz, CSp)
  endif

  call log_param(param_file, "initialize_sponges_file", &
                 "SPONGE_DAMPING_FILE", damping_file, &
                 "The name of the file with the sponge damping rates.")
  call log_param(param_file, "initialize_sponges_file", &
                 "SPONGE_STATE_FILE", state_file, &
                 "The name of the file with the state to damp toward.", &
                 default=damping_file)
  call log_param(param_file, "initialize_sponges_file", &
                 "SPONGE_PTEMP_VAR", potemp_var, &
                 "The name of the potential temperature variable in \n"//&
                 "SPONGE_STATE_FILE.", default="PTEMP")
  call log_param(param_file, "initialize_sponges_file", &
                 "SPONGE_SALT_VAR", salin_var, &
                 "The name of the salinity variable in \n"//&
                 "SPONGE_STATE_FILE.", default="SALT")
  call log_param(param_file, "initialize_sponges_file", &
                 "SPONGE_ETA_VAR", eta_var, &
                 "The name of the interface height variable in \n"//&
                 "SPONGE_STATE_FILE.", default="ETA")
  call log_param(param_file, "initialize_sponges_file", &
                 "SPONGE_IDAMP_VAR", Idamp_var, &
                 "The name of the inverse damping rate variable in \n"//&
                 "SPONGE_DAMPING_FILE.", default="IDAMP")

end subroutine initialize_sponges_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_Open_Bdry_Conds(OBC, tv, G, param_file, advect_tracer_CSp)
  type(ocean_OBC_type),  pointer    :: OBC
  type(thermo_var_ptrs), intent(in) :: tv
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  type(advect_tracer_CS), pointer   :: advect_tracer_CSp
!   This subroutine sets the properties of flow at open boundary conditions.
! This particular example is for the DOME inflow describe in Legg et al. 2006.

! Arguments: OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (out)     tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

  logical :: any_OBC        ! Set to true if any points in this subdomain use
                            ! open boundary conditions.
  logical, pointer, dimension(:,:) :: &
    OBC_mask_u => NULL(), & ! These arrays are true at zonal or meridional
    OBC_mask_v => NULL()    ! velocity points that have prescribed open boundary
                            ! conditions.
  real, pointer, dimension(:,:,:) :: &
    OBC_T_u => NULL(), &    ! These arrays should be allocated and set to
    OBC_T_v => NULL(), &    ! specify the values of T and S that should come
    OBC_S_u => NULL(), &    ! in through u- and v- points through the open
    OBC_S_v => NULL()       ! boundary conditions, in C and psu.
  logical :: apply_OBC_u = .false., apply_OBC_v = .false.
  ! The following variables are used to set the target temperature and salinity.
  real :: T0(SZK_(G)), S0(SZK_(G))
  real :: pres(SZK_(G))      ! An array of the reference pressure in Pa.
  real :: drho_dT(SZK_(G))   ! Derivative of density with temperature in kg m-3 K-1.                              !
  real :: drho_dS(SZK_(G))   ! Derivative of density with salinity in kg m-3 PSU-1.                             !
  real :: rho_guess(SZK_(G)) ! Potential density at T0 & S0 in kg m-3.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz, yhalo
  integer :: Isdq, Iedq, Jsdq, Jedq

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq
  yhalo = G%jsc-G%jsd

  call read_param(param_file,"APPLY_OBC_U",apply_OBC_u)
  call read_param(param_file,"APPLY_OBC_V",apply_OBC_v)

  if (apply_OBC_u) then
    ! Determine where u points are applied.
    allocate(OBC_mask_u(Isdq:Iedq,jsd:jed)) ; OBC_mask_u(:,:) = .false.
    any_OBC = .false.
    do j=jsd,jed ; do I=Isdq,Iedq
    ! if (SOME_TEST_FOR_U_OPEN_BCS) then
    !   OBC_mask_u(I,j) = .true. ; any_OBC = .true.
    ! endif
    enddo ; enddo
    if (.not.any_OBC) then
      ! This processor has no u points at which open boundary conditions are
      ! to be applied.
      apply_OBC_u = .false.
      deallocate(OBC_mask_u)
    endif
  endif
  if (apply_OBC_v) then
    ! Determine where v points are applied.
    allocate(OBC_mask_v(isd:ied,Jsdq:Jedq)) ; OBC_mask_v(:,:) = .false.
    any_OBC = .false.
    do J=Jsdq,Jedq ; do i=isd,ied
    ! if (SOME_TEST_FOR_V_OPEN_BCS) then
    !   OBC_mask_v(i,J) = .true. ; any_OBC = .true.
    ! endif
    enddo ; enddo
    if (.not.any_OBC) then
      ! This processor has no v points at which open boundary conditions are
      ! to be applied.
      apply_OBC_v = .false.
      deallocate(OBC_mask_v)
    endif
  endif

  if (.not.(apply_OBC_u .or. apply_OBC_v)) return

  if (.not.associated(OBC)) allocate(OBC)

  if (apply_OBC_u) then
    OBC%apply_OBC_u = .true.
    OBC%OBC_mask_u => OBC_mask_u
    allocate(OBC%u(Isdq:Iedq,jsd:jed,nz)) ; OBC%u(:,:,:) = 0.0
    allocate(OBC%uh(Isdq:Iedq,jsd:jed,nz)) ; OBC%uh(:,:,:) = 0.0
    allocate(OBC%OBC_kind_u(Isdq:Iedq,jsd:jed)) ; OBC%OBC_kind_u(:,:) = OBC_NONE
    do j=jsd,jed ; do I=Isdq,Iedq
      if (OBC%OBC_mask_u(I,j)) OBC%OBC_kind_u(I,j) = OBC_SIMPLE
    enddo ; enddo
  endif
  if (apply_OBC_v) then
    OBC%apply_OBC_v = .true.
    OBC%OBC_mask_v => OBC_mask_v
    allocate(OBC%v(isd:ied,Jsdq:Jedq,nz)) ; OBC%v(:,:,:) = 0.0
    allocate(OBC%vh(isd:ied,Jsdq:Jedq,nz)) ; OBC%vh(:,:,:) = 0.0
    allocate(OBC%OBC_kind_v(isd:ied,Jsdq:Jedq)) ; OBC%OBC_kind_v(:,:) = OBC_NONE
    do J=Jsdq,Jedq ; do i=isd,ied
      if (OBC%OBC_mask_v(i,J)) OBC%OBC_kind_v(i,J) = OBC_SIMPLE
    enddo ; enddo
  endif

  if (apply_OBC_v) then
    do k=1,nz ; do J=Jsd,Jed ; do i=isd,ied
      if (OBC_mask_v(i,J)) then
        ! An appropriate expression for the meridional inflow velocities and
        ! transports should go here.
        OBC%vh(i,J,k) = 0.0 * G%m_to_H ; OBC%v(i,J,k) = 0.0
      else
        OBC%vh(i,J,k) = 0.0 ; OBC%v(i,J,k) = 0.0
      endif
    enddo ; enddo ; enddo
  endif

  if (apply_OBC_u) then
    do k=1,nz ; do j=jsd,jed ; do I=Isdq,Iedq
      if (OBC_mask_u(I,j)) then
        ! An appropriate expression for the zonal inflow velocities and
        ! transports should go here.
        OBC%uh(I,j,k) = 0.0 * G%m_to_H ; OBC%u(I,j,k) = 0.0
      else
        OBC%uh(I,j,k) = 0.0 ; OBC%u(I,j,k) = 0.0
      endif
    enddo ; enddo ; enddo
  endif

  !   The inflow values of temperature and salinity also need to be set here if
  ! these variables are used.  The following code is just a naive example.
  if (apply_OBC_u .or. apply_OBC_v) then
    if (associated(tv%S)) then
      ! In this example, all S inflows have values of 35 psu.
      call add_tracer_OBC_values("S", advect_tracer_CSp, OBC_inflow=35.0)
    endif
    if (associated(tv%T)) then
      ! In this example, the T values are set to be consistent with the layer
      ! target density and a salinity of 35 psu.  This code is taken from
      !  initialize_temp_sal.
      pres(:) = tv%P_Ref ; S0(:) = 35.0 ; T0(1) = 25.0
      call calculate_density(T0(1),S0(1),pres(1),rho_guess(1),1,1,tv%eqn_of_state)
      call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,1,tv%eqn_of_state)

      do k=1,nz ; T0(k) = T0(1) + (G%Rlay(k)-rho_guess(1)) / drho_dT(1) ; enddo
      do itt=1,6
        call calculate_density(T0,S0,pres,rho_guess,1,nz,tv%eqn_of_state)
        call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,nz,tv%eqn_of_state)
        do k=1,nz ; T0(k) = T0(k) + (G%Rlay(k)-rho_guess(k)) / drho_dT(k) ; enddo
      enddo

      if (apply_OBC_u) then
        allocate(OBC_T_u(Isdq:Iedq,jsd:jed,nz))
        do k=1,nz ; do j=jsd,jed ; do I=Isdq,Iedq
          OBC_T_u(I,j,k) = T0(k)
        enddo ; enddo ; enddo
      endif
      if (apply_OBC_v) then
        allocate(OBC_T_v(isd:ied,Jsdq:Jedq,nz))
        do k=1,nz ; do J=Jsdq,Jedq ; do i=isd,ied
          OBC_T_v(i,J,k) = T0(k)
        enddo ; enddo ; enddo
      endif
      call add_tracer_OBC_values("T", advect_tracer_CSp, OBC_in_u=OBC_T_u, &
                                                         OBC_in_v=OBC_T_v)
    endif
  endif

end subroutine set_Open_Bdry_Conds
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_Flather_Bdry_Conds(OBC, tv, G, param_file, advect_tracer_CSp)
  type(ocean_OBC_type),  pointer    :: OBC
  type(thermo_var_ptrs), intent(in) :: tv
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  type(advect_tracer_CS), pointer   :: advect_tracer_CSp
!   This subroutine sets the initial definitions of the characteristic open boundary
!   conditions. Written by Mehmet Ilicak

! Arguments: OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (out)     tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

  logical :: any_OBC        ! Set to true if any points in this subdomain use
                            ! open boundary conditions.

!!!!!!!!!!!!!!!!!!!!!! Mehmet !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical :: apply_OBC_u_flather_east = .false., apply_OBC_u_flather_west = .false.
  logical :: apply_OBC_v_flather_north = .false., apply_OBC_v_flather_south = .false.  
  integer :: isd_global, jsd_global
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz, yhalo
  integer :: Isdq, Iedq, Jsdq, Jedq
  integer :: east_boundary, west_boundary, north_boundary, south_boundary

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq
  yhalo = G%jsc-G%jsd 
  
!!!!!!!!!!!!!!!!! Mehmet !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  isd_global = G%isd_global
  jsd_global = G%jsd_global
  call read_param(param_file,"APPLY_OBC_U_FLATHER_EAST",apply_OBC_u_flather_east)
  call read_param(param_file,"APPLY_OBC_U_FLATHER_WEST",apply_OBC_u_flather_west)
  call read_param(param_file,"APPLY_OBC_V_FLATHER_NORTH",apply_OBC_v_flather_north)
  call read_param(param_file,"APPLY_OBC_V_FLATHER_SOUTH",apply_OBC_v_flather_south)

  if (.not.(apply_OBC_u_flather_east  .or. apply_OBC_u_flather_west .or. &            
            apply_OBC_v_flather_north .or. apply_OBC_v_flather_south)) return
  
  if (.not.associated(OBC)) allocate(OBC)
       
  OBC%apply_OBC_u_flather_east = apply_OBC_u_flather_east
  OBC%apply_OBC_u_flather_west = apply_OBC_u_flather_west 
  OBC%apply_OBC_v_flather_north = apply_OBC_v_flather_north 
  OBC%apply_OBC_v_flather_south = apply_OBC_v_flather_south

  if (G%symmetric) then
    east_boundary = G%domain%nxtot+G%domain%nx_halo
    west_boundary = G%domain%nx_halo
    north_boundary = G%domain%nytot+G%domain%ny_halo
    south_boundary = G%domain%ny_halo
  else
    east_boundary = G%domain%nxtot+G%domain%nx_halo-1
    west_boundary = G%domain%nx_halo+1
    north_boundary = G%domain%nytot+G%domain%ny_halo-1
    south_boundary = G%domain%ny_halo+1
  endif

  if (.not.associated(OBC%OBC_mask_u)) then
    allocate(OBC%OBC_mask_u(Isdq:Iedq,jsd:jed)) ; OBC%OBC_mask_u(:,:) = .false.
  endif
  if (.not.associated(OBC%OBC_kind_u)) then
    allocate(OBC%OBC_kind_u(Isdq:Iedq,jsd:jed)) ; OBC%OBC_kind_u(:,:) = OBC_NONE
  endif
  if (.not.associated(OBC%OBC_mask_v)) then
    allocate(OBC%OBC_mask_v(isd:ied,Jsdq:Jedq)) ; OBC%OBC_mask_v(:,:) = .false.
  endif
  if (.not.associated(OBC%OBC_kind_v)) then
    allocate(OBC%OBC_kind_v(isd:ied,Jsdq:Jedq)) ; OBC%OBC_kind_v(:,:) = OBC_NONE
  endif

  if (apply_OBC_u_flather_east) then
    ! Determine where u points are applied at east side 
    do j=jsd,jed ; do I=Isdq,Iedq
      if ((I+isd_global-isd) .eq. east_boundary) then !eastern side
        OBC%OBC_mask_u(I,j) = .true.
        OBC%OBC_kind_u(I,j) = OBC_FLATHER_E
        if ((i+1>isd) .and. (i+1<ied) .and. (J>Jsdq) .and. (J<Jedq)) then
          OBC%OBC_mask_v(i+1,J) = .true.
          if (OBC%OBC_kind_v(i+1,J) == OBC_NONE) OBC%OBC_kind_v(i+1,J) = OBC_FLATHER_E
        endif
        if ((i+1>isd) .and. (i+1<ied) .and. (J-1>Jsdq) .and. (J-1<Jedq)) then
          OBC%OBC_mask_v(i+1,J-1) = .true.
          if (OBC%OBC_kind_v(i+1,J-1) == OBC_NONE) OBC%OBC_kind_v(i+1,J-1) = OBC_FLATHER_E
        endif
      endif
    enddo ; enddo
  endif

  if (apply_OBC_u_flather_west) then
    ! Determine where u points are applied at west side 
    do j=jsd,jed ; do I=Isdq,Iedq
      if ((I+isd_global-isd) .eq. west_boundary) then !western side
        OBC%OBC_mask_u(I,j) = .true.
        OBC%OBC_kind_u(I,j) = OBC_FLATHER_W
        if ((i>isd) .and. (i<ied) .and. (J>Jsdq) .and. (J<Jedq)) then
          OBC%OBC_mask_v(i,J) = .true.
          if (OBC%OBC_kind_v(i,J) == OBC_NONE) OBC%OBC_kind_v(i,J) = OBC_FLATHER_W
        endif
        if ((i>isd) .and. (i<ied) .and. (J-1>Jsdq) .and. (J-1<Jedq)) then
          OBC%OBC_mask_v(i,J-1) = .true.
          if (OBC%OBC_kind_v(i,J-1) == OBC_NONE) OBC%OBC_kind_v(i,J-1) = OBC_FLATHER_W
        endif
      endif
    enddo ; enddo
  endif

  if (apply_OBC_v_flather_north) then
    ! Determine where v points are applied at north side 
    do J=Jsdq,Jedq ; do i=isd,ied
      if ((J+jsd_global-jsd) .eq. north_boundary) then         !northern side
        OBC%OBC_mask_v(i,J) = .true.
        OBC%OBC_kind_v(i,J) = OBC_FLATHER_N
        if ((I>Isdq) .and. (I<Iedq) .and. (j+1>jsd) .and. (j+1<jed)) then
          OBC%OBC_mask_u(I,j+1) = .true.
          if (OBC%OBC_kind_u(I,j+1) == OBC_NONE) OBC%OBC_kind_u(I,j+1) = OBC_FLATHER_N
        endif
        if ((I-1>Isdq) .and. (I-1<Iedq) .and. (j+1>jsd) .and. (j+1<jed)) then
          OBC%OBC_mask_u(I-1,j+1) = .true.
          if (OBC%OBC_kind_u(I-1,j+1) == OBC_NONE) OBC%OBC_kind_u(I-1,j+1) = OBC_FLATHER_N
        endif
     endif
    enddo ; enddo
  endif
  
  if (apply_OBC_v_flather_south) then
    ! Determine where v points are applied at south side 
    do J=Jsdq,Jedq ; do i=isd,ied
      if ((J+jsd_global-jsd) .eq. south_boundary) then         !southern side
        OBC%OBC_mask_v(i,J) = .true.
        OBC%OBC_kind_v(i,J) = OBC_FLATHER_S
        if ((I>Isdq) .and. (I<Iedq) .and. (j>jsd) .and. (j<jed)) then
          OBC%OBC_mask_u(I,j) = .true.
          if (OBC%OBC_kind_u(I,j) == OBC_NONE) OBC%OBC_kind_u(I,j) = OBC_FLATHER_S
        endif
        if ((I-1>Isdq) .and. (I-1<Iedq) .and. (j>jsd) .and. (j<jed)) then
          OBC%OBC_mask_u(I-1,j) = .true.
          if (OBC%OBC_kind_u(I-1,j) == OBC_NONE) OBC%OBC_kind_u(I-1,j) = OBC_FLATHER_S
        endif
      endif
    enddo ; enddo
  endif

  !   If there are no OBC points on this PE, there is no reason to keep the OBC
  ! type, and it could be deallocated.


  ! Define radiation coefficients r[xy]_old_[uvh] as needed.  For now, there are
  ! no radiation conditions applied to the thicknesses, since the thicknesses
  ! might not be physically motivated.  Instead, sponges should be used to
  ! enforce the near-boundary layer structure.
  if (apply_OBC_u_flather_west .or. apply_OBC_u_flather_east) then
    allocate(OBC%rx_old_u(Isdq:Iedq,jsd:jed,nz)) ; OBC%rx_old_u(:,:,:) = 0.0
 !   allocate(OBC%rx_old_h(Isd:Ied,jsd:jed,nz))   ; OBC%rx_old_h(:,:,:) = 0.0
  endif
  if (apply_OBC_v_flather_south .or. apply_OBC_v_flather_north) then
    allocate(OBC%ry_old_v(isd:ied,Jsdq:Jedq,nz)) ; OBC%ry_old_v(:,:,:) = 0.0
 !   allocate(OBC%ry_old_h(isd:ied,Jsd:Jed,nz))   ; OBC%ry_old_h(:,:,:) = 0.0
  endif

end subroutine set_Flather_Bdry_Conds   
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_ref_profile(G, PF, Compress_CSp)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in) :: PF
  type(Compress_CS), pointer   :: Compress_CSp
!   This subroutine specifies a reference profile to use to compensate for the
! compressibility of seawater.  This profile should be both representative and
! smooth.  A fairly uniform (with depth) distribution of points should be
! selected.
!   This subroutine should only be used with a nonlinear equation of state.                             
! Arguments: G - The ocean's grid structure.
!  (in)      PF - A structure indicating the open file to parse for model
!                 parameter values.
!  (in/out)  Compress_CSp - A pointer to the control structure for the
!                           CompressComp module.

  character(len=200) :: in_dir, filename   ! Strings for directory/file.
  character(len=200) :: Compress_file      ! The name of the file that
    ! the model will generate or reread.
  character(len=200) :: Ref_file, Reference_file ! The name of a file with
    ! spatially varying reference temperature, salinities, and depths.
  character(len=80)  :: Temp_name, Salt_name, Z_name  ! The names of the
    ! temperature, salinity, and depth variables in Reference_file.
  character(len=40)  :: mod = "set_ref_profile" ! This subroutine's name.
  real :: smooth_len   ! The length scale in m over which the reference
                       ! compressibility file is smoothed where the ocean is
                       ! 4000 m deep (the smoothing scales as the ocean depth
                       ! squared).
  logical :: ref_comp_3d ! If true, the reference compressibilities vary both
                         ! horizontally and vertically.  Otherwise, they vary
                         ! only in the vertical.
  integer, parameter :: PROF_LEN = 27  ! The length of a hard-wired reference
                         ! profile of temperature, salinity, and height.
  real, dimension(PROF_LEN) :: &
    T_ref, &      ! A hard-wired reference potential temperature profile, in C.
    S_ref, &      ! A hard-wired reference salinity profile, in psu.
    e_ref         ! A hard-wired reference height profile, in m.

  call GOLD_mesg("  GOLD_initialization.F90, set_ref_profile: "// &
                "setting reference compressibility", 5)

  ref_comp_3d = .true. ; call read_param(PF,"REF_COMPRESS_3D",ref_comp_3d)

  in_dir = "."; call read_param(PF,"INPUTDIR",in_dir) ; in_dir = slasher(in_dir)
  Ref_file = "" ; call read_param(PF,"REFERENCE_COMPRESS_FILE",Ref_file)
  Reference_file = trim(in_dir)//trim(Ref_file)
  Temp_name = "PTEMP" ; Salt_name = "SALT" ; Z_name = "ZT"
  call read_param(PF,"REF_COMPRESS_FILE_TEMP",Temp_name)
  call read_param(PF,"REF_COMPRESS_FILE_SALT",Salt_name)
  call read_param(PF,"REF_COMPRESS_FILE_DEPTH",Z_name)

  if (ref_comp_3d) then
    filename = "GOLD_Compress.nc" ; call read_param(PF,"COMPRESS_FILE",filename)
    Compress_file = trim(in_dir)//trim(filename)

    smooth_len = 1.0e6; call read_param(PF,"REF_COMPRESS_SMOOTH_LEN",smooth_len)

    call read_compress(Reference_file, Temp_name, Salt_name, Z_name, &
                       smooth_len,4000.0,.true.,Compress_file, G, Compress_CSp)
  else
    if (len_trim(Ref_file) > 0.0) then
      call GOLD_error(FATAL,"  GOLD_initialize, set_ref_profile: "//&
                "1-D reference profiles from files not implemented yet.")
    ! call read_compress_1d(Reference_file, Temp_name, Salt_name, Z_name, &
    !                       G, Compress_CSp)
    else
      ! This particular example is for a reference profile in the eastern
      ! North Atlantic, somewhere near 20 W, 30 N.  PROF_LEN must be 27.
      T_ref(:) = 3.81
      S_ref(:) = 34.98
      e_ref = (/    0,  -100,  -350,  -600,  -800, -1000, -1200, -1400, -1600, -1800, &
                -2000, -2200, -2400, -2600, -2800, -3000, -3200, -3400, -3600, -3800, &
                -4000, -4200, -4400, -4600, -4800, -5000, -5200/)
      call register_compress(T_ref, S_ref, e_ref, PROF_LEN, G, Compress_CSp)
    endif
  endif

  call log_param(PF, mod, "REF_COMPRESS_3D", ref_comp_3d, &
                 "If true, compensate for a 3-d reference compressiblity \n"//&
                 "read from REFERENCE_COMPRESS_FILE with variable names \n"//&
                 "given by REF_COMPRESS_{TEMP,SALT,DEPTH}.", default=.true.)

  call log_param(PF, mod, "REFERENCE_COMPRESS_FILE", Ref_file, &
                 "The file from which to read reference compressiblities, \n"//&
                 "or blank to use constant temperatures and salnities.", &
                 default="")
  if (ref_comp_3d) then
    call log_param(PF, mod, "INPUTDIR/REFERENCE_COMPRESS_FILE", Reference_file)
    call log_param(PF, mod, "REF_COMPRESS_FILE_TEMP", Temp_name, &
                 "The reference compressibility temperature variable name.", &
                 default="PTEMP")
    call log_param(PF, mod, "REF_COMPRESS_FILE_SALT", Salt_name, &
                 "The reference compressibility salinity variable name.", &
                 default="SALT")
    call log_param(PF, mod, "REF_COMPRESS_FILE_DEPTH", Z_name, &
                 "The reference compressibility depth variable name.", &
                 default="ZT")
  endif

end subroutine set_ref_profile
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine reset_face_lengths_global_1deg(G, param_file)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in) :: param_file
!   This subroutine sets the open face lengths at selected points to restrict
! passages to their observed widths.

! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

  character(len=256) :: mesg    ! Message for error messages.
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq

  do j=jsd,jed ; do I=Isdq,Iedq
    if ((abs(G%geolatu(I,j)-35.5) < 0.5) .and. (G%geolonu(I,j) < -4.5) .and. &
        (G%geolonu(I,j) > -6.5)) &
      G%dy_u(I,j) = G%umask(I,j)*12000.0   ! Gibraltar

    if ((abs(G%geolatu(I,j)-12.5) < 0.5) .and. (abs(G%geolonu(I,j)-43.0) < 0.5)) &
      G%dy_u(I,j) = G%umask(I,j)*10000.0   ! Red Sea

    if ((abs(G%geolatu(i,j)-40.5) < 0.5) .and. (abs(G%geolonu(i,j)-26.0) < 0.5)) &
      G%dy_u(i,j) = G%umask(i,j)*5000.0   ! Dardanelles

    if ((abs(G%geolatu(I,j)-80.84) < 0.2) .and. (abs(G%geolonu(I,j)+64.9) < 0.8)) &
      G%dy_u(I,j) = G%umask(I,j)*38000.0   ! Smith Sound in Canadian Arch - tripolar region

    if ((abs(G%geolatu(I,j)-41.5) < 0.5) .and. (abs(G%geolonu(I,j)+220.0) < 0.5)) &
      G%dy_u(I,j) = G%umask(I,j)*35000.0   ! Tsugaru strait at 140.0e

    if ((abs(G%geolatu(I,j)-45.5) < 0.5) .and. (abs(G%geolonu(I,j)+217.5) < 0.9)) &
      G%dy_u(I,j) = G%umask(I,j)*15000.0   ! Betw Hokkaido and Sakhalin at 217&218 = 142e

    ! Any u-face lengths should be added above this point.

    if (G%dy_u(I,j) > G%DYu(I,j)) then
      write(mesg,'("dy_u of ",ES11.4," exceeds unrestricted width of ",ES11.4,&
                   &" by ",ES11.4," at lon/lat of ", ES11.4, ES11.4)') &
                   G%dy_u(I,j), G%DYu(I,j), G%dy_u(I,j)-G%DYu(I,j), &
                   G%geolonu(I,j), G%geolatu(I,j)
      call GOLD_error(FATAL,"reset_face_lengths_global_1deg "//mesg)
    endif
    G%dxdy_u(I,j) = G%DXu(I,j)*G%dy_u(I,j)
    G%Idxdy_u(I,j) = 0.0
    if (G%dxdy_u(I,j) > 0.0) G%Idxdy_u(I,j) = G%umask(I,j) / G%dxdy_u(I,j)
  enddo ; enddo

  do J=Jsdq,Jedq ; do i=isd,ied
    if ((abs(G%geolatv(i,J)-41.0) < 0.5) .and. (abs(G%geolonv(i,J)-28.5) < 0.5)) &
      G%dx_v(i,J) = G%vmask(i,J)*2500.0   ! Bosporus - should be 1000.0 m wide.

    if ((abs(G%geolatv(i,J)-13.0) < 0.5) .and. (abs(G%geolonv(i,J)-42.5) < 0.5)) &
      G%dx_v(i,J) = G%vmask(i,J)*10000.0   ! Red Sea

    if ((abs(G%geolatv(i,J)+2.8) < 0.8) .and. (abs(G%geolonv(i,J)+241.5) < 0.5)) &
      G%dx_v(i,J) = G%vmask(i,J)*40000.0   ! Makassar Straits at 241.5 W = 118.5 E

    if ((abs(G%geolatv(i,J)-0.56) < 0.5) .and. (abs(G%geolonv(i,J)+240.5) < 0.5)) &
      G%dx_v(i,J) = G%vmask(i,J)*80000.0   ! entry to Makassar Straits at 240.5 W = 119.5 E

    if ((abs(G%geolatv(i,J)-0.19) < 0.5) .and. (abs(G%geolonv(i,J)+230.5) < 0.5)) &
      G%dx_v(i,J) = G%vmask(i,J)*25000.0   ! Channel betw N Guinea and Halmahara 230.5 W = 129.5 E

    if ((abs(G%geolatv(i,J)-0.19) < 0.5) .and. (abs(G%geolonv(i,J)+229.5) < 0.5)) &
      G%dx_v(i,J) = G%vmask(i,J)*25000.0   ! Channel betw N Guinea and Halmahara 229.5 W = 130.5 E

    if ((abs(G%geolatv(i,J)-0.0) < 0.25) .and. (abs(G%geolonv(i,J)+228.5) < 0.5)) &
      G%dx_v(i,J) = G%vmask(i,J)*25000.0   ! Channel betw N Guinea and Halmahara 228.5 W = 131.5 E

    if ((abs(G%geolatv(i,J)+8.5) < 0.5) .and. (abs(G%geolonv(i,J)+244.5) < 0.5)) &
      G%dx_v(i,J) = G%vmask(i,J)*20000.0   ! Lombok Straits at 244.5 W = 115.5 E

    if ((abs(G%geolatv(i,J)+8.5) < 0.5) .and. (abs(G%geolonv(i,J)+235.5) < 0.5)) &
      G%dx_v(i,J) = G%vmask(i,J)*20000.0   ! Timor Straits at 235.5 W = 124.5 E

    if ((abs(G%geolatv(i,J)-76.8) < 0.06) .and. (abs(G%geolonv(i,J)+88.7) < 0.5)) &
      G%dx_v(i,J) = G%vmask(i,J)*8400.0   ! Jones Sound in Canadian Arch - tripolar region

    if ((abs(G%geolatv(i,J)-52.5) < 0.5) .and. (abs(G%geolonv(i,J)+218.5) < 0.5)) &
      G%dx_v(i,J) = G%vmask(i,J)*2500.0   ! Russia and Sakhalin Straits at 218.5 W = 141.5 E

    ! Any v-face lengths should be added above this point.

    if (G%dx_v(i,J) > G%DXv(i,J)) then
      write(mesg,'("dx_v of ",ES11.4," exceeds unrestricted width of ",ES11.4,&
                   &" by ",ES11.4, " at lon/lat of ", ES11.4, ES11.4)') &
                   G%dx_v(i,J), G%DXv(i,J), G%dx_v(i,J)-G%DXv(i,J), &
                   G%geolonv(i,J), G%geolatv(i,J)

      call GOLD_error(FATAL,"reset_face_lengths_global_1deg "//mesg)
    endif
    G%dxdy_v(i,J) = G%DYv(i,J)*G%dx_v(i,J)
    G%Idxdy_v(i,J) = 0.0
    if (G%dxdy_v(i,J) > 0.0) G%Idxdy_v(i,J) = G%vmask(i,J) / G%dxdy_v(i,J)
  enddo ; enddo

end subroutine reset_face_lengths_global_1deg
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_velocity_depth_max(G)
  type(ocean_grid_type), intent(inout) :: G
  ! This subroutine sets the 4 bottom depths at velocity points to be the
  ! maximum of the adjacent depths.
  integer :: i, j

  do I=G%isd,G%ied-1 ; do j=G%jsd,G%jed
    G%Dblock_u(I,j) = G%umask(I,j)*max(G%D(i,j),G%D(i+1,j))
    G%Dopen_u(I,j) = G%Dblock_u(I,j)
  enddo ; enddo
  do i=G%isd,G%ied ; do J=G%jsd,G%jed-1
    G%Dblock_v(I,J) = G%vmask(i,J)*max(G%D(i,j),G%D(i,j+1))
    G%Dopen_v(I,J) = G%Dblock_v(I,J)
  enddo ; enddo
end subroutine set_velocity_depth_max
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_velocity_depth_min(G)
  type(ocean_grid_type), intent(inout) :: G
  ! This subroutine sets the 4 bottom depths at velocity points to be the
  ! minimum of the adjacent depths.
  integer :: i, j

  do I=G%isd,G%ied-1 ; do j=G%jsd,G%jed
    G%Dblock_u(I,j) = G%umask(I,j)*min(G%D(i,j),G%D(i+1,j))
    G%Dopen_u(I,j) = G%Dblock_u(I,j)
  enddo ; enddo
  do i=G%isd,G%ied ; do J=G%jsd,G%jed-1
    G%Dblock_v(I,J) = G%vmask(i,J)*min(G%D(i,j),G%D(i,j+1))
    G%Dopen_v(I,J) = G%Dblock_v(I,J)
  enddo ; enddo
end subroutine set_velocity_depth_min
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine write_ocean_geometry_file(G, param_file, directory)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in) :: param_file
  character(len=*),      intent(in) :: directory
!   This subroutine writes out a file containing all of the ocean geometry
! and grid data uses by the GOLD ocean model.
! Arguments: G - The ocean's grid structure.  Effectively intent in.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      directory - The directory into which to place the file.
  character(len=120) :: filepath
  integer, parameter :: nFlds=23
  type(vardesc) :: vars(nFlds)
  type(fieldtype) :: fields(nFlds)
  integer :: unit
  integer :: file_threading
  integer :: nFlds_used
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq
  logical :: multiple_files = .false.
  real :: out_h(SZI_(G),SZJ_(G))
  real :: out_u(SZIQ_(G),SZJ_(G))
  real :: out_v(SZI_(G),SZJQ_(G))
  real :: out_q(SZIQ_(G),SZJQ_(G))
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq

!   vardesc is a structure defined in GOLD_io.F90.  The elements of
! this structure, in order, are:
! (1) the variable name for the NetCDF file
! (2) the variable's long name
! (3) a character indicating the  horizontal grid, which may be '1' (column),
!     'h', 'q', 'u', or 'v', for the corresponding C-grid variable
! (4) a character indicating the vertical grid, which may be 'L' (layer),
!     'i' (interface), or '1' (no vertical location)
! (5) a character indicating the time levels of the field, which may be
!    's' (snap-shot), 'a' (average between snapshots), 'm' (monthly average),
!     or '1' (no time variation) between snapshots
! (6) the variable's units
! (7) a character indicating the size in memory to write, which may be
!     'd' (8-byte) or 'f' (4-byte).
  vars(1) = vardesc("geolatb","latitude at q points",'q','1','1',"degree",'d')
  vars(2) = vardesc("geolonb","longitude at q points",'q','1','1',"degree",'d')
  vars(3) = vardesc("geolat", "latitude at h points", 'h','1','1',"degree",'d')
  vars(4) = vardesc("geolon","longitude at h points",'h','1','1',"degree",'d')
  vars(5) = vardesc("D","Basin Depth",'h','1','1',"meter",'d')
  vars(6) = vardesc("f","Coriolis Parameter",'q','1','1',"second-1", 'd')
  vars(7) = vardesc("dxv","Zonal grid spacing at v points",'v','1','1',"m",'d')
  vars(8) = vardesc("dyu","Meridional grid spacing at u points",'u','1','1',"m",'d')
  vars(9) = vardesc("dxu","Zonal grid spacing at u points",'u','1','1',"m",'d')
  vars(10)= vardesc("dyv","Meridional grid spacing at v points",'v','1','1',"m",'d')
  vars(11)= vardesc("dxh","Zonal grid spacing at h points",'h','1','1',"m",'d')
  vars(12)= vardesc("dyh","Meridional grid spacing at h points",'h','1','1',"m",'d')
  vars(13)= vardesc("dxq","Zonal grid spacing at q points",'q','1','1',"m",'d')
  vars(14)= vardesc("dyq","Meridional grid spacing at q points",'q','1','1',"m",'d')
  vars(15)= vardesc("Ah","Area of h cells",'h','1','1',"m2",'d')
  vars(16)= vardesc("Aq","Area of q cells",'q','1','1',"m2",'d')

  vars(17)= vardesc("dxvo","Open zonal grid spacing at v points",'v','1','1',"m",'d')
  vars(18)= vardesc("dyuo","Open meridional grid spacing at u points",'u','1','1',"m",'d')
  vars(19)= vardesc("wet", "land or ocean?", 'h','1','1',"none",'d')

  vars(20) = vardesc("Dblock_u","Blocked depth at u points",'u','1','1',"meter",'d')
  vars(21) = vardesc("Dopen_u","Open depth at u points",'u','1','1',"meter",'d')
  vars(22) = vardesc("Dblock_v","Blocked depth at v points",'v','1','1',"meter",'d')
  vars(23) = vardesc("Dopen_v","Open depth at v points",'v','1','1',"meter",'d')

  nFlds_used = 19 ; if (G%bathymetry_at_vel) nFlds_used = 23

  filepath = trim(directory) // "ocean_geometry"

  out_h(:,:) = 0.0
  out_u(:,:) = 0.0
  out_v(:,:) = 0.0
  out_q(:,:) = 0.0

  call read_param(param_file,"PARALLEL_RESTARTFILES",multiple_files)
  file_threading = SINGLE_FILE
  if (multiple_files) file_threading = MULTIPLE

  call create_file(unit, trim(filepath), vars, nFlds_used, G, fields, file_threading)

  do J=Jsq,Jeq; do I=Isq,Ieq; out_q(I,J) = G%geolatq(I,J); enddo; enddo
  call write_field(unit, fields(1), G%Domain%mpp_domain, out_q)
  do J=Jsq,Jeq; do I=Isq,Ieq; out_q(I,J) = G%geolonq(I,J); enddo; enddo
  call write_field(unit, fields(2), G%Domain%mpp_domain, out_q)
  call write_field(unit, fields(3), G%Domain%mpp_domain, G%geolath)
  call write_field(unit, fields(4), G%Domain%mpp_domain, G%geolonh)

  call write_field(unit, fields(5), G%Domain%mpp_domain, G%D)
  call write_field(unit, fields(6), G%Domain%mpp_domain, G%f)

  do J=Jsq,Jeq; do i=is,ie; out_v(i,J) = G%DXv(i,J); enddo; enddo
  call write_field(unit, fields(7), G%Domain%mpp_domain, out_v)
  do j=js,je; do I=Isq,Ieq; out_u(I,j) = G%DYu(I,j); enddo; enddo
  call write_field(unit, fields(8), G%Domain%mpp_domain, out_u)

  do J=Jsq,Jeq; do i=is,ie; out_u(i,J) = G%DXu(i,J); enddo; enddo
  call write_field(unit, fields(9), G%Domain%mpp_domain, out_u)
  do j=js,je; do I=Isq,Ieq; out_v(I,j) = G%DYv(I,j); enddo; enddo
  call write_field(unit, fields(10), G%Domain%mpp_domain, out_v)

  do J=Jsq,Jeq; do i=is,ie; out_h(i,J) = G%DXh(i,J); enddo; enddo
  call write_field(unit, fields(11), G%Domain%mpp_domain, out_h)
  do j=js,je; do I=Isq,Ieq; out_h(I,j) = G%DYh(I,j); enddo; enddo
  call write_field(unit, fields(12), G%Domain%mpp_domain, out_h)

  do J=Jsq,Jeq; do i=is,ie; out_q(i,J) = G%DXq(i,J); enddo; enddo
  call write_field(unit, fields(13), G%Domain%mpp_domain, out_q)
  do j=js,je; do I=Isq,Ieq; out_q(I,j) = G%DYq(I,j); enddo; enddo
  call write_field(unit, fields(14), G%Domain%mpp_domain, out_q)

  do j=js,je; do i=is,ie; out_h(i,j) = G%DXDYh(i,j); enddo; enddo
  call write_field(unit, fields(15), G%Domain%mpp_domain, out_h)
  do j=js,je; do i=is,ie; out_q(i,j) = G%DXDYq(i,j); enddo; enddo
  call write_field(unit, fields(16), G%Domain%mpp_domain, out_q)

!  do J=Jsq,Jeq; do i=is,ie; out_v(i,J) = G%dx_v(i,J); enddo; enddo
  call write_field(unit, fields(17), G%Domain%mpp_domain, G%dx_v)
!  do j=js,je; do I=Isq,Ieq; out_u(I,j) = G%dy_u(I,j); enddo; enddo
  call write_field(unit, fields(18), G%Domain%mpp_domain, G%dy_u)
  call write_field(unit, fields(19), G%Domain%mpp_domain, G%hmask)

  if (G%bathymetry_at_vel) then
    call write_field(unit, fields(20), G%Domain%mpp_domain, G%Dblock_u)
    call write_field(unit, fields(21), G%Domain%mpp_domain, G%Dopen_u)
    call write_field(unit, fields(22), G%Domain%mpp_domain, G%Dblock_v)
    call write_field(unit, fields(23), G%Domain%mpp_domain, G%Dopen_v)
  endif

  call close_file(unit)

end subroutine write_ocean_geometry_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine write_vertgrid_file(G, param_file, directory)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in) :: param_file
  character(len=*),      intent(in) :: directory
!   This subroutine writes out a file containing any available data related
! to the vertical grid used by the GOLD ocean model.
! Arguments: G - The ocean's grid structure.  Effectively intent in.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      directory - The directory into which to place the file.
  character(len=120) :: filepath
  type(vardesc) :: vars(2)
  type(fieldtype) :: fields(2)
  integer :: unit

  filepath = trim(directory) // trim("Vertical_coordinate")

  vars(1) = vardesc("R","Target Potential Density",'1','L','1',"kilogram meter-3", 'd')
  vars(2) = vardesc("g","Reduced gravity",'1','L','1',"meter second-2", 'd')

  call create_file(unit, trim(filepath), vars, 2, G, fields, SINGLE_FILE)

  call write_field(unit, fields(1), G%Rlay)
  call write_field(unit, fields(2), G%g_prime)

  call close_file(unit)

end subroutine write_vertgrid_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine Get_GOLD_Input(param_file, dirs, check_params)
  use GOLD_io, only : open_namelist_file, check_nml_error
  type(param_file_type), optional, intent(out) :: param_file
  type(directories),     optional, intent(out) :: dirs
  logical,               optional, intent(in)  :: check_params

!    See if the run is to be started from saved conditions, and get  !
!  the names of the I/O directories and initialization file.  This   !
!  subroutine also calls the subroutine that allows run-time changes !
!  in parameters.                                                    !
!    Read from a namelist instead of getting strings from the command line as
!  in the C-version of GOLD.
  integer, parameter :: npf = 5 ! Maximum number of parameter files
  character(len=120) :: &
    parameter_filename(npf) = ' ', & ! List of files containing parameters.
    output_directory = ' ', &   ! Directory to use to write the model output.
    restart_input_dir = ' ', &  ! Directory for reading restart and input files.
    restart_output_dir = ' ', & ! Directory into which to write restart files.
    input_filename  = ' '       ! A string that indicates the input files or how
                                ! the run segment should be started.
  integer :: unit, io, ierr

  namelist /GOLD_input_nml/ output_directory, input_filename, parameter_filename, &
                           restart_input_dir, restart_output_dir

  if (file_exists('input.nml')) then
    unit = open_namelist_file(file='input.nml')
  else
    call GOLD_error(FATAL,'Required namelist file input.nml does not exist.')
  endif

  ierr=1 ; do while (ierr /= 0)
    read(unit, nml=GOLD_input_nml, iostat=io, end=10)
    ierr = check_nml_error(io, 'GOLD_input_nml')
  enddo
10 call close_file(unit)
  if (present(dirs)) then
    dirs%output_directory = slasher(output_directory)
    dirs%restart_output_dir = slasher(restart_output_dir)
    dirs%restart_input_dir = slasher(restart_input_dir)
    dirs%input_filename = input_filename
  endif

  if (present(param_file)) then ; do io = 1, npf
    if (len_trim(trim(parameter_filename(io))) > 0) &
      call open_param_file(trim(parameter_filename(io)), param_file, check_params)
  enddo ; endif

end subroutine Get_GOLD_Input
! -----------------------------------------------------------------------------

end module GOLD_initialization
