module GOLD
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
!*            The Generalized Ocean Layer Dynamics model               *
!*                               GOLD                                  *
!*                                                                     *
!*  By Robert Hallberg and Alistair Adcroft                            *
!*                                                                     *
!*  With software contributions from:                                  *
!*    Whit Anderson, Brian Arbic, Will Cooke, Anand Gnanadesikan,      *
!*    Matthew Harrison, Mehmet Ilicak, Laura Jackson, Jasmine John,    *
!*    Harper Simmons, and Niki Zadeh                                   *
!*                                                                     *
!*  GOLD ice-shelf code by Daniel Goldberg, Robert Hallberg            *
!*    Chris Little, and Olga Sergienko                                 *
!*  Original C-language HIM code by R. Hallberg, translated into F90   *
!*    by Will Cooke, 2003-2004.                                        *
!*                                                                     *
!*  This file worked on 1992 - 2011.                                   *
!*                                                                     *
!*    This program (GOLD) simulates the ocean by numerically solving   *
!*  the hydrostatic primitive equations in generalized Lagrangian      *
!*  vertical coordinates, often following isopycnals in the ocean's    *
!*  interior, and general orthogonal horizontal coordinates.  These    *
!*  equations are horizontally discretized on an Arakawa C-grid.       *
!*  There are a range of options for the physical parameterizations,   *
!*  from those most appropriate to highly idealized models for studies *
!*  of geophysical fluid dynamics to a rich suite of processes         *
!*  appropriate for realistic ocean simulations.  The thermodynamic    *
!*  options range from an adiabatic model with fixed density layers    *
!*  to a model with temperature and salinity as state variables and    *
!*  using a full nonlinear equation of state.  GOLD was originally     *
!*  developed starting from the F90 version of the Hallberg Isopycnal  *
!*  Model (HIM) about 2005, and has been greatly extended beyond what  *
!*  was available at that time.  GOLD has also benefitted tremendously *
!*  from the FMS infrastucture, which it utilizes, and from extensive  *
!*  comparisons with simulations using the MOM4.1 ocean model.         *
!*                                                                     *
!*    When run is isopycnal-coordinate mode, the uppermost few layers  *
!*  are often used to describe a bulk mixed layer, including the       *
!*  effects of penetrating shortwave radiation.  Either a split-       *
!*  explicit time stepping scheme or a non-split scheme may be used    *
!*  for the dynamics, while the time stepping may be split (and use    *
!*  different numbers of steps to cover the same interval) for the     *
!*  forcing, the thermodynamics, and for the dynamics.  Most of the    *
!*  numerics are second order accurate in space.  GOLD can run with an *
!*  absurdly thin minimum layer thickness. A variety of non-isopycnal  *
!*  vertical coordinate options are under development, but all exploit *
!*  the advantages of a Lagrangian vertical coordinate, as discussed   *
!*  in detail by Adcroft and Hallberg (Ocean Modelling, 2006).         *
!*                                                                     *
!*    Details of the numerics and physical parameterizations are       *
!*  provided in the appropriate source files.  All of the available    *
!*  options are selected by the settings in the input file, GOLD_input.*
!*                                                                     *
!*   The present file (GOLD.F90) contains the main time stepping loops.*
!*  One time integration option for the dynamics uses a split explicit *
!*  time stepping scheme to rapidly step the barotropic pressure and   *
!*  velocity fields. The barotropic velocities are averaged over the   *
!*  baroclinic time step before they are used to advect thickness      *
!*  and determine the baroclinic accelerations. At the end of every    *
!*  time step, the free surface height perturbation is determining     *
!*  by adding up the layer thicknesses; this perturbation is used to   *
!*  drive the free surface heights from the barotropic calculation     *
!*  and from the sum of the layer thicknesses toward each other over   *
!*  subsequent time steps. The barotropic and baroclinic velocities    *
!*  are synchronized as part of the vertical viscosity algorithm and   *
!*  be recalculating. the barotropic velocities from the baroclinic    *
!*  velocities each time step. This scheme is described in Hallberg,   *
!*  1997, J. Comp. Phys. 135, 54-65.                                   *
!*                                                                     *
!*    The other time integration option uses a non-split time stepping *
!*  scheme based on the 3-step third order Runge-Kutta scheme          *
!*  described in Matsuno, 1966, J. Met. Soc. Japan,  44, 85-88.        *
!*                                                                     *
!*    There are a range of closure options available.  Horizontal      *
!*  velocities are subject to a combination of horizontal biharmonic   *
!*  and Laplacian friction (based on a stress tensor formalism) and a  *
!*  vertical Fickian viscosity (perhaps using the kinematic viscosity  *
!*  of water).  The horizontal viscosities may be constant, spatially  *
!*  varying or may be dynamically calculated using Smagorinsky's       *
!*  approach.  A diapycnal diffusion of density and thermodynamic      *
!*  quantities is also allowed, but not required, as is horizontal     *
!*  diffusion of interface heights (akin to the Gent-McWilliams        *
!*  closure of geopotential coordinate models).  The diapycnal mixing  *
!*  may use a fixed diffusivity or it may use the shear Richardson     *
!*  number dependent closure described in Hallberg (MWR, 2000).        *
!*  When there is diapycnal diffusion, it applies to momentum as well. *
!*  As this is in addition to the vertical viscosity, the vertical     *
!*  Prandtl always exceeds 1.                                          *
!*                                                                     *
!*    GOLD has a number of noteworthy debugging capabilities.          *
!*  Excessively large velocities are truncated and GOLD will stop      *
!*  itself after a number of such instances to keep the model from     *
!*  crashing altogether, and the model state is output with a          *
!*  reported time of 9.9e9.  This is useful in diagnosing failures,    *
!*  or (by accepting some truncations) it may be useful for getting    *
!*  the model past the adjustment from an ill-balanced initial         *
!*  condition.  In addition, all of the accelerations in the columns   *
!*  with excessively large velocities may be directed to a text file.  *
!*  Parallelization errors may be diagnosed with the CHECK_PARALLEL    *
!*  option, whereby ostensibly identical model incarnations are run    *
!*  simultaneously on one and multiple processors and any differences  *
!*  are reported.                                                      *
!*                                                                     *
!*    About 35 other files of source code and 6 header files must      *
!*  be used to make the model work.  A makefile is included to         *
!*  make compiling easy.  Type "make GOLD" to compile.  Some run       *
!*  time input is required, but this is prompted for.  It may be       *
!*  convenient to direct a small file containing the needed inform-    *
!*  ation to standard input, and direct the output to another file.    *
!*                                                                     *
!*    The ~35 source files contain the following subroutines:          *
!*  GOLD/src/core/GOLD.F90:                                            *
!*    step_GOLD steps GOLD over a specified interval of time.          *
!*    GOLD_initialize calls initialize and does other initialization   *
!*      that does not warrant user modification.                       *
!*    set_restart_fields is used to specify those fields that are      *
!*      written to and read from the restart file.                     *
!*    calculate_surface_state determines the surface (mixed layer)     *
!*      properties of the current model state and packages pointers    *
!*      to these fields into an exported structure.                    *
!*  GOLD/config_src/solo_driver/GOLD_driver.F90:                       *
!*    main is where GOLD starts.  Inside of main are the calls that    *
!*      set up the run, step the model, and orchestrate output and     *
!*      normal termination of the run.                                 *
!*                                                                     *
!*                                                                     *
!*     THE FOLLOWING FILES ARE WHERE INITIAL CONDITIONS, FORCING,      *
!*     AND DOMAIN PROPERTIES ARE PRINCIPALLY SPECIFIED.                *
!*                                                                     *
!*  GOLD/src/initialization/GOLD_initialization.F90:                   *
!*    GOLD_initialize does just that to all of the fields that are     *
!*      needed to specify the initial conditions of the model.         *
!*      GOLD_initialize calls a number of other subroutines in         *
!*      GOLD_initialization.F90, each of which initializes a single    *
!*      field (or a few closely related fields) that are indicated by  *
!*      the subroutine name.                                           *
!*       (initialize_output.c has been supplanted by diag_manager.F90, *
!*       which parses the diag table, these last two routines are in   *
!*       GOLD_initialization.F90)                                      *
!*    Get_GOLD_Input gets 5 controlling inputs from a namelist, setting*
!*      the directories for I/O and the parameter specification file.  *
!*    write_grid_file writes out a file describing the model grid.     *
!*  GOLD/config_src/solo_driver/GOLD_surface_forcing.F90:              *
!*    set_forcing sets the current values of surface forcing fields.   *
!*    wind_forcing sets the current surface wind stresses.             *
!*    buoyancy_forcing sets the current surface heat, fresh water,     *
!*      buoyancy or other appropriate tracer fluxes.                   *
!*    set_forcing_output sets up the output of any forcing fields.     *
!*    average_forcing accumulates time averages of indicated forcing   *
!*      fields.                                                        *
!*    register_forcing_restarts is used to specify the forcing-related *
!*      fields that are written to and read from the restart file.     *
!*  GOLD/src/initialization/GOLD_grid_init.F90:                        *
!*    set_metrics calculates the horizontal grid spacings and related  *
!*      metric fields, along with the gridpoint locations.             *
!*    initialize_masks initializes the land masks.                     *
!*                                                                     *
!*     THE FOLLOWING FILES CONTAIN THE PRINCIPAL DYNAMIC ROUTINES.     *
!*                                                                     *
!*  GOLD/src/core/GOLD_CoriolisAdv.F90:                                *
!*    CorAdCalc calculates the Coriolis and advective accelerations.   *
!*  GOLD/src/core/GOLD_PressureForce.F90:                              *
!*    PressureForce calculates the pressure acceleration.              *
!*  GOLD/src/core/GOLD_CompressComp.F90:                               *
!*    register_compress is used to specify the reference profile of    *
!*      potential temperature and salinity that is used to compensate  *
!*      for compressibility.                                           *
!*    uncompress_e_rho makes internally consistent changes to profiles *
!*      of interface height and density to offset compressibility and  *
!*      to minimize the non-solenoidal pressure gradient term.         *
!*  GOLD/src/core/GOLD_continuity.F90, GOLD_continuity_PPM.F90,        *
!*  or GOLD_continuity_Hallberg95.F90:                                 *
!*    continuity time steps the layer thicknesses.                     *
!*  GOLD/src/core/GOLD_barotropic.F90:                                 *
!*    btstep time steps the linearized barotropic equations for use    *
!*      with the split explicit time stepping scheme.                  *
!*    btcalc calculates the barotropic velocities from the layer       *
!*      velocities.                                                    *
!*    barotropic_init initializes several split-related variables and  *
!*      calculates several static quantities for use by btstep.        *
!*    register_barotropic_restarts indicates those time splitting-     *
!*      related fields that are to be in the restart file.             *
!*                                                                     *
!*  GOLD/src/parameterizations/lateral/GOLD_hor_visc.F90:              *
!*    horizontal_viscosity calculates the convergence of momentum      *
!*      due to Laplacian or biharmonic horizontal viscosity.           *
!*    set_up_hor_visc calculates combinations of metric coefficients   *
!*      and other static quantities used in horizontal_viscosity.      *
!*  GOLD/src/parameterizations/vertical/GOLD_vert_friction.F90:        *
!*    vertvisc changes the velocity due to vertical viscosity,         *
!*      including application of a surface stress and bottom drag.     *
!*  GOLD/src/parameterizations/vertical/GOLD_set_viscosity.F90:        *
!*    set_viscous_BBL determines the bottom boundary layer thickness   *
!*      and viscosity according to a linear or quadratic drag law.     *
!*  GOLD/src/parameterizations/lateral/GOLD_thickness_diffuse.F90:     *
!*    thickness_diffuse moves fluid adiabatically to horizontally      *
!*      diffuse interface heights.                                     *
!*                                                                     *
!*                                                                     *
!*     THE FOLLOWING FILES CONTAIN THE THERMODYNAMIC ROUTINES.         *
!*                                                                     *
!*  GOLD/src/parameterizations/vertical/GOLD_diabatic_driver.F90:      *
!*    diabatic orchestrates the calculation of vertical advection and  *
!*      diffusion of momentum and tracers due to diapycnal mixing and  *
!*      mixed layer (or other diabatic) processes.  mixedlayer,        *
!*      Calculate_Entrainment, apply_sponge, and any user-specified    *
!*      tracer column physics routines are all called by diabatic.     *
!*  GOLD/src/parameterizations/vertical/GOLD_diabatic_entrain.F90:     *
!*    Calculate_Entrainment calculates the diapycnal mass fluxes due   *
!*      to interior diapycnal mixing processes, which may include a    *
!*      Richardson number dependent entrainment.                       *
!*    Calculate_Rino_flux estimates the Richardson number dependent    *
!*      entrainment in the absence of interactions between layers,     *
!*      from which the full interacting entrainment can be found.      *
!*    Estimate_u_h estimates what the velocities at thickness points   *
!*      will be after entrainment.                                     *
!*  GOLD/src/parameterizations/vertical/GOLD_mixed_layer.F90:          *
!*    mixed_layer implements a bulk mixed layer, including entrainment *
!*      and detrainment, related advection of dynamically active       *
!*      tracers, and buffer layer splitting.  The bulk mixed layer     *
!*      may consist of several layers.                                 *
!*  GOLD/src/parameterizations/vertical/GOLD_sponge.F90:               *
!*    apply_sponge damps fields back to reference profiles.            *
!*    initialize_sponge stores the damping rates and allocates the     *
!*      memory for the reference profiles.                             *
!*    set_up_sponge_field registers reference profiles and associates  *
!*      them with the fields to be damped.                             *
!*                                                                     *
!*  GOLD/src/tracer/GOLD_tracer.F90:                                   *
!*    register_tracer is called to indicate a field that is to be      *
!*      advected by advect_tracer and diffused by tracer_hordiff       *
!*    advect_tracer does along-isopycnal advection of tracer fields.   *
!*    tracer_hordiff diffuses tracers along isopycnals.                *
!*    register_tracer_init_fn registers a user-specified tracer        *
!*      initialization subroutine.                                     *
!*    call_tracer_init_fns calls any user-specified tracer initial-    *
!*      ization subroutines that have been registered.                 *
!*    register_tracer_column_fn registers a user-specified tracer      *
!*      column processes subroutine.                                   *
!*    call_tracer_column_fns calls any user-specified tracer column    *
!*      processes subroutines that have been registered.               *
!*                                                                     *
!* In GOLD/src/equation_of_state/                                      *
!*  GOLD/src/equation_of_state/GOLD_EOS.F90,                           *
!*  GOLD_EOS_linear.F90, or GOLD_EOS_UNESCO.F90:                       *
!*    calculate_density calculates a list of densities at given        *
!*      potential temperatures, salinities and pressures.              *
!*    calculate_density_derivs calculates a list of the partial        *
!*      derivatives with temperature and salinity at the given         *
!*      potential temperatures, salinities and pressures.              *
!*    calculate_compress calculates a list of the compressibilities    *
!*      (partial derivatives of density with pressure) at the given    *
!*      potential temperatures, salinities and pressures.              *
!*    calculate_2_densities calculates a list of the densities at two  *
!*      specified reference pressures at the given potential           *
!*      temperatures and salinities.                                   *
!*  GOLD/src/equation_of_state/GOLD_fit_compressibility.F90:           *
!*    fit_compressibility determines the best fit of compressibility   *
!*      with pressure using a fixed 5-coefficient functional form,     *
!*      based on a provided reference profile of potential temperature *
!*      and salinity with depth.  This fit is used in PressureForce.   *
!*                                                                     *
!*                                                                     *
!*     THE FOLLOWING FILES CONTAIN INFRASTRUCTURAL ROUTINES.           *
!*                                                                     *
!*  GOLD/src/framework/GOLD_restart.F90:                               *
!*    save_restart saves a restart file (or multiple files if they     *
!*      would otherwise be too large).                                 *
!*    register_restart_field is called to specify a field that is to   *
!*      written to and read from the restart file.                     *
!*    restore_state reads the model state from restart or other files. *
!*    query_initialized indicates whether a specific field or all      *
!*      restart fields have been read from the restart files.          *
!*  GOLD/src/framework/GOLD_parser.F90:                                *
!*    GOLD_parser parses a parameter specification file to enable      *
!*      parameters to be set at run time.                              *
!*  GOLD/src/framework/GOLD_domains.F90:                               *
!*    pass_var passes a 2-D or 3-D variable to neighboring processors  *
!*      applies corresponding boundary conditions.                     *
!*    pass_vector passes a 2-D or 3-D pair of vector components or     *
!*      scalars to neighboring processors.                             *
!*    pass_var_start and pass_var_complete provide non-blocking halo   *
!*      updates akin to pass_var.                                      *
!*    pass_vector_start and pass_vector_complete provide non-blocking  *
!*      halo updates akin to pass_vector.                              *
!*    GOLD_domains_init initializes the computational domain.          *
!*    chksum sums the bits in an array and writes out the total.       *
!*  GOLD/src/framework/GOLD_diag_mediator.F90:                         *
!*    axes_info initiates the output axes and stores groupings of      *
!*      them in the ocean grid.                                        *
!*    post_data uses the time-weighting and time_end from              *
!*      enable_averaging in a call to send_data.                       *
!*    enable_averaging enables averaging for a time interval.          *
!*    disable_averaging disables the accumulation of averages.         *
!*    query_averaging_enabled indicates whether averaging is           *
!*      currently enabled.                                             *
!*   (GOLD_field_output.c has been replaced by the FMS diag_manager.)  *
!*  GOLD/src/framework/GOLD_io.F90: (Input/Output utility subroutines.)*
!*    create_file creates a new file, set up structures that are       *
!*      needed for subsequent output, and write the coordinates.       *
!*    reopen_file reopens an existing file for writing and set up      *
!*      structures that are needed for subsequent output.              *
!*  (The remaining functions from GOLD_io.c have been replaced by the  *
!*   FMS mpp_io, fms_io, and diag_manager.)                            *
!*  (GOLD_parallel.c has been replaced by the FMS mpp calls.)          *
!*                                                                     *
!*                                                                     *
!*     THE FOLLOWING FILES CONTAIN PURELY DIAGNOSTIC ROUTINES.         *
!*                                                                     *
!*  GOLD/src/diagnostics/GOLD_diagnostics.F90:                         *
!*    calculate_diagnostic_fields is used to calculate several         *
!*      diagnostic fields that are not naturally calculated elsewhere. *
!*    register_time_deriv is used to register the information needed   *
!*      for diagnostically calculating a time derivative.              *
!*    calculate_derivs calculates any registered time derivatives.     *
!*  GOLD/src/diagnostics/GOLD_sum_output.F90:                          *
!*    write_energy writes the layer energies and masses and other      *
!*      spatially integrated quantities and monitors CPU time use.     *
!*    depth_list_setup generates a list of the volumes of fluid below  *
!*      various depths.                                                *
!*  GOLD/src/diagnostics/PointAccel.c has not yet been translated.     *
!*    write_u_accel writes a long list of zonal accelerations and      *
!*      related quantities for one column out to a file.  This is      *
!*      typically called for diagnostic purposes from vertvisc when a  *
!*      zonal velocity exceeds the specified threshold.                *
!*    write_v_accel writes a long list of meridional accelerations and *
!*      related quantities for one column out to a file.  This is      *
!*      typically called for diagnostic purposes from vertvisc when a  *
!*      meridional velocity exceeds the specified threshold.           *
!*                                                                     *
!*                                                                     *
!*    In addition there are 4 include files:                           *
!*  GOLD/config_src/dynamic/GOLD_memory.h sets macros related to       *
!*    memory allocation, although in dynamic memory mode, the actual   *
!*    values are set at run time.                                      *
!*  GOLD/src/framework/GOLD_memory_macros.h contians a number of       *
!*    macros to enable the use of static or dynamic memory allocation. *
!*  GOLD/src/core/GOLD_grid_macros.h contains the descriptions for a   *
!*    number of metric terms.                                          *
!*  GOLD/src/parameterizations/lateral/GOLD_hor_visc.h contains the    *
!*    descriptions for a number of metric-related fields that are only *
!*    used in GOLD_hor_visc.F90 to calculate horizontal viscosity.     *
!*                                                                     *
!*                                                                     *
!*    Most simulations can be set up by modifying only the files       *
!*  GOLD_input, GOLD_initialization.F90, and GOLD_surface_forcing.F90. *
!*  In addition, the diag_table (GOLD_diag_table) will commonly be     *
!*  modified to tailor the output to the needs of the question at      *
!*  hand.  The FMS utility mkmf works with a file called path_names    *
!*  to build an appropriate makefile, and path_names should be edited  *
!*  to reflect the actual location of the desired source code.         *
!*                                                                     *
!*  Macros written all in capital letters are defined in GOLD_memory.h.*
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, f                                     *
!*    j+1  > o > o >   At ^:  v, PFv, CAv, vh, diffv, tauy, vbt, vhtr  *
!*    j    x ^ x ^ x   At >:  u, PFu, CAu, uh, diffu, taux, ubt, uhtr  *
!*    j    > o > o >   At o:  h, D, eta, T, S, tr, actflux             *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1                                                  *
!*           i  i+1                                                    *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_variables, only : directories, vertvisc_type, ocean_OBC_type
use GOLD_variables, only : BT_cont_type, alloc_bt_cont_type, dealloc_bt_cont_type
use GOLD_variables, only : &
  forcing, &      ! A structure containing pointers to the forcing fields
                  ! which may be used to drive GOLD.  All fluxes are
                  ! positive downward.
  surface, &      ! A structure containing pointers to various fields which
                  ! may be used describe the surface state of GOLD, and
                  ! which will be returned to the calling program
  thermo_var_ptrs, & ! A structure containing pointers to an assortment of
                  ! thermodynamic fields that may be available, including
                  ! potential temperature, salinity and mixed layer density.
  ocean_internal_state  ! A structure containing pointers to most of the above.

use GOLD_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use GOLD_cpu_clock, only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
use GOLD_cpu_clock, only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use GOLD_diag_mediator, only : diag_mediator_init, enable_averaging
use GOLD_diag_mediator, only : disable_averaging, post_data, safe_alloc_ptr
use GOLD_diag_mediator, only : register_diag_field, register_static_field
use GOLD_diag_mediator, only : set_diag_mediator_grid, diag_ptrs
use GOLD_domains, only : GOLD_domains_init, pass_var, pass_vector
use GOLD_domains, only : pass_var_start, pass_var_complete
use GOLD_domains, only : pass_vector_start, pass_vector_complete
use GOLD_domains, only : To_South, To_West, To_All, CGRID_NE, SCALAR_PAIR
use GOLD_checksums, only : GOLD_checksums_init, hchksum, uchksum, vchksum
use GOLD_error_handler, only : GOLD_error, GOLD_mesg, FATAL, WARNING, is_root_pe
use GOLD_error_handler, only : GOLD_set_verbosity
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_io, only : GOLD_io_init, vardesc
use GOLD_obsolete_params, only : find_obsolete_params
use GOLD_restart, only : register_restart_field, query_initialized, save_restart
use GOLD_restart, only : restart_init, GOLD_restart_CS
use GOLD_time_manager, only : time_type, set_time, time_type_to_real, operator(+)
use GOLD_time_manager, only : operator(-), operator(>), operator(*), operator(/)
use GOLD_initialization, only : GOLD_initialize, Get_GOLD_Input
use GOLD_initialization, only : GOLD_initialization_struct

use GOLD_barotropic, only : barotropic_init, btstep, btcalc, bt_mass_source
use GOLD_barotropic, only : register_barotropic_restarts, set_dtbt, barotropic_CS
use GOLD_continuity, only : continuity, continuity_init, continuity_CS
use GOLD_CoriolisAdv, only : CorAdCalc, CoriolisAdv_init, CoriolisAdv_CS
use GOLD_diabatic_driver, only : diabatic, diabatic_driver_init, diabatic_CS
use GOLD_diagnostics, only : calculate_diagnostic_fields, GOLD_diagnostics_init
use GOLD_diagnostics, only : diagnostics_CS
use GOLD_diag_to_Z, only : calculate_Z_diag_fields, calculate_Z_transport
use GOLD_diag_to_Z, only : GOLD_diag_to_Z_init, register_Z_tracer, diag_to_Z_CS
use GOLD_diag_to_Z, only : GOLD_diag_to_Z_end
use GOLD_EOS, only : select_eqn_of_state
use GOLD_error_checking, only : check_redundant
use GOLD_grid, only : GOLD_grid_init, ocean_grid_type, get_thickness_units
use GOLD_grid, only : get_flux_units, get_tr_flux_units
use GOLD_hor_visc, only : horizontal_viscosity, hor_visc_init, hor_visc_CS
use GOLD_lateral_mixing_coeffs, only : calc_slope_function, VarMix_init
use GOLD_lateral_mixing_coeffs, only : calc_resoln_function, VarMix_CS
use GOLD_interface_heights, only : find_eta
use GOLD_MEKE, only : MEKE_init, MEKE_alloc_register_restart, step_forward_MEKE, MEKE_CS
use GOLD_MEKE_types, only : MEKE_type
use GOLD_mixed_layer_restrat, only : mixedlayer_restrat, mixedlayer_restrat_init, mixedlayer_restrat_CS
use GOLD_open_boundary, only : Radiation_Open_Bdry_Conds, open_boundary_init
use GOLD_open_boundary, only : open_boundary_CS
use GOLD_PressureForce, only : PressureForce, PressureForce_init, PressureForce_CS
use GOLD_thickness_diffuse, only : thickness_diffuse, thickness_diffuse_init, thickness_diffuse_CS
use GOLD_tidal_forcing, only : tidal_forcing_init, tidal_forcing_CS
use GOLD_tracer, only : advect_tracer, register_tracer, add_tracer_diagnostics
use GOLD_tracer, only : add_tracer_2d_diagnostics, tracer_hordiff
use GOLD_tracer, only : advect_tracer_init, advect_tracer_diag_init, advect_tracer_CS
use GOLD_tracer_flow_control, only : call_tracer_register, tracer_flow_control_CS
use GOLD_tracer_flow_control, only : tracer_flow_control_init, call_tracer_surface_state
use GOLD_vert_friction, only : vertvisc, vertvisc_coef, vertvisc_remnant
use GOLD_vert_friction, only : vertvisc_limit_vel, vertvisc_init, vertvisc_CS
use GOLD_set_visc, only : set_viscous_BBL, set_viscous_ML, set_visc_init, set_visc_CS

implicit none ; private

#include <GOLD_memory.h>

public initialize_GOLD, step_GOLD, GOLD_end, calculate_surface_state


type, public :: GOLD_control_struct
  real PTR_, dimension(NXMEMQP_,NYMEM_,NZ_,C2_) :: &
    u         ! Zonal velocity, in m s-1.
  real PTR_, dimension(NXMEM_,NYMEMQP_,NZ_,C2_) :: &
    v         ! Meridional velocity, in m s-1.
  real PTR_, dimension(NXMEM_,NYMEM_,NZ_,C2_) :: &
    h         ! Layer thickness, in m.
  real PTR_, dimension(NXMEM_,NYMEM_,C2_) :: &
    eta       ! Instantaneous free surface height, in m.
  real PTR_, dimension(NXMEM_,NYMEM_,NZ_) :: &
    T, &      ! Potential temperature in C.
    S, &      ! Salinity in PSU.
    Rml       ! The mixed and buffer layer potential densities in kg m-3.  Rml has
              ! nkml+nkbl active layers and is only used with the bulk mixed layer.
  real PTR_, dimension(NXMEMQP_,NYMEM_,NZ_) :: &
    uh, &     ! uh = u * h * dy at u grid points in m3 s-1.
    CAu, &    ! CAu = f*v - u.grad(u) in m s-2.
    PFu, &    ! PFu = -dM/dx, in m s-2.
    diffu, &  ! Zonal acceleration due to convergence of the along-isopycnal
              ! stress tensor, in m s-2.
    visc_rem_u, & ! Both the fraction of the zonal momentum originally in a
              ! layer that remains after a time-step of viscosity, and the
              ! fraction of a time-step's worth of a barotropic acceleration
              ! that a layer experiences after viscosity is applied.
              ! Nondimensional between 0 (at the bottom) and 1 (far above).
    uhtr      ! Accumlated zonal thickness fluxes used to advect tracers, in m3.
  real PTR_, dimension(NXMEM_,NYMEMQP_,NZ_) :: &
    vh, &     ! vh = v * h * dx at v grid points in m3 s-1.
    CAv, &    ! CAv = -f*u - u.grad(v) in m s-2.
    PFv, &    ! PFv = -dM/dy, in m s-2.
    diffv, &  ! Meridional acceleration due to convergence of the
              ! along-isopycnal stress tensor, in m s-2.
    visc_rem_v, & ! Both the fraction of the meridional momentum originally in
              ! a layer that remains after a time-step of viscosity, and the
              ! fraction of a time-step's worth of a barotropic acceleration
              ! that a layer experiences after viscosity is applied.
              ! Nondimensional between 0 (at the bottom) and 1 (far above).
    vhtr      ! Accumlated meridional thickness fluxes used to advect tracers, in m3.
  real PTR_, dimension(NXMEM_,NYMEM_) :: &
    ave_ssh, &! The time-averaged sea surface height in m.
    eta_PF    ! The instantaneous SSH used in calculating PFu and PFv, in m.
  real PTR_, dimension(NXMEMQP_,NYMEM_) :: uhbt
  real PTR_, dimension(NXMEM_,NYMEMQP_) :: vhbt
    ! uhbt and vhbt are the average volume or mass fluxes determined by the
    ! barotropic solver in m3 s-1 or kg s-1.  uhbt and vhbt should (roughly?) 
    ! equal the verticals sum of uh and vh, respectively.
  real PTR_, dimension(NXMEMQP_,NYMEM_) :: uhbt_in
  real PTR_, dimension(NXMEM_,NYMEMQP_) :: vhbt_in
    ! uhbt_in and vhbt_in are the vertically summed transports from based on
    ! the final thicknessses and velocities from the previous dynamics time
    ! step, both in units of m3 s-1 or kg s-1.
! The following 6 variables are only used with the split time stepping scheme.
  real PTR_, dimension(NXMEM_,NYMEM_,NZ_) :: pbce
      ! pbce times eta gives the baroclinic pressure anomaly in each layer due
      ! to free surface height anomalies.  pbce has units of m2 H-1 s-2.
  real PTR_, dimension(NXMEMQP_,NYMEM_,NZ_) :: u_accel_bt
  real PTR_, dimension(NXMEM_,NYMEMQP_,NZ_) :: v_accel_bt
    ! u_accel_bt and v_accel_bt are layers' accelerations due to
    ! the difference between the accelerations from the barotropic calculation
    ! baroclinic accelerations that were fed into the barotropic
    ! calculation, in m s-2.
  real PTR_, dimension(NXMEMQP_,NYMEM_,NZ_) :: u_av
  real PTR_, dimension(NXMEM_,NYMEMQP_,NZ_) :: v_av
    ! u_av and v_av are the layer velocities with the vertical mean replaced by
    ! the time-mean barotropic velocity over a baroclinic timestep, in m s-1.
  real PTR_, dimension(NXMEM_,NYMEM_,NZ_)  :: h_av
    ! the difference between the accelerations from the barotropic calculation
    ! baroclinic accelerations that were fed into the barotropic
    ! calculation, in m s-2.

  type(ocean_grid_type) :: grid ! A structure containing metrics and grid info.
  type(thermo_var_ptrs) :: tv ! A structure containing pointers to an assortment
                              ! of thermodynamic fields that may be available.
  type(diag_ptrs) :: diag     ! A structure containing pointers to
                              ! diagnostic fields that might be calculated
                              ! and shared between modules.
  type(vertvisc_type) :: visc ! A structure containing vertical viscosities,
                              ! bottom drag viscosities, and related fields.
  type(BT_cont_type), pointer :: BT_cont => NULL()
                              ! A structure with elements that describe the
                              ! effective summed open face areas as a function
                              ! of barotropic flow.
  type(MEKE_type), pointer :: MEKE => NULL()  ! A structure containing fields
                             ! relateded to the Mesoscale Eddy Kinetic Energy.

  logical :: split           ! If true, use the split time stepping scheme.
  logical :: adiabatic       ! If true, there are no diapycnal mass fluxes, and
                             ! the subroutine calls to calculate and apply such
                             ! diapycnal fluxes are eliminated.
  logical :: use_temperature ! If true, temperature and salinity are used as
                             ! state variables.
  logical :: use_frazil      ! If true, water freezes if it gets too cold, and
                             ! the accumulated heat deficit is returned in the
                             ! surface state.
  logical :: bound_salinity  ! If true, salt is added to keep the salinity above
                             ! a minimum value, and the deficit is reported.
  logical :: bulkmixedlayer  ! If true, a refined bulk mixed layer is used with
                             ! nkml sublayers and nkbl buffer layer.
  logical :: thickness_diffuse ! If true, interfaces are diffused with a
                             ! coefficient of KHTH.
  logical :: thickness_diffuse_first ! If true, diffuse thickness before dynamics.
  logical :: mixedlayer_restrat ! If true, a density-gradient dependent
                             ! restratifying flow is imposed in the mixed layer.
  logical :: debug           ! If true, write verbose checksums for debugging purposes.
  logical :: debug_truncations  ! If true, make sure that all diagnostics that
                             ! could be useful for debugging any truncations are
                             ! calculated.

  real    :: dt              ! The (baroclinic) dynamics time step, in s.
  real    :: dt_therm        ! The thermodynamics time step, in s.
  real    :: be              ! A nondimensional number from 0.5 to 1 that controls
                             ! the backward weighting of the time stepping scheme.
  real    :: begw            ! A nondimensional number from 0 to 1 that controls
                             ! the extent to which the treatment of gravity waves
                             ! is forward-backward (0) or simulated backward
                             ! Euler (1).  0 is almost always used.
  type(time_type) :: Z_diag_interval  !   The amount of time between calls to
                             ! calculate Z-space diagnostics.
  type(time_type) :: Z_diag_time  ! The next time at which Z-space diagnostics
                             ! should be calculated.
  real    :: Hmix            ! The diagnostic mixed layer thickness in m when
                             ! the bulk mixed layer is not used.
  real    :: C_p             !   The heat capacity of seawater, in J K-1 kg-1.
  logical :: calc_bbl ! If true, the BBL viscosity and thickness need to be
                      ! calculated. This only applies with BOTTOMDRAGLAW true.
  real :: bbl_calc_time_interval ! The amount of time to use in diagnostics of
                                 ! the BBL properties.

  integer :: ntrunc           ! The number of times the velocity has been
                              ! truncated since the last call to write_energy.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(ocean_OBC_type), pointer :: OBC => NULL() ! A pointer to an open boundary
                             ! condition type that specifies whether, where, and
                             ! what open boundary conditions are used.  If no
                             ! open BCs are used, this pointer stays nullified.
  real :: rel_time = 0.0     ! Relative time in s since the start
                             ! of the current execution.

  ! This is to allow the previous, velocity-based coupling with between the
  ! baroclinic and barotropic modes.
  logical :: interp_p_surf     ! If true, linearly interpolate the surface
                               ! pressure over the coupling time step, using 
                               ! the specified value at the end of the coupling
                               ! step. False by default.
  real    :: smooth_ssh_passes  ! If greater than 0, apply this number of 
                               ! spatial smoothing passes to the sea surface
                               ! heights that are reported back to the calling
                               ! program by GOLD.  The default is 0.
  logical :: p_surf_prev_set   ! If true, p_surf_prev has been properly set from
                               ! a previous time-step or the ocean restart file.
                               ! This is only valid when interp_p_surf is true.
  logical :: flux_BT_coupling  ! If true, use volume fluxes, not velocities,
                               ! to couple the baroclinic and barotropic modes.
  logical :: readjust_BT_trans ! If true, readjust the barotropic transport of
                               ! the input velocities to agree with CS%uhbt_in
                               ! and CS%vhbt_in after the diabatic step.
  logical :: BT_include_udhdt  ! If true, estimate the sum of u dh/dt and v dh/dt
                               ! and provide them to the barotropic solver.
  real    :: dtbt_reset_period ! The time interval in seconds between dynamic
                               ! recalculation of the barotropic time step.  If
                               ! this is negative, it is never calculated, and
                               ! if it is 0, it is calculated every step.
  logical :: calc_dtbt         ! If true, calculate the barotropic time-step
                               ! dynamically.
  logical :: readjust_velocity ! A flag that varies with time that determines
                               ! whether the velocities currently need to be
                               ! readjusted to agree with CS%uhbt_in and
                               ! CS%vhbt_in.  This is only used if 
                               ! CS%readjust_BT_trans or BT_include_udhdt are true.
  logical :: check_bad_surface_vals ! If true, scans the surface state for
                                    ! ridiculous values
  real    :: bad_val_ssh_max   ! Maximum SSH before triggering bad value message
  real    :: bad_val_sst_max   ! Maximum SST before triggering bad value message
  real    :: bad_val_sst_min   ! Minimum SST before triggering bad value message
  real    :: bad_val_sss_max   ! Maximum SSS before triggering bad value message

  real, pointer, dimension(:,:) :: &
    p_surf_prev, &  ! The value of the surface pressure at the end of the
                    ! previous call to step_GOLD, in Pa.
    p_surf_begin, & ! The values of the surface pressure at the beginning and
    p_surf_end      ! end of a call to step_GOLD_dyn_..., in Pa.

  ! Arrays that can be used to store advective and diffusive tracer fluxes.
  real, pointer, dimension(:,:,:) :: &
    T_adx => NULL(), T_ady => NULL(), T_diffx => NULL(), T_diffy => NULL(), &
    S_adx => NULL(), S_ady => NULL(), S_diffx => NULL(), S_diffy => NULL()
  ! Arrays that can be used to store vertically integrated advective and
  ! diffusive tracer fluxes.
  real, pointer, dimension(:,:) :: &
    T_adx_2d => NULL(), T_ady_2d => NULL(), T_diffx_2d => NULL(), T_diffy_2d => NULL(), &
    S_adx_2d => NULL(), S_ady_2d => NULL(), S_diffx_2d => NULL(), S_diffy_2d => NULL(), &
    SST_sq => NULL()

! The following are the ids of various diagnostics.
  integer :: id_u = -1, id_v = -1, id_h = -1, id_uh = -1, id_vh = -1
  integer :: id_uav = -1, id_vav = -1
  integer :: id_T = -1, id_S = -1, id_Rml = -1, id_ssh = -1, id_fraz = -1
  integer :: id_salt_deficit = -1, id_Heat_PmE = -1, id_intern_heat = -1
  integer :: id_du_adj = -1, id_dv_adj = -1, id_du_adj2 = -1, id_dv_adj2 = -1
  integer :: id_h_dudt = -1, id_h_dvdt = -1
  integer :: id_sst = -1, id_sst_sq = -1, id_sss = -1, id_ssu = -1, id_ssv = -1

  integer :: id_PFu = -1, id_PFv = -1, id_CAu = -1, id_CAv = -1
  integer :: id_u_BT_accel = -1, id_v_BT_accel = -1
  integer :: id_Tadx = -1, id_Tady = -1, id_Tdiffx = -1, id_Tdiffy = -1
  integer :: id_Sadx = -1, id_Sady = -1, id_Sdiffx = -1, id_Sdiffy = -1
  integer :: id_Tadx_2d = -1, id_Tady_2d = -1, id_Tdiffx_2d = -1, id_Tdiffy_2d = -1
  integer :: id_Sadx_2d = -1, id_Sady_2d = -1, id_Sdiffx_2d = -1, id_Sdiffy_2d = -1
  integer :: id_u_predia = -1, id_v_predia = -1, id_h_predia = -1
  integer :: id_T_predia = -1, id_S_predia = -1, id_e_predia = -1
! The remainder of the structure is pointers to child subroutines' control strings.
  type(hor_visc_CS), pointer :: hor_visc_CSp => NULL()
  type(continuity_CS), pointer :: continuity_CSp => NULL()
  type(CoriolisAdv_CS), pointer :: CoriolisAdv_CSp => NULL()
  type(PressureForce_CS), pointer :: PressureForce_CSp => NULL()
  type(barotropic_CS), pointer :: barotropic_CSp => NULL()
  type(vertvisc_CS), pointer :: vertvisc_CSp => NULL()
  type(set_visc_CS), pointer :: set_visc_CSp => NULL()
  type(diabatic_CS), pointer :: diabatic_CSp => NULL()
  type(thickness_diffuse_CS), pointer :: thickness_diffuse_CSp => NULL()
  type(mixedlayer_restrat_CS), pointer :: mixedlayer_restrat_CSp => NULL()
  type(open_boundary_CS), pointer :: open_boundary_CSp => NULL()
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL()
  type(MEKE_CS),  pointer :: MEKE_CSp => NULL()
  type(VarMix_CS),  pointer :: VarMix => NULL()
  type(advect_tracer_CS), pointer :: tracer_CSp => NULL()
  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp => NULL()
  type(diagnostics_CS), pointer :: diagnostics_CSp => NULL()
  type(diag_to_Z_CS), pointer :: diag_to_Z_CSp => NULL()
  type(GOLD_restart_CS),  pointer :: restart_CSp => NULL()
end type GOLD_control_struct

integer :: id_clock_ocean, id_clock_dynamics, id_clock_thermo
integer :: id_clock_tracer, id_clock_diabatic, id_clock_pass
integer :: id_clock_Cor, id_clock_pres, id_clock_continuity, id_clock_vertvisc
integer :: id_clock_horvisc, id_clock_thick_diff, id_clock_ml_restrat
integer :: id_clock_mom_update, id_clock_diagnostics, id_clock_Z_diag
integer :: id_clock_btstep, id_clock_btcalc, id_clock_btforce
integer :: id_clock_init, id_clock_GOLD_init, id_clock_pass_init

!#######################################################################

contains

function step_GOLD(fluxes, state, Time_start, time_interval, CS)
  ! ==================================================================
  !   This subroutine time steps GOLD.
  ! ==================================================================

  integer                              :: step_GOLD
  type(forcing), intent(inout)         :: fluxes
  type(surface), intent(inout)         :: state
  type(time_type), intent(in)          :: Time_start
  real,           intent(in)           :: time_interval
  type(GOLD_control_struct), pointer   :: CS
! Arguments: fluxes - An intent in structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (out)     state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (in)      Time_start - The starting time of a run segment, as a time type.
!  (in)      time_interval - The interval of time over which to integrate in s.
!  (in)      CS - The control structure returned by a previous call to
!                 initialize_GOLD.
  type(ocean_grid_type), pointer :: grid ! A pointer to a structure containing
                                  ! metrics and related information.
  integer, save :: nt = 1 ! The running number of iterations.
  integer :: ntstep ! The number of time steps between tracer updates
                    ! or diabatic forcing.
  integer :: n_max  ! The number of steps to take in this call.
  integer :: m = 1  ! The current time level (1, 2, or 3).
  integer :: mp     ! The previous value of m.
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, n
  integer :: isd, ied, jsd, jed
  real :: dt        ! The baroclinic time step in s.
  real :: dtth      ! The time step used for thickness diffusion, in s.
  real :: dtnt      ! The elapsed time since updating the tracers and applying
                    ! diabatic processes, in s.
  real :: dtbt_reset_time ! The value of CS%rel_time when DTBT was last
                    ! calculated, in s.
  real :: wt_end, wt_beg
  real, dimension(SZI_(CS%grid),SZJ_(CS%grid)) :: &
    eta_av, &       ! The average sea surface height or column mass over
                    ! a time step, in m or kg m-2.
    ssh             ! The sea surface height based on eta_av, in m.
  real, allocatable, dimension(:,:) :: &
    frazil_ave, &   ! The average heat frazil heat flux
                    ! required to keep the temperature above freezing, in W m-2.
    salt_deficit_ave, &  ! The average salt flux required to keep the
                    ! salinity above 0.01 PSU, in gSalt m-2 s-1.
    Heat_PmE_ave, & !   The average effective heat flux into the ocean due to
                    ! the exchange of water with other components, times the
                    ! heat capacity of water, in W m-2.   
    intern_heat_ave !   The average heat flux into the ocean from geothermal or
                    ! other internal heat sources, in W m-2.   
  real, pointer, dimension(:,:,:,:) :: &
    u, &                     ! u : Zonal velocity, in m s-1.
    v, &                     ! v : Meridional velocity, in m s-1.
    h                        ! h : Layer thickness, in m.
  real, dimension(SZI_(CS%grid),SZJ_(CS%grid),SZK_(CS%grid)+1) :: eta_predia
  real :: tot_wt_ssh, Itot_wt_ssh, I_time_int
  type(time_type) :: Time_local
  integer :: pid_tau, pid_ustar, pid_psurf, pid_u, pid_h
  integer :: pid_Rml, pid_T, pid_S

  grid => CS%grid
  is = grid%isc ; ie = grid%iec ; js = grid%jsc ; je = grid%jec ; nz = grid%ke
  Isq = grid%Iscq ; Ieq = grid%Iecq ; Jsq = grid%Jscq ; Jeq = grid%Jecq
  isd = grid%isd ; ied = grid%ied ; jsd = grid%jsd ; jed = grid%jed
  u => CS%u ; v => CS%v ; h => CS%h

  call cpu_clock_begin(id_clock_ocean)
 !   First determine the time step that is consistent with this call.
 ! It is anticipated that the time step will almost always coincide
 ! with dt.  In addition, ntstep is determined, subject to the constraint
 ! that ntstep cannot exceed n_max.
  if (time_interval <= CS%dt) then
    n_max = 1
  else
    n_max = ceiling(time_interval/CS%dt - 0.001)
  endif

  dt = real(time_interval) / real(n_max)
  dtnt = 0.0
  ntstep = MAX(1,MIN(n_max,floor(CS%dt_therm/dt + 0.001)))

  CS%calc_bbl = .true.
  if (.not.ASSOCIATED(fluxes%p_surf)) CS%interp_p_surf = .false.

  call cpu_clock_begin(id_clock_pass)
  if (grid%nonblocking_updates) then
    pid_tau = pass_vector_start(fluxes%taux, fluxes%tauy, grid%Domain)
    if (ASSOCIATED(fluxes%ustar)) &
      pid_ustar = pass_var_start(fluxes%ustar(:,:), grid%Domain)
    if (ASSOCIATED(fluxes%p_surf)) &
      pid_psurf = pass_var_start(fluxes%p_surf(:,:), grid%Domain)
  else
    call pass_vector(fluxes%taux, fluxes%tauy, grid%Domain)
    if (ASSOCIATED(fluxes%ustar)) call pass_var(fluxes%ustar(:,:), grid%Domain)
    if (ASSOCIATED(fluxes%p_surf)) call pass_var(fluxes%p_surf(:,:), grid%Domain)
  endif
  call cpu_clock_end(id_clock_pass)
  if (ASSOCIATED(CS%tv%frazil)) CS%tv%frazil(:,:) = 0.0
  if (ASSOCIATED(CS%tv%salt_deficit)) CS%tv%salt_deficit(:,:) = 0.0   
  if (ASSOCIATED(CS%tv%TempxPmE)) CS%tv%TempxPmE(:,:) = 0.0
  if (ASSOCIATED(CS%tv%internal_heat)) CS%tv%internal_heat(:,:) = 0.0

  CS%rel_time = 0.0

  tot_wt_ssh = 0.0
  do j=js,je ; do i=is,ie ; CS%ave_ssh(i,j) = 0.0 ; enddo ; enddo

  if (CS%interp_p_surf) then
    if (.not.ASSOCIATED(CS%p_surf_end)) allocate(CS%p_surf_end(isd:ied,jsd:jed))
    if (.not.ASSOCIATED(CS%p_surf_begin)) allocate(CS%p_surf_begin(isd:ied,jsd:jed))
    
    if (.not.CS%p_surf_prev_set) then
      do j=jsd,jed ; do i=isd,ied
        CS%p_surf_prev(i,j) = fluxes%p_surf(i,j)
      enddo ; enddo
      CS%p_surf_prev_set = .true.
    endif
  else
    CS%p_surf_end  => fluxes%p_surf
  endif

  mp = 1 ; if (CS%split) mp = MOD(nt+1,2) + 1
  if (CS%debug) then
    call GOLD_state_chksum("Before steps ", u(:,:,:,mp), v(:,:,:,mp), &
                          h(:,:,:,mp), CS%uh, CS%vh, grid)
    call check_redundant("Before steps mp ", u(:,:,:,mp), v(:,:,:,mp), grid)
  endif

  if (associated(CS%VarMix)) then
    call enable_averaging(time_interval, Time_start+set_time(int(time_interval)), &
                          CS%diag)
    call calc_resoln_function(h(:,:,:,mp), CS%tv, grid, CS%VarMix)
    call disable_averaging(CS%diag)
  endif

  if (grid%nonblocking_updates) then
    call cpu_clock_begin(id_clock_pass)
    call pass_vector_complete(pid_tau, fluxes%taux, fluxes%tauy, grid%Domain)
    if (ASSOCIATED(fluxes%ustar)) &
      call pass_var_complete(pid_ustar, fluxes%ustar(:,:), grid%Domain)
    if (ASSOCIATED(fluxes%p_surf)) &
      call pass_var_complete(pid_psurf, fluxes%p_surf(:,:), grid%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  do n=1,n_max
    nt = nt + 1
    ! Set the universally visible time to the middle of the time step
    CS%Time = Time_start + set_time(int(floor(CS%rel_time+0.5*dt+0.5)))
    CS%rel_time = CS%rel_time + dt
    ! Set the local time to the end of the time step.
    Time_local = Time_start + set_time(int(floor(CS%rel_time+0.5)))

    call cpu_clock_begin(id_clock_dynamics)
    call disable_averaging(CS%diag)

    if (CS%thickness_diffuse .and. CS%thickness_diffuse_first) then
      if (MOD(n-1,ntstep) == 0) then
        dtth = dt*min(ntstep,n_max-n+1)
        call enable_averaging(dtth,Time_local+set_time(int(floor(dtth-dt+0.5))), CS%diag)
        mp = 1 ; if (CS%split) mp = MOD(nt,2) + 1
        call cpu_clock_begin(id_clock_thick_diff)
        if (associated(CS%VarMix)) &
          call calc_slope_function(h(:,:,:,mp), CS%tv, grid, CS%VarMix)
        call thickness_diffuse(h(:,:,:,mp), CS%uhtr, CS%vhtr, CS%tv, dtth, grid, &
                               CS%MEKE, CS%VarMix, CS%thickness_diffuse_CSp)
        call cpu_clock_end(id_clock_thick_diff)
        call cpu_clock_begin(id_clock_pass)
        call pass_var(h(:,:,:,mp), grid%Domain)
        call cpu_clock_end(id_clock_pass)
        call disable_averaging(CS%diag)
      endif
    endif

    if (CS%interp_p_surf) then
      wt_end = real(n) / real(n_max)
      wt_beg = real(n-1) / real(n_max)
      do j=jsd,jed ; do i=isd,ied
        CS%p_surf_end(i,j) = wt_end * fluxes%p_surf(i,j) + &
                        (1.0-wt_end) * CS%p_surf_prev(i,j)
        CS%p_surf_begin(i,j) = wt_beg * fluxes%p_surf(i,j) + &
                        (1.0-wt_beg) * CS%p_surf_prev(i,j)
      enddo ; enddo
    endif

    if (CS%calc_bbl) &
      CS%bbl_calc_time_interval = dt*real(1+MIN(ntstep-MOD(n,ntstep),n_max-n))

    if (CS%split) then !--------------------------- start SPLIT
!   This section uses a predictor corrector scheme, that is somewhere
! (determined by be) between the forward-backward (be=0.5) scheme and
! the backward Euler scheme (be=1.0) to time step the dynamic equations.
      mp = MOD(nt,2) + 1
      m  = 3 - mp

      CS%calc_dtbt = .false.
      if ((CS%dtbt_reset_period >= 0.0) .and. &
          ((n==1) .or. (CS%dtbt_reset_period == 0.0) .or. &
           (CS%rel_time >= dtbt_reset_time + 0.999*CS%dtbt_reset_period))) then
        CS%calc_dtbt = .true.
        dtbt_reset_time = CS%rel_time
      endif

      call step_GOLD_dyn_split(u(:,:,:,mp), v(:,:,:,mp), h(:,:,:,mp), &
                    CS%eta(:,:,mp), CS%uhbt_in, CS%vhbt_in, Time_local, dt, &
                    fluxes, CS%p_surf_begin, CS%p_surf_end, dtnt, dt*ntstep, &
                    CS%uh, CS%vh, CS%u_av, CS%v_av, CS%h_av, eta_av, &
                    u(:,:,:,m), v(:,:,:,m), h(:,:,:,m), CS%eta(:,:,m), grid, CS)

    else ! --------------------------------------------------- not SPLIT

      call step_GOLD_dyn_unsplit_RK3(u(:,:,:,1), v(:,:,:,1), h(:,:,:,1), &
                    Time_local, dt, fluxes, CS%p_surf_begin, CS%p_surf_end, &
                    CS%uh, CS%vh, eta_av, grid, CS)

    endif ! -------------------------------------------------- end SPLIT

    call disable_averaging(CS%diag)
    call cpu_clock_end(id_clock_dynamics)

    dtnt = dtnt + dt
    if ((MOD(n,ntstep) == 0) .or. (n==n_max)) then
      if (CS%debug) then
        call uchksum(u(:,:,:,m),"Pre-advection u",grid,haloshift=2)
        call vchksum(v(:,:,:,m),"Pre-advection v",grid,haloshift=2)
        call hchksum(h(:,:,:,m),"Pre-advection h",grid,haloshift=1)
        call uchksum(CS%uhtr,"Pre-advection uh",grid,haloshift=0)
        call vchksum(CS%vhtr,"Pre-advection vh",grid,haloshift=0)
      ! call GOLD_state_chksum("Pre-advection ", u(:,:,:,m), v(:,:,:,m), &
      !                       h(:,:,:,m), CS%uhtr, CS%vhtr, grid, haloshift=1)
          if (associated(CS%tv%T)) call hchksum(CS%tv%T, "Pre-advection T",grid,haloshift=1)
          if (associated(CS%tv%S)) call hchksum(CS%tv%S, "Pre-advection S",grid,haloshift=1)
          if (associated(CS%tv%Rml)) call hchksum(CS%tv%Rml, "Pre-advection Rml",grid,haloshift=0)
          if (associated(CS%tv%frazil)) call hchksum(CS%tv%frazil, "Pre-advection frazil",grid,haloshift=0)
          if (associated(CS%tv%salt_deficit)) call hchksum(CS%tv%salt_deficit, "Pre-advection salt deficit",grid,haloshift=0)
      ! call GOLD_thermo_chksum("Pre-advection ", CS%tv, grid)
        call check_redundant("Pre-advection ", u(:,:,:,m), v(:,:,:,m), grid)
      endif

      call cpu_clock_begin(id_clock_thermo)
      call enable_averaging(dtnt,Time_local, CS%diag)

      call cpu_clock_begin(id_clock_tracer)
      call advect_tracer(h(:,:,:,m), CS%uhtr, CS%vhtr, CS%OBC, dtnt, grid, &
                         CS%tracer_CSp)
      call tracer_hordiff(h(:,:,:,m), dtnt, CS%MEKE, CS%VarMix, grid, CS%tracer_CSp, &
                          CS%tv)
      call cpu_clock_end(id_clock_tracer)

      call cpu_clock_begin(id_clock_Z_diag)
      call calculate_Z_transport(CS%uhtr, CS%vhtr, h(:,:,:,m), dtnt, grid, &
                                 CS%diag_to_Z_CSp)
      call cpu_clock_end(id_clock_Z_diag)

      if (CS%id_u_predia > 0) call post_data(CS%id_u_predia, u(:,:,:,m), CS%diag)
      if (CS%id_v_predia > 0) call post_data(CS%id_v_predia, v(:,:,:,m), CS%diag)
      if (CS%id_h_predia > 0) call post_data(CS%id_h_predia, h(:,:,:,m), CS%diag)
      if (CS%id_T_predia > 0) call post_data(CS%id_T_predia, CS%tv%T, CS%diag)
      if (CS%id_S_predia > 0) call post_data(CS%id_S_predia, CS%tv%S, CS%diag)
      if (CS%id_e_predia > 0) then
        call find_eta(h(:,:,:,m), CS%tv, grid%g_Earth, grid, eta_predia)
        call post_data(CS%id_e_predia, eta_predia, CS%diag)
      endif

      if (.not.CS%adiabatic) then
        if (CS%debug) then
          call uchksum(u(:,:,:,m),"Pre-diabatic u",grid,haloshift=2)
          call vchksum(v(:,:,:,m),"Pre-diabatic v",grid,haloshift=2)
          call hchksum(h(:,:,:,m),"Pre-diabatic h",grid,haloshift=1)
          call uchksum(CS%uhtr,"Pre-diabatic uh",grid,haloshift=0)
          call vchksum(CS%vhtr,"Pre-diabatic vh",grid,haloshift=0)
        ! call GOLD_state_chksum("Pre-diabatic ", u(:,:,:,m), v(:,:,:,m), &
        !                       h(:,:,:,m), CS%uhtr, CS%vhtr, grid)
          call GOLD_thermo_chksum("Pre-diabatic ", CS%tv, grid,haloshift=0)
          call check_redundant("Pre-diabatic ", u(:,:,:,m), v(:,:,:,m), grid)
        endif

        if (CS%readjust_BT_trans .and. .not.CS%BT_include_udhdt) then
          call find_total_transport(u(:,:,:,m), v(:,:,:,m), h(:,:,:,m), &
                                    CS%uhbt_in, CS%vhbt_in, dt, grid, CS)
          CS%readjust_velocity = .true.
        endif

        call cpu_clock_begin(id_clock_diabatic)
        call diabatic(u(:,:,:,m),v(:,:,:,m),h(:,:,:,m),CS%tv,fluxes,CS%visc,dtnt, &
                      grid, CS%diabatic_CSp)
        call cpu_clock_end(id_clock_diabatic)

        call cpu_clock_begin(id_clock_pass)
        if (grid%nonblocking_updates) then        
          pid_u = pass_vector_start(u(:,:,:,m), v(:,:,:,m), grid%Domain)
          if (CS%bulkmixedlayer .and. (.not.associated(CS%tv%eqn_of_state))) &
            pid_Rml = pass_var_start(CS%tv%Rml, grid%Domain)
          if (CS%use_temperature) then
            pid_T = pass_var_start(CS%tv%T, grid%Domain)
            pid_S = pass_var_start(CS%tv%S, grid%Domain)
          endif
          pid_h = pass_var_start(h(:,:,:,m), grid%Domain)

          call pass_vector_complete(pid_u, u(:,:,:,m), v(:,:,:,m), grid%Domain)
          if (CS%bulkmixedlayer .and. (.not.associated(CS%tv%eqn_of_state))) &
            call pass_var_complete(pid_Rml, CS%tv%Rml, grid%Domain)
          if (CS%use_temperature) then
            call pass_var_complete(pid_T, CS%tv%T, grid%Domain)
            call pass_var_complete(pid_S, CS%tv%S, grid%Domain)
          endif
          call pass_var_complete(pid_h, h(:,:,:,m), grid%Domain)
        else 
          call pass_vector(u(:,:,:,m), v(:,:,:,m), grid%Domain)
          if (CS%bulkmixedlayer .and. (.not.associated(CS%tv%eqn_of_state))) &
            call pass_var(CS%tv%Rml, grid%Domain, complete=.false.)
          if (CS%use_temperature) then
            call pass_var(CS%tv%T, grid%Domain, complete=.false.)
            call pass_var(CS%tv%S, grid%Domain, complete=.false.)
          endif
          call pass_var(h(:,:,:,m), grid%Domain)
        endif
        call cpu_clock_end(id_clock_pass)

        if (CS%debug) then
          call uchksum(u(:,:,:,m),"Post-diabatic u",grid,haloshift=2)
          call vchksum(v(:,:,:,m),"Post-diabatic v",grid,haloshift=2)
          call hchksum(h(:,:,:,m),"Post-diabatic h",grid,haloshift=1)
          call uchksum(CS%uhtr,"Post-diabatic uh",grid,haloshift=0)
          call vchksum(CS%vhtr,"Post-diabatic vh",grid,haloshift=0)
        ! call GOLD_state_chksum("Post-diabatic ", u(:,:,:,m), v(:,:,:,m), &
        !                       h(:,:,:,m), CS%uhtr, CS%vhtr, grid, haloshift=1)
          if (associated(CS%tv%T)) call hchksum(CS%tv%T, "Post-diabatic T",grid,haloshift=1)
          if (associated(CS%tv%S)) call hchksum(CS%tv%S, "Post-diabatic S",grid,haloshift=1)
          if (associated(CS%tv%Rml)) call hchksum(CS%tv%Rml, "Post-diabatic Rml",grid,haloshift=0)
          if (associated(CS%tv%frazil)) call hchksum(CS%tv%frazil, "Post-diabatic frazil",grid,haloshift=0)
          if (associated(CS%tv%salt_deficit)) call hchksum(CS%tv%salt_deficit, "Post-diabatic salt deficit",grid,haloshift=0)
        ! call GOLD_thermo_chksum("Post-diabatic ", CS%tv, grid)
          call check_redundant("Post-diabatic ", u(:,:,:,m), v(:,:,:,m), grid)
        endif
      endif                                                  ! ADIABATIC
      CS%uhtr(:,:,:) = 0.0
      CS%vhtr(:,:,:) = 0.0
      call cpu_clock_end(id_clock_thermo)

      call cpu_clock_begin(id_clock_diagnostics)
      if (CS%split) then
        call calculate_diagnostic_fields(u(:,:,:,m),v(:,:,:,m),h(:,:,:,m), &
                 CS%uh, CS%vh, m, CS%tv, dtnt, grid, CS%diagnostics_CSp, &
                 CS%eta(:,:,m)) ! Could be eta? RWH
      else
        call calculate_diagnostic_fields(u(:,:,:,m),v(:,:,:,m),h(:,:,:,m), &
                 CS%uh, CS%vh, m, CS%tv, dtnt, grid, CS%diagnostics_CSp)
      endif
      call cpu_clock_end(id_clock_diagnostics)

      if (CS%id_T > 0) call post_data(CS%id_T, CS%tv%T, CS%diag)
      if (CS%id_S > 0) call post_data(CS%id_S, CS%tv%S, CS%diag)
      if (CS%id_Rml > 0) call post_data(CS%id_Rml, CS%tv%Rml, CS%diag)

      if (CS%id_Tadx > 0) call post_data(CS%id_Tadx, CS%T_adx, CS%diag)
      if (CS%id_Tady > 0) call post_data(CS%id_Tady, CS%T_ady, CS%diag)
      if (CS%id_Tdiffx > 0) call post_data(CS%id_Tdiffx, CS%T_diffx, CS%diag)
      if (CS%id_Tdiffy > 0) call post_data(CS%id_Tdiffy, CS%T_diffy, CS%diag)

      if (CS%id_Sadx > 0) call post_data(CS%id_Sadx, CS%S_adx, CS%diag)
      if (CS%id_Sady > 0) call post_data(CS%id_Sady, CS%S_ady, CS%diag)
      if (CS%id_Sdiffx > 0) call post_data(CS%id_Sdiffx, CS%S_diffx, CS%diag)
      if (CS%id_Sdiffy > 0) call post_data(CS%id_Sdiffy, CS%S_diffy, CS%diag)

      if (CS%id_Tadx_2d > 0) call post_data(CS%id_Tadx_2d, CS%T_adx_2d, CS%diag)
      if (CS%id_Tady_2d > 0) call post_data(CS%id_Tady_2d, CS%T_ady_2d, CS%diag)
      if (CS%id_Tdiffx_2d > 0) call post_data(CS%id_Tdiffx_2d, CS%T_diffx_2d, CS%diag)
      if (CS%id_Tdiffy_2d > 0) call post_data(CS%id_Tdiffy_2d, CS%T_diffy_2d, CS%diag)

      if (CS%id_Sadx_2d > 0) call post_data(CS%id_Sadx_2d, CS%S_adx_2d, CS%diag)
      if (CS%id_Sady_2d > 0) call post_data(CS%id_Sady_2d, CS%S_ady_2d, CS%diag)
      if (CS%id_Sdiffx_2d > 0) call post_data(CS%id_Sdiffx_2d, CS%S_diffx_2d, CS%diag)
      if (CS%id_Sdiffy_2d > 0) call post_data(CS%id_Sdiffy_2d, CS%S_diffy_2d, CS%diag)

      call disable_averaging(CS%diag)

      call cpu_clock_begin(id_clock_Z_diag)

      if (Time_local + set_time(int(0.5*CS%dt_therm)) > CS%Z_diag_time) then
        call enable_averaging(time_type_to_real(CS%Z_diag_interval), &
                              CS%Z_diag_time, CS%diag)
        call calculate_Z_diag_fields(u(:,:,:,m),v(:,:,:,m),h(:,:,:,m), dtnt, &
                                     grid, CS%diag_to_Z_CSp)
        CS%Z_diag_time = CS%Z_diag_time + CS%Z_diag_interval
        call disable_averaging(CS%diag)
      endif
      call cpu_clock_end(id_clock_Z_diag)

      dtnt = 0.0
      CS%calc_bbl = .true.
    else  ! It is not time to do thermodynamics.
      if (.not.CS%BT_include_udhdt) CS%readjust_velocity = .false.
    endif

    call enable_averaging(dt,Time_local, CS%diag)
    if (CS%id_u > 0) call post_data(CS%id_u, u(:,:,:,m), CS%diag)
    if (CS%id_v > 0) call post_data(CS%id_v, v(:,:,:,m), CS%diag)
    if (CS%id_h > 0) call post_data(CS%id_h, h(:,:,:,m), CS%diag)
    tot_wt_ssh = tot_wt_ssh + dt

    call find_eta(h(:,:,:,1), CS%tv, grid%g_Earth, grid, ssh, eta_av)
    do j=js,je ; do i=is,ie
      CS%ave_ssh(i,j) = CS%ave_ssh(i,j) + dt*ssh(i,j)
    enddo ; enddo
    call disable_averaging(CS%diag)

  enddo ! End of n loop

  Itot_wt_ssh = 1.0/tot_wt_ssh
  do j=js,je ; do i=is,ie
    CS%ave_ssh(i,j) = CS%ave_ssh(i,j)*Itot_wt_ssh
  enddo ; enddo

  call enable_averaging(dt*n_max,Time_local, CS%diag)
    I_time_int = 1.0/(dt*n_max)
    if (CS%id_ssh > 0) call post_data(CS%id_ssh, CS%ave_ssh, CS%diag)
    if (ASSOCIATED(CS%tv%frazil) .and. (CS%id_fraz > 0)) then
      allocate(frazil_ave(grid%isd:grid%ied,grid%jsd:grid%jed))
      do j=js,je ; do i=is,ie
        frazil_ave(i,j) = CS%tv%frazil(i,j) * I_time_int
      enddo ; enddo
      call post_data(CS%id_fraz, frazil_ave, CS%diag)
      deallocate(frazil_ave)
    endif
    if (ASSOCIATED(CS%tv%salt_deficit) .and. (CS%id_salt_deficit > 0)) then
      allocate(salt_deficit_ave(grid%isd:grid%ied,grid%jsd:grid%jed))
      do j=js,je ; do i=is,ie
        salt_deficit_ave(i,j) = CS%tv%salt_deficit(i,j) * I_time_int
      enddo ; enddo
      call post_data(CS%id_salt_deficit, salt_deficit_ave, CS%diag)
      deallocate(salt_deficit_ave)
    endif
    if (ASSOCIATED(CS%tv%TempxPmE) .and. (CS%id_Heat_PmE > 0)) then
      allocate(Heat_PmE_ave(grid%isd:grid%ied,grid%jsd:grid%jed))
      do j=js,je ; do i=is,ie
        Heat_PmE_ave(i,j) = CS%tv%TempxPmE(i,j) * (CS%C_p * I_time_int)
      enddo ; enddo
      call post_data(CS%id_Heat_PmE, Heat_PmE_ave, CS%diag)      
      deallocate(Heat_PmE_ave)
    endif
    if (ASSOCIATED(CS%tv%internal_heat) .and. (CS%id_intern_heat > 0)) then
      allocate(intern_heat_ave(grid%isd:grid%ied,grid%jsd:grid%jed))
      do j=js,je ; do i=is,ie
        intern_heat_ave(i,j) = CS%tv%internal_heat(i,j) * (CS%C_p * I_time_int)
      enddo ; enddo
      call post_data(CS%id_intern_heat, intern_heat_ave, CS%diag)      
      deallocate(intern_heat_ave)
    endif
  call disable_averaging(CS%diag)

  call calculate_surface_state(state, u(:,:,:,m), v(:,:,:,m), h(:,:,:,m), &
                               CS%ave_ssh, grid, CS, fluxes%p_surf_full)

  call enable_averaging(dt*n_max,Time_local, CS%diag)
    if (CS%id_sst > 0) call post_data(CS%id_sst, state%SST, CS%diag)
    if (CS%id_sst_sq > 0) then
      do j=js,je ; do i=is,ie
        CS%SST_sq(i,j) = state%SST(i,j)*state%SST(i,j)
      enddo ; enddo
      call post_data(CS%id_sst_sq, CS%SST_sq, CS%diag)
    endif
    if (CS%id_sss > 0) call post_data(CS%id_sss, state%SSS, CS%diag)
    if (CS%id_ssu > 0) call post_data(CS%id_ssu, state%u, CS%diag)
    if (CS%id_ssv > 0) call post_data(CS%id_ssv, state%v, CS%diag)
  call disable_averaging(CS%diag)

  if (CS%interp_p_surf) then ; do j=jsd,jed ; do i=isd,ied
    CS%p_surf_prev(i,j) = fluxes%p_surf(i,j)
  enddo ; enddo ; endif

  call cpu_clock_end(id_clock_ocean)
  step_GOLD = m

end function step_GOLD

subroutine step_GOLD_dyn_split(u_in, v_in, h_in, eta_in, uhbt_in, vhbt_in, &
                 Time_local, dt, fluxes, p_surf_begin, p_surf_end, &
                 dt_since_flux, dt_therm, uh, vh, u_av, v_av, h_av, eta_av, &
                 u_out, v_out, h_out, eta_out, G, CS)
  real, dimension(NXMEMQ_,NYMEM_,NZ_), target, intent(in) :: u_in
  real, dimension(NXMEM_,NYMEMQ_,NZ_), target, intent(in) :: v_in
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(in)    :: h_in
  real, dimension(NXMEM_,NYMEM_),      intent(inout) :: eta_in
  real, dimension(NXMEMQ_,NYMEM_),     intent(inout) :: uhbt_in
  real, dimension(NXMEM_,NYMEMQ_),     intent(inout) :: vhbt_in
  type(time_type),                     intent(in)    :: Time_local
  real,                                intent(in)    :: dt
  type(forcing),                       intent(in)    :: fluxes
  real, dimension(:,:),                pointer       :: p_surf_begin, p_surf_end
  real,                                intent(in)    :: dt_since_flux, dt_therm
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(inout) :: uh
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(inout) :: vh
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(inout) :: u_av
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(inout) :: v_av
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(inout) :: h_av
  real, dimension(NXMEM_,NYMEM_),      intent(out)   :: eta_av
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(out)   :: u_out
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(out)   :: v_out
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(out)   :: h_out
  real, dimension(NXMEM_,NYMEM_),      intent(out)   :: eta_out
  type(ocean_grid_type),               intent(inout) :: G
  type(GOLD_control_struct),           pointer       :: CS
! Arguments: u_in - The input zonal velocity, in m s-1.
!  (in)      v_in - The input meridional velocity, in m s-1.
!  (in)      h_in - The input layer thicknesses, in m or kg m-2, depending on
!                   whether the Boussinesq approximation is made.
!  (in)      eta_in - The input free surface height or column mass, in m or
!                     kg m-2.  (Intent inout for halo updates.)
!  (in)      Time_local - The model time at the end of the time step.
!  (in)      dt - The time step in s.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      p_surf_begin - A pointer (perhaps NULL) to the surface pressure
!                     at the beginning of this dynamic step, in Pa.
!  (in)      p_surf_end - A pointer (perhaps NULL) to the surface pressure
!                     at the end of this dynamic step, in Pa.
!  (in)      dt_since_flux - The elapsed time since fluxes were applied, in s.
!  (in)      dt_therm - The thermodynamic time step, in s.
!  (inout)   uh - The zonal volume or mass transport, in m3 s-1 or kg s-1.
!  (inout)   vh - The meridional volume or mass transport, in m3 s-1 or kg s-1.
!  (inout)   u_av - The zonal velocity time-averaged over a time step, in m s-1.
!  (inout)   v_av - The meridional velocity time-averaged over a time step, in m s-1.
!  (inout)   h_av - The layer thickness time-averaged over a time step, in m or
!                   kg m-2.
!  (out)     eta_av - The free surface height or column mass time-averaged
!                     over a time step, in m or kg m-2.
!  (out)     u_out - The output zonal velocity, in m s-1.
!  (out)     v_out - The output meridional velocity, in m s-1.
!  (out)     h_out - The output layer thicknesses, in m or kg m-2.
!  (out)     eta_out - The output free surface height or column mass, in m or
!                      kg m-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure set up by initialize_GOLD.

  real :: dt_pred   ! The time step for the predictor part of the baroclinic
                    ! time stepping.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_tmp
    ! A temporary estimated thickness, in m.
  real, dimension(SZIQ_(G),SZJ_(G),SZK_(G)) :: uh_tmp
  real, dimension(SZI_(G),SZJQ_(G),SZK_(G)) :: vh_tmp
    ! Temporary transports, in m3 s-1 or kg s-1.
  real, dimension(SZIQ_(G),SZJ_(G)) :: u_dhdt
  real, dimension(SZI_(G),SZJQ_(G)) :: v_dhdt
    ! Vertically summed transport tendencies in m3 s-2 or kg s-2.

  real, dimension(SZIQ_(G),SZJ_(G),SZK_(G)) :: u_bc_accel
  real, dimension(SZI_(G),SZJQ_(G),SZK_(G)) :: v_bc_accel
    ! u_bc_accel and v_bc_accel are the summed baroclinic accelerations of each
    ! layer calculated by the non-barotropic part of the model, both in m s-2.
  real, dimension(SZIQ_(G),SZJ_(G),SZK_(G)) :: uh_in
  real, dimension(SZI_(G),SZJQ_(G),SZK_(G)) :: vh_in
    ! uh_in and vh_in are the zonal or meridional mass transports that would be
    ! obtained using the velocities u_in and v_in, both in m3 s-1 or kg s-1.
  real, dimension(SZIQ_(G),SZJ_(G)) :: uhbt_out
  real, dimension(SZI_(G),SZJQ_(G)) :: vhbt_out
    ! uhbt_out and vhbt_out are the vertically summed transports from the
    ! barotropic solver based on its final velocities, both in m3 s-1 or kg s-1.
  real, dimension(SZIQ_(G),SZJ_(G),SZK_(G)), target :: u_adj
  real, dimension(SZI_(G),SZJQ_(G),SZK_(G)), target :: v_adj
    ! u_adj and v_adj are the zonal or meridional velocities after u_in and v_in
    ! have been barotropically adjusted so the resulting transports match
    ! uhbt_out and vhbt_out, both in m s-1.
  real :: Pa_to_eta ! A factor that converts pressures to the units of eta.
  real, pointer, dimension(:,:,:) :: u_init, v_init  ! Pointers to u_in and v_in
                                                     ! or u_adj and v_adj.
  real, pointer, dimension(:,:)   :: p_surf => NULL(), eta_PF_start => NULL()
  real :: Idt
  logical :: dyn_p_surf
  integer :: pid_Ray, pid_bbl_h, pid_kv_bbl, pid_eta_PF, pid_eta_in, pid_visc
  integer :: pid_h, pid_u, pid_u_av, pid_uh, pid_uhbt_in
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  Idt = 1.0 / dt

  if (CS%debug) then
    call GOLD_state_chksum("Start predictor ", u_in, v_in, h_in, uh, vh, G)
    call check_redundant("Start predictor u_in ", u_in, v_in, G)
    call check_redundant("Start predictor uh ", uh, vh, G)
  endif

  dyn_p_surf = CS%interp_p_surf .and. associated(p_surf_begin) .and. &
               associated(p_surf_end)
  if (dyn_p_surf) then
    p_surf => p_surf_end
    call safe_alloc_ptr(eta_PF_start,G%isd,G%ied,G%jsd,G%jed)
    eta_PF_start(:,:) = 0.0
  else
    p_surf => fluxes%p_surf
  endif

  if (CS%calc_bbl) then
    ! Calculate the BBL properties and store them inside CS%visc (u_in,h_in).
    call cpu_clock_begin(id_clock_vertvisc)
    call enable_averaging(CS%bbl_calc_time_interval, &
                          Time_local-set_time(int(dt)), CS%diag)
    call set_viscous_BBL(u_in, v_in, h_in, CS%tv, CS%visc, G, CS%set_visc_CSp)
    call disable_averaging(CS%diag)
    call cpu_clock_end(id_clock_vertvisc)

    call cpu_clock_begin(id_clock_pass)
    if (G%nonblocking_updates) then   
      if (associated(CS%visc%Ray_u) .and. associated(CS%visc%Ray_v)) &
        pid_Ray = pass_vector_start(CS%visc%Ray_u, CS%visc%Ray_v, G%Domain, &
                       To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(CS%visc%bbl_thick_u) .and. associated(CS%visc%bbl_thick_v)) &
        pid_bbl_h = pass_vector_start(CS%visc%bbl_thick_u, CS%visc%bbl_thick_v, &
                       G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(CS%visc%kv_bbl_u) .and. associated(CS%visc%kv_bbl_v)) &
        pid_kv_bbl = pass_vector_start(CS%visc%kv_bbl_u, CS%visc%kv_bbl_v, &
                       G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
      ! CS%calc_bbl will be set to .false. when the message passing is complete.
    else
      if (associated(CS%visc%Ray_u) .and. associated(CS%visc%Ray_v)) &
        call pass_vector(CS%visc%Ray_u, CS%visc%Ray_v, G%Domain, &
                       To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(CS%visc%kv_bbl_u) .and. associated(CS%visc%kv_bbl_v)) then
        call pass_vector(CS%visc%bbl_thick_u, CS%visc%bbl_thick_v, G%Domain, &
                       To_All+SCALAR_PAIR, CGRID_NE, complete=.false.)
        call pass_vector(CS%visc%kv_bbl_u, CS%visc%kv_bbl_v, G%Domain, &
                       To_All+SCALAR_PAIR, CGRID_NE)
      endif
      CS%calc_bbl = .false.
    endif
    call cpu_clock_end(id_clock_pass)
  endif

! PFu = d/dx M(h_in,T,S)
! pbce = dM/deta
  if (CS%begw == 0.0) call enable_averaging(dt, Time_local, CS%diag)
  call cpu_clock_begin(id_clock_pres)
  call PressureForce(h_in, CS%tv, CS%PFu, CS%PFv, G, CS%PressureForce_CSp, &
                     p_surf, CS%pbce, CS%eta_PF)
  if (dyn_p_surf) then
    if (G%Boussinesq) then
      Pa_to_eta = 1.0 / (G%Rho0*G%g_Earth)
    else
      Pa_to_eta = 1.0 / G%H_to_Pa
    endif
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_PF_start(i,j) = CS%eta_PF(i,j) - Pa_to_eta * &
                          (p_surf_begin(i,j) - p_surf_end(i,j))
    enddo ; enddo
  endif
  call cpu_clock_end(id_clock_pres)
  call disable_averaging(CS%diag)

  if (G%nonblocking_updates) then
    call cpu_clock_begin(id_clock_pass)
    pid_eta_PF = pass_var_start(CS%eta_PF(:,:), G%Domain)
    pid_eta_in = pass_var_start(eta_in(:,:), G%Domain)
    if (CS%readjust_velocity) &
      pid_uhbt_in = pass_vector_start(uhbt_in, vhbt_in, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

! CAu = -(f+zeta_av)/h_av vh + d/dx KE_av
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(u_av, v_av, h_av, uh, vh, CS%CAu, CS%CAv, G,CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)

! u_bc_accel = CAu + PFu + diffu(u[n-1])
  call cpu_clock_begin(id_clock_btforce)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u_bc_accel(I,j,k) = (CS%Cau(I,j,k) + CS%PFu(I,j,k)) + CS%diffu(I,j,k)
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v_bc_accel(i,J,k) = (CS%Cav(i,J,k) + CS%PFv(i,J,k)) + CS%diffv(i,J,k)
  enddo ; enddo ; enddo
  call cpu_clock_end(id_clock_btforce)

  if (CS%debug) then
    call check_redundant("pre-btstep CS%Ca ", CS%Cau, CS%Cav, G)
    call check_redundant("pre-btstep CS%PF ", CS%PFu, CS%PFv, G)
    call check_redundant("pre-btstep CS%diff ", CS%diffu, CS%diffv, G)
    call check_redundant("pre-btstep u_bc_accel ", u_bc_accel, v_bc_accel, G)
  endif

  if (G%nonblocking_updates) then
    call cpu_clock_begin(id_clock_pass)
    if (CS%calc_bbl) then
      if (associated(CS%visc%Ray_u) .and. associated(CS%visc%Ray_v)) &
        call pass_vector_complete(pid_Ray, CS%visc%Ray_u, CS%visc%Ray_v, G%Domain, &
                         To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(CS%visc%bbl_thick_u) .and. associated(CS%visc%bbl_thick_v)) &
        call pass_vector_complete(pid_bbl_h, CS%visc%bbl_thick_u, CS%visc%bbl_thick_v, &
                       G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(CS%visc%kv_bbl_u) .and. associated(CS%visc%kv_bbl_v)) &
        call pass_vector_complete(pid_kv_bbl, CS%visc%kv_bbl_u, CS%visc%kv_bbl_v, &
                       G%Domain, To_All+SCALAR_PAIR, CGRID_NE)

      ! CS%calc_bbl is set to .false. now that the message passing is completed.
      CS%calc_bbl = .false.
    endif
    call pass_var_complete(pid_eta_PF, CS%eta_PF(:,:), G%Domain)
    call pass_var_complete(pid_eta_in, eta_in(:,:), G%Domain)
    if (CS%readjust_velocity) &
      call pass_vector_complete(pid_uhbt_in, uhbt_in, vhbt_in, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  call cpu_clock_begin(id_clock_vertvisc)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u_out(i,j,k) = G%umask(i,j) * (u_in(i,j,k) + dt * u_bc_accel(I,j,k))
  enddo ; enddo ;  enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v_out(i,j,k) = G%vmask(i,j) * (v_in(i,j,k) + dt * v_bc_accel(i,J,k))
  enddo ; enddo ;  enddo
  call enable_averaging(dt, Time_local, CS%diag)
  call set_viscous_ML(u_in, v_in, h_in, CS%tv, fluxes, CS%visc, dt, G, &
                      CS%set_visc_CSp)
  call disable_averaging(CS%diag)

  call vertvisc_coef(u_out, v_out, h_in, fluxes, CS%visc, dt, G, CS%vertvisc_CSp)
  call vertvisc_remnant(CS%visc, CS%visc_rem_u, CS%visc_rem_v, dt, G, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)

  call cpu_clock_begin(id_clock_pass)
  if (G%nonblocking_updates) then
    pid_visc = pass_vector_start(CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                                 To_All+SCALAR_PAIR, CGRID_NE)
  else
    call pass_var(CS%eta_PF(:,:), G%Domain, complete=.false.)
    call pass_var(eta_in(:,:), G%Domain)
    if (CS%readjust_velocity) call pass_vector(uhbt_in, vhbt_in, G%Domain)
    call pass_vector(CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE)
  endif
  call cpu_clock_end(id_clock_pass)

  call cpu_clock_begin(id_clock_btcalc)
  ! Calculate the relative layer weights for determining barotropic quantities.
  call btcalc(h_in, G, CS%barotropic_CSp)
  call bt_mass_source(h_in, eta_in, fluxes, .true., dt_therm, dt_since_flux, &
                      G, CS%barotropic_CSp)
  call cpu_clock_end(id_clock_btcalc)

  if (G%nonblocking_updates) then
    call cpu_clock_begin(id_clock_pass)
    call pass_vector_complete(pid_visc, CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE)
    call cpu_clock_end(id_clock_pass)
  endif

! u_accel_bt = layer accelerations due to barotropic solver
  if (CS%flux_BT_coupling) then
    call cpu_clock_begin(id_clock_continuity)
    if (CS%readjust_velocity) then
      ! Adjust the input velocites so that their transports match uhbt_out & vhbt_out.
      call continuity(u_in, v_in, h_in, h_out, uh_in, vh_in, dt, G, &
                      CS%continuity_CSp, uhbt_in, vhbt_in, CS%OBC, &
                      CS%visc_rem_u, CS%visc_rem_v, u_adj, v_adj, &
                      BT_cont=CS%BT_cont)
      u_init => u_adj ; v_init => v_adj
      if (ASSOCIATED(CS%diag%du_adj2)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
        CS%diag%du_adj2(I,j,k) = u_adj(I,j,k) - u_in(I,j,k)
      enddo ; enddo ; enddo ; endif
      if (ASSOCIATED(CS%diag%dv_adj2)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
        CS%diag%dv_adj2(i,J,k) = v_adj(i,J,k) - v_in(i,J,k)
      enddo ; enddo ; enddo ; endif
    else
      call continuity(u_in, v_in, h_in, h_out, uh_in, vh_in, dt, G, CS%continuity_CSp, &
                      OBC=CS%OBC, BT_cont=CS%BT_cont)
      u_init => u_in ; v_init => v_in
    endif
    call cpu_clock_end(id_clock_continuity)
    call cpu_clock_begin(id_clock_btstep)
    if (CS%calc_dtbt) call set_dtbt(G, CS%barotropic_CSp, eta_in, CS%pbce, &
                                    CS%BT_cont)
    call btstep(.true., uh_in, vh_in, eta_in, dt, u_bc_accel, v_bc_accel, &
                fluxes, CS%pbce, CS%eta_PF, uh, vh, CS%u_accel_bt, &
                CS%v_accel_bt, eta_out, CS%uhbt, CS%vhbt, G, &
                CS%barotropic_CSp, CS%visc_rem_u, CS%visc_rem_v, &
                uhbt_out = uhbt_out, vhbt_out = vhbt_out, OBC = CS%OBC, &
                BT_cont = CS%BT_cont, eta_PF_start = eta_PF_start)
    call cpu_clock_end(id_clock_btstep)
  else
    u_init => u_in ; v_init => v_in
    call cpu_clock_begin(id_clock_btstep)
    if (CS%calc_dtbt) call set_dtbt(G, CS%barotropic_CSp, eta_in, CS%pbce)
    call btstep(.false., u_in, v_in, eta_in, dt, u_bc_accel, v_bc_accel, &
                fluxes, CS%pbce, CS%eta_PF, u_av, v_av, CS%u_accel_bt, &
                CS%v_accel_bt, eta_out, CS%uhbt, CS%vhbt, G, CS%barotropic_CSp,&
                CS%visc_rem_u, CS%visc_rem_v, OBC=CS%OBC, &
                eta_PF_start=eta_PF_start)
    call cpu_clock_end(id_clock_btstep)
  endif

! u_out = u_in + dt_pred*( u_bc_accel + u_accel_bt )
  dt_pred = dt * CS%be
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v_out(i,J,k) = G%vmask(i,J) * (v_init(i,J,k) + dt_pred * &
                    (v_bc_accel(i,J,k) + CS%v_accel_bt(i,J,k)))
  enddo ; enddo ;  enddo

  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u_out(i,j,k) = G%umask(i,j) * (u_init(i,j,k) + dt_pred  * &
                    (u_bc_accel(I,j,k) + CS%u_accel_bt(I,j,k)))
  enddo ; enddo ;  enddo
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) then
    call uchksum(u_out,"Predictor 1 u",G,haloshift=0)
    call vchksum(v_out,"Predictor 1 v",G,haloshift=0)
    call hchksum(G%H_to_kg_m2*h_in,"Predictor 1 h",G,haloshift=1)
    call uchksum(G%H_to_kg_m2*uh,"Predictor 1 uh",G,haloshift=2)
    call vchksum(G%H_to_kg_m2*vh,"Predictor 1 vh",G,haloshift=2)
!   call GOLD_state_chksum("Predictor 1", u_out, v_out, h_in, uh, vh, G, haloshift=1)
    call GOLD_accel_chksum("Predictor accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv, &
             CS%diffu, CS%diffv, G, CS%pbce, CS%u_accel_bt, CS%v_accel_bt)
    call GOLD_state_chksum("Predictor 1 init", u_init, v_init, h_in, uh, vh, G, haloshift=2)
    call check_redundant("Predictor 1 u_out", u_out, v_out, G)
    call check_redundant("Predictor 1 uh", uh, vh, G)
  endif

! u_out <- u_out + dt_pred d/dz visc d/dz u_out
! u_av  <- u_av  + dt_pred d/dz visc d/dz u_av
  call cpu_clock_begin(id_clock_vertvisc)
  call vertvisc_coef(u_out, v_out, h_in, fluxes, CS%visc, dt_pred, G, CS%vertvisc_CSp)
  call vertvisc(u_out, v_out, h_in, fluxes, CS%visc, dt_pred, CS%OBC, G, &
                CS%vertvisc_CSp)
  if (G%nonblocking_updates) then
    call cpu_clock_end(id_clock_vertvisc) ; call cpu_clock_begin(id_clock_pass)
    pid_u = pass_vector_start(u_out, v_out, G%Domain)
    call cpu_clock_end(id_clock_pass) ; call cpu_clock_begin(id_clock_vertvisc)
  endif
  call vertvisc_remnant(CS%visc, CS%visc_rem_u, CS%visc_rem_v, dt_pred, G, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)

  call cpu_clock_begin(id_clock_pass)
  call pass_vector(CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                   To_All+SCALAR_PAIR, CGRID_NE)
  if (G%nonblocking_updates) then
    call pass_vector_complete(pid_u, u_out, v_out, G%Domain)
  else
    call pass_vector(u_out, v_out, G%Domain)
  endif
  call cpu_clock_end(id_clock_pass)

! uh = u_av * h_in
! h_out = h_in + dt * div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(u_out, v_out, h_in, h_out, uh, vh, dt, G, CS%continuity_CSp, &
                  CS%uhbt, CS%vhbt, CS%OBC, &
                  CS%visc_rem_u, CS%visc_rem_v, u_av, v_av, BT_cont=CS%BT_cont)
  call cpu_clock_end(id_clock_continuity)

  call cpu_clock_begin(id_clock_pass)
  call pass_var(h_out, G%Domain)
  if (G%nonblocking_updates) then
    pid_u_av = pass_vector_start(u_av, v_av, G%Domain)
    pid_uh = pass_vector_start(uh(:,:,:), vh(:,:,:), G%Domain)
  else
    call pass_vector(u_av, v_av, G%Domain, complete=.false.)
    call pass_vector(uh(:,:,:), vh(:,:,:), G%Domain)
  endif
  call cpu_clock_end(id_clock_pass)

  if (CS%BT_include_udhdt) then
    call cpu_clock_begin(id_clock_continuity)
    ! Estimate u dh_dt for driving the barotropic solver.
    call continuity(u_av, v_av, h_out, h_tmp, uh_tmp, vh_tmp, dt, G, &
                    CS%continuity_CSp, OBC=CS%OBC)
    do j=js,je ; do I=Isq,Ieq ; u_dhdt(I,j) = 0.0 ; enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie ; v_dhdt(i,J) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js,je ; do I=is-1,ie
      u_dhdt(I,j) = u_dhdt(I,j) + (uh_tmp(I,j,k) - uh(I,j,k)) * Idt
    enddo ; enddo ; enddo
    do k=1,nz ; do J=js-1,je ; do i=is,ie
      v_dhdt(i,J) = v_dhdt(i,J) + (vh_tmp(i,J,k) - vh(i,J,k)) * Idt
    enddo ; enddo ; enddo
    call cpu_clock_end(id_clock_continuity)
  endif

! h_av = (h_in + h_out)/2
  do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
    h_av(i,j,k) = 0.5*(h_in(i,j,k) + h_out(i,j,k))
  enddo ; enddo ; enddo

! The correction phase of the time step starts here.
  call enable_averaging(dt, Time_local, CS%diag)

!   Calculate a revised estimate of the free-surface height correction to be
! used in the next call to btstep.  This call is at this point so that
! h_out can be changed if CS%begw /= 0.
! eta_cor = ...                 (hidden inside CS%barotropic_CSp)
  call cpu_clock_begin(id_clock_btcalc)
  call bt_mass_source(h_out, eta_out, fluxes, .false., dt_therm, &
                      dt_since_flux+dt, G, CS%barotropic_CSp)
  call cpu_clock_end(id_clock_btcalc)

  if (CS%begw /= 0.0) then
    ! h_out <- (1-begw)*h_in + begw*h_out
    ! Back up h_out to the value it would have had after a time-step of
    ! begw*dt.  h_out is not used again until recalculated by continuity.
    do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
      h_out(i,j,k) = (1.0-CS%begw)*h_in(i,j,k) + CS%begw*h_out(i,j,k)
    enddo ; enddo ; enddo

! PFu = d/dx M(h_out,T,S)
! pbce = dM/deta
    call cpu_clock_begin(id_clock_pres)
    call PressureForce(h_out, CS%tv, CS%PFu, CS%PFv, G, &
                       CS%PressureForce_CSp, p_surf, CS%pbce, CS%eta_PF)
    call cpu_clock_end(id_clock_pres)
    call cpu_clock_begin(id_clock_pass)
    call pass_var(CS%eta_PF(:,:), G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  if (G%nonblocking_updates) then
    call cpu_clock_begin(id_clock_pass)
    call pass_vector_complete(pid_u_av, u_av, v_av, G%Domain)
    call pass_vector_complete(pid_uh, uh(:,:,:), vh(:,:,:), G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  if (CS%debug) then
    call GOLD_state_chksum("Predictor ", u_out, v_out, h_out, uh, vh, G)
    call uchksum(u_av,"Predictor avg u",G,haloshift=1)
    call vchksum(v_av,"Predictor avg v",G,haloshift=1)
    call hchksum(G%H_to_kg_m2*h_av,"Predictor avg h",G,haloshift=0)
  ! call GOLD_state_chksum("Predictor avg ", u_av, v_av,  h_av,uh, vh, G)
    call check_redundant("Predictor u_out ", u_out, v_out, G)
    call check_redundant("Predictor u_out ", uh, vh, G)
  endif

! diffu = horizontal viscosity terms (u_av)
  call cpu_clock_begin(id_clock_horvisc)
  call horizontal_viscosity(u_av, v_av, h_av, CS%diffu, CS%diffv, &
                            CS%MEKE, CS%Varmix, G, CS%hor_visc_CSp, OBC=CS%OBC)
  call cpu_clock_end(id_clock_horvisc)

! CAu = -(f+zeta_av)/h_av vh + d/dx KE_av
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(u_av, v_av, h_av, uh, vh, CS%CAu, CS%CAv, G,CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)

! Calculate the momentum forcing terms for the barotropic equations.

! u_bc_accel = CAu + PFu + diffu(u[n-1])
  call cpu_clock_begin(id_clock_btforce)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u_bc_accel(I,j,k) = (CS%Cau(I,j,k) + CS%PFu(I,j,k)) + CS%diffu(I,j,k)
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v_bc_accel(i,J,k) = (CS%Cav(i,J,k) + CS%PFv(i,J,k)) + CS%diffv(i,J,k)
  enddo ; enddo ; enddo
  call cpu_clock_end(id_clock_btforce)

  if (CS%debug) then
    call check_redundant("corr pre-btstep CS%Ca ", CS%Cau, CS%Cav, G)
    call check_redundant("corr pre-btstep CS%PF ", CS%PFu, CS%PFv, G)
    call check_redundant("corr pre-btstep CS%diff ", CS%diffu, CS%diffv, G)
    call check_redundant("corr pre-btstep u_bc_accel ", u_bc_accel, v_bc_accel, G)
  endif

! u_accel_bt = layer accelerations due to barotropic solver
! pbce = dM/deta
  call cpu_clock_begin(id_clock_btstep)
  if (CS%flux_BT_coupling .and. CS%BT_include_udhdt) then
    call btstep(.true., uh_in, vh_in, eta_in, dt, u_bc_accel, v_bc_accel, &
                fluxes, CS%pbce, CS%eta_PF, uh, vh, CS%u_accel_bt, &
                CS%v_accel_bt, eta_out, CS%uhbt, CS%vhbt, G, &
                CS%barotropic_CSp, CS%visc_rem_u, CS%visc_rem_v, etaav=eta_av, &
                uhbt_out = uhbt_out, vhbt_out = vhbt_out, OBC=CS%OBC, &
                BT_cont = CS%BT_cont, eta_PF_start = eta_PF_start, &
                sum_u_dhdt=u_dhdt, sum_v_dhdt=v_dhdt)
  elseif (CS%flux_BT_coupling) then
    call btstep(.true., uh_in, vh_in, eta_in, dt, u_bc_accel, v_bc_accel, &
                fluxes, CS%pbce, CS%eta_PF, uh, vh, CS%u_accel_bt, &
                CS%v_accel_bt, eta_out, CS%uhbt, CS%vhbt, G, &
                CS%barotropic_CSp, CS%visc_rem_u, CS%visc_rem_v, etaav=eta_av, &
                uhbt_out = uhbt_out, vhbt_out = vhbt_out, OBC=CS%OBC, &
                BT_cont = CS%BT_cont, eta_PF_start = eta_PF_start)
  else
    call btstep(.false., u_in, v_in, eta_in, dt, u_bc_accel, v_bc_accel, &
                fluxes, CS%pbce, CS%eta_PF, u_av, v_av, CS%u_accel_bt, &
                CS%v_accel_bt, eta_out, CS%uhbt, CS%vhbt, G, &
                CS%barotropic_CSp, CS%visc_rem_u, CS%visc_rem_v, &
                etaav=eta_av, OBC=CS%OBC, eta_PF_start=eta_PF_start)
  endif
  call cpu_clock_end(id_clock_btstep)

  if (CS%debug) then
    call check_redundant("u_accel_bt ", CS%u_accel_bt, CS%v_accel_bt, G)
  endif

! u_out = u_in + dt*( u_bc_accel + u_accel_bt )
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u_out(i,j,k) = G%umask(i,j) * (u_init(i,j,k) + dt * &
                    (u_bc_accel(I,j,k) + CS%u_accel_bt(I,j,k)))
  enddo ; enddo ; enddo
  if (ASSOCIATED(CS%diag%PFu_tot)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      CS%diag%PFu_tot(i,j,k) = CS%PFu(i,j,k)
  enddo ; enddo ; enddo ; endif
  if (ASSOCIATED(CS%diag%CAu_tot)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      CS%diag%CAu_tot(i,j,k) = CS%CAu(i,j,k) !+ CS%u_accel_bt(i,j) - CS%diag%PFu_bt(i,j)
  enddo ; enddo ; enddo ; endif

  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v_out(i,j,k) = G%vmask(i,j) * (v_init(i,j,k) + dt * &
                    (v_bc_accel(i,J,k) + CS%v_accel_bt(i,J,k)))
  enddo ; enddo ; enddo
  if (ASSOCIATED(CS%diag%PFv_tot)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    CS%diag%PFv_tot(i,j,k) = CS%PFv(i,j,k)
  enddo ; enddo ; enddo ; endif
  if (ASSOCIATED(CS%diag%CAv_tot)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    CS%diag%CAv_tot(i,j,k) = CS%CAv(i,j,k) !+ CS%v_accel_bt(i,j) - CS%diag%PFv_bt(i,j)
  enddo ; enddo ; enddo ; endif
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) then
    call uchksum(u_out,"Corrector 1 u",G,haloshift=0)
    call vchksum(v_out,"Corrector 1 v",G,haloshift=0)
    call hchksum(G%H_to_kg_m2*h_in,"Corrector 1 h",G,haloshift=2)
    call uchksum(G%H_to_kg_m2*uh,"Corrector 1 uh",G,haloshift=2)
    call vchksum(G%H_to_kg_m2*vh,"Corrector 1 vh",G,haloshift=2)
  ! call GOLD_state_chksum("Corrector 1", u_out, v_out, h_in, uh, vh, G, haloshift=1)
    call GOLD_accel_chksum("Corrector accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv, &
             CS%diffu, CS%diffv, G, CS%pbce, CS%u_accel_bt, CS%v_accel_bt)
  endif

! u_out <- u_out + dt d/dz visc d/dz u_out
! u_av  <- u_av  + dt d/dz visc d/dz u_av
  call cpu_clock_begin(id_clock_vertvisc)
  call vertvisc_coef(u_out, v_out, h_in, fluxes, CS%visc, dt, G, CS%vertvisc_CSp)
  call vertvisc(u_out, v_out, h_in, fluxes, CS%visc, dt, CS%OBC, G, &
                CS%vertvisc_CSp)
  if (G%nonblocking_updates) then
    call cpu_clock_end(id_clock_vertvisc) ; call cpu_clock_begin(id_clock_pass)
    pid_u = pass_vector_start(u_out, v_out, G%Domain)
    call cpu_clock_end(id_clock_pass) ; call cpu_clock_begin(id_clock_vertvisc)
  endif
  call vertvisc_remnant(CS%visc, CS%visc_rem_u, CS%visc_rem_v, dt, G, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)

  call cpu_clock_begin(id_clock_pass)
  call pass_vector(CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                   To_All+SCALAR_PAIR, CGRID_NE)
  if (G%nonblocking_updates) then
    call pass_vector_complete(pid_u, u_out, v_out, G%Domain)
  else
    call pass_vector(u_out, v_out, G%Domain)
  endif
  call cpu_clock_end(id_clock_pass)

! uh = u_av * h_in
! h_out = h_in + dt * div . uh
  if (CS%flux_BT_coupling) then
    ! u_av and v_av adjusted so their mass transports match uhbt and vhbt.
    ! Also, determine the values of u_out and v_out so that their transports
    ! that agree with uhbt_out and vhbt_out.
    if (ASSOCIATED(CS%diag%du_adj)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      CS%diag%du_adj(I,j,k) = u_out(I,j,k)
    enddo ; enddo ; enddo ; endif
    if (ASSOCIATED(CS%diag%dv_adj)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      CS%diag%dv_adj(i,J,k) = v_out(i,J,k)
    enddo ; enddo ; enddo ; endif
    if (CS%BT_include_udhdt) then
      do j=js,je ; do I=is-1,ie ; uhbt_in(I,j) = uhbt_out(I,j) ; enddo ; enddo
      do J=js-1,je ; do i=is,ie ; vhbt_in(i,J) = vhbt_out(i,J) ; enddo ; enddo
      CS%readjust_velocity = .true.
    endif
    call cpu_clock_begin(id_clock_continuity)
    call continuity(u_out, v_out, h_in, h_out, uh, vh, dt, G, &
                    CS%continuity_CSp, CS%uhbt, CS%vhbt, CS%OBC, &
                    CS%visc_rem_u, CS%visc_rem_v, u_av, v_av, &
                    uhbt_out, vhbt_out, u_out, v_out)
    call cpu_clock_end(id_clock_continuity)
    if (G%nonblocking_updates) then
      call cpu_clock_begin(id_clock_pass)
      pid_h = pass_var_start(h_out, G%Domain)
      call cpu_clock_end(id_clock_pass)
    endif
    if (ASSOCIATED(CS%diag%du_adj)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      CS%diag%du_adj(I,j,k) = u_out(I,j,k) - CS%diag%du_adj(I,j,k)
    enddo ; enddo ; enddo ; endif
    if (ASSOCIATED(CS%diag%dv_adj)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      CS%diag%dv_adj(i,J,k) = v_out(i,J,k) - CS%diag%dv_adj(i,J,k)
    enddo ; enddo ; enddo ; endif

    call cpu_clock_begin(id_clock_vertvisc)
    call vertvisc_limit_vel(u_out, v_out, h_in, fluxes, CS%visc, dt, G, CS%vertvisc_CSp)
    if (G%nonblocking_updates) then
      call cpu_clock_end(id_clock_vertvisc) ; call cpu_clock_begin(id_clock_pass)
      pid_u = pass_vector_start(u_out, v_out, G%Domain)
      call cpu_clock_end(id_clock_pass) ; call cpu_clock_begin(id_clock_vertvisc)
    endif
    call vertvisc_limit_vel(u_av, v_av, h_in, fluxes, CS%visc, dt, G, CS%vertvisc_CSp)
    call cpu_clock_end(id_clock_vertvisc)

    call cpu_clock_begin(id_clock_pass)
    if (G%nonblocking_updates) then
      call pass_var_complete(pid_h, h_out, G%Domain)
      call pass_vector_complete(pid_u, u_out, v_out, G%Domain)
    else
      call pass_var(h_out, G%Domain)
      call pass_vector(u_out, v_out, G%Domain, complete=.false.)
    endif
    call cpu_clock_end(id_clock_pass)
  else
    ! u_av and v_av adjusted so their mass transports match uhbt and vhbt.
    call cpu_clock_begin(id_clock_continuity)
    call continuity(u_out, v_out, h_in, h_out, uh, vh, dt, G, &
                    CS%continuity_CSp, CS%uhbt, CS%vhbt, CS%OBC, &
                    CS%visc_rem_u, CS%visc_rem_v, u_av, v_av)
    call cpu_clock_end(id_clock_continuity)
    call cpu_clock_begin(id_clock_pass)
    call pass_var(h_out, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  call cpu_clock_begin(id_clock_pass)
  if (G%nonblocking_updates) then
    pid_uh = pass_vector_start(uh(:,:,:), vh(:,:,:), G%Domain)
    pid_u_av = pass_vector_start(u_av, v_av, G%Domain)
  else
    call pass_vector(u_av, v_av, G%Domain, complete=.false.)
    call pass_vector(uh(:,:,:), vh(:,:,:), G%Domain)
  endif
  call cpu_clock_end(id_clock_pass)

  if (associated(CS%OBC)) then
    call Radiation_Open_Bdry_Conds(CS%OBC, u_out, u_in, v_out, v_in, &
                                   h_out, h_in, G, CS%open_boundary_CSp)
  endif

! h_av = (h_in + h_out)/2
  do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
    h_av(i,j,k) = 0.5*(h_in(i,j,k) + h_out(i,j,k))
  enddo ; enddo ; enddo

  if (G%nonblocking_updates) then
    call cpu_clock_begin(id_clock_pass)
    call pass_vector_complete(pid_uh, uh(:,:,:), vh(:,:,:), G%Domain)
    call pass_vector_complete(pid_u_av, u_av, v_av, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  do k=1,nz ; do j=js-2,je+2 ; do I=Isq-2,Ieq+2
    CS%uhtr(I,j,k) = CS%uhtr(I,j,k) + uh(I,j,k)*dt
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq-2,Jeq+2 ; do i=is-2,ie+2
    CS%vhtr(i,J,k) = CS%vhtr(i,J,k) + vh(i,J,k)*dt
  enddo ; enddo ; enddo

  if (CS%thickness_diffuse .and. .not.CS%thickness_diffuse_first) then
    call cpu_clock_begin(id_clock_thick_diff)
    if (associated(CS%VarMix)) &
      call calc_slope_function(h_out, CS%tv, G, CS%VarMix)
    call thickness_diffuse(h_out, CS%uhtr, CS%vhtr, CS%tv, dt, G, &
                           CS%MEKE, CS%VarMix, CS%thickness_diffuse_CSp)
    call cpu_clock_end(id_clock_thick_diff)
    call cpu_clock_begin(id_clock_pass)
    call pass_var(h_out, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  if (CS%mixedlayer_restrat) then
    call cpu_clock_begin(id_clock_ml_restrat)
    call mixedlayer_restrat(h_out, CS%uhtr,CS%vhtr,CS%tv, fluxes, dt, &
                            G, CS%mixedlayer_restrat_CSp)
    call cpu_clock_end(id_clock_ml_restrat)
    call cpu_clock_begin(id_clock_pass)
    call pass_var(h_out, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  if (associated(CS%MEKE)) then
    call step_forward_MEKE(CS%MEKE, h_out, CS%visc, dt, G, CS%MEKE_CSp)
  endif

!   The time-averaged free surface height has already been set by the last
!  call to btstep.

!   Here various terms used in to update the momentum equations are
! offered for averaging.
  if (CS%id_PFu > 0) call post_data(CS%id_PFu, CS%diag%PFu_tot, CS%diag)
  if (CS%id_PFv > 0) call post_data(CS%id_PFv, CS%diag%PFv_tot, CS%diag)
  if (CS%id_CAu > 0) call post_data(CS%id_CAu, CS%diag%CAu_tot, CS%diag)
  if (CS%id_CAv > 0) call post_data(CS%id_CAv, CS%diag%CAv_tot, CS%diag)

!   Here the thickness fluxes are offered for averaging.
  if (CS%id_uh > 0) call post_data(CS%id_uh, uh, CS%diag)
  if (CS%id_vh > 0) call post_data(CS%id_vh, vh, CS%diag)
  if (CS%id_uav > 0) call post_data(CS%id_uav, u_av, CS%diag)
  if (CS%id_vav > 0) call post_data(CS%id_vav, v_av, CS%diag)
  if (CS%id_u_BT_accel > 0) call post_data(CS%id_u_BT_accel, CS%u_accel_bt, CS%diag)
  if (CS%id_v_BT_accel > 0) call post_data(CS%id_v_BT_accel, CS%v_accel_bt, CS%diag)
  if (CS%id_du_adj > 0) call post_data(CS%id_du_adj, CS%diag%du_adj, CS%diag)
  if (CS%id_dv_adj > 0) call post_data(CS%id_dv_adj, CS%diag%dv_adj, CS%diag)
  if (CS%id_du_adj2 > 0) call post_data(CS%id_du_adj2, CS%diag%du_adj2, CS%diag)
  if (CS%id_dv_adj2 > 0) call post_data(CS%id_dv_adj2, CS%diag%dv_adj2, CS%diag)
  if (CS%BT_include_udhdt) then
    if (CS%id_h_dudt > 0) call post_data(CS%id_h_dudt, u_dhdt, CS%diag)
    if (CS%id_h_dvdt > 0) call post_data(CS%id_h_dvdt, v_dhdt, CS%diag)
  endif
  if (CS%debug) then
    call GOLD_state_chksum("Corrector ", u_out, v_out, h_out, uh, vh, G)
    call uchksum(u_av,"Corrector avg u",G,haloshift=1)
    call vchksum(v_av,"Corrector avg v",G,haloshift=1)
    call hchksum(G%H_to_kg_m2*h_av,"Corrector avg h",G,haloshift=1)
 !  call GOLD_state_chksum("Corrector avg ", u_av, v_av, h_av, uh, vh, G)
  endif

end subroutine step_GOLD_dyn_split

subroutine step_GOLD_dyn_unsplit_RK3(u, v, h, Time_local, dt, fluxes, &
                  p_surf_begin, p_surf_end, uh, vh, eta_av, G, CS)
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(inout) :: u
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(inout) :: v
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(inout) :: h
  type(time_type),                     intent(in)    :: Time_local
  real,                                intent(in)    :: dt
  type(forcing),                       intent(in)    :: fluxes
  real, dimension(:,:),                pointer       :: p_surf_begin, p_surf_end
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(inout) :: uh
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(inout) :: vh
  real, dimension(NXMEM_,NYMEM_),      intent(out)   :: eta_av
  type(ocean_grid_type),               intent(inout) :: G
  type(GOLD_control_struct),            pointer      :: CS
! Arguments: u - The input and output zonal velocity, in m s-1.
!  (in)      v - The input and output meridional velocity, in m s-1.
!  (in)      h - The input and output layer thicknesses, in m or kg m-2,
!                depending on whether the Boussinesq approximation is made.
!  (in)      Time_local - The model time at the end of the time step.
!  (in)      dt - The time step in s.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      p_surf_begin - A pointer (perhaps NULL) to the surface pressure
!                     at the beginning of this dynamic step, in Pa.
!  (in)      p_surf_end - A pointer (perhaps NULL) to the surface pressure
!                     at the end of this dynamic step, in Pa.
!  (inout)   uh - The zonal volume or mass transport, in m3 s-1 or kg s-1.
!  (inout)   vh - The meridional volume or mass transport, in m3 s-1 or kg s-1.
!  (out)     eta_av - The time-mean free surface height or column mass, in m or
!                     kg m-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure set up by initialize_GOLD.

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_av, hp
  real, dimension(SZIQ_(G),SZJ_(G),SZK_(G)) :: up, upp
  real, dimension(SZI_(G),SZJQ_(G),SZK_(G)) :: vp, vpp
  real, dimension(:,:), pointer :: p_surf
  real :: dt_pred   ! The time step for the predictor part of the baroclinic
                    ! time stepping.
  logical :: dyn_p_surf
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  dt_pred = dt / 3.0

  h_av(:,:,:) = 0; hp(:,:,:) = 0
  up(:,:,:) = 0; upp(:,:,:) = 0
  vp(:,:,:) = 0; vpp(:,:,:) = 0

  dyn_p_surf = CS%interp_p_surf .and. associated(p_surf_begin) .and. &
               associated(p_surf_end)
  if (dyn_p_surf) then
    call safe_alloc_ptr(p_surf,G%isd,G%ied,G%jsd,G%jed) ; p_surf(:,:) = 0.0
  else
    p_surf => fluxes%p_surf
  endif

! Matsuno's third order accurate three step scheme is used to step
! all of the fields except h.  h is stepped separately.

  if (CS%debug) then
    call GOLD_state_chksum("Start First Predictor ", u, v, h, uh, vh, G)
  endif

! diffu = horizontal viscosity terms (u,h)
  call enable_averaging(dt,Time_local, CS%diag)
  call cpu_clock_begin(id_clock_horvisc)
  call horizontal_viscosity(u, v, h, CS%diffu, CS%diffv, CS%MEKE, CS%Varmix, &
                            G, CS%hor_visc_CSp)
  call cpu_clock_end(id_clock_horvisc)
  call disable_averaging(CS%diag)

! uh = u*h
! hp = h + dt/2 div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(u, v, h, hp, uh, vh, dt*0.5, G, CS%continuity_CSp, OBC=CS%OBC)
  call cpu_clock_end(id_clock_continuity)
  call cpu_clock_begin(id_clock_pass)
  call pass_var(hp, G%Domain)
  call pass_vector(uh, vh, G%Domain)
  call cpu_clock_end(id_clock_pass)

  call enable_averaging(0.5*dt,Time_local-set_time(int(0.5*dt)), CS%diag)
!   Here the thickness fluxes are offered for averaging.
  if (CS%id_uh > 0) call post_data(CS%id_uh, uh, CS%diag)
  if (CS%id_vh > 0) call post_data(CS%id_vh, vh, CS%diag)
  call disable_averaging(CS%diag)

! h_av = (h + hp)/2
! u = u + dt diffu
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      h_av(i,j,k) = (h(i,j,k) + hp(i,j,k)) * 0.5
    enddo ; enddo
    do j=js,je ; do I=Isq,Ieq
      u(i,j,k) = u(i,j,k) + dt * CS%diffu(i,j,k) * G%umask(i,j)
    enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie
      v(i,j,k) = v(i,j,k) + dt * CS%diffv(i,j,k) * G%vmask(i,j)
    enddo ; enddo
    do j=js-2,je+2 ; do I=Isq-2,Ieq+2
      CS%uhtr(i,j,k) = CS%uhtr(i,j,k) + 0.5*dt*uh(i,j,k)
    enddo ; enddo
    do J=Jsq-2,Jeq+2 ; do i=is-2,ie+2
      CS%vhtr(i,j,k) = CS%vhtr(i,j,k) + 0.5*dt*vh(i,j,k)
    enddo ; enddo
  enddo
  call cpu_clock_end(id_clock_mom_update)
  call cpu_clock_begin(id_clock_pass)
  call pass_vector(u, v, G%Domain)
  call cpu_clock_end(id_clock_pass)

! CAu = -(f+zeta)/h_av vh + d/dx KE
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(u, v, h_av, uh, vh, CS%CAu, CS%CAv, G, CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)

! PFu = d/dx M(h_av,T,S)
  call cpu_clock_begin(id_clock_pres)
  if (dyn_p_surf) then ; do j=js-2,je+2 ; do i=is-2,ie+2
    p_surf(i,j) = 0.75*p_surf_begin(i,j) + 0.25*p_surf_end(i,j)
  enddo ; enddo ; endif
  call PressureForce(h_av, CS%tv, CS%PFu, CS%PFv, G, &
                     CS%PressureForce_CSp, p_surf)
  call cpu_clock_end(id_clock_pres)

! up = u + dt_pred * (PFu + CAu)
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    up(i,j,k) = G%umask(i,j) * (u(i,j,k) + dt_pred * &
                               (CS%PFu(i,j,k) + CS%CAu(i,j,k)))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    vp(i,j,k) = G%vmask(i,j) * (v(i,j,k) + dt_pred * &
                               (CS%PFv(i,j,k) + CS%CAv(i,j,k)))
  enddo ; enddo ; enddo
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) then
    call GOLD_state_chksum("Predictor 1", up, vp, h_av, uh, vh, G)
    call GOLD_accel_chksum("Predictor 1 accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv,&
                          CS%diffu, CS%diffv, G)
  endif

! CS%visc contains viscosity and BBL thickness (u_in,h_in)
  if (CS%calc_bbl) then
    call enable_averaging(CS%bbl_calc_time_interval, &
                          Time_local-set_time(int(dt)), CS%diag)
    call set_viscous_BBL(u, v, h_av, CS%tv, CS%visc, G, CS%set_visc_CSp)
    call cpu_clock_begin(id_clock_pass)
    if (associated(CS%visc%Ray_u) .and. associated(CS%visc%Ray_v)) &
      call pass_vector(CS%visc%Ray_u, CS%visc%Ray_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE)
    if (associated(CS%visc%kv_bbl_u) .and. associated(CS%visc%kv_bbl_v)) then
      call pass_vector(CS%visc%bbl_thick_u, CS%visc%bbl_thick_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE, complete=.false.)
      call pass_vector(CS%visc%kv_bbl_u, CS%visc%kv_bbl_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE)
    endif
    call cpu_clock_end(id_clock_pass)
    call disable_averaging(CS%diag)
    CS%calc_bbl = .false.
  endif

 ! up <- up + dt/2 d/dz visc d/dz up
  call cpu_clock_begin(id_clock_vertvisc)
  call enable_averaging(dt, Time_local, CS%diag)
  call set_viscous_ML(u, v, h_av, CS%tv, fluxes, CS%visc, dt*0.5, G, &
                      CS%set_visc_CSp)
  call disable_averaging(CS%diag)
  call vertvisc_coef(up, vp, h_av, fluxes, CS%visc, dt*0.5, G, CS%vertvisc_CSp)
  call vertvisc(up, vp, h_av, fluxes, CS%visc, dt*0.5, CS%OBC, G, &
                CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)
  call cpu_clock_begin(id_clock_pass)
  call pass_vector(up, vp, G%Domain)
  call cpu_clock_end(id_clock_pass)

! uh = up * hp
! h_av = hp + dt/2 div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(up, vp, hp, h_av, uh, vh, &
                  (0.5*dt), G, CS%continuity_CSp, OBC=CS%OBC)
  call cpu_clock_end(id_clock_continuity)
  call cpu_clock_begin(id_clock_pass)
  call pass_var(h_av, G%Domain)
  call pass_vector(uh, vh, G%Domain)
  call cpu_clock_end(id_clock_pass)

! h_av <- (hp + h_av)/2
  do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
    h_av(i,j,k) = (hp(i,j,k) + h_av(i,j,k)) * 0.5
  enddo ; enddo ; enddo

! CAu = -(f+zeta(up))/h_av vh + d/dx KE(up)
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(up, vp, h_av, uh, vh, CS%CAu, CS%CAv, &
                 G, CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)

! PFu = d/dx M(h_av,T,S)
  call cpu_clock_begin(id_clock_pres)
  if (dyn_p_surf) then ; do j=js-2,je+2 ; do i=is-2,ie+2
    p_surf(i,j) = 0.25*p_surf_begin(i,j) + 0.75*p_surf_end(i,j)
  enddo ; enddo ; endif
  call PressureForce(h_av, CS%tv, CS%PFu, CS%PFv, G, &
                     CS%PressureForce_CSp, p_surf)
  call cpu_clock_end(id_clock_pres)

! upp = u + dt/2 * ( PFu + CAu )
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    upp(i,j,k) = G%umask(i,j) * (u(i,j,k) + dt * 0.5 * &
            (CS%PFu(i,j,k) + CS%CAu(i,j,k)))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    vpp(i,j,k) = G%vmask(i,j) * (v(i,j,k) + dt * 0.5 * &
            (CS%PFv(i,j,k) + CS%CAv(i,j,k)))
  enddo ; enddo ; enddo
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) then
    call GOLD_state_chksum("Predictor 2", upp, vpp, h_av, uh, vh, G)
    call GOLD_accel_chksum("Predictor 2 accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv,&
                          CS%diffu, CS%diffv, G)
  endif

! upp <- upp + dt/2 d/dz visc d/dz upp
  call cpu_clock_begin(id_clock_vertvisc)
  call vertvisc_coef(upp, vpp, hp, fluxes, CS%visc, dt*0.5, G, CS%vertvisc_CSp)
  call vertvisc(upp, vpp, hp, fluxes, CS%visc, dt*0.5, CS%OBC, G, &
                CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)
  call cpu_clock_begin(id_clock_pass)
  call pass_vector(upp, vpp, G%Domain)
  call cpu_clock_end(id_clock_pass)

! uh = upp * hp
! h = hp + dt/2 div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(upp, vpp, hp, h, uh, vh, &
                  (dt*0.5), G, CS%continuity_CSp, OBC=CS%OBC)
  call cpu_clock_end(id_clock_continuity)
  call cpu_clock_begin(id_clock_pass)
  call pass_var(h, G%Domain)
  call pass_vector(uh, vh, G%Domain)
  call cpu_clock_end(id_clock_pass)

  call enable_averaging(0.5*dt, Time_local, CS%diag)
!   Here the thickness fluxes are offered for averaging.
  if (CS%id_uh > 0) call post_data(CS%id_uh, uh, CS%diag)
  if (CS%id_vh > 0) call post_data(CS%id_vh, vh, CS%diag)
  call disable_averaging(CS%diag)

! h_av = (h + hp)/2
  do k=1,nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      h_av(i,j,k) = 0.5*(h(i,j,k) + hp(i,j,k))
    enddo ; enddo
    do j=js-2,je+2 ; do I=Isq-2,Ieq+2
      CS%uhtr(i,j,k) = CS%uhtr(i,j,k) + 0.5*dt*uh(i,j,k)
    enddo ; enddo
    do J=Jsq-2,Jeq+2 ; do i=is-2,ie+2
      CS%vhtr(i,j,k) = CS%vhtr(i,j,k) + 0.5*dt*vh(i,j,k)
    enddo ; enddo
  enddo

  call enable_averaging(dt,Time_local, CS%diag)

! CAu = -(f+zeta(upp))/h_av vh + d/dx KE(upp)
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(upp, vpp, h_av, uh, vh, CS%CAu, CS%CAv, &
                 G, CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)

! PFu = d/dx M(h_av,T,S)
  call cpu_clock_begin(id_clock_pres)
  call PressureForce(h_av, CS%tv, CS%PFu, CS%PFv, G, &
                     CS%PressureForce_CSp, p_surf)
  call cpu_clock_end(id_clock_pres)

! u = u + dt * ( PFu + CAu )
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u(i,j,k) = G%umask(i,j) * (u(i,j,k) + dt * &
            (CS%PFu(i,j,k) + CS%CAu(i,j,k)))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v(i,j,k) = G%vmask(i,j) * (v(i,j,k) + dt * &
            (CS%PFv(i,j,k) + CS%CAv(i,j,k)))
  enddo ; enddo ; enddo

! u <- u + dt d/dz visc d/dz u
  call cpu_clock_begin(id_clock_vertvisc)
  call vertvisc_coef(u, v, h_av, fluxes, CS%visc, dt, G, CS%vertvisc_CSp)
  call vertvisc(u, v, h_av, fluxes, CS%visc, dt, CS%OBC, G, &
                CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)
  call cpu_clock_begin(id_clock_pass)
  call pass_vector(u, v, G%Domain)
  call cpu_clock_end(id_clock_pass)

  if (CS%thickness_diffuse .and. .not.CS%thickness_diffuse_first) then
    call cpu_clock_begin(id_clock_thick_diff)
    if (associated(CS%VarMix)) &
      call calc_slope_function(h, CS%tv, G, CS%VarMix)
    call thickness_diffuse(h, CS%uhtr, CS%vhtr, CS%tv, dt, G, &
                           CS%MEKE, CS%VarMix, CS%thickness_diffuse_CSp)
    call cpu_clock_end(id_clock_thick_diff)
    call cpu_clock_begin(id_clock_pass)
    call pass_var(h, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  if (CS%mixedlayer_restrat) then
    call cpu_clock_begin(id_clock_ml_restrat)
    call mixedlayer_restrat(h, CS%uhtr ,CS%vhtr, CS%tv, fluxes, dt, &
                            G, CS%mixedlayer_restrat_CSp)
    call cpu_clock_end(id_clock_ml_restrat)
    call cpu_clock_begin(id_clock_pass)
    call pass_var(h, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  if (associated(CS%MEKE)) then
    call step_forward_MEKE(CS%MEKE, h, CS%visc, dt, G, CS%MEKE_CSp)
  endif

  if (CS%debug) then
    call GOLD_state_chksum("Corrector", u, v, h, uh, vh, G)
    call GOLD_accel_chksum("Corrector accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv, &
                          CS%diffu, CS%diffv, G)
  endif

  if (G%Boussinesq) then
    do j=js,je ; do i=is,ie ; eta_av(i,j) = -G%D(i,j) ; enddo ; enddo
  else
    do j=js,je ; do i=is,ie ; eta_av(i,j) = 0.0 ; enddo ; enddo
  endif
  do k=1,nz ; do j=js,je ; do i=is,ie
    eta_av(i,j) = eta_av(i,j) + h_av(i,j,k)
  enddo ; enddo ; enddo 
  
  if (dyn_p_surf) deallocate(p_surf)

!   Here various terms used in to update the momentum equations are
! offered for averaging.
  if (CS%id_PFu > 0) call post_data(CS%id_PFu, CS%PFu, CS%diag)
  if (CS%id_PFv > 0) call post_data(CS%id_PFv, CS%PFv, CS%diag)
  if (CS%id_CAu > 0) call post_data(CS%id_CAu, CS%CAu, CS%diag)
  if (CS%id_CAv > 0) call post_data(CS%id_CAv, CS%CAv, CS%diag)

end subroutine step_GOLD_dyn_unsplit_RK3

subroutine find_total_transport(u_in, v_in, h_in, uh_tot, vh_tot, dt, G, CS)
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(in)    :: u_in
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(in)    :: v_in
  real,                                intent(in)    :: dt
  real, dimension(NXMEMQ_,NYMEM_),     intent(out)   :: uh_tot
  real, dimension(NXMEM_,NYMEMQ_),     intent(out)   :: vh_tot
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(in)    :: h_in
  type(ocean_grid_type),               intent(inout) :: G
  type(GOLD_control_struct),           pointer       :: CS
!   This subroutine determines the vertically summed transport based on input
! velocities and thicknesses.  The individual layers' transports are not retained.

! Arguments: u_in - The input zonal velocity, in m s-1. (Intent in.)
!  (in)      v_in - The input meridional velocity, in m s-1.
!  (in)      h_in - The input layer thicknesses, in m or kg m-2, depending on
!                   whether the Boussinesq approximation is made.
!  (out)     uh_tot - The vertically summed zonal and meridional volume or mass
!  (out)     vh_tot -  transports, in m3 s-1 or kg s-1.
!  (in)      dt - The time step in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure set up by initialize_GOLD.

  ! Temporary arrays to contain layer thickness fluxes in m3 s-1 or kg s-1.
  real, dimension(SZIQ_(G),SZJ_(G),SZK_(G)) :: uh_temp 
  real, dimension(SZI_(G),SZJQ_(G),SZK_(G)) :: vh_temp
  ! A temporary array to contain layer projected thicknesses in m or kg m-2.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G))  :: h_temp
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call cpu_clock_begin(id_clock_continuity)
  call continuity(u_in, v_in, h_in, h_temp, uh_temp, vh_temp, dt, G, &
                  CS%continuity_CSp, OBC=CS%OBC)
  call cpu_clock_end(id_clock_continuity)

  do j=js,je ; do I=is-1,ie ; uh_tot(I,j) = uh_temp(I,j,1) ; enddo ; enddo
  do k=2,nz ; do j=js,je ; do I=is-1,ie
    uh_tot(I,j) = uh_tot(I,j) + uh_temp(I,j,k)
  enddo ; enddo ; enddo
  do J=js-1,je ; do i=is,ie ; vh_tot(i,J) = vh_temp(i,J,1) ; enddo ; enddo
  do k=2,nz ; do J=js-1,je ; do i=is,ie
    vh_tot(i,J) = vh_tot(i,J) + vh_temp(i,J,k)
  enddo ; enddo ; enddo

end subroutine find_total_transport

!#######################################################################
subroutine initialize_GOLD(Time, param_file, dirs, CS, Time_in)
  type(time_type), target, intent(inout) :: Time
  type(param_file_type), intent(out)     :: param_file
  type(directories), intent(out)         :: dirs
  type(GOLD_control_struct), pointer      :: CS
  type(time_type), optional, intent(in)  :: Time_in
! Arguments: Time - The model time, set in this routine.
!  (out)     param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      dirs - A structure containing several relevant directory paths.
!  (out)     CS - A pointer set in this routine to the GOLD control structure.
!  (in)      Time_in - An optional time passed to GOLD_initialize to use when
!                      the model is not being started from a restart file.
  type(ocean_grid_type), pointer :: grid ! A pointer to a structure containing
                                  ! metrics and related information.
  type(diag_ptrs), pointer :: diag
  character(len=4), parameter :: vers_num = 'v2.0'
  character(len=128) :: version = '$Id: GOLD.F90,v 13.0.2.32.2.108 2011/10/12 14:12:46 Alistair.Adcroft Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: Isdq, Iedq, Jsdq, Jedq
  real    :: dtbt
  real    :: Z_diag_int      ! The minimum interval between calculations of
                             ! Depth-space diagnostic quantities, in s.
  real, pointer, dimension(:,:,:,:) :: &
    u, &                     ! u : Zonal velocity, in m s-1.
    v, &                     ! v : Meridional velocity, in m s-1.
    h                        ! h : Layer thickness, in m or kg m-2.
  real, allocatable, dimension(:,:,:) :: e ! The interface heights in m.
  type(GOLD_restart_CS),  pointer :: restart_CSp_tmp => NULL()
  integer :: nkml, nkbl, ntstep, verbosity
  real    :: SSH_max_def, SSS_max_def, SST_max_def, SST_min_def ! Defaults. 
  logical :: new_sim
  logical :: use_BT_cont_type
  logical :: use_geothermal  ! If true, apply geothermal heating.
  logical :: use_EOS         ! If true, density is calculated from T & S using
                             ! an equation of state.
  logical :: use_tides       ! If true, tidal momentum forcing is used.
  logical :: save_IC         ! If true, save the initial conditions.
  character(len=80) :: IC_file ! A file into which the initial conditions are
                             ! written in a new run if save_IC is true.
  type(vardesc) :: vd
  type(time_type) :: Start_time
  type(GOLD_initialization_struct) :: init_CS
  type(ocean_internal_state) :: GOLD_internal_state

  if (associated(CS)) then
    call GOLD_error(WARNING, "initialize_GOLD called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)
  grid => CS%grid ; CS%Time => Time
  diag => CS%diag

  id_clock_init = cpu_clock_id('Ocean Initialization', grain=CLOCK_SUBCOMPONENT)
  call cpu_clock_begin(id_clock_init)

  Start_time = Time ; if (present(Time_in)) Start_time = Time_in

  call Get_GOLD_Input(param_file, dirs)
  verbosity=2; call read_param(param_file,"VERBOSITY",verbosity)
  call GOLD_set_verbosity(verbosity)
  call find_obsolete_params(param_file)
  call GOLD_domains_init(grid%domain, param_file)
  call GOLD_checksums_init(param_file)

  call GOLD_io_init(param_file)
  call diag_mediator_init(param_file)
  call GOLD_grid_init(grid, param_file)
  is = grid%isc ; ie = grid%iec ; js = grid%jsc ; je = grid%jec ; nz = grid%ke
  isd = grid%isd ; ied = grid%ied ; jsd = grid%jsd ; jed = grid%jed
  Isdq = grid%Isdq ; Iedq = grid%Iedq ; Jsdq = grid%Jsdq ; Jedq = grid%Jedq
  call set_diag_mediator_grid(grid, diag)

! Initialize the local control structure.  The model will fail if SPLIT or
! TEMPERATURE is not actively defined or undefined, so they have no defaults.
  CS%adiabatic = .false. ; CS%bulkmixedlayer = .false. ; CS%debug = .false.
  CS%thickness_diffuse = .false. ; CS%mixedlayer_restrat = .false.
  CS%Hmix = 1.0 ; CS%begw = 0.0 ; Z_diag_int = 0.0 ; CS%dtbt_reset_period = -1.0
  CS%flux_BT_coupling = .true. ; CS%readjust_BT_trans = .false.
  CS%thickness_diffuse_first = .false. ; CS%debug_truncations = .false.
  use_BT_cont_type = .false. ; CS%interp_p_surf = .false.
  CS%BT_include_udhdt = .false.
  CS%smooth_ssh_passes = 0.0
  call read_param(param_file,"SPLIT",CS%split,.true.)
  call read_param(param_file,"TEMPERATURE",CS%use_temperature,.true.)
  use_EOS = CS%use_temperature
  call read_param(param_file,"NONLINEAR_EOS",use_EOS)
  if (use_EOS /= CS%use_temperature) then
    if (is_root_pe()) call GOLD_error(WARNING, &
      "GOLD: NONLINEAR_EOS is a depricated option.  Instead define " // &
      "USE_EOS to use an equation of state to calculate density.")
  endif
  call read_param(param_file,"USE_EOS",use_EOS)
  if (use_EOS .and. .not.CS%use_temperature) call GOLD_error(FATAL, &
         "GOLD: TEMPERATURE must be defined to use USE_EOS")
  call read_param(param_file,"ADIABATIC",CS%adiabatic)
  call read_param(param_file,"BULKMIXEDLAYER",CS%bulkmixedlayer)
  if (CS%adiabatic .and. CS%bulkmixedlayer) call GOLD_error(FATAL, &
         "GOLD: ADIABATIC and BULKMIXEDLAYER can not both be defined.")
  call read_param(param_file,"THICKNESSDIFFUSE",CS%thickness_diffuse)
  call read_param(param_file,"THICKNESSDIFFUSE_FIRST",CS%thickness_diffuse_first)
  if (CS%bulkmixedlayer) &
    call read_param(param_file,"MIXEDLAYER_RESTRAT",CS%mixedlayer_restrat)
  call read_param(param_file,"DEBUG",CS%debug)
  call read_param(param_file,"DEBUG_TRUNCATIONS",CS%debug_truncations)
  call read_param(param_file,"DT",CS%dt,.true.)
  CS%dt_therm = CS%dt ; call read_param(param_file,"DT_THERM",CS%dt_therm)
  call read_param(param_file,"BE",CS%be,.true.)
  call read_param(param_file,"BEGW",CS%begw)
  call read_param(param_file,"HMIX",CS%Hmix)
  call read_param(param_file,"MIN_Z_DIAG_INTERVAL",Z_diag_int)
  call read_param(param_file,"FLUX_BT_COUPLING",CS%flux_BT_coupling)
  call read_param(param_file,"INTERPOLATE_P_SURF",CS%interp_p_surf)
  call read_param(param_file,"SSH_SMOOTHING_PASSES",CS%smooth_ssh_passes)
  if (CS%split .and. CS%flux_BT_coupling) then
    call read_param(param_file,"USE_BT_CONT_TYPE",use_BT_cont_type)
    call read_param(param_file,"BT_INCLUDE_UDHDT", CS%BT_include_udhdt)
  endif
  if (CS%split) then
    dtbt = 0.0 ; call read_param(param_file, "DTBT", dtbt)
    if (dtbt <= 0.0) &
      call read_param(param_file, "DTBT_RESET_PERIOD", CS%dtbt_reset_period)
    call read_param(param_file, "READJUST_BT_TRANS", CS%readjust_BT_trans)
  endif
  if (.not.(CS%split .and. CS%flux_BT_coupling) .or. CS%adiabatic) &
    CS%readjust_BT_trans = .false.
  if (use_EOS) call read_param(param_file,"P_REF",CS%tv%P_Ref,.true.)
  if (CS%use_temperature) call read_param(param_file,"C_P",CS%C_p,.true.)

  CS%check_bad_surface_vals=.FALSE.
  call read_param(param_file,"CHECK_BAD_SURFACE_VALS",CS%check_bad_surface_vals)
  if (CS%check_bad_surface_vals) then
    SSH_max_def = 20.0 ; SSS_max_def = 45.0
    SST_max_def = 45.0 ; SST_min_def = -2.1
    CS%bad_val_ssh_max = SSH_max_def ; CS%bad_val_sss_max = SSS_max_def
    CS%bad_val_sst_max = SST_max_def ; CS%bad_val_sst_min = SST_min_def
    call read_param(param_file,"BAD_VAL_SSH_MAX",CS%bad_val_ssh_max)
    call read_param(param_file,"BAD_VAL_SSS_MAX",CS%bad_val_sss_max)
    call read_param(param_file,"BAD_VAL_SST_MAX",CS%bad_val_sst_max)
    call read_param(param_file,"BAD_VAL_SST_MIN",CS%bad_val_sst_min)
  endif

  ! Allocate the auxiliary non-symmetric domain for debugging purposes.
  if (CS%debug) &
    call GOLD_domains_init(grid%Domain_aux, param_file, symmetric=.false.)

  CS%Z_diag_interval = set_time(int((CS%dt_therm) * &
       max(1,floor(0.01 + Z_diag_int/(CS%dt_therm)))))

  CS%use_frazil = .false. ; CS%bound_salinity = .false.
  if (CS%use_temperature) then
    call read_param(param_file,"FRAZIL",CS%use_frazil)
    call read_param(param_file,"BOUND_SALINITY",CS%bound_salinity)
  endif
  use_tides = .false.  ; call read_param(param_file,"TIDES",use_tides)
  if (CS%bulkmixedlayer) then
    nkml = 1 ; call read_param(param_file,"NKML",nkml)
    nkbl = 1 ; call read_param(param_file,"NKBL",nkbl)
  endif

  ! Write all relevant parameters to the model log.
  call log_version(param_file, "GOLD", version, tagname, "")
  call log_param(param_file, "GOLD", "VERBOSITY", verbosity,  &
                 "Integer controlling level of messaging\n" // &
                 "\t0 = Only FATAL messages\n" // &
                 "\t2 = Only FATAL, WARNING, NOTE [default]\n" // &
                 "\t9 = All)")
  call log_param(param_file, "GOLD", "SPLIT", CS%split, &
                 "Use the split time stepping if true.")
  call log_param(param_file, "GOLD", "TEMPERATURE", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state \n"//&
                 "variables.")
  call log_param(param_file, "GOLD", "USE_EOS", use_EOS, &
                 "If true,  density is calculated from temperature and \n"//&
                 "salinity with an equation of state.  TEMPERATURE must \n"//&
                 "be true if USE_EOS is.")
  call log_param(param_file, "GOLD", "ADIABATIC", CS%adiabatic, &
                 "There are no diapycnal mass fluxes if ADIABATIC is \n"//&
                 "true. This assumes that KD = KDML = 0.0 and that \n"//&
                 "there is no buoyancy forcing, but makes the model \n"//&
                 "faster by eliminating subroutine calls.", default=.false.)
  call log_param(param_file, "GOLD", "BULKMIXEDLAYER", CS%bulkmixedlayer, &
                 "If true, use a Kraus-Turner-like bulk mixed layer \n"//&
                 "with transitional buffer layers.  Layers 1 through  \n"//&
                 "NKML+NKBL have variable densities. There must be at \n"//&
                 "least NKML+NKBL+1 layers if BULKMIXEDLAYER is true.", &
                 default=.false.)
  call log_param(param_file, "GOLD", "THICKNESSDIFFUSE", CS%thickness_diffuse, &
                 "If true, interfaces or isopycnal surfaces are diffused, \n"//&
                 "depending on the value of FULL_THICKNESSDIFFUSE.", &
                 default=.false.)
  call log_param(param_file, "GOLD", "THICKNESSDIFFUSE_FIRST", &
                                     CS%thickness_diffuse_first, &
                 "If true, do thickness diffusion before dynamics.\n"//&
                 "This is only used if THICKNESSDIFFUSE is true.", &
                 default=.false.)
  call log_param(param_file, "GOLD", "MIXEDLAYER_RESTRAT",CS%mixedlayer_restrat, &
                 "If true, a density-gradient dependent re-stratifying \n"//&
                 "flow is imposed in the mixed layer. \n"//&
                 "This is only used if BULKMIXEDLAYER is true.", default=.false.)
  call log_param(param_file, "GOLD", "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", default=.false.)
  call log_param(param_file, "GOLD", "DEBUG_TRUNCATIONS", CS%debug_truncations, &
                 "If true, calculate all diagnostics that are useful for \n"//&
                 "debugging truncations.", default=.false.)
  call log_param(param_file, "GOLD", "DT", CS%dt, &
                 "The (baroclinic) dynamics time step.  The time-step that \n"//&
                 "is actually used will be an integer fraction of the \n"//&
                 "forcing time-step (DT_FORCING in ocean-only mode or the \n"//&
                 "coupling timestep in coupled mode.)", units="s")
  call log_param(param_file, "GOLD", "DT_THERM", CS%dt_therm, &
                 "The thermodynamic and tracer advection time step. \n"//&
                 "Ideally DT_THERM should be an integer multiple of DT \n"//&
                 "and less than the forcing or coupling time-step. \n"//&
                 "By default DT_THERM is set to DT.", units="s", default=CS%dt)
  call log_param(param_file, "GOLD", "BE", CS%be, &
                 "If SPLIT is true, BE determines the relative weighting \n"//&
                 "of a  2nd-order Runga-Kutta baroclinic time stepping \n"//&
                 "scheme (0.5) and a backward Euler scheme (1) that is \n"//&
                 "used for the Coriolis and inertial terms.  BE may be \n"//&
                 "from 0.5 to 1, but instability may occur near 0.5.", &
                 units="nondim")
  call log_param(param_file, "GOLD", "BEGW", CS%begw, &
                 "If SPILT is true, BEGW is a number from 0 to 1 that \n"//&
                 "controls the extent to which the treatment of gravity \n"//&
                 "waves is forward-backward (0) or simulated backward \n"//&
                 "Euler (1).  0 is almost always used.", units="nondim", &
                 default=0.0)
  call log_param(param_file, "GOLD", "HMIX", CS%Hmix, &
                 "If BULKMIXEDLAYER is false, HMIX is the depth over \n"//&
                 "which to average to find surface properties like SST \n"//&
                 "and SSS, and over which the vertical viscosity and \n"//&
                 "diapycnal diffusivity are elevated.  HMIX is only used \n"//&
                 "directly if BULKMIXEDLAYER is false, but provides a \n"//&
                 "default value for other variables if BULKMIXEDLAYER is \n"//&
                 "true.", units="m", default=1.0)
  call log_param(param_file, "GOLD", "MIN_Z_DIAG_INTERVAL", Z_diag_int, &
                 "The minimum amount of time in seconds between \n"//&
                 "calculations of depth-space diagnostics. Making this \n"//&
                 "larger than DT_THERM reduces the  performance penalty \n"//&
                 "of regridding to depth online.", units="s", default=0.0)
  call log_param(param_file, "GOLD", "FLUX_BT_COUPLING", CS%flux_BT_coupling, &
                 "If true, use mass fluxes to ensure consistency between \n"//&
                 "the baroclinic and barotropic modes. This is only used \n"//&
                 "if SPLIT is true.", default=.true.)
  call log_param(param_file, "GOLD", "INTERPOLATE_P_SURF", CS%interp_p_surf, &
                 "If true, linearly interpolate the surface pressure \n"//&
                 "over the coupling time step, using the specified value \n"//&
                 "at the end of the step.", default=.false.)
  call log_param(param_file, "GOLD", "SSH_SMOOTHING_PASSES", CS%smooth_ssh_passes, &
                 "The number of Laplacian smoothing passes to apply to the \n"//&
                 "the sea surface height that is reported to the sea-ice.", &
                 units="nondim", default=0.0)
  if (CS%split .and. CS%flux_BT_coupling) then
    call log_param(param_file, "GOLD", "USE_BT_CONT_TYPE", use_BT_cont_type, &
                 "If true, use a structure with elements that describe \n"//&
                 "effective face areas from the summed continuity solver \n"//&
                 "as a function the barotropic flow in coupling between \n"//&
                 "the barotropic and baroclinic flow.  This is only used \n"//&
                 "if SPLIT and FLUX_BT_COUPLING are true. \n", default=.false.)
    call log_param(param_file, "GOLD", "BT_INCLUDE_UDHDT", CS%BT_include_udhdt, &
                 "If true, include the barotropic transport tendancies \n"//&
                 "from sum(u dhdt) and sum(v dhdt) in the barotropic \n"//&
                 "solver.", default=.false.)
  endif
  if (CS%split) then
    call log_param(param_file, "GOLD", "READJUST_BT_TRANS", CS%readjust_BT_trans, &
                 "If true, make a barotropic adjustment to the layer \n"//&
                 "velocities after the thermodynamic part of the step \n"//&
                 "to ensure that the interaction between the thermodynamics \n"//&
                 "and the continuity solver do not change the barotropic \n"//&
                 "transport.  This is only used if FLUX_BT_COUPLING and \n"//&
                 "SPLIT are true.", default=.false.)
    call log_param(param_file, "GOLD", "DTBT_RESET_PERIOD", CS%dtbt_reset_period, &
                 "The period between recalculations of DTBT (if DTBT <= 0). \n"//&
                 "If DTBT_RESET_PERIOD is negative, DTBT is set based \n"//&
                 "only on information available at initialization.  If \n"//&
                 "dynamic, DTBT will be set at least every forcing time \n"//&
                 "step, and if 0, every dynamics time step.  This is \n"//&
                 "only used if SPLIT is true.", units="s", default=-1.0)
  endif
  if (CS%use_temperature) then
    call log_param(param_file, "GOLD", "FRAZIL", CS%use_frazil, &
                 "If true, water freezes if it gets too cold, and the \n"//&
                 "the accumulated heat deficit is returned in the \n"//&
                 "surface state.  This is only used if TEMPERATURE is true.", &
                 default=.false.)
    call log_param(param_file, "GOLD", "BOUND_SALINITY", CS%bound_salinity, &
                 "If true, limit salinity to being positive. (The sea-ice \n"//&
                 "model may ask for more salt than is available and \n"//&
                 "drive the salinity negative otherwise.)", default=.false.)
    call log_param(param_file, "GOLD", "C_P", CS%C_p, &
                 "The heat capacity of sea water, approximated as a \n"//&
                 "constant. This is only used if TEMPERATURE is true.", &
                 units="J kg-1 K-1")
  endif
  call log_param(param_file, "GOLD", "TIDES", use_tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)
  if (CS%bulkmixedlayer) then
    call log_param(param_file, "GOLD", "NKML", nkml, &
                 "The number of sublayers within the mixed layer if \n"//&
                 "BULKMIXEDLAYER is true.", units="nondim", default=1)
    call log_param(param_file, "GOLD", "NKBL", nkbl, &
                 "The number of layers that are used as variable density \n"//&
                 "buffer layers if BULKMIXEDLAYER is true.", units="nondim", &
                 default=1)
  endif
  if (use_EOS) call log_param(param_file, "GOLD", "P_REF", CS%tv%P_Ref, &
                 "The pressure that is used for calculating the coordinate \n"//&
                 "density.  (1 Pa = 1e4 dbar, so 2e7 is commonly used.) \n"//&
                 "This is only used if USE_EOS and TEMPERATURE are true.", &
                 units="Pa")

  call log_param(param_file, "GOLD", "CHECK_BAD_SURFACE_VALS", &
                                     CS%check_bad_surface_vals, &
                 "If true, check the surface state for ridiculous values.", &
                 default=.false.)
  if (CS%check_bad_surface_vals) then
    call log_param(param_file, "GOLD", "BAD_VAL_SSH_MAX", CS%bad_val_ssh_max, &
                 "The value of SSH above which a bad value message is \n"//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", units="m", &
                 default=SSH_max_def)
    call log_param(param_file, "GOLD", "BAD_VAL_SSS_MAX", CS%bad_val_sss_max, &
                 "The value of SSS above which a bad value message is \n"//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", units="PSU", &
                 default=SSS_max_def)
    call log_param(param_file, "GOLD", "BAD_VAL_SST_MAX", CS%bad_val_sst_max, &
                 "The value of SST above which a bad value message is \n"//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", &
                 units="deg C", default=SST_max_def)
    call log_param(param_file, "GOLD", "BAD_VAL_SST_MIN", CS%bad_val_sst_min, &
                 "The value of SST below which a bad value message is \n"//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", &
                 units="deg C", default=SST_min_def)
  endif

  call GOLD_timing_init(CS)

  call advect_tracer_init(param_file, CS%tracer_CSp)

! Allocate and initialize space for the primary GOLD variables.
  ALLOC(CS%u(Isdq:Iedq,jsd:jed,nz,2)) ; CS%u(:,:,:,:) = 0.0
  ALLOC(CS%v(isd:ied,Jsdq:Jedq,nz,2)) ; CS%v(:,:,:,:) = 0.0
  ALLOC(CS%h(isd:ied,jsd:jed,nz,2))   ; CS%h(:,:,:,:) = grid%Angstrom
  u => CS%u ; v => CS%v ; h => CS%h
  ALLOC(CS%uh(Isdq:Iedq,jsd:jed,nz))  ; CS%uh(:,:,:) = 0.0
  ALLOC(CS%vh(isd:ied,Jsdq:Jedq,nz))  ; CS%vh(:,:,:) = 0.0
  ALLOC(CS%diffu(Isdq:Iedq,jsd:jed,nz)) ; CS%diffu(:,:,:) = 0.0
  ALLOC(CS%diffv(isd:ied,Jsdq:Jedq,nz)) ; CS%diffv(:,:,:) = 0.0
  ALLOC(CS%CAu(Isdq:Iedq,jsd:jed,nz)) ; CS%CAu(:,:,:) = 0.0
  ALLOC(CS%CAv(isd:ied,Jsdq:Jedq,nz)) ; CS%CAv(:,:,:) = 0.0
  ALLOC(CS%PFu(Isdq:Iedq,jsd:jed,nz)) ; CS%PFu(:,:,:) = 0.0
  ALLOC(CS%PFv(isd:ied,Jsdq:Jedq,nz)) ; CS%PFv(:,:,:) = 0.0
  if (CS%use_temperature) then
    ALLOC(CS%T(isd:ied,jsd:jed,nz))   ; CS%T(:,:,:) = 0.0
    ALLOC(CS%S(isd:ied,jsd:jed,nz))   ; CS%S(:,:,:) = 0.0
    CS%tv%T => CS%T ; CS%tv%S => CS%S
    call register_tracer(CS%tv%T, "T", param_file, CS%tracer_CSp)
    call register_tracer(CS%tv%S, "S", param_file, CS%tracer_CSp)
  endif
  if (CS%use_frazil) then
    allocate(CS%tv%frazil(isd:ied,jsd:jed)) ; CS%tv%frazil(:,:) = 0.0
  endif
  if (CS%bound_salinity) then
    allocate(CS%tv%salt_deficit(isd:ied,jsd:jed)) ; CS%tv%salt_deficit(:,:)=0.0
  endif
  
  if (CS%bulkmixedlayer) then
    CS%tv%nk_Rml = nkml+nkbl
    ALLOC(CS%Rml(isd:ied,jsd:jed,nz)) ; CS%Rml(:,:,:) = 0.0
    CS%tv%Rml => CS%Rml
    allocate(CS%tv%Hml(isd:ied,jsd:jed)) ; CS%tv%Hml(:,:) = 0.0
    if (.not.CS%use_temperature) &
      call register_tracer(CS%tv%Rml, "Rml", param_file, CS%tracer_CSp)
  endif

  ALLOC(CS%uhtr(Isdq:Iedq,jsd:jed,nz)) ; CS%uhtr(:,:,:) = 0.0
  ALLOC(CS%vhtr(isd:ied,Jsdq:Jedq,nz)) ; CS%vhtr(:,:,:) = 0.0

  GOLD_internal_state%u => u ; GOLD_internal_state%v => v ; GOLD_internal_state%h =>h
  GOLD_internal_state%uh => CS%uh ; GOLD_internal_state%vh => CS%vh
  GOLD_internal_state%diffu => CS%diffu ; GOLD_internal_state%diffv => CS%diffv
  GOLD_internal_state%PFu => CS%PFu ; GOLD_internal_state%PFv => CS%PFv
  GOLD_internal_state%CAu => CS%CAu ; GOLD_internal_state%CAv => CS%CAv
  if (CS%use_temperature) then
    GOLD_internal_state%T => CS%T ; GOLD_internal_state%S => CS%S
  endif
  if (CS%split) then
    ALLOC(CS%eta(isd:ied,jsd:jed,2))    ; CS%eta(:,:,:) = 0.0
    ALLOC(CS%uhbt(Isdq:Iedq,jsd:jed))   ; CS%uhbt(:,:) = 0.0
    ALLOC(CS%vhbt(isd:ied,Jsdq:Jedq))   ; CS%vhbt(:,:) = 0.0
    ALLOC(CS%uhbt_in(Isdq:Iedq,jsd:jed)) ; CS%uhbt_in(:,:) = 0.0
    ALLOC(CS%vhbt_in(isd:ied,Jsdq:Jedq)) ; CS%vhbt_in(:,:) = 0.0
    ALLOC(CS%u_av(Isdq:Iedq,jsd:jed,nz)); CS%u_av(:,:,:) = 0.0
    ALLOC(CS%v_av(isd:ied,Jsdq:Jedq,nz)); CS%v_av(:,:,:) = 0.0
    ALLOC(CS%h_av(isd:ied,jsd:jed,nz))  ; CS%h_av(:,:,:) = grid%Angstrom
    ALLOC(CS%eta_PF(isd:ied,jsd:jed))   ; CS%eta_PF(:,:) = 0.0
    ALLOC(CS%pbce(isd:ied,jsd:jed,nz))  ; CS%pbce(:,:,:) = 0.0
    ALLOC(CS%u_accel_bt(Isdq:Iedq,jsd:jed,nz)) ; CS%u_accel_bt(:,:,:) = 0.0
    ALLOC(CS%v_accel_bt(isd:ied,Jsdq:Jedq,nz)) ; CS%v_accel_bt(:,:,:) = 0.0
    ALLOC(CS%visc_rem_u(Isdq:Iedq,jsd:jed,nz)) ; CS%visc_rem_u(:,:,:) = 0.0
    ALLOC(CS%visc_rem_v(isd:ied,Jsdq:Jedq,nz)) ; CS%visc_rem_v(:,:,:) = 0.0

    GOLD_internal_state%u_accel_bt => CS%u_accel_bt
    GOLD_internal_state%v_accel_bt => CS%v_accel_bt
    GOLD_internal_state%pbce => CS%pbce
    GOLD_internal_state%eta => CS%eta
    GOLD_internal_state%u_av => CS%u_av
    GOLD_internal_state%v_av => CS%v_av
  endif
  if (CS%interp_p_surf) then
    allocate(CS%p_surf_prev(isd:ied,jsd:jed)) ; CS%p_surf_prev(:,:) = 0.0
  endif

  ALLOC(CS%ave_ssh(isd:ied,jsd:jed)) ; CS%ave_ssh(:,:) = 0.0

! Use the Wright equation of state by default, unless otherwise specified
! Note: this line and the following block ought to be in a separate
! initialization routine for tv.
  if (use_EOS) call select_eqn_of_state(param_file,CS%tv%eqn_of_state)
  if (CS%use_temperature) then
    allocate(CS%tv%TempxPmE(isd:ied,jsd:jed))
    CS%tv%TempxPmE(:,:) = 0.0
    use_geothermal = .false.
    call read_param(param_file, "DO_GEOTHERMAL", use_geothermal)
    call log_param(param_file, "GOLD", "DO_GEOTHERMAL", use_geothermal, &
                 "If true, apply geothermal heating.", default=.false.)
    if (use_geothermal) then
      allocate(CS%tv%internal_heat(isd:ied,jsd:jed))
      CS%tv%internal_heat(:,:) = 0.0
    endif
  endif
  if (use_BT_cont_type) call alloc_BT_cont_type(CS%BT_cont, grid)

!   Set the fields that are needed for bitwise identical restarting
! the time stepping scheme.
  call restart_init(param_file, CS%restart_CSp)
  call set_restart_fields(grid, param_file, CS)
!   This subroutine calls user-specified tracer registration routines.
! Additional calls can be added to GOLD_tracer_flow_control.F90.
  call call_tracer_register(grid, param_file, CS%tracer_flow_CSp, &
                            diag, CS%tracer_CSp, CS%restart_CSp)
  call MEKE_alloc_register_restart(grid, param_file, CS%MEKE, CS%restart_CSp)

!   Initialize all of the relevant fields.
  if (associated(CS%tracer_CSp)) &
    init_CS%advect_tracer_CSp => CS%tracer_CSp

  call cpu_clock_begin(id_clock_GOLD_init)
  call GOLD_initialize(u(:,:,:,1), v(:,:,:,1), h(:,:,:,1), CS%tv, Time, &
                       grid, param_file, dirs, CS%restart_CSp, init_CS, Time_in)
  call cpu_clock_end(id_clock_GOLD_init)

  if (associated(init_CS%advect_tracer_CSp)) &
    CS%tracer_CSp => init_CS%advect_tracer_CSp
  if (associated(init_CS%OBC)) then
    CS%OBC => init_CS%OBC
    call open_boundary_init(Time, grid, param_file, diag, CS%open_boundary_CSp)
  endif
! if (associated(init_CS%sponge_CSp)) CS%sponge_CSp => init_CS%sponge_CSp

  if (use_tides) call tidal_forcing_init(Time, grid, param_file, CS%tides_CSp)

  call continuity_init(Time, grid, param_file, diag, CS%continuity_CSp)
  call CoriolisAdv_init(Time, grid, param_file, diag, CS%CoriolisAdv_CSp)
  call PressureForce_init(Time, grid, param_file, diag, CS%PressureForce_CSp, &
                          init_CS%compress_CSp, CS%tides_CSp)

  call hor_visc_init(Time, grid, param_file, diag, CS%hor_visc_CSp)
  call vertvisc_init(GOLD_internal_state, Time, grid, param_file, diag, dirs, &
                     CS%ntrunc, CS%vertvisc_CSp)
  call set_visc_init(Time, grid, param_file, diag, CS%visc, CS%set_visc_CSp)
  call thickness_diffuse_init(Time, grid, param_file, diag, CS%thickness_diffuse_CSp)
  if (CS%mixedlayer_restrat) &
    call mixedlayer_restrat_init(Time, grid, param_file, diag, CS%mixedlayer_restrat_CSp)
  call MEKE_init(Time, grid, param_file, diag, CS%MEKE_CSp, CS%MEKE)
  call VarMix_init(Time, grid, param_file, diag, CS%VarMix)
  call GOLD_diagnostics_init(GOLD_internal_state, Time, grid, param_file, &
                             diag, CS%diagnostics_CSp)
  call GOLD_diag_to_Z_init(Time, grid, param_file, diag, CS%diag_to_Z_CSp)
  CS%Z_diag_time = Start_time + CS%Z_diag_interval * (1 + &
    ((Time + set_time(int(CS%dt_therm))) - Start_time) / CS%Z_diag_interval)

  if (.not.CS%adiabatic) then
    call diabatic_driver_init(Time, grid, param_file, diag, CS%diabatic_CSp, &
                     CS%tracer_flow_CSp, init_CS%sponge_CSp, CS%diag_to_Z_CSp)
  endif

  call register_diags(Time, grid, CS)

  call advect_tracer_diag_init(Time, grid, diag, CS%tracer_CSp)
  if (CS%use_temperature) then
    ! If needed T_adx, etc., would have been allocated in register_diags.
    call add_tracer_diagnostics("T", CS%tracer_CSp, CS%T_adx, CS%T_ady, &
                                CS%T_diffx, CS%T_diffy)
    call add_tracer_diagnostics("S", CS%tracer_CSp, CS%S_adx, CS%S_ady, &
                                CS%S_diffx, CS%S_diffy)
    call add_tracer_2d_diagnostics("T", CS%tracer_CSp, CS%T_adx_2d, CS%T_ady_2d, &
                                   CS%T_diffx_2d, CS%T_diffy_2d)
    call add_tracer_2d_diagnostics("S", CS%tracer_CSp, CS%S_adx_2d, CS%S_ady_2d, &
                                    CS%S_diffx_2d, CS%S_diffy_2d)
    call register_Z_tracer(CS%tv%T, "temp_z", "Potential Temperature", "degC", Time, &
                           grid, CS%diag_to_Z_CSp)
    call register_Z_tracer(CS%tv%S, "salt_z", "Salinity", "PSU", Time, &
                           grid, CS%diag_to_Z_CSp)
  endif

  ! This subroutine initializes any tracer packages.
  new_sim = ((dirs%input_filename(1:1) == 'n') .and. &
             (LEN_TRIM(dirs%input_filename) == 1))
  call tracer_flow_control_init(.not.new_sim, Time, grid, h(:,:,:,1), CS%OBC, &
           CS%tracer_flow_CSp, init_CS%sponge_CSp, CS%diag_to_Z_CSp)

  call cpu_clock_begin(id_clock_pass_init)
  call pass_vector(u(:,:,:,1),v(:,:,:,1), grid%Domain)
  call pass_var(h(:,:,:,1), grid%Domain)

  if (CS%bulkmixedlayer .and. .not.use_eos) then
    call pass_var(CS%tv%Rml, grid%Domain)
  endif
  if (CS%use_temperature) then
    call pass_var(CS%tv%T, grid%Domain)
    call pass_var(CS%tv%S, grid%Domain)
  endif
  call cpu_clock_end(id_clock_pass_init)

  if (CS%split) then
    if (.not. query_initialized(CS%eta(:,:,1),"sfc",CS%restart_CSp))  then
      ! Estimate eta based on the layer thicknesses - h.  With the Boussinesq
      ! approximation, eta is the free surface height anomaly, while without it
      ! eta is the mass of ocean per unit area.  eta always has the same
      ! dimensions as h, either m or kg m-3.  
      !   CS%eta(:,:,1) = 0.0 already from initialization.
      if (grid%Boussinesq) then
        do j=js,je ; do i=is,ie ; CS%eta(i,j,1) = -grid%D(i,j) ; enddo ; enddo
      endif
      do k=1,nz ; do j=js,je ; do i=is,ie
        CS%eta(i,j,1) = CS%eta(i,j,1) + h(i,j,k,1)
      enddo ; enddo ; enddo
    endif  

    call barotropic_init(u(:,:,:,1), v(:,:,:,1), h(:,:,:,1), CS%eta(:,:,1), Time, grid, &
                         param_file, diag, CS%barotropic_CSp, CS%restart_CSp, &
                         CS%tides_CSp)

    if (.not. query_initialized(CS%diffu,"diffu",CS%restart_CSp) .or. &
        .not. query_initialized(CS%diffv,"diffv",CS%restart_CSp)) &
      call horizontal_viscosity(u(:,:,:,1), v(:,:,:,1), h(:,:,:,1), &
                                CS%diffu, CS%diffv, CS%MEKE, CS%VarMix, &
                                grid, CS%hor_visc_CSp)
    if (.not. query_initialized(CS%u_av,"u2", CS%restart_CSp) .or. &
        .not. query_initialized(CS%u_av,"v2", CS%restart_CSp)) then
      CS%u_av(:,:,:) = u(:,:,:,1)
      CS%v_av(:,:,:) = v(:,:,:,1)
    endif
  ! This call is just here to initialize uh and vh.
    if (.not. query_initialized(CS%uh,"uh",CS%restart_CSp) .or. &
        .not. query_initialized(CS%vh,"vh",CS%restart_CSp)) then
      call continuity(u(:,:,:,1), v(:,:,:,1), h(:,:,:,1), h(:,:,:,2), CS%uh, CS%vh, &
                      CS%dt, grid, CS%continuity_CSp, OBC=CS%OBC)
      call cpu_clock_begin(id_clock_pass_init)
      call pass_var(h(:,:,:,2), grid%Domain)
      call cpu_clock_end(id_clock_pass_init)
      CS%h_av(:,:,:) = 0.5*(h(:,:,:,1) + h(:,:,:,2))
    else
      if (.not. query_initialized(CS%h_av,"h2",CS%restart_CSp)) &
        CS%h_av(:,:,:) = h(:,:,:,1)
    endif

    !   Determine whether there is a barotropic transport that is to be used
    ! to adjust the layers' velocities.
    CS%readjust_velocity = .false.
    if (CS%readjust_BT_trans .or. CS%BT_include_udhdt) then
      if (query_initialized(CS%uhbt_in,"uhbt_in",CS%restart_CSp) .and. &
          query_initialized(CS%vhbt_in,"vhbt_in",CS%restart_CSp)) then
        CS%readjust_velocity = .true.
        call cpu_clock_begin(id_clock_pass_init)
        call pass_vector(CS%uhbt_in, CS%vhbt_in, grid%Domain)
        call cpu_clock_end(id_clock_pass_init)
      endif
    endif

    call cpu_clock_begin(id_clock_pass_init)
    call pass_vector(CS%u_av,CS%v_av, grid%Domain)
    call pass_var(CS%h_av, grid%Domain)
    call pass_vector(CS%uh, CS%vh, grid%Domain)
    call cpu_clock_end(id_clock_pass_init)
  endif

  call write_static_fields(grid, CS%diag)
  call enable_averaging(0.0,Time, CS%diag)

!  call calculate_diagnostic_fields(u(:,:,:,1),v(:,:,:,1),h(:,:,:,1), &
!            uh, vh, 1, CS%tv, 0.0, grid, CS%diagnostics_CSp)

!  if (CS%id_u > 0) call post_data(CS%id_u, CS%u(:,:,:,1), CS%diag)
!  if (CS%id_v > 0) call post_data(CS%id_v, CS%v(:,:,:,1), CS%diag)
!  if (CS%id_h > 0) call post_data(CS%id_h, CS%h(:,:,:,1), CS%diag)
!  if (CS%id_T > 0) call post_data(CS%id_T, CS%tv%T, CS%diag)
!  if (CS%id_S > 0) call post_data(CS%id_S, CS%tv%S, CS%diag)
!  if (CS%id_Rml > 0) call post_data(CS%id_Rml, CS%tv%Rml, CS%diag)

  call disable_averaging(CS%diag)

  if (CS%use_frazil) then
    if (.not.query_initialized(CS%tv%frazil,"frazil",CS%restart_CSp)) &
      CS%tv%frazil(:,:) = 0.0
  endif

  if (CS%interp_p_surf) CS%p_surf_prev_set = &
    query_initialized(CS%p_surf_prev,"p_surf_prev",CS%restart_CSp)

  if (.not.query_initialized(CS%ave_ssh,"ave_ssh",CS%restart_CSp)) then
    if (CS%split) then
      call find_eta(h(:,:,:,1), CS%tv, grid%g_Earth, grid, CS%ave_ssh, CS%eta(:,:,1))
    else
      call find_eta(h(:,:,:,1), CS%tv, grid%g_Earth, grid, CS%ave_ssh)
    endif
  endif

  save_IC = .false. ; call read_param(param_file,"SAVE_INITIAL_CONDS",save_IC)
  IC_file = "GOLD_IC" ; call read_param(param_file,"IC_OUTPUT_FILE",IC_file)
  if (save_IC .and. .not.((dirs%input_filename(1:1) == 'r') .and. &
                          (LEN_TRIM(dirs%input_filename) == 1))) then
    allocate(restart_CSp_tmp)
    restart_CSp_tmp = CS%restart_CSp
    allocate(e(SZI_(grid),SZJ_(grid),SZK_(grid)+1))
    call find_eta(h(:,:,:,1), CS%tv, grid%g_Earth, grid, e)
    vd = vardesc("eta","Interface heights",'h','i','1',"meter", 'd')
    call register_restart_field(e, e, vd, .true., restart_CSp_tmp)
    
    call save_restart(dirs%output_directory, Time, 1, grid, &
                      restart_CSp_tmp, filename=IC_file)
    deallocate(e)
    deallocate(restart_CSp_tmp)
  endif

!  call calculate_surface_state(state, u(:,:,:,1), v(:,:,:,1), h(:,:,:,1), &
!                               CS%ave_ssh, grid, CS)

  call cpu_clock_end(id_clock_init)

end subroutine initialize_GOLD


!#######################################################################

subroutine register_diags(Time, G, CS)
  type(time_type),           intent(in)    :: Time
  type(ocean_grid_type),     intent(inout) :: G
  type(GOLD_control_struct), intent(inout) :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure set up by initialize_GOLD.
  character(len=48) :: thickness_units, flux_units, T_flux_units, S_flux_units
  integer :: isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq, nz
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq

  thickness_units = get_thickness_units(G)
  flux_units = get_flux_units(G)
  T_flux_units = get_tr_flux_units(G, "Celsius")
  S_flux_units = get_tr_flux_units(G, "PSU")

  CS%id_u = register_diag_field('ocean_model', 'u', G%axesuL, Time, &
      'Zonal velocity', 'meter second-1')
  CS%id_v = register_diag_field('ocean_model', 'v', G%axesvL, Time, &
      'Meridional velocity', 'meter second-1')
  CS%id_h = register_diag_field('ocean_model', 'h', G%axeshL, Time, &
      'Layer Thickness', thickness_units)
  CS%id_Rml = register_diag_field('ocean_model', 'Rml', G%axeshL, Time, &
      'Mixed Layer Density', 'kg meter-3')
  CS%id_ssh = register_diag_field('ocean_model', 'SSH', G%axesh1, Time, &
      'Sea Surface Height', 'meter')
  CS%id_ssu = register_diag_field('ocean_model', 'SSU', G%axesu1, Time, &
      'Sea Surface Zonal Velocity', 'meter second-1')
  CS%id_ssv = register_diag_field('ocean_model', 'SSV', G%axesv1, Time, &
      'Sea Surface Meridional Velocity', 'meter second-1')
  if (CS%use_temperature) then
    CS%id_T = register_diag_field('ocean_model', 'temp', G%axeshL, Time, &
        'Potential Temperature', 'Celsius')
    CS%id_S = register_diag_field('ocean_model', 'salt', G%axeshL, Time, &
        'Salinity', 'PSU')
    CS%id_sst = register_diag_field('ocean_model', 'SST', G%axesh1, Time, &
        'Sea Surface Temperature', 'Celsius')
    CS%id_sst_sq = register_diag_field('ocean_model', 'SST_sq', G%axesh1, Time, &
        'Sea Surface Temperature Squared', 'Celsius**2')    
    CS%id_sss = register_diag_field('ocean_model', 'SSS', G%axesh1, Time, &
        'Sea Surface Salinity', 'PSU')
    if (CS%id_sst_sq > 0) call safe_alloc_ptr(CS%SST_sq,isd,ied,jsd,jed)    
  endif
  if (CS%use_temperature .and. CS%use_frazil) then
    CS%id_fraz = register_diag_field('ocean_model', 'frazil', G%axesh1, Time, &
          'Heat sink from frazil formation', 'Watt meter-2')
  endif

  CS%id_salt_deficit = register_diag_field('ocean_model', 'salt_deficit', G%axesh1, Time, &
         'Salt sink in ocean due to ice flux', 'g Salt meter-2 s-1')
  CS%id_Heat_PmE = register_diag_field('ocean_model', 'Heat_PmE', G%axesh1, Time, &
         'Heat flux into ocean from mass flux into ocean', 'Watt meter-2')
  CS%id_intern_heat = register_diag_field('ocean_model', 'internal_heat', G%axesh1, Time, &
         'Heat flux into ocean from geothermal or other internal sources', &
         'Watt meter-2')

  CS%id_CAu = register_diag_field('ocean_model', 'CAu', G%axesuL, Time, &
      'Zonal Coriolis and Advective Acceleration', 'meter second-2')
  CS%id_CAv = register_diag_field('ocean_model', 'CAv', G%axesvL, Time, &
      'Meridional Coriolis and Advective Acceleration', 'meter second-2')
  CS%id_PFu = register_diag_field('ocean_model', 'PFu', G%axesuL, Time, &
      'Zonal Pressure Force Acceleration', 'meter second-2')
  CS%id_PFv = register_diag_field('ocean_model', 'PFv', G%axesvL, Time, &
      'Meridional Pressure Force Acceleration', 'meter second-2')
  if (CS%split) then
    if (CS%id_PFu > 0) call safe_alloc_ptr(CS%diag%PFu_tot,Isdq,Iedq,jsd,jed,nz)
    if (CS%id_PFv > 0) call safe_alloc_ptr(CS%diag%PFv_tot,isd,ied,Jsdq,Jedq,nz)
    if (CS%id_CAu > 0) call safe_alloc_ptr(CS%diag%CAu_tot,Isdq,Iedq,jsd,jed,nz)
    if (CS%id_CAv > 0) call safe_alloc_ptr(CS%diag%CAv_tot,isd,ied,Jsdq,Jedq,nz)
  endif
  CS%id_u_BT_accel = register_diag_field('ocean_model', 'u_BT_accel', G%axesul, Time, &
    'Barotropic Anomaly Zonal Acceleration', 'meter second-1')
  CS%id_v_BT_accel = register_diag_field('ocean_model', 'v_BT_accel', G%axesvl, Time, &
    'Barotropic Anomaly Meridional Acceleration', 'meter second-1')

  CS%id_Tadx = register_diag_field('ocean_model', 'T_adx', G%axesul, Time, &
      'Advective Zonal Flux of Potential Temperature', T_flux_units)
  CS%id_Tady = register_diag_field('ocean_model', 'T_ady', G%axesvl, Time, &
      'Advective Meridional Flux of Potential Temperature', T_flux_units)
  CS%id_Tdiffx = register_diag_field('ocean_model', 'T_diffx', G%axesul, Time, &
      'Diffusive Zonal Flux of Potential Temperature', T_flux_units)
  CS%id_Tdiffy = register_diag_field('ocean_model', 'T_diffy', G%axesvl, Time, &
      'Diffusive Meridional Flux of Potential Temperature', T_flux_units)
  if (CS%id_Tadx > 0)   call safe_alloc_ptr(CS%T_adx,Isdq,Iedq,jsd,jed,nz)
  if (CS%id_Tady > 0)   call safe_alloc_ptr(CS%T_ady,isd,ied,Jsdq,Jedq,nz)
  if (CS%id_Tdiffx > 0) call safe_alloc_ptr(CS%T_diffx,Isdq,Iedq,jsd,jed,nz)
  if (CS%id_Tdiffy > 0) call safe_alloc_ptr(CS%T_diffy,isd,ied,Jsdq,Jedq,nz)

  CS%id_Sadx = register_diag_field('ocean_model', 'S_adx', G%axesul, Time, &
      'Advective Zonal Flux of Salinity', S_flux_units)
  CS%id_Sady = register_diag_field('ocean_model', 'S_ady', G%axesvl, Time, &
      'Advective Meridional Flux of Salinity', S_flux_units)
  CS%id_Sdiffx = register_diag_field('ocean_model', 'S_diffx', G%axesul, Time, &
      'Diffusive Zonal Flux of Salinity', S_flux_units)
  CS%id_Sdiffy = register_diag_field('ocean_model', 'S_diffy', G%axesvl, Time, &
      'Diffusive Meridional Flux of Salinity', S_flux_units)
  if (CS%id_Sadx > 0)   call safe_alloc_ptr(CS%S_adx,Isdq,Iedq,jsd,jed,nz)
  if (CS%id_Sady > 0)   call safe_alloc_ptr(CS%S_ady,isd,ied,Jsdq,Jedq,nz)
  if (CS%id_Sdiffx > 0) call safe_alloc_ptr(CS%S_diffx,Isdq,Iedq,jsd,jed,nz)
  if (CS%id_Sdiffy > 0) call safe_alloc_ptr(CS%S_diffy,isd,ied,Jsdq,Jedq,nz)

  CS%id_Tadx_2d = register_diag_field('ocean_model', 'T_adx_2d', G%axesu1, Time, &
      'Vertically Integrated Advective Zonal Flux of Potential Temperature', T_flux_units)
  CS%id_Tady_2d = register_diag_field('ocean_model', 'T_ady_2d', G%axesv1, Time, &
      'Vertically Integrated Advective Meridional Flux of Potential Temperature', T_flux_units)
  CS%id_Tdiffx_2d = register_diag_field('ocean_model', 'T_diffx_2d', G%axesu1, Time, &
      'Vertically Integrated Diffusive Zonal Flux of Potential Temperature', T_flux_units)
  CS%id_Tdiffy_2d = register_diag_field('ocean_model', 'T_diffy_2d', G%axesv1, Time, &
      'Vertically Integrated Diffusive Meridional Flux of Potential Temperature', T_flux_units)
  if (CS%id_Tadx_2d > 0)   call safe_alloc_ptr(CS%T_adx_2d,Isdq,Iedq,jsd,jed)
  if (CS%id_Tady_2d > 0)   call safe_alloc_ptr(CS%T_ady_2d,isd,ied,Jsdq,Jedq)
  if (CS%id_Tdiffx_2d > 0) call safe_alloc_ptr(CS%T_diffx_2d,Isdq,Iedq,jsd,jed)
  if (CS%id_Tdiffy_2d > 0) call safe_alloc_ptr(CS%T_diffy_2d,isd,ied,Jsdq,Jedq)

  CS%id_Sadx_2d = register_diag_field('ocean_model', 'S_adx_2d', G%axesu1, Time, &
      'Vertically Integrated Advective Zonal Flux of Salinity', S_flux_units)
  CS%id_Sady_2d = register_diag_field('ocean_model', 'S_ady_2d', G%axesv1, Time, &
      'Vertically Integrated Advective Meridional Flux of Salinity', S_flux_units)
  CS%id_Sdiffx_2d = register_diag_field('ocean_model', 'S_diffx_2d', G%axesu1, Time, &
      'Vertically Integrated Diffusive Zonal Flux of Salinity', S_flux_units)
  CS%id_Sdiffy_2d = register_diag_field('ocean_model', 'S_diffy_2d', G%axesv1, Time, &
      'Vertically Integrated Diffusive Meridional Flux of Salinity', S_flux_units)
  if (CS%id_Sadx_2d > 0)   call safe_alloc_ptr(CS%S_adx_2d,Isdq,Iedq,jsd,jed)
  if (CS%id_Sady_2d > 0)   call safe_alloc_ptr(CS%S_ady_2d,isd,ied,Jsdq,Jedq)
  if (CS%id_Sdiffx_2d > 0) call safe_alloc_ptr(CS%S_diffx_2d,Isdq,Iedq,jsd,jed)
  if (CS%id_Sdiffy_2d > 0) call safe_alloc_ptr(CS%S_diffy_2d,isd,ied,Jsdq,Jedq)

  CS%id_uh = register_diag_field('ocean_model', 'uh', G%axesul, Time, &
      'Zonal Thickness Flux', flux_units)
  CS%id_vh = register_diag_field('ocean_model', 'vh', G%axesvl, Time, &
      'Meridional Thickness Flux', flux_units)
  CS%id_uav = register_diag_field('ocean_model', 'uav', G%axesul, Time, &
      'Barotropic-step Averaged Zonal Velocity', 'meter second-1')
  CS%id_vav = register_diag_field('ocean_model', 'vav', G%axesvl, Time, &
      'Barotropic-step Averaged Meridional Velocity', 'meter second-1')

  if (CS%debug_truncations) then
    if (CS%split .and. CS%flux_BT_coupling) then
      call safe_alloc_ptr(CS%diag%du_adj,Isdq,Iedq,jsd,jed,nz)
      call safe_alloc_ptr(CS%diag%dv_adj,isd,ied,Jsdq,Jedq,nz)
    endif
    if (CS%split .and. CS%flux_BT_coupling .and. CS%readjust_BT_trans) then
      call safe_alloc_ptr(CS%diag%du_adj2,Isdq,Iedq,jsd,jed,nz)
      call safe_alloc_ptr(CS%diag%dv_adj2,isd,ied,Jsdq,Jedq,nz)
    endif
    call safe_alloc_ptr(CS%diag%du_dt_visc,Isdq,Iedq,jsd,jed,nz)
    call safe_alloc_ptr(CS%diag%dv_dt_visc,isd,ied,Jsdq,Jedq,nz)
    if (.not.CS%adiabatic) then
      call safe_alloc_ptr(CS%diag%du_dt_dia,Isdq,Iedq,jsd,jed,nz)
      call safe_alloc_ptr(CS%diag%dv_dt_dia,isd,ied,Jsdq,Jedq,nz)
    endif
  endif

  if (CS%split .and. CS%flux_BT_coupling) then
    CS%id_du_adj = register_diag_field('ocean_model', 'du_adj', G%axesuL, Time, &
        'Zonal velocity Adjustment 1', 'meter second-1')
    CS%id_dv_adj = register_diag_field('ocean_model', 'dv_adj', G%axesvL, Time, &
        'Meridional velocity Adjustment 1', 'meter second-1')
    if (CS%id_du_adj > 0) call safe_alloc_ptr(CS%diag%du_adj,Isdq,Iedq,jsd,jed,nz)
    if (CS%id_dv_adj > 0) call safe_alloc_ptr(CS%diag%dv_adj,isd,ied,Jsdq,Jedq,nz)
    if (CS%readjust_BT_trans) then
      CS%id_du_adj2 = register_diag_field('ocean_model', 'du_adj2', G%axesuL, Time, &
          'Zonal velocity Adjustment 2', 'meter second-1')
      CS%id_dv_adj2 = register_diag_field('ocean_model', 'dv_adj2', G%axesvL, Time, &
          'Meridional velocity Adjustment 2', 'meter second-1')
      if (CS%id_du_adj2 > 0) call safe_alloc_ptr(CS%diag%du_adj2,Isdq,Iedq,jsd,jed,nz)
      if (CS%id_dv_adj2 > 0) call safe_alloc_ptr(CS%diag%dv_adj2,isd,ied,Jsdq,Jedq,nz)
    endif
    if (CS%BT_include_udhdt) then
      CS%id_h_dudt = register_diag_field('ocean_model', 'BT_u_dhdt', G%axesu1, Time, &
          'Barotropic zonal transport tendancy from sum(h du_dt)', 'meter3 second-2')
      CS%id_h_dvdt = register_diag_field('ocean_model', 'BT_v_dhdt', G%axesv1, Time, &
          'Barotropic meridional transport tendancy from sum(h du_dt)', 'meter3 second-2')
    endif
  endif


  CS%id_u_predia = register_diag_field('ocean_model', 'u_predia', G%axesuL, Time, &
      'Zonal velocity', 'meter second-1')
  CS%id_v_predia = register_diag_field('ocean_model', 'v_predia', G%axesvL, Time, &
      'Meridional velocity', 'meter second-1')
  CS%id_h_predia = register_diag_field('ocean_model', 'h_predia', G%axeshL, Time, &
      'Layer Thickness', thickness_units)
  CS%id_e_predia = register_diag_field('ocean_model', 'e_predia', G%axeshi, Time, &
      'Interface Heights', 'meter')
  if (CS%use_temperature) then
    CS%id_T_predia = register_diag_field('ocean_model', 'temp_predia', G%axeshL, Time, &
        'Potential Temperature', 'Celsius')
    CS%id_S_predia = register_diag_field('ocean_model', 'salt_predia', G%axeshL, Time, &
        'Salinity', 'PSU')
  endif


end subroutine register_diags

subroutine GOLD_timing_init(CS)
  type(GOLD_control_struct), intent(in)    :: CS
! Arguments: CS - The control structure set up by initialize_GOLD.
  ! This subroutine sets up clock IDs for timing various subroutines.

  id_clock_ocean = cpu_clock_id('Ocean', grain=CLOCK_COMPONENT)
  id_clock_dynamics = cpu_clock_id('Ocean dynamics', grain=CLOCK_SUBCOMPONENT)
  id_clock_thermo = cpu_clock_id('Ocean thermodynamics and tracers', grain=CLOCK_SUBCOMPONENT)
  id_clock_tracer = cpu_clock_id('(Ocean tracer advection)', grain=CLOCK_MODULE_DRIVER)
  if (.not.CS%adiabatic) &
    id_clock_diabatic = cpu_clock_id('(Ocean diabatic driver)', grain=CLOCK_MODULE_DRIVER)

  id_clock_Cor = cpu_clock_id('(Ocean Coriolis & mom advection)', grain=CLOCK_MODULE)
  id_clock_continuity = cpu_clock_id('(Ocean continuity equation)', grain=CLOCK_MODULE)
  id_clock_pres = cpu_clock_id('(Ocean pressure force)', grain=CLOCK_MODULE)
  id_clock_vertvisc = cpu_clock_id('(Ocean vertical viscosity)', grain=CLOCK_MODULE)
  id_clock_horvisc = cpu_clock_id('(Ocean horizontal viscosity)', grain=CLOCK_MODULE)
  id_clock_mom_update = cpu_clock_id('(Ocean momentum increments)', grain=CLOCK_MODULE)
  id_clock_pass = cpu_clock_id('(Ocean message passing)', grain=CLOCK_MODULE)
  id_clock_GOLD_init = cpu_clock_id('(Ocean GOLD_initialize)', grain=CLOCK_MODULE)
  id_clock_pass_init = cpu_clock_id('(Ocean init message passing)', grain=CLOCK_ROUTINE)
  if (CS%split) then
    id_clock_btcalc = cpu_clock_id('(Ocean barotropic mode calc)', grain=CLOCK_MODULE)
    id_clock_btstep = cpu_clock_id('(Ocean barotropic mode stepping)', grain=CLOCK_MODULE)
    id_clock_btforce = cpu_clock_id('(Ocean barotropic forcing calc)', grain=CLOCK_MODULE)
  endif
  if (CS%thickness_diffuse) &
    id_clock_thick_diff = cpu_clock_id('(Ocean thickness diffusion)', grain=CLOCK_MODULE)
  if (CS%mixedlayer_restrat) &
    id_clock_ml_restrat = cpu_clock_id('(Ocean mixed layer restrat)', grain=CLOCK_MODULE)
  id_clock_diagnostics = cpu_clock_id('(Ocean collective diagnostics)', grain=CLOCK_MODULE)
  id_clock_Z_diag = cpu_clock_id('(Ocean Z-space diagnostics)', grain=CLOCK_MODULE)

end subroutine GOLD_timing_init

subroutine write_static_fields(G, diag)
  type(ocean_grid_type),   intent(in) :: G
  type(diag_ptrs), target, intent(in) :: diag
!   This subroutine offers the static fields in the ocean grid type
! for output via the diag_manager.
! Arguments: G - The ocean's grid structure.  Effectively intent in.
!  (in)      diag - A structure containing pointers to common diagnostic fields.

  ! The out_X arrays are needed because some of the elements of the grid
  ! type may be reduced rank macros.
  real :: out_h(SZI_(G),SZJ_(G))
  integer :: id, i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  out_h(:,:) = 0.0

  id = register_static_field('ocean_model', 'geolat', G%axesh1, &
        'Latitude of tracer (h) points', 'degrees_N')
  if (id > 0) call post_data(id, G%geolath, diag, .true.)

  id = register_static_field('ocean_model', 'geolon', G%axesh1, &
        'Longitude of tracer (h) points', 'degrees_E')
  if (id > 0) call post_data(id, G%geolonh, diag, .true.)

  id = register_static_field('ocean_model', 'geolat_c', G%axesq1, &
        'Latitude of corner (q) points', 'degrees_N')
  if (id > 0) call post_data(id, G%geolatq, diag, .true.)

  id = register_static_field('ocean_model', 'geolon_c', G%axesq1, &
        'Longitude of corner (q) points', 'degrees_E')
  if (id > 0) call post_data(id, G%geolonq, diag, .true.)

  id = register_static_field('ocean_model', 'geolat_v', G%axesv1, &
        'Latitude of meridional velocity (v) points', 'degrees_N')
  if (id > 0) call post_data(id, G%geolatv, diag, .true.)

  id = register_static_field('ocean_model', 'geolon_v', G%axesv1, &
        'Longitude of meridional velocity (v) points', 'degrees_E')
  if (id > 0) call post_data(id, G%geolonv, diag, .true.)

  id = register_static_field('ocean_model', 'geolat_u', G%axesu1, &
        'Latitude of zonal velocity (u) points', 'degrees_N')
  if (id > 0) call post_data(id, G%geolatu, diag, .true.)

  id = register_static_field('ocean_model', 'geolon_u', G%axesu1, &
        'Longitude of zonal velocity (u) points', 'degrees_E')
  if (id > 0) call post_data(id, G%geolonu, diag, .true.)

  id = register_static_field('ocean_model', 'area_t', G%axesh1, &
        'Surface area of tracer (h) cells', 'degrees_E')
  if (id > 0) then
    do j=js,je ; do i=is,ie ; out_h(i,j) = G%DXDYh(i,j) ; enddo ; enddo
    call post_data(id, out_h, diag, .true.)
  endif

  id = register_static_field('ocean_model', 'depth_ocean', G%axesh1, &
        'Depth of the ocean at tracer points', 'm', &
        standard_name='sea_floor_depth_below_geoid')
  if (id > 0) call post_data(id, G%D, diag, .true.)

  id = register_static_field('ocean_model', 'wet', G%axesh1, &
        '0 if land, 1 if ocean at tracer points', 'none')
  if (id > 0) call post_data(id, G%hmask, diag, .true.)

  id = register_static_field('ocean_model', 'wet_c', G%axesq1, &
        '0 if land, 1 if ocean at corner (q) points', 'none')
  if (id > 0) call post_data(id, G%qmask, diag, .true.)

  id = register_static_field('ocean_model', 'wet_u', G%axesu1, &
        '0 if land, 1 if ocean at zonal velocity (u) points', 'none')
  if (id > 0) call post_data(id, G%umask, diag, .true.)

  id = register_static_field('ocean_model', 'wet_v', G%axesv1, &
        '0 if land, 1 if ocean at meridional velocity (v) points', 'none')
  if (id > 0) call post_data(id, G%vmask, diag, .true.)

  id = register_static_field('ocean_model', 'Coriolis', G%axesq1, &
        'Coriolis parameter at corner (q) points', 's-1')
  if (id > 0) call post_data(id, G%f, diag, .true.)

end subroutine write_static_fields

subroutine set_restart_fields(grid, param_file, CS)
  type(ocean_grid_type),     intent(in) :: grid
  type(param_file_type),     intent(in) :: param_file
  type(GOLD_control_struct), intent(in) :: CS
!   Set the fields that are needed for bitwise identical restarting
! the time stepping scheme.  In addition to those specified here
! directly, there may be fields related to the forcing or to the
! barotropic solver that are needed; these are specified in sub-
! routines that are called from this one.
!   This routine should be altered if there are any changes to the
! time stepping scheme.  The CHECK_RESTART facility may be used to
! confirm that all needed restart fields have been included.
!
! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      CS - The control structure set up by initialize_GOLD.
  type(vardesc) :: vd
  character(len=48) :: thickness_units, flux_units

  thickness_units = get_thickness_units(grid)
  flux_units = get_flux_units(grid)

  vd = vardesc("u","Zonal velocity",'u','L','s',"meter second-1", 'd')
  call register_restart_field(CS%u(:,:,:,1), CS%u(:,:,:,2), vd, .true., CS%restart_CSp)

  vd = vardesc("v","Meridional velocity",'v','L','s',"meter second-1", 'd')
  call register_restart_field(CS%v(:,:,:,1), CS%v(:,:,:,2), vd, .true., CS%restart_CSp)

  vd = vardesc("h","Layer Thickness",'h','L','s',thickness_units, 'd')
  call register_restart_field(CS%h(:,:,:,1), CS%h(:,:,:,2), vd, .true., CS%restart_CSp)

  if (CS%split) then
   ! if (G%Boussinesq) then
      vd = vardesc("sfc","Free surface Height",'h','1','s',thickness_units, 'd')
   ! else
   !   vd(1) = vardesc("ocean_mass?","Ocean column mass",'h','1','s',"kg meter-2", 'd')
   ! endif
    call register_restart_field(CS%eta(:,:,1), CS%eta(:,:,2), vd, .false., CS%restart_CSp)

    vd = vardesc("u2","Auxiliary Zonal velocity",'u','L','s',"meter second-1", 'd')
    call register_restart_field(CS%u_av, CS%u_av, vd, .false., CS%restart_CSp)

    vd = vardesc("v2","Auxiliary Meridional velocity",'v','L','s',"meter second-1", 'd')
    call register_restart_field(CS%v_av, CS%v_av, vd, .false., CS%restart_CSp)

    vd = vardesc("h2","Auxiliary Layer Thickness",'h','L','s',thickness_units, 'd')
    call register_restart_field(CS%h_av, CS%h_av, vd, .false., CS%restart_CSp)

    vd = vardesc("uh","Zonal thickness flux",'u','L','s',flux_units, 'd')
    call register_restart_field(CS%uh, CS%uh, vd, .false., CS%restart_CSp)

    vd = vardesc("vh","Meridional thickness flux",'v','L','s',flux_units, 'd')
    call register_restart_field(CS%vh, CS%vh, vd, .false., CS%restart_CSp)

    vd = vardesc("diffu","Zonal horizontal viscous acceleration",'u','L','s', &
                 "meter second-2", 'd')
    call register_restart_field(CS%diffu, CS%diffu, vd, .false., CS%restart_CSp)

    vd = vardesc("diffv","Meridional horizontal viscous acceleration",'v','L','s',&
                 "meter second-2", 'd')
    call register_restart_field(CS%diffv, CS%diffv, vd, .false., CS%restart_CSp)

    call register_barotropic_restarts(grid, param_file, CS%barotropic_CSp, &
                                      CS%restart_CSp)

    if (CS%readjust_bt_trans) then
      vd = vardesc("uhbt_in","Final instantaneous barotropic zonal thickness flux",'u','1','s',flux_units, 'd')
      call register_restart_field(CS%uhbt_in, CS%uhbt_in, vd, .false., CS%restart_CSp)

      vd = vardesc("vhbt_in","Final instantaneous barotropic meridional thickness flux",'v','1','s',flux_units, 'd')
      call register_restart_field(CS%vhbt_in, CS%vhbt_in, vd, .false., CS%restart_CSp)
    endif
  endif

  if (CS%bulkmixedlayer .and. .not.associated(CS%tv%eqn_of_state)) then
    vd = vardesc("Rml","Mixed Layer Potential Density",'h','L','s',"kg meter-3", 'd')
    call register_restart_field(CS%tv%Rml, CS%tv%Rml, vd, .true., CS%restart_CSp)
  endif

  if (CS%use_temperature) then
    vd = vardesc("Temp","Potential Temperature",'h','L','s',"degC", 'd')
    call register_restart_field(CS%tv%T, CS%tv%T, vd, .true., CS%restart_CSp)

    vd = vardesc("Salt","Salinity",'h','L','s',"PSU", 'd')
    call register_restart_field(CS%tv%S, CS%tv%S, vd, .true., CS%restart_CSp)
  endif

  if (CS%use_frazil) then
    vd = vardesc("frazil","Frazil heat flux into ocean",'h','1','s',"J m-2", 'd')
    call register_restart_field(CS%tv%frazil, CS%tv%frazil, vd, .false., CS%restart_CSp)
  endif

  if (CS%interp_p_surf) then
    vd = vardesc("p_surf_prev","Previous ocean surface pressure",'h','1','s',"Pa", 'd')
    call register_restart_field(CS%p_surf_prev, CS%p_surf_prev, vd, .false., CS%restart_CSp)
  endif

  vd = vardesc("ave_ssh","Time average sea surface height",'h','1','s',"meter", 'd')
  call register_restart_field(CS%ave_ssh, CS%ave_ssh, vd, .false., CS%restart_CSp)

end subroutine set_restart_fields


subroutine calculate_surface_state(state, u, v, h, ssh, G, CS, p_atm)
  type(surface),                               intent(inout) :: state
  real, target, dimension(NXMEMQ_,NYMEM_,NZ_), intent(in)    :: u
  real, target, dimension(NXMEM_,NYMEMQ_,NZ_), intent(in)    :: v
  real, target, dimension(NXMEM_,NYMEM_,NZ_),  intent(in)    :: h
  real, target, dimension(NXMEM_,NYMEM_),      intent(inout) :: ssh
  type(ocean_grid_type),                       intent(inout) :: G
  type(GOLD_control_struct),                   intent(in)    :: CS
  real, optional, pointer, dimension(:,:)                    :: p_atm
!   This subroutine sets the surface (return) properties of the ocean
! model by setting the appropriate pointers in state.  Unused fields
! are set to NULL.
!
! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m.
!  (in/out)  ssh - Time mean sea surface hieght, in m.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure set up by initialize_GOLD.
!  (in)      p_atm - The atmospheric pressure, in Pa.
!  (out)     state - A structure containing fields that describe the
!                    surface state of the ocean.
  real :: depth(SZ1_(u))    ! The distance from the surface, in m.
  real :: depth_ml          ! The depth over which to average to
                            ! determine mixed layer properties, m.
  real :: dh                ! The thickness of a layer that is within the
                            ! mixed layer, in m.
  real :: mass              ! The mass per unit area of a layer, in kg m-2.
  real :: IgR0
  integer :: i, j, k, is, ie, js, je, nz, num_pnts, num_errs
  integer :: isd, ied, jsd, jed
  character(128) :: msg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  state%sea_lev => ssh

  if (present(p_atm)) then ; if (ASSOCIATED(p_atm)) then
    IgR0 = 1.0 / (G%Rho0 * G%g_Earth)
    do j=js,je ; do i=is,ie
      ssh(i,j) = ssh(i,j) + p_atm(i,j) * IgR0
    enddo ; enddo
  endif ; endif

  if (CS%smooth_ssh_passes > 0.0) then
    call smooth_SSH(ssh, G, CS%smooth_ssh_passes)
  endif

  if (CS%bulkmixedlayer) then
    if (CS%use_temperature) then
      state%SST => CS%tv%T(:,:,1)
      state%SSS => CS%tv%S(:,:,1)
      nullify(state%Rml)
    else
      nullify(state%SST)
      nullify(state%SSS)
      state%Rml => CS%tv%Rml(:,:,1)
    endif
    if (associated(CS%tv%Hml)) state%Hml => CS%tv%Hml
  else
    if (CS%use_temperature) then
      if (.not.associated(state%SST)) then
        allocate(state%SST(isd:ied,jsd:jed)) ; state%SST(:,:) = 0.0
      endif
      if (.not.associated(state%SSS)) then
        allocate(state%SSS(isd:ied,jsd:jed)) ; state%SSS(:,:) = 0.0
      endif
      nullify(state%Rml)
    else
      if (.not.associated(state%Rml)) then
        allocate(state%Rml(isd:ied,jsd:jed)) ; state%Rml(:,:) = 0.0
      endif
      nullify(state%SST) ; nullify(state%SSS)
    endif
    if (.not.associated(state%Hml)) allocate(state%Hml(isd:ied,jsd:jed))

    depth_ml = CS%Hmix
  !   Determine the mean properties of the uppermost depth_ml fluid.
    do j=js,je
      do i=is,ie
        depth(i) = 0.0
        if (CS%use_temperature) then
          state%SST(i,j) = 0.0 ; state%SSS(i,j) = 0.0
        else
          state%Rml(i,j) = 0.0
        endif
      enddo

      do k=1,nz ; do i=is,ie
        if (depth(i) + h(i,j,k) < depth_ml) then
          dh = h(i,j,k)
        elseif (depth(i) < depth_ml) then
          dh = depth_ml - depth(i)
        else
          dh = 0.0
        endif
        if (CS%use_temperature) then
          state%SST(i,j) = state%SST(i,j) + dh * CS%tv%T(i,j,k)
          state%SSS(i,j) = state%SSS(i,j) + dh * CS%tv%S(i,j,k)
        else
          state%Rml(i,j) = state%Rml(i,j) + dh * G%Rlay(k)
        endif
        depth(i) = depth(i) + dh
      enddo ; enddo
  ! Calculate the average properties of the mixed layer depth.
      do i=is,ie
        if (depth(i) < G%H_subroundoff) depth(i) = G%H_subroundoff
        if (CS%use_temperature) then
          state%SST(i,j) = state%SST(i,j) / depth(i)
          state%SSS(i,j) = state%SSS(i,j) / depth(i)
        else
          state%Rml(i,j) = state%Rml(i,j) / depth(i)
        endif
        state%Hml(:,:) = depth(i)
      enddo
    enddo ! end of j loop
  endif                                             ! end BULKMIXEDLAYER

  state%u => u(:,:,1)
  state%v => v(:,:,1)
  state%frazil => CS%tv%frazil
  state%TempxPmE => CS%tv%TempxPmE
  state%internal_heat => CS%tv%internal_heat

  if (associated(state%salt_deficit) .and. associated(CS%tv%salt_deficit)) then
    do j=js,je ; do i=is,ie
      ! Convert from gSalt to kgSalt
      state%salt_deficit(i,j) = 1000.0 * CS%tv%salt_deficit(i,j)
    enddo ; enddo
  endif

  ! Allocate structures for ocean_mass, ocean_heat, and ocean_salt.  This could
  ! be wrapped in a run-time flag to disable it for economy, since the 3-d
  ! sums are not negligible.
  if (.not.associated(state%ocean_mass)) then
    allocate(state%ocean_mass(isd:ied,jsd:jed)) ; state%ocean_mass(:,:) = 0.0
  endif
  if (CS%use_temperature) then
    if (.not.associated(state%ocean_heat)) then
      allocate(state%ocean_heat(isd:ied,jsd:jed)) ; state%ocean_heat(:,:) = 0.0
    endif
    if (.not.associated(state%ocean_salt)) then
      allocate(state%ocean_salt(isd:ied,jsd:jed)) ; state%ocean_salt(:,:) = 0.0
    endif
  endif

  if (associated(state%ocean_mass) .and. associated(state%ocean_heat) .and. &
      associated(state%ocean_salt)) then
    do j=js,je ; do i=is,ie
      state%ocean_mass(i,j) = 0.0
      state%ocean_heat(i,j) = 0.0 ; state%ocean_salt(i,j) = 0.0
    enddo ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      mass = G%H_to_kg_m2*h(i,j,k)
      state%ocean_mass(i,j) = state%ocean_mass(i,j) + mass
      state%ocean_heat(i,j) = state%ocean_heat(i,j) + mass*CS%tv%T(i,j,k)
      state%ocean_salt(i,j) = state%ocean_salt(i,j) + &
                              mass * (1.0e-3*CS%tv%S(i,j,k))
    enddo ; enddo ; enddo
  else
    if (associated(state%ocean_mass)) then
      do j=js,je ; do i=is,ie ; state%ocean_mass(i,j) = 0.0 ; enddo ; enddo
      do k=1,nz ; do j=js,je ; do i=is,ie
        state%ocean_mass(i,j) = state%ocean_mass(i,j) + G%H_to_kg_m2*h(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (associated(state%ocean_heat)) then
      do j=js,je ; do i=is,ie ; state%ocean_heat(i,j) = 0.0 ; enddo ; enddo
      do k=1,nz ; do j=js,je ; do i=is,ie
        mass = G%H_to_kg_m2*h(i,j,k)
        state%ocean_heat(i,j) = state%ocean_heat(i,j) + mass*CS%tv%T(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (associated(state%ocean_salt)) then
      do j=js,je ; do i=is,ie ; state%ocean_salt(i,j) = 0.0 ; enddo ; enddo
      do k=1,nz ; do j=js,je ; do i=is,ie
        mass = G%H_to_kg_m2*h(i,j,k)
        state%ocean_salt(i,j) = state%ocean_salt(i,j) + &
                                mass * (1.0e-3*CS%tv%S(i,j,k))
      enddo ; enddo ; enddo
    endif
  endif

  if (associated(CS%visc%taux_shelf)) state%taux_shelf => CS%visc%taux_shelf
  if (associated(CS%visc%tauy_shelf)) state%tauy_shelf => CS%visc%tauy_shelf

  if (associated(CS%tracer_flow_CSp)) then
    if (.not.associated(state%tr_fields)) allocate(state%tr_fields)
    call call_tracer_surface_state(state, h, G, CS%tracer_flow_CSp)
  endif

  if (CS%check_bad_surface_vals) then
    num_errs=0 ! count number of errors
    num_pnts=0 ! count number of errors
    do j=js,je; do i=is,ie
      if (num_errs>99) exit ! If things get really bad, stop filling up the tty
      k=0 ! Num errors at this point
      if (G%hmask(i,j)>0.) then
        if (state%sea_lev(i,j)<=-G%D(i,j)) then
          k=k+1
          write(msg(1:128),'(2(a,i4,x),2(a,f8.3,x),a,es)') &
             'Sea level < bathymetry at i=',i,'j=',j,'x=',G%geolonh(i,j),'y=',G%geolath(i,j),'SSH=',state%sea_lev(i,j)
          call GOLD_error(WARNING, trim(msg))
        endif
        if (state%sea_lev(i,j)>=CS%bad_val_ssh_max) then
          k=k+1
          write(msg(1:128),'(2(a,i4,x),2(a,f8.3,x),a,es)') &
             'Very high sea level at i=',i,'j=',j,'x=',G%geolonh(i,j),'y=',G%geolath(i,j),'SSH=',state%sea_lev(i,j)
          call GOLD_error(WARNING, trim(msg))
        endif
        if (CS%use_temperature) then
          if (state%SSS(i,j)<0.) then
            k=k+1
            write(msg(1:128),'(2(a,i4,x),2(a,f8.3,x),a,es)') &
               'Negative salinity at i=',i,'j=',j,'x=',G%geolonh(i,j),'y=',G%geolath(i,j),'SSS=',state%SSS(i,j)
            call GOLD_error(WARNING, trim(msg))
          elseif (state%SSS(i,j)>=CS%bad_val_sss_max) then
            k=k+1
            write(msg(1:128),'(2(a,i4,x),2(a,f8.3,x),a,es)') &
               'Very high salinity at i=',i,'j=',j,'x=',G%geolonh(i,j),'y=',G%geolath(i,j),'SSS=',state%SSS(i,j)
            call GOLD_error(WARNING, trim(msg))
          endif
          if (state%SST(i,j)<CS%bad_val_sst_min) then
            k=k+1
            write(msg(1:128),'(2(a,i4,x),2(a,f8.3,x),a,es)') &
               'Very cold SST at i=',i,'j=',j,'x=',G%geolonh(i,j),'y=',G%geolath(i,j),'SST=',state%SST(i,j)
            call GOLD_error(WARNING, trim(msg))
          elseif (state%SST(i,j)>=CS%bad_val_sst_max) then
            k=k+1
            write(msg(1:128),'(2(a,i4,x),2(a,f8.3,x),a,es)') &
               'Very hot SST at i=',i,'j=',j,'x=',G%geolonh(i,j),'y=',G%geolath(i,j),'SST=',state%SST(i,j)
            call GOLD_error(WARNING, trim(msg))
          endif
        endif ! use_temperature
        num_pnts=num_pnts+min(1,k)
        num_errs=num_errs+k
      endif ! hmask
    enddo; enddo
    if (num_errs>0) then
      write(msg(1:128),'(3(a,i4,x))') 'There were',num_errs,'errors involving',num_pnts,'points'
      call GOLD_error(FATAL, trim(msg))
    endif
  endif

end subroutine calculate_surface_state

subroutine smooth_SSH(ssh, G, smooth_passes)
  real, dimension(NXMEM_,NYMEM_),      intent(inout) :: ssh
  type(ocean_grid_type),               intent(inout) :: G
  real,                                intent(in)    :: smooth_passes
!   This subroutine applies a number of 2-D smoothing passes, each of which
! applies a nominal filter with the following weights:
!         1/8
!    1/8  1/2  1/8
!         1/8  
! Arguments: ssh - Time mean sea surface hieght, in m.  (Intent inout.)
!  (in)      G - The ocean's grid structure.
!  (in)      smooth_passes - the number of smoothing passes to apply.
  real, dimension(SZIQ_(G), SZJ_(G)) :: flux_x, area_x
  real, dimension(SZI_(G), SZJQ_(G)) :: flux_y, area_y
  real :: wt
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed, isl, iel, jsl, jel, halo
  integer :: pass, tot_pass

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (smooth_passes <= 0.0) return

  tot_pass = ceiling(smooth_passes)
  wt = 0.125 * smooth_passes / real(tot_pass)
  halo = -1

  do j=jsd+1,jed-1 ; do I=isd,ied-1
    area_x(I,j) = min(G%dy_u(I,j)*G%dxu(I,j), G%dxdyh(i,j), G%dxdyh(i+1,j))
  enddo ; enddo

  do J=jsd,jed-1 ; do i=isd+1,ied-1
    area_y(i,J) = min(G%dx_v(i,J)*G%dyv(i,J), G%dxdyh(i,j), G%dxdyh(i,j+1))
  enddo ; enddo

  do pass=1,tot_pass
    if (halo < 0) then
      call pass_var(ssh, G%domain)
      halo = min(is-isd-1, ied-ie-1, js-jsd-1, jed-je-1, tot_pass-pass)
    endif
    isl = is-halo ; iel = ie+halo ; jsl = js-halo ; jel = je+halo

    do j=jsl,jel ; do I=isl-1,iel
      flux_x(I,j) =  (wt * area_x(I,j)) * (ssh(i,j) - ssh(i+1,j))
    enddo ; enddo

    do J=jsl-1,jel ; do i=isl,iel
      flux_y(i,J) =  (wt * area_y(i,J)) * (ssh(i,j) - ssh(i,j+1))
    enddo ; enddo
  
    do j=jsl,jel ; do i=isl,iel
      ssh(i,j) = ssh(i,j) + ((flux_x(I-1,j) - flux_x(I,j)) + &
                             (flux_y(i,J-1) - flux_y(i,J))) * G%Idxdyh(i,j)
    enddo ; enddo

    halo = halo - 1
  enddo

end subroutine smooth_SSH

subroutine GOLD_state_chksum(mesg, u, v, h, uh, vh, G, haloshift)
  character(len=*),                    intent(in) :: mesg
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(in) :: u
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(in) :: v
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(in) :: h
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(in) :: uh
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(in) :: vh
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
  integer :: is, ie, js, je, nz, hs
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  hs=1; if (present(haloshift)) hs=haloshift
  call uchksum(u, mesg//" u",G,haloshift=hs)
  call vchksum(v, mesg//" v",G,haloshift=hs)
  call hchksum(G%H_to_kg_m2*h, mesg//" h",G,haloshift=hs)
  call uchksum(G%H_to_kg_m2*uh, mesg//" uh",G,haloshift=hs)
  call vchksum(G%H_to_kg_m2*vh, mesg//" vh",G,haloshift=hs)
end subroutine GOLD_state_chksum

subroutine GOLD_thermo_chksum(mesg, tv, G, haloshift)
  character(len=*),         intent(in) :: mesg
  type(thermo_var_ptrs),    intent(in) :: tv
  type(ocean_grid_type),    intent(in) :: G
  integer, optional,        intent(in) :: haloshift
!   This subroutine writes out chksums for the model's thermodynamic state
! variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      tv - A structure containing pointers to any thermodynamic
!                 fields that are in use.
!  (in)      G - The ocean's grid structure.
  integer :: is, ie, js, je, nz, hs
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  hs=1; if (present(haloshift)) hs=haloshift

  if (associated(tv%T)) call hchksum(tv%T, mesg//" T",G,haloshift=hs)
  if (associated(tv%S)) call hchksum(tv%S, mesg//" S",G,haloshift=hs)
  if (associated(tv%Rml)) call hchksum(tv%Rml, mesg//" Rml",G,haloshift=hs)
  if (associated(tv%frazil)) call hchksum(tv%frazil, mesg//" frazil",G,haloshift=hs)
  if (associated(tv%salt_deficit)) call hchksum(tv%salt_deficit, mesg//" salt deficit",G,haloshift=hs)

end subroutine GOLD_thermo_chksum

subroutine GOLD_accel_chksum(mesg, CAu, CAv, PFu, PFv, diffu, diffv, G, pbce, &
                            u_accel_bt, v_accel_bt)
  character(len=*),                    intent(in) :: mesg
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(in) :: CAu
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(in) :: CAv
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(in) :: PFu
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(in) :: PFv
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(in) :: diffu
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(in) :: diffv
  type(ocean_grid_type),               intent(in) :: G
  real, dimension(NXMEM_,NYMEM_,NZ_),  optional, intent(in) :: pbce
  real, dimension(NXMEMQ_,NYMEM_,NZ_), optional, intent(in) :: u_accel_bt
  real, dimension(NXMEM_,NYMEMQ_,NZ_), optional, intent(in) :: v_accel_bt
!   This subroutine writes out chksums for the model's accelerations.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      CAu - Zonal acceleration due to Coriolis and momentum
!                  advection terms, in m s-2.
!  (in)      CAv - Meridional acceleration due to Coriolis and
!                  momentum advection terms, in m s-2.
!  (in)      PFu - Zonal acceleration due to pressure gradients
!                  (equal to -dM/dx) in m s-2.
!  (in)      PFv - Meridional acceleration due to pressure
!                  gradients (equal to -dM/dy) in m s-2.
!  (in)      diffu - Zonal acceleration due to convergence of the
!                    along-isopycnal stress tensor, in m s-2.
!  (in)      diffv - Meridional acceleration due to convergence of
!                    the along-isopycnal stress tensor, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      pbce - the baroclinic pressure anomaly in each layer
!                   due to free surface height anomalies, in m s-2.
!                   pbce points to a space with nz layers or NULL.
!  (in)      u_accel_bt - The zonal acceleration from terms in the barotropic
!                         solver, in m s-2.
!  (in)      v_accel_bt - The meridional acceleration from terms in the
!                         barotropic solver, in m s-2.
  integer :: is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  call uchksum(CAu, mesg//" CAu",G,haloshift=0)
  call vchksum(CAv, mesg//" CAv",G,haloshift=0)
  call uchksum(PFu, mesg//" PFu",G,haloshift=0)
  call vchksum(PFv, mesg//" PFv",G,haloshift=0)
  call uchksum(diffu, mesg//" diffu",G,haloshift=0)
  call vchksum(diffv, mesg//" diffv",G,haloshift=0)
  if (present(pbce)) &
    call hchksum(G%kg_m2_to_H*pbce, mesg//" pbce",G,haloshift=0)
  if (present(u_accel_bt)) &
    call uchksum(u_accel_bt, mesg//" u_accel_bt",G,haloshift=0)
  if (present(v_accel_bt)) &
    call vchksum(v_accel_bt, mesg//" v_accel_bt",G,haloshift=0)
end subroutine GOLD_accel_chksum

subroutine GOLD_end(CS)
  type(GOLD_control_struct), pointer      :: CS

  DEALLOC(CS%u) ; DEALLOC(CS%v) ; DEALLOC(CS%h)
  DEALLOC(CS%uh) ; DEALLOC(CS%vh)
  DEALLOC(CS%diffu) ; DEALLOC(CS%diffv)
  DEALLOC(CS%CAu) ; DEALLOC(CS%CAv)
  DEALLOC(CS%PFu) ; DEALLOC(CS%PFv)
  if (CS%use_temperature) then
    DEALLOC(CS%tv%T) ; DEALLOC(CS%tv%S)
  endif
  if (associated(CS%tv%frazil)) deallocate(CS%tv%frazil)
  if (associated(CS%tv%salt_deficit)) deallocate(CS%tv%salt_deficit)  
  if (CS%bulkmixedlayer) then ; DEALLOC(CS%tv%Rml) ; endif
  if (associated(CS%tv%Hml)) deallocate(CS%tv%Hml)

  DEALLOC(CS%uhtr) ; DEALLOC(CS%vhtr)
  if (CS%split) then
    DEALLOC(CS%eta) ; DEALLOC(CS%uhbt) ; DEALLOC(CS%vhbt)
    DEALLOC(CS%uhbt_in) ; DEALLOC(CS%vhbt_in)
    DEALLOC(CS%h_av) ; DEALLOC(CS%u_av) ; DEALLOC(CS%v_av)
    DEALLOC(CS%eta_PF) ; DEALLOC(CS%pbce)
    DEALLOC(CS%u_accel_bt) ; DEALLOC(CS%v_accel_bt)
    DEALLOC(CS%visc_rem_u) ; DEALLOC(CS%visc_rem_v)
    call dealloc_BT_cont_type(CS%BT_cont)
  endif
  DEALLOC(CS%ave_ssh)

end subroutine GOLD_end

end module GOLD
