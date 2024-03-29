module GOLD_tracer_flow_control
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
!*  By Will Cooke, April 2003                                          *
!*                                                                     *
!*    This module contains two subroutines into which calls to other   *
!*  tracer initialization (call_tracer_init_fns) and column physics    *
!*  routines (call_tracer_column_fns) can be inserted.                 *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_diag_mediator, only : time_type, diag_ptrs
use GOLD_diag_to_Z, only : diag_to_Z_CS
use GOLD_error_handler, only : GOLD_error, FATAL, WARNING
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_grid, only : ocean_grid_type
use GOLD_restart, only : GOLD_restart_CS
use GOLD_sponge, only : sponge_CS
use GOLD_tracer, only : advect_tracer_CS
use GOLD_variables, only : forcing, surface, ocean_OBC_type, thermo_var_ptrs
use GOLD_variables, only : optics_type

#include <GOLD_memory.h>

! Add references to other user-provide tracer modules here.
use USER_tracer_example, only : tracer_column_physics, USER_initialize_tracer, USER_tracer_stock
use USER_tracer_example, only : USER_register_tracer_example, USER_tracer_surface_state
use USER_tracer_example, only : USER_tracer_example_end, USER_tracer_example_CS
use DOME_tracer, only : register_DOME_tracer, initialize_DOME_tracer
use DOME_tracer, only : DOME_tracer_column_physics, DOME_tracer_surface_state
use DOME_tracer, only : DOME_tracer_end, DOME_tracer_CS
use ideal_age_example, only : register_ideal_age_tracer, initialize_ideal_age_tracer
use ideal_age_example, only : ideal_age_tracer_column_physics, ideal_age_tracer_surface_state
use ideal_age_example, only : ideal_age_stock, ideal_age_example_end, ideal_age_tracer_CS
use GOLD_OCMIP2_CFC, only : register_OCMIP2_CFC, initialize_OCMIP2_CFC
use GOLD_OCMIP2_CFC, only : OCMIP2_CFC_column_physics, OCMIP2_CFC_surface_state
use GOLD_OCMIP2_CFC, only : OCMIP2_CFC_stock, OCMIP2_CFC_end, OCMIP2_CFC_CS
use oil_tracer, only : register_oil_tracer, initialize_oil_tracer
use oil_tracer, only : oil_tracer_column_physics, oil_tracer_surface_state
use oil_tracer, only : oil_stock, oil_tracer_end, oil_tracer_CS
use advection_test_tracer, only : register_advection_test_tracer, initialize_advection_test_tracer
use advection_test_tracer, only : advection_test_tracer_column_physics, advection_test_tracer_surface_state
use advection_test_tracer, only : advection_test_tracer_end, advection_test_tracer_CS
#ifdef _USE_TOPAZ
use GOLD_OCEAN_TOPAZ_MOD, only : register_TOPAZ, initialize_TOPAZ
use GOLD_OCEAN_TOPAZ_MOD, only : TOPAZ_column_physics, TOPAZ_surface_state
use GOLD_OCEAN_TOPAZ_MOD, only : TOPAZ_end, TOPAZ_CS, get_chl_from_TOPAZ
#endif
#ifdef _USE_GENERIC_TRACER
use GOLD_generic_tracer, only : register_GOLD_generic_tracer, initialize_GOLD_generic_tracer
use GOLD_generic_tracer, only : GOLD_generic_tracer_column_physics, GOLD_generic_tracer_surface_state
use GOLD_generic_tracer, only : end_GOLD_generic_tracer, GOLD_generic_tracer_get
use GOLD_generic_tracer, only : GOLD_generic_tracer_stock, GOLD_generic_tracer_min_max, GOLD_generic_tracer_CS
#endif

implicit none ; private

public call_tracer_register, tracer_flow_control_init, call_tracer_set_forcing
public call_tracer_column_fns, call_tracer_surface_state, call_tracer_stocks
public get_chl_from_model

type, public :: tracer_flow_control_CS ; private
  logical :: use_USER_tracer_example = .false.
  logical :: use_DOME_tracer = .false.
  logical :: use_ideal_age = .false.
  logical :: use_oil = .false.
  logical :: use_advection_test_tracer = .false.
  logical :: use_OCMIP2_CFC = .false.
  logical :: use_TOPAZ = .false.
  logical :: use_GOLD_generic_tracer = .false.
  type(USER_tracer_example_CS), pointer :: USER_tracer_example_CSp => NULL()
  type(DOME_tracer_CS), pointer :: DOME_tracer_CSp => NULL()
  type(ideal_age_tracer_CS), pointer :: ideal_age_tracer_CSp => NULL()
  type(oil_tracer_CS), pointer :: oil_tracer_CSp => NULL()
  type(advection_test_tracer_CS), pointer :: advection_test_tracer_CSp => NULL()
  type(OCMIP2_CFC_CS), pointer :: OCMIP2_CFC_CSp => NULL()
#ifdef _USE_TOPAZ
  type(TOPAZ_CS), pointer :: TOPAZ_CSp => NULL()
#endif
#ifdef _USE_GENERIC_TRACER
  type(GOLD_generic_tracer_CS), pointer :: GOLD_generic_tracer_CSp => NULL()
#endif
end type tracer_flow_control_CS

contains

! The following 5 subroutines and associated definitions provide the
! machinery to register and call the subroutines that initialize
! tracers and apply vertical column processes to tracers.

subroutine call_tracer_register(G, param_file, CS, diag, tr_adv_CSp, restart_CS)
  type(ocean_grid_type),        intent(in) :: G
  type(param_file_type),        intent(in) :: param_file
  type(tracer_flow_control_CS), pointer    :: CS
  type(diag_ptrs), target,      intent(in) :: diag
  type(advect_tracer_CS),       pointer    :: tr_adv_CSp
  type(GOLD_restart_CS),        pointer    :: restart_CS
! Argument:  G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  tr_adv_CSp - A pointer that is set to point to the control structure
!                  for the tracer advection and diffusion module.
!  (in)      restart_CS - A pointer to the restart control structure.

  character(len=128) :: version = '$Id: GOLD_tracer_flow_control.F90,v 13.0.2.3.2.20 2011/10/12 14:18:40 Alistair.Adcroft Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "GOLD_tracer_flow_control" ! This module's name.

  if (associated(CS)) then
    call GOLD_error(WARNING, "call_tracer_register called with an associated "// &
                            "control structure.")
    return
  else ; allocate(CS) ; endif

  CS%use_USER_tracer_example = .false. ; CS%use_ideal_age = .false.
  CS%use_DOME_tracer = .false. ; CS%use_OCMIP2_CFC = .false.
  CS%use_TOPAZ = .false.; CS%use_GOLD_generic_tracer = .false.
  CS%use_oil = .false.; CS%use_advection_test_tracer = .false.
  call read_param(param_file,"USE_USER_TRACER_EXAMPLE",CS%use_USER_tracer_example)
  call read_param(param_file,"USE_DOME_TRACER",CS%use_DOME_tracer)
  call read_param(param_file,"USE_IDEAL_AGE_TRACER",CS%use_ideal_age)
  call read_param(param_file,"USE_OIL_TRACER",CS%use_oil)
  call read_param(param_file,"USE_ADVECTION_TEST_TRACER",CS%use_advection_test_tracer)
  call read_param(param_file,"USE_OCMIP2_CFC",CS%use_OCMIP2_CFC)
  call read_param(param_file,"USE_TOPAZ",CS%use_TOPAZ)
  call read_param(param_file,"USE_generic_tracer",CS%use_GOLD_generic_tracer)
#ifndef _USE_TOPAZ
  if (CS%use_TOPAZ) call GOLD_error(FATAL, "call_tracer_register: use_TOPAZ=.true. BUT not compiled")
#endif
#ifndef _USE_GENERIC_TRACER
  if (CS%use_GOLD_generic_tracer) call GOLD_error(FATAL, &
       "call_tracer_register: use_GOLD_generic_tracer=.true. BUT not compiled")
#endif

  ! Write all relevant parameters to the model log.
  call log_version(param_file, mod, version, tagname, "")
  call log_param(param_file, mod, "USE_USER_TRACER_EXAMPLE", &
                                CS%use_USER_tracer_example, &
                 "If true, use the USER_tracer_example tracer package.", &
                 default=.false.)
  call log_param(param_file, mod, "USE_DOME_TRACER", CS%use_DOME_tracer, &
                 "If true, use the DOME_tracer tracer package.", &
                 default=.false.)
  call log_param(param_file, mod, "USE_IDEAL_AGE_TRACER", CS%use_ideal_age, &
                 "If true, use the ideal_age_example tracer package.", &
                 default=.false.)
  call log_param(param_file, mod, "USE_OIL_TRACER", CS%use_oil, &
                 "If true, use the oil_tracer tracer package.", &
                 default=.false.)
  call log_param(param_file, mod, "USE_ADVECTION_TEST_TRACER", CS%use_advection_test_tracer, &
                 "If true, use the advection_test_tracer tracer package.", &
                 default=.false.)
  call log_param(param_file, mod, "USE_OCMIP2_CFC", CS%use_OCMIP2_CFC, &
                 "If true, use the GOLD_OCMIP2_CFC tracer package.", &
                 default=.false.)
  call log_param(param_file, mod, "USE_TOPAZ", CS%use_TOPAZ, &
                 "If true and _USE_TOPAZ is defined as a preprocessor \n"//&
                 "macro, use the GOLD-specific TOPAZ tracer package.", &
                 default=.false.)
  call log_param(param_file, mod, "USE_generic_tracer", &
                                CS%use_GOLD_generic_tracer, &
                 "If true and _USE_GENERIC_TRACER is defined as a \n"//&
                 "preprocessor macro, use the GOLD_generic_tracer packages.", &
                 default=.false.)

!    Add other user-provided calls to register tracers for restarting here. Each
!  tracer package registration call returns a logical false if it cannot be run
!  for some reason.  This then overrides the run-time selection from above.
  if (CS%use_USER_tracer_example) CS%use_USER_tracer_example = &
    USER_register_tracer_example(G, param_file, CS%USER_tracer_example_CSp, &
                                 diag, tr_adv_CSp, restart_CS)
  if (CS%use_DOME_tracer) CS%use_DOME_tracer = &
    register_DOME_tracer(G, param_file, CS%DOME_tracer_CSp, &
                         diag, tr_adv_CSp, restart_CS)
  if (CS%use_ideal_age) CS%use_ideal_age = &
    register_ideal_age_tracer(G, param_file,  CS%ideal_age_tracer_CSp, &
                              diag, tr_adv_CSp, restart_CS)
  if (CS%use_oil) CS%use_oil = &
    register_oil_tracer(G, param_file,  CS%oil_tracer_CSp, &
                              diag, tr_adv_CSp, restart_CS)
  if (CS%use_advection_test_tracer) CS%use_advection_test_tracer = &
    register_advection_test_tracer(G, param_file, CS%advection_test_tracer_CSp, &
                         diag, tr_adv_CSp, restart_CS)
  if (CS%use_OCMIP2_CFC) CS%use_OCMIP2_CFC = &
    register_OCMIP2_CFC(G, param_file,  CS%OCMIP2_CFC_CSp, &
                        diag, tr_adv_CSp, restart_CS)
#ifdef _USE_TOPAZ
  if (CS%use_TOPAZ) CS%use_TOPAZ = &
    register_TOPAZ(G, param_file,  CS%TOPAZ_CSp, diag, tr_adv_CSp, restart_CS)
#endif
#ifdef _USE_GENERIC_TRACER
  if (CS%use_GOLD_generic_tracer) CS%use_GOLD_generic_tracer = &
    register_GOLD_generic_tracer(G, param_file,  CS%GOLD_generic_tracer_CSp, &
                        diag, tr_adv_CSp, restart_CS)
#endif

end subroutine call_tracer_register

subroutine tracer_flow_control_init(restart, day, G, h, OBC, CS, sponge_CSp, &
                                    diag_to_Z_CSp)
  logical,                            intent(in) :: restart
  type(time_type), target,            intent(in) :: day
  type(ocean_grid_type),              intent(in) :: G
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in) :: h
  type(ocean_OBC_type),               pointer    :: OBC
  type(tracer_flow_control_CS),       pointer    :: CS
  type(sponge_CS),                    pointer    :: sponge_CSp
  type(diag_to_Z_CS),                 pointer    :: diag_to_Z_CSp
!   This subroutine calls all registered tracer initialization
! subroutines.

! Arguments: restart - 1 if the fields have already been read from
!                     a restart file.
!  (in)      day - Time of the start of the run.
!  (in)      G - The ocean's grid structure.
!  (in)      h - Layer thickness, in m (Boussinesq) or kg m-2 (non-Boussinesq).
!  (in)      OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in)      CS - The control structure returned by a previous call to
!                 call_tracer_register.
!  (in/out)  sponge_CSp - A pointer to the control structure for the sponges, if
!                         they are in use.  Otherwise this may be unassociated.
!  (in/out)  diag_to_Z_Csp - A pointer to the control structure for diagnostics
!                            in depth space.
  if (.not. associated(CS)) call GOLD_error(FATAL, "tracer_flow_control_init: "// &
         "Module must be initialized via call_tracer_register before it is used.")

!  Add other user-provided calls here.
  if (CS%use_USER_tracer_example) &
    call USER_initialize_tracer(restart, day, G, h, OBC, CS%USER_tracer_example_CSp, &
                                sponge_CSp, diag_to_Z_CSp)
  if (CS%use_DOME_tracer) &
    call initialize_DOME_tracer(restart, day, G, h, OBC, CS%DOME_tracer_CSp, &
                                sponge_CSp, diag_to_Z_CSp)
  if (CS%use_ideal_age) &
    call initialize_ideal_age_tracer(restart, day, G, h, OBC, CS%ideal_age_tracer_CSp, &
                                     sponge_CSp, diag_to_Z_CSp)
  if (CS%use_oil) &
    call initialize_oil_tracer(restart, day, G, h, OBC, CS%oil_tracer_CSp, &
                                     sponge_CSp, diag_to_Z_CSp)
  if (CS%use_advection_test_tracer) &
    call initialize_advection_test_tracer(restart, day, G, h, OBC, CS%advection_test_tracer_CSp, &
                                sponge_CSp, diag_to_Z_CSp)
  if (CS%use_OCMIP2_CFC) &
    call initialize_OCMIP2_CFC(restart, day, G, h, OBC, CS%OCMIP2_CFC_CSp, &
                                sponge_CSp, diag_to_Z_CSp)
#ifdef _USE_TOPAZ
  if (CS%use_TOPAZ) &
    call initialize_TOPAZ(restart, day, G, h, OBC, CS%TOPAZ_CSp, &
                                sponge_CSp, diag_to_Z_CSp)
#endif
#ifdef _USE_GENERIC_TRACER
  if (CS%use_GOLD_generic_tracer) &
    call initialize_GOLD_generic_tracer(restart, day, G, h, OBC, CS%GOLD_generic_tracer_CSp, &
                                sponge_CSp, diag_to_Z_CSp)
#endif

end subroutine tracer_flow_control_init

subroutine get_chl_from_model(Chl_array, G, CS)
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(out) :: Chl_array
  type(ocean_grid_type),              intent(in)  :: G
  type(tracer_flow_control_CS),       pointer     :: CS
! Arguments: Chl_array - The array into which the model's Chlorophyll-A
!                        concentrations in mg m-3 are to be read.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 call_tracer_register.

#ifdef _USE_TOPAZ
  if (CS%use_TOPAZ) then
    call get_chl_from_TOPAZ(Chl_array, CS%TOPAZ_CSp)
#ifdef _USE_GENERIC_TRACER
  elseif (CS%use_GOLD_generic_tracer) then
    call GOLD_generic_tracer_get('chl','field',Chl_array, CS%GOLD_generic_tracer_CSp)
#endif
  else
    call GOLD_error(FATAL, "get_chl_from_model was called in a configuration "// &
             "that is unable to provide a sensible model-based value.\n"// &
             "CS%use_GOLD_generic_tracer and CS%use_TOPAZ are false and "// &
             "no other options are currently viable.")
  endif  
#elif _USE_GENERIC_TRACER
  if (CS%use_GOLD_generic_tracer) then
    call GOLD_generic_tracer_get('chl','field',Chl_array, CS%GOLD_generic_tracer_CSp)
  else
    call GOLD_error(FATAL, "get_chl_from_model was called in a configuration "// &
             "that is unable to provide a sensible model-based value.\n"// &
             "CS%use_GOLD_generic_tracer is false and no other viable options are on.")
  endif  
#else
  call GOLD_error(FATAL, "get_chl_from_model was called in a configuration "// &
           "that is unable to provide a sensible model-based value.\n"// &
           "_USE_TOPAZ and _USE_GENERIC_TRACER are undefined and no other options "//&
           "are currently viable.")
#endif

end subroutine get_chl_from_model

subroutine call_tracer_set_forcing(state, fluxes, day_start, day_interval, G, CS)

  type(surface),                intent(inout) :: state
  type(forcing),                intent(inout) :: fluxes
  type(time_type),              intent(in)    :: day_start
  type(time_type),              intent(in)    :: day_interval
  type(ocean_grid_type),        intent(in)    :: G
  type(tracer_flow_control_CS), pointer       :: CS
!   This subroutine calls the individual tracer modules' subroutines to
! specify or read quantities related to their surface forcing.
! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (out)     fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      day_start - Start time of the fluxes.
!  (in)      day_interval - Length of time over which these fluxes
!                           will be applied.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 call_tracer_register.

  if (.not. associated(CS)) call GOLD_error(FATAL, "call_tracer_set_forcing"// &
         "Module must be initialized via call_tracer_register before it is used.")
!  if (CS%use_ideal_age) &
!    call ideal_age_tracer_set_forcing(state, fluxes, day_start, day_interval, &
!                                      G, CS%ideal_age_tracer_CSp)

end subroutine call_tracer_set_forcing

subroutine call_tracer_column_fns(h_old, h_new, ea, eb, fluxes, dt, G, tv, optics, CS)
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in) :: h_old, h_new, ea, eb
  type(forcing),                      intent(in) :: fluxes
  real,                               intent(in) :: dt
  type(ocean_grid_type),              intent(in) :: G
  type(thermo_var_ptrs),              intent(in) :: tv
  type(optics_type),                  pointer    :: optics
  type(tracer_flow_control_CS),       pointer    :: CS
!   This subroutine calls all registered tracer column physics
! subroutines.

! Arguments: h_old -  Layer thickness before entrainment, in m (Boussinesq)
!                     or kg m-2 (non-Boussinesq).
!  (in)      h_new -  Layer thickness after entrainment, in m or kg m-2.
!  (in)      ea - an array to which the amount of fluid entrained
!                 from the layer above during this call will be
!                 added, in m or kg m-2, the same as h_old.
!  (in)      eb - an array to which the amount of fluid entrained
!                 from the layer below during this call will be
!                 added, in m or kg m-2, the same as h_old.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.
!  (in)      tv - The structure containing thermodynamic variables.
!  (in)      optics - The structure containing optical properties.
!  (in)      CS - The control structure returned by a previous call to
!                 call_tracer_register.

  if (.not. associated(CS)) call GOLD_error(FATAL, "call_tracer_column_fns: "// &
         "Module must be initialized via call_tracer_register before it is used.")
! Add calls to tracer column functions here.
  if (CS%use_USER_tracer_example) &
    call tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                               G, CS%USER_tracer_example_CSp)
  if (CS%use_DOME_tracer) &
    call DOME_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                    G, CS%DOME_tracer_CSp)
  if (CS%use_ideal_age) &
    call ideal_age_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                         G, CS%ideal_age_tracer_CSp)
  if (CS%use_oil) &
    call oil_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                         G, CS%oil_tracer_CSp, tv)
  if (CS%use_advection_test_tracer) &
    call advection_test_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                    G, CS%advection_test_tracer_CSp)
  if (CS%use_OCMIP2_CFC) &
    call OCMIP2_CFC_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                   G, CS%OCMIP2_CFC_CSp)
#ifdef _USE_TOPAZ
  if (CS%use_TOPAZ) &
    call TOPAZ_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                   G, CS%TOPAZ_CSp, tv, optics)
#endif
#ifdef _USE_GENERIC_TRACER
  if (CS%use_GOLD_generic_tracer) &
    call GOLD_generic_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                   G, CS%GOLD_generic_tracer_CSp, tv, optics)
#endif


end subroutine call_tracer_column_fns


subroutine call_tracer_stocks(h, stock_values, G, CS, stock_names, stock_units, &
                              num_stocks, stock_index, got_min_max,global_min,  global_max,xgmin, ygmin, zgmin, xgmax, ygmax, zgmax)
  real, dimension(NXMEM_,NYMEM_,NZ_),       intent(in)  :: h
  real, dimension(:),                       intent(out) :: stock_values
  type(ocean_grid_type),                    intent(in)  :: G
  type(tracer_flow_control_CS),             pointer     :: CS
  character(len=*), dimension(:), optional, intent(out) :: stock_names
  character(len=*), dimension(:), optional, intent(out) :: stock_units
  integer,                        optional, intent(out) :: num_stocks
  integer,                        optional, intent(in)  :: stock_index
  logical,  dimension(:),         optional, intent(inout) :: got_min_max
  real, dimension(:),             optional, intent(out) :: global_min,  global_max
  real, dimension(:),             optional, intent(out) :: xgmin, ygmin, zgmin, xgmax, ygmax, zgmax
!   This subroutine calls all registered tracer packages to enable them to
! add to the surface state returned to the coupler. These routines are optional.

! Arguments: h - Layer thickness, in m (Boussinesq) or kg m-2 (non-Boussinesq).
!  (out)     stock_values - The integrated amounts of a tracer on the current
!                           PE, usually in kg x concentration.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 call_tracer_register.
!  (out,opt) stock_names - Diagnostic names to use for each stock.
!  (out,opt) stock_units - Units to use in the metadata for each stock.
!  (out,opt) num_stocks - The number of tracer stocks being returned.
!  (in,opt)  stock_index - The integer stock index from stocks_constans_mod of
!                          the stock to be returned.  If this is present and
!                          greater than 0, only a single stock can be returned.
  character(len=200), dimension(MAX_FIELDS) :: names, units
  character(len=200) :: set_pkg_name
  real, dimension(MAX_FIELDS) :: values
  integer :: max_ns, ns_tot, ns, index, pkg, max_pkgs, nn

  if (.not. associated(CS)) call GOLD_error(FATAL, "call_tracer_stocks: "// &
       "Module must be initialized via call_tracer_register before it is used.")

  index = -1 ; if (present(stock_index)) index = stock_index
  ns_tot = 0
  max_ns = size(stock_values)
  if (present(stock_names)) max_ns = min(max_ns,size(stock_names))
  if (present(stock_units)) max_ns = min(max_ns,size(stock_units))

!  Add other user-provided calls here.
  if (CS%use_USER_tracer_example) then
    ns = USER_tracer_stock(h, values, G, CS%USER_tracer_example_CSp, &
                           names, units, stock_index)
    call store_stocks("tracer_example", ns, names, units, values, index, stock_values, &
                       set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  endif
! if (CS%use_DOME_tracer) then
!   ns = DOME_tracer_stock(h, values, G, CS%DOME_tracer_CSp, &
!                          names, units, stock_index)
!   call store_stocks("DOME_tracer", ns, names, units, values, index, stock_values, &
!                      set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
! endif
  if (CS%use_ideal_age) then
    ns = ideal_age_stock(h, values, G, CS%ideal_age_tracer_CSp, &
                         names, units, stock_index)
    call store_stocks("ideal_age_example", ns, names, units, values, index, &
           stock_values, set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  endif
  if (CS%use_oil) then
    ns = oil_stock(h, values, G, CS%oil_tracer_CSp, &
                         names, units, stock_index)
    call store_stocks("oil_tracer", ns, names, units, values, index, &
           stock_values, set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  endif
  if (CS%use_OCMIP2_CFC) then
    ns = OCMIP2_CFC_stock(h, values, G, CS%OCMIP2_CFC_CSp, names, units, stock_index)
    call store_stocks("GOLD_OCMIP2_CFC", ns, names, units, values, index, stock_values, &
                       set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  endif
#ifdef _USE_TOPAZ
! if (CS%use_TOPAZ) then
!   ns = TOPAZ_stock(h, values, G, CS%TOPAZ_CSp, names, units, stock_index)
!   call store_stocks("GOLD_OCEAN_TOPAZ", ns, names, units, values, index, stock_values, &
!                      set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
! endif
#endif
#ifdef _USE_GENERIC_TRACER
  if (CS%use_GOLD_generic_tracer) then
    ns = GOLD_generic_tracer_stock(h, values, G, CS%GOLD_generic_tracer_CSp, &
                                   names, units, stock_index)
    call store_stocks("GOLD_generic_tracer", ns, names, units, values, index, stock_values, &
                       set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
    nn=ns_tot-ns+1
    nn=GOLD_generic_tracer_min_max(nn, got_min_max, global_min,  global_max, xgmin, ygmin, zgmin, xgmax, ygmax, zgmax ,&
                                     G, CS%GOLD_generic_tracer_CSp,names, units)   

  endif
#endif

  if (ns_tot == 0) stock_values(1) = 0.0

  if (present(num_stocks)) num_stocks = ns_tot

end subroutine call_tracer_stocks
  
subroutine store_stocks(pkg_name, ns, names, units, values, index, stock_values, &
                        set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  character(len=*),                         intent(in)    :: pkg_name
  integer,                                  intent(in)    :: ns
  character(len=*), dimension(:),           intent(in)    :: names, units
  real, dimension(:),                       intent(in)    :: values
  integer,                                  intent(in)    :: index
  real, dimension(:),                       intent(inout) :: stock_values
  character(len=*),                         intent(inout) :: set_pkg_name
  integer,                                  intent(in)    :: max_ns
  integer,                                  intent(inout) :: ns_tot
  character(len=*), dimension(:), optional, intent(inout) :: stock_names, stock_units

! This routine stores the stocks and does error handling for call_tracer_stocks.
  character(len=16) :: ind_text, ns_text, max_text
  integer :: n

  if ((index > 0) .and. (ns > 0)) then
    write(ind_text,'(i8)') index
    if (ns > 1) then
      call GOLD_error(FATAL,"Tracer package "//trim(pkg_name)//&
          " is not permitted to return more than one value when queried"//&
          " for specific stock index "//trim(adjustl(ind_text))//".")
    elseif (ns+ns_tot > 1) then
      call GOLD_error(FATAL,"Tracer packages "//trim(pkg_name)//" and "//&
          trim(set_pkg_name)//" both attempted to set values for"//&
          " specific stock index "//trim(adjustl(ind_text))//".")
    else
      set_pkg_name = pkg_name
    endif
  endif

  if (ns_tot+ns > max_ns) then
    write(ns_text,'(i8)') ns_tot+ns ; write(max_text,'(i8)') max_ns 
    call GOLD_error(FATAL,"Attempted to return more tracer stock values (at least "//&
      trim(adjustl(ns_text))//") than the size "//trim(adjustl(max_text))//&
      "of the smallest value, name, or units array.")
  endif

  do n=1,ns
    stock_values(ns_tot+n) = values(n)
    if (present(stock_names)) stock_names(ns_tot+n) = names(n)
    if (present(stock_units)) stock_units(ns_tot+n) = units(n)
  enddo
  ns_tot = ns_tot + ns

end subroutine store_stocks

subroutine call_tracer_surface_state(state, h, G, CS)
  type(surface),                      intent(inout) :: state
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in) :: h
  type(ocean_grid_type),              intent(in) :: G
  type(tracer_flow_control_CS),       pointer    :: CS
!   This subroutine calls all registered tracer packages to enable them to
! add to the surface state returned to the coupler. These routines are optional.

! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (in)      h - Layer thickness, in m (Boussinesq) or kg m-2 (non-Boussinesq).
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 call_tracer_register.

  if (.not. associated(CS)) call GOLD_error(FATAL, "call_tracer_surface_state: "// &
         "Module must be initialized via call_tracer_register before it is used.")

!  Add other user-provided calls here.
  if (CS%use_USER_tracer_example) &
    call USER_tracer_surface_state(state, h, G, CS%USER_tracer_example_CSp)
  if (CS%use_DOME_tracer) &
    call DOME_tracer_surface_state(state, h, G, CS%DOME_tracer_CSp)
  if (CS%use_ideal_age) &
    call ideal_age_tracer_surface_state(state, h, G, CS%ideal_age_tracer_CSp)
  if (CS%use_oil) &
    call oil_tracer_surface_state(state, h, G, CS%oil_tracer_CSp)
  if (CS%use_advection_test_tracer) &
    call advection_test_tracer_surface_state(state, h, G, CS%advection_test_tracer_CSp)
  if (CS%use_OCMIP2_CFC) &
    call OCMIP2_CFC_surface_state(state, h, G, CS%OCMIP2_CFC_CSp)
#ifdef _USE_TOPAZ
  if (CS%use_TOPAZ) &
    call TOPAZ_surface_state(state, h, G, CS%TOPAZ_CSp)
#endif
#ifdef _USE_GENERIC_TRACER
  if (CS%use_GOLD_generic_tracer) &
    call GOLD_generic_tracer_surface_state(state, h, G, CS%GOLD_generic_tracer_CSp)
#endif

end subroutine call_tracer_surface_state

subroutine tracer_flow_control_end(CS)
  type(tracer_flow_control_CS), pointer :: CS

  if (CS%use_USER_tracer_example) &
    call USER_tracer_example_end(CS%USER_tracer_example_CSp)
  if (CS%use_DOME_tracer) call DOME_tracer_end(CS%DOME_tracer_CSp)
  if (CS%use_ideal_age) call ideal_age_example_end(CS%ideal_age_tracer_CSp)
  if (CS%use_oil) call oil_tracer_end(CS%oil_tracer_CSp)
  if (CS%use_advection_test_tracer) call advection_test_tracer_end(CS%advection_test_tracer_CSp)
  if (CS%use_OCMIP2_CFC) call OCMIP2_CFC_end(CS%OCMIP2_CFC_CSp)
#ifdef _USE_TOPAZ
  if (CS%use_TOPAZ) call TOPAZ_end(CS%TOPAZ_CSp)
#endif
#ifdef _USE_GENERIC_TRACER
  if (CS%use_GOLD_generic_tracer) call end_GOLD_generic_tracer(CS%GOLD_generic_tracer_CSp)
#endif

  if (associated(CS)) deallocate(CS)
end subroutine tracer_flow_control_end

end module GOLD_tracer_flow_control
