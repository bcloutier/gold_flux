module GOLD_diag_mediator

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

use GOLD_error_handler, only : GOLD_error, FATAL, is_root_pe
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_file_parser, only : lowercase
use GOLD_grid, only : ocean_grid_type
use GOLD_io, only : vardesc
use GOLD_time_manager, only : time_type

use diag_manager_mod, only : diag_manager_init, diag_manager_end
use diag_manager_mod, only : send_data, diag_axis_init
use diag_manager_mod, only : register_diag_field_fms=>register_diag_field
use diag_manager_mod, only : register_static_field

implicit none ; private

#include <GOLD_memory.h>

public set_axes_info, post_data, safe_alloc_ptr, register_diag_field, time_type
public enable_averaging, disable_averaging, query_averaging_enabled
public diag_mediator_init, diag_mediator_end, set_diag_mediator_grid
public diag_mediator_close_registration
public diag_axis_init, ocean_register_diag, register_static_field

interface safe_alloc_ptr
  module procedure safe_alloc_ptr_3d_2arg, safe_alloc_ptr_2d_2arg
  module procedure safe_alloc_ptr_3d, safe_alloc_ptr_2d, safe_alloc_ptr_1d
end interface safe_alloc_ptr

interface post_data
  module procedure post_data_3d, post_data_2d
end interface post_data

!   The following data type contains pointers to diagnostic fields that might
! be shared between modules, and also to the variables that control the handling
! of model output.
type, public :: diag_ptrs
! Each of the following fields has nz+1 levels.
  real, pointer :: diapyc_vel(:,:,:) => NULL()! The net diapycnal velocity,
                                              ! in m s-1.  Has nz+1 layers.
! Each of the following fields has nz layers.
  real, pointer :: du_dt_visc(:,:,:) => NULL()! Accelerations due to vertical
  real, pointer :: dv_dt_visc(:,:,:) => NULL()! viscosity, in m s-2.
  real, pointer :: du_dt_dia(:,:,:) => NULL()! Accelerations due to diapycnal
  real, pointer :: dv_dt_dia(:,:,:) => NULL()! mixing, in m s-2.
  real, pointer :: du_adj(:,:,:) => NULL()   ! Velocity changes due to the split
  real, pointer :: dv_adj(:,:,:) => NULL()   ! adjustment at the start of a step
                                             ! in m s-1.
  real, pointer :: du_adj2(:,:,:) => NULL()  ! A second set of velocity changes
  real, pointer :: dv_adj2(:,:,:) => NULL()  ! due to the split adjustment at
                                             ! the start of a step in m s-1.

  real, pointer :: diffu(:,:,:) => NULL()    ! Accelerations due to along iso-
  real, pointer :: diffv(:,:,:) => NULL()    ! pycnal viscosity, in m s-2.
  real, pointer :: CAu(:,:,:) => NULL()      ! Coriolis and momentum advection
  real, pointer :: CAv(:,:,:) => NULL()      ! accelerations, in m s-2.
  real, pointer :: PFu(:,:,:) => NULL()      ! Accelerations due to pressure
  real, pointer :: PFv(:,:,:) => NULL()      ! forces, in m s-2.
  real, pointer :: gradKEu(:,:,:) => NULL()  ! gradKEu = - d/dx(u2), in m s-2.
  real, pointer :: gradKEv(:,:,:) => NULL()  ! gradKEv = - d/dy(u2), in m s-2.
  real, pointer :: rv_x_v(:,:,:) => NULL()   ! rv_x_v = rv * v at u, in m s-2.
  real, pointer :: rv_x_u(:,:,:) => NULL()   ! rv_x_u = rv * u at v, in m s-2.
  real, pointer :: uhGM(:,:,:) => NULL()     ! Thickness diffusion induced
  real, pointer :: vhGM(:,:,:) => NULL()     ! volume fluxes in m3 s-1.
  real, pointer :: rv(:,:,:) => NULL()       ! Relative vorticity in s-1.
  real, pointer :: q(:,:,:) => NULL()        ! Potential vorticity, s-1 m-1.
  real, pointer :: PFu_tot(:,:,:) => NULL()  ! Accelerations due to both baro-
  real, pointer :: PFv_tot(:,:,:) => NULL()  ! clinic and barotropic pressure
                                             ! gradients in m s-2.
  real, pointer :: CAu_tot(:,:,:) => NULL()  ! Accelerations due to both baro-
  real, pointer :: CAv_tot(:,:,:) => NULL()  ! clinic and barotropic Coriolis
                                             ! forces in m s-2.

  real, pointer :: Ah_h(:,:,:) => NULL()     ! Biharmonic viscosity at h or q
  real, pointer :: Ah_q(:,:,:) => NULL()     ! points in m4 s-1.
  real, pointer :: Kh_h(:,:,:) => NULL()     ! Laplacian viscosity at h or q
  real, pointer :: Kh_q(:,:,:) => NULL()     ! points in m2 s-1.
  real, pointer :: Kd(:,:,:) => NULL()       ! Diapycnal diffusivity in m2 s-1.
  real, pointer :: PFu_bc(:,:,:) => NULL()   ! Accelerations due to pressure
  real, pointer :: PFv_bc(:,:,:) => NULL()   ! gradients deriving from density
                                             ! gradients within layers, m s-2.
  real, pointer :: eta(:,:) => NULL()        ! SSH, m
  real, pointer :: bott_press(:,:) => NULL() ! bottom pressure, Pa
 
! Each of the following fields has 1 layer.
  real, pointer :: PFu_bt(:,:) => NULL()     ! Barotropic pressure gradient
  real, pointer :: PFv_bt(:,:) => NULL()     ! accelerations, in m s-2.
  real, pointer :: Coru_bt(:,:) => NULL()    ! Barotropic Coriolis accel-
  real, pointer :: Corv_bt(:,:) => NULL()    ! erations, in m s-2.
  real, pointer :: Nonlnu_bt(:,:) => NULL()  ! Barotropic nonlinear accel-
  real, pointer :: Nonlnv_bt(:,:) => NULL()  ! erations, in m s-2.
  real, pointer :: ubt_flux(:,:) => NULL()   ! Barotropic mass fluxes across
  real, pointer :: vbt_flux(:,:) => NULL()   ! cell faces, in m3 s-1.

  real, pointer :: ML_depth(:,:) => NULL()   ! The mixed layer depth in m.
! These are terms in the mixed layer TKE budget, all in m3 s-2.
  real, pointer :: TKE_wind(:,:) => NULL()   ! The wind source of TKE.
  real, pointer :: TKE_RiBulk(:,:) => NULL() ! The resolved KE source of TKE.
  real, pointer :: TKE_conv(:,:) => NULL()   ! The convective source of TKE.
  real, pointer :: TKE_pen_SW(:,:) => NULL() ! The TKE sink required to mix
                                             ! penetrating shortwave heating.
  real, pointer :: TKE_mech_decay(:,:) => NULL() ! The decay of mechanical TKE.
  real, pointer :: TKE_conv_decay(:,:) => NULL() ! The decay of convective TKE.
  real, pointer :: TKE_mixing(:,:) => NULL() ! The work done by TKE to deepen
                                             ! the mixed layer.
  real, pointer :: TKE_conv_s2(:,:) => NULL()! The convective source of TKE due to
                                             ! to mixing in sigma2.
  real, pointer :: PE_detrain(:,:) => NULL() ! The spurious source of potential
                                             ! energy due to mixed layer
                                             ! detrainment, W m-2.
  real, pointer :: PE_detrain2(:,:) => NULL()! The spurious source of potential
                                             ! energy due to mixed layer only
                                             ! detrainment, W m-2.

! The following are a number of estimates of the thickness fluxes, in m3 s-1.
  real, pointer :: uh_min(:,:,:) => NULL()
  real, pointer :: uh_max(:,:,:) => NULL()
  real, pointer :: uh_lay(:,:,:) => NULL()
  real, pointer :: uh_cent(:,:,:) => NULL()

  real, pointer :: vh_min(:,:,:) => NULL()
  real, pointer :: vh_max(:,:,:) => NULL()
  real, pointer :: vh_lay(:,:,:) => NULL()
  real, pointer :: vh_cent(:,:,:) => NULL()

! The following fields are used for the output of the data.
  integer :: is, ie, js, je
  integer :: isd, ied, jsd, jed
  real :: time_int              ! The time interval in s for any fields
                                ! that are offered for averaging.
  type(time_type) :: time_end   ! The end time of the valid
                                ! interval for any offered field.
  logical :: ave_enabled = .false. ! .true. if averaging is enabled.
end type diag_ptrs

integer :: doc_unit = -1

contains

subroutine set_axes_info(latq, lath, lonq, lonh, G, param_file, set_vertical)
  real, intent(in) :: latq(:), lath(:), lonq(:), lonh(:)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
  logical, optional,     intent(in)    :: set_vertical
! Arguments: latq - The latitude of q points in the entire domain.
!  (in)      lath - The latitude of h points in the entire domain.
!  (in)      lonq - The longitude of q points in the entire domain.
!  (in)      lonh - The longitude of h points in the entire domain.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in,opt)  set_vertical - If true (or missing), set up the vertical axes.
  integer :: id_xq, id_yq, id_zl, id_zi, id_xh, id_yh, k, nz
  real :: zlev(SZK_(G)), zinter(SZK_(G)+1)
  logical :: set_vert, Cartesian_grid
  character(len=80) :: grid_config, units_temp
  character(len=128) :: version = '$Id: GOLD_diag_mediator.F90,v 13.0.2.6.2.21 2011/07/21 17:12:42 Robert.Hallberg Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod  = "GOLD_diag_mediator" ! This module's name.
  nz = G%ke

  set_vert = .true. ; if (present(set_vertical)) set_vert = set_vertical

  G%x_axis_units = "degrees_E"
  G%y_axis_units = "degrees_N"
  call read_param(param_file,"GRID_CONFIG",grid_config,.true.)
  if (index(lowercase(trim(grid_config)),"cartesian") > 0) then
    ! This is a cartesian grid, and may have different axis units.
    Cartesian_grid = .true.
    units_temp = 'd' ; call read_param(param_file,"AXIS_UNITS",units_temp)
    if (units_temp(1:1) == 'k') then
      G%x_axis_units = "kilometers" ; G%y_axis_units = "kilometers"
    elseif (units_temp(1:1) == 'm') then
      G%x_axis_units = "meters" ; G%y_axis_units = "meters"
    endif
  else
    Cartesian_grid = .false.
  endif
  
  do k=1,nz ; zlev(k) = G%Rlay(k) ; enddo
  zinter(1) = 1.5*G%Rlay(1) - 0.5*G%Rlay(2)
  do k=2,nz ; zinter(k) = 0.5*(G%Rlay(k) + G%Rlay(k-1)) ; enddo
  zinter(nz+1) = 1.5*G%Rlay(nz) - 0.5*G%Rlay(nz-1)

!  do i=1,nz ; zlev(i) = real(i) ; enddo
!  do i=1,nz+1 ; zinter(i) = real(i) - 0.5 ; enddo
  if(G%symmetric) then 
    id_xq = diag_axis_init('xq', lonq(G%isgq:G%iegq), G%x_axis_units, 'x', &
              'q point nominal longitude', Domain2=G%Domain%mpp_domain)
    id_yq = diag_axis_init('yq', latq(G%jsgq:G%jegq), G%y_axis_units, 'y', &
              'q point nominal latitude', Domain2=G%Domain%mpp_domain)
  else
    id_xq = diag_axis_init('xq', lonq(G%isg:G%ieg), G%x_axis_units, 'x', &
              'q point nominal longitude', Domain2=G%Domain%mpp_domain)
    id_yq = diag_axis_init('yq', latq(G%jsg:G%jeg), G%y_axis_units, 'y', &
              'q point nominal latitude', Domain2=G%Domain%mpp_domain)
  endif           
  id_xh = diag_axis_init('xh', lonh(G%isg:G%ieg), G%x_axis_units, 'x', &
              'h point nominal longitude', Domain2=G%Domain%mpp_domain)
  id_yh = diag_axis_init('yh', lath(G%jsg:G%jeg), G%y_axis_units, 'y', &
              'h point nominal latitude', Domain2=G%Domain%mpp_domain)

  if (set_vert) then
    id_zl = diag_axis_init('zl', zlev, 'layer', 'z', 'cell depth')
    id_zi = diag_axis_init('zi', zinter, 'interface', 'z', 'cell interface depth')
  else
    id_zl = -1 ; id_zi = -1
  endif

  ! Vertical axes for the interfaces and layers.
  G%axeszi(1) = id_zi ; G%axeszL(1) = id_zL

  ! Axis groupings for the model layers.
  G%axeshL(:) = (/ id_xh, id_yh, id_zL /)
  G%axesqL(:) = (/ id_xq, id_yq, id_zL /)
  G%axesuL(:) = (/ id_xq, id_yh, id_zL /)
  G%axesvL(:) = (/ id_xh, id_yq, id_zL /)

  ! Axis groupings for the model interfaces.
  G%axeshi(:) = (/ id_xh, id_yh, id_zi /)
  G%axesui(:) = (/ id_xq, id_yh, id_zi /)
  G%axesvi(:) = (/ id_xh, id_yq, id_zi /)
  G%axesqi(:) = (/ id_xq, id_yq, id_zi /)

  ! Axis groupings for 2-D arrays.
  G%axesh1(:) = (/ id_xh, id_yh /)
  G%axesq1(:) = (/ id_xq, id_yq /)
  G%axesu1(:) = (/ id_xq, id_yh /)
  G%axesv1(:) = (/ id_xh, id_yq /)

  call log_version(param_file, mod, version, tagname)
  call log_param(param_file, mod, "GRID_CONFIG", grid_config, &
                 "The method for defining the horizontal grid.  Valid \n"//&
                 "entries include:\n"//&
                 "\t file - read the grid from GRID_FILE \n"//&
                 "\t mosaic - read the grid from a mosaic grid file \n"//&
                 "\t cartesian - a Cartesian grid \n"//&
                 "\t spherical - a spherical grid \n"//&
                 "\t mercator  - a Mercator grid")
  if (Cartesian_grid) &
    call log_param(param_file, mod, "AXIS_UNITS", G%x_axis_units, &
                 "The units for the x- and y- axis labels.  AXIS_UNITS \n"//&
                 "should be defined as 'k' for km, 'm' for m, or omitted \n"//&
                 "for degrees of latitude and longitude (the default). \n"//&
                 "Except on a Cartesian grid, only degrees are currently \n"//&
                 "implemented.")
 
end subroutine set_axes_info

subroutine set_diag_mediator_grid(G, diag)
  type(ocean_grid_type), intent(inout) :: G
  type(diag_ptrs),       intent(inout) :: diag
! Arguments: G - The ocean's grid structure.
!  (inout)   diag - A structure containing pointers to common diagnostic fields
!                   and some control information for diagnostics.
  diag%is = G%isc ; diag%ie = G%iec ; diag%js = G%jsc ; diag%je = G%jec
  diag%isd = G%isd ; diag%ied = G%ied ; diag%jsd = G%jsd ; diag%jed = G%jed
end subroutine set_diag_mediator_grid

subroutine post_data_2d(diag_field_id, field, diag, is_static)
  integer, intent(in)         :: diag_field_id
  real, intent(in)            :: field(:,:)
  type(diag_ptrs), intent(in) :: diag
  logical, intent(in), optional :: is_static
! Arguments: diag_field_id - the id for an output variable returned by a
!                            previous call to register_diag_field.
!  (in)      field - The 2-d array being offered for output or averaging.
!  (inout)   diag - A structure containing pointers to common diagnostic fields
!                   and some control information for diagnostics.
!  (in)      static - If true, this is a static field that is always offered.
  logical :: used, is_stat
  integer :: ishift, jshift

  is_stat = .false. ; if (present(is_static)) is_stat = is_static 

  ishift = 0 ; jshift = 0
  if ( size(field,1) == diag%ied-diag%isd +1 ) then
    ishift = 0
  elseif ( size(field,1) == diag%ied-diag%isd +2 ) then
    ishift = 1
!  else
!    call GOLD_error(FATAL,"post_data_2d: peculiar size in i-direction")
  endif
  if ( size(field,2) == diag%jed-diag%jsd +1 ) then
    jshift = 0
  elseif ( size(field,2) == diag%jed-diag%jsd +2 ) then
    jshift = 1
!  else
!    call GOLD_error(FATAL,"post_data_2d: peculiar size in j-direction")
  endif

  if (is_stat) then
    used = send_data(diag_field_id, field, &
                     is_in = diag%is-ishift, js_in = diag%js-jshift, &
                     ie_in = diag%ie, je_in = diag%je)
  elseif (diag%ave_enabled) then
    used = send_data(diag_field_id, field, diag%time_end, &
                     is_in = diag%is-ishift, js_in = diag%js-jshift, &
                     ie_in = diag%ie, je_in = diag%je, &
                     weight=diag%time_int)
  endif

end subroutine post_data_2d

subroutine post_data_3d(diag_field_id, field, diag, is_static)
  integer, intent(in)         :: diag_field_id
  real, intent(in)            :: field(:,:,:)
  type(diag_ptrs), intent(in) :: diag
  logical, intent(in), optional :: is_static
! Arguments: diag_field_id - the id for an output variable returned by a
!                            previous call to register_diag_field.
!  (in)      field - The 3-d array being offered for output or averaging.
!  (inout)   diag - A structure containing pointers to common diagnostic fields
!                   and some control information for diagnostics.
!  (in)      static - If true, this is a static field that is always offered.
  logical :: used  ! The return value of send_data is not used for anything.
  logical :: is_stat
  integer :: ishift, jshift
  is_stat = .false. ; if (present(is_static)) is_stat = is_static 
  
  ishift = 0 ; jshift = 0
  if ( size(field,1) == diag%ied-diag%isd +1 ) then
    ishift = 0
  elseif ( size(field,1) == diag%ied-diag%isd +2 ) then
    ishift = 1
!  else
!    call GOLD_error(FATAL,"post_data_3d: peculiar size in i-direction")
  endif
  if ( size(field,2) == diag%jed-diag%jsd +1 ) then
    jshift = 0
  elseif ( size(field,2) == diag%jed-diag%jsd +2 ) then
    jshift = 1  
!  else
!    call GOLD_error(FATAL,"post_data_3d: peculiar size in j-direction")
  endif

  if (is_stat) then
    used = send_data(diag_field_id, field, &
                     is_in = diag%is-ishift, js_in = diag%js-jshift, &
                     ie_in = diag%ie, je_in = diag%je)
  elseif (diag%ave_enabled) then
    used = send_data(diag_field_id, field, diag%time_end, &
                     is_in = diag%is-ishift, js_in = diag%js-jshift, &
                     ie_in = diag%ie, je_in = diag%je, &
                     weight=diag%time_int)
  endif

end subroutine post_data_3d


subroutine enable_averaging(time_int_in, time_end_in, diag)
  real, intent(in) :: time_int_in
  type(time_type), intent(in) :: time_end_in
  type(diag_ptrs), intent(inout) :: diag
! This subroutine enables the accumulation of time averages over the
! specified time interval.

! Arguments: time_int_in - the time interval in s over which any
!                          values that are offered are valid.
!  (in)      time_end_in - the end time in s of the valid interval.
!  (inout)   diag - A structure containing pointers to common diagnostic fields
!                   and some control information for diagnostics.
!  if (num_file==0) return
  diag%time_int = time_int_in
  diag%time_end = time_end_in
  diag%ave_enabled = .true.
end subroutine enable_averaging

! Call this subroutine to avoid averaging any offered fields.
subroutine disable_averaging(diag)
  type(diag_ptrs), intent(inout) :: diag
! Argument: diag - A structure containing pointers to common diagnostic fields
!                  and some control information for diagnostics.

  diag%time_int = 0.0
  diag%ave_enabled = .false.

end subroutine disable_averaging

! Call this subroutine to determine whether the averaging is
! currently enabled.  .true. is returned if it is.
function query_averaging_enabled(diag, time_int, time_end)
  type(diag_ptrs),           intent(in)  :: diag
  real,            optional, intent(out) :: time_int
  type(time_type), optional, intent(out) :: time_end
  logical :: query_averaging_enabled
! Arguments: diag - A structure containing pointers to common diagnostic fields
!                   and some control information for diagnostics.
!  (out,opt) time_int - The current setting of diag%time_int, in s.
!  (out,opt) time_end - The current setting of diag%time_end.

  if (present(time_int)) time_int = diag%time_int
  if (present(time_end)) time_end = diag%time_end
  query_averaging_enabled = diag%ave_enabled
end function query_averaging_enabled

subroutine safe_alloc_ptr_1d(ptr, i1, i2)
  real, pointer :: ptr(:)
  integer, intent(in) :: i1
  integer, optional, intent(in) :: i2
  if (.not.ASSOCIATED(ptr)) then
    if (present(i2)) then
      allocate(ptr(i1:i2))
    else
      allocate(ptr(i1))
    endif
    ptr(:) = 0.0
  endif
end subroutine safe_alloc_ptr_1d

subroutine safe_alloc_ptr_2d_2arg(ptr, ni, nj)
  real, pointer :: ptr(:,:)
  integer, intent(in) :: ni, nj
  if (.not.ASSOCIATED(ptr)) then
    allocate(ptr(ni,nj))
    ptr(:,:) = 0.0
  endif
end subroutine safe_alloc_ptr_2d_2arg

subroutine safe_alloc_ptr_3d_2arg(ptr, ni, nj, nk)
  real, pointer :: ptr(:,:,:)
  integer, intent(in) :: ni, nj, nk
  if (.not.ASSOCIATED(ptr)) then
    allocate(ptr(ni,nj,nk))
    ptr(:,:,:) = 0.0
  endif
end subroutine safe_alloc_ptr_3d_2arg

subroutine safe_alloc_ptr_2d(ptr, is, ie, js, je)
  real, pointer :: ptr(:,:)
  integer, intent(in) :: is, ie, js, je
  if (.not.ASSOCIATED(ptr)) then
    allocate(ptr(is:ie,js:je))
    ptr(:,:) = 0.0
  endif
end subroutine safe_alloc_ptr_2d

subroutine safe_alloc_ptr_3d(ptr, is, ie, js, je, nk)
  real, pointer :: ptr(:,:,:)
  integer, intent(in) :: is, ie, js, je, nk
  if (.not.ASSOCIATED(ptr)) then
    allocate(ptr(is:ie,js:je,nk))
    ptr(:,:,:) = 0.0
  endif
end subroutine safe_alloc_ptr_3d

function register_diag_field(module_name, field_name, axes, init_time, &
             long_name, units, missing_value, range, mask_variant, standard_name, &
             verbose, do_not_log, err_msg, interp_method, tile_count)
  integer :: register_diag_field
  character(len=*), intent(in) :: module_name, field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), optional, intent(in) :: long_name, units, standard_name
  real,             optional, intent(in) :: missing_value, range(2)
  logical,          optional, intent(in) :: mask_variant, verbose, do_not_log
  character(len=*), optional, intent(out):: err_msg
  character(len=*), optional, intent(in) :: interp_method
  integer,          optional, intent(in) :: tile_count
! Output:    An integer handle for a diagnostic array.
! Arguments: module_name - The name of this module, usually "ocean_model" or "ice_shelf_model".
!  (in)      field_name - The name of the diagnostic field.
!  (in)      axes - A set of up to 3 integers that indicates the axes for this field.
!  (in)      init_time - The time at which a field is first available?
!  (in,opt)  long_name - The long name of a field.
!  (in,opt)  units - The units of a field.
!  (in,opt)  standard_name - The standardized name associated with a field. (Not yet used in GOLD.)
!  (in,opt)  missing_value - A value that indicates missing values.
!  (in,opt)  range - The valid range of a variable. (Not used in GOLD.)
!  (in,opt)  mask_variant - If true a logical mask must be provided with post_data calls.  (Not used in GOLD.)
!  (in,opt)  verbose - If true, FMS is verbosed. (Not used in GOLD.)
!  (in,opt)  do_not_log - If true, do not log something. (Not used in GOLD.)
!  (out,opt) err_msg - An character string into which an error message might be placed. (Not used in GOLD.)
!  (in,opt)  interp_method - No clue. (Not used in GOLD.)
!  (in,opt)  tile_count - No clue. (Not used in GOLD.)
  character(len=240) :: mesg

  register_diag_field = register_diag_field_fms(module_name, field_name, axes, &
         init_time, long_name=long_name, units=units, missing_value=missing_value, &
         range=range, mask_variant=mask_variant, standard_name=standard_name, &
         verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
         interp_method=interp_method, tile_count=tile_count)

  if (is_root_pe() .and. doc_unit > 0) then
    if (register_diag_field > 0) then
      mesg = '"'//trim(module_name)//'", "'//trim(field_name)//'"  [Used]'
    else
      mesg = '"'//trim(module_name)//'", "'//trim(field_name)//'"  [Unused]'
    endif
    write(doc_unit, '(a)') trim(mesg)
    if (present(long_name)) call describe_option("long_name", long_name)
    if (present(units)) call describe_option("units", units)
    if (present(standard_name)) call describe_option("standard_name", standard_name)
  endif

end function register_diag_field

subroutine describe_option(opt_name, value)
  character(len=*), intent(in) :: opt_name, value

  character(len=240) :: mesg
  integer :: start_ind = 1, end_ind, len_ind

  len_ind = len_trim(value)

  mesg = "    ! "//trim(opt_name)//": "//trim(value)
  write(doc_unit, '(a)') trim(mesg)
end subroutine describe_option

function ocean_register_diag(var_desc, grid, day)
  integer :: ocean_register_diag
  type(vardesc), intent(in) :: var_desc
  type(ocean_grid_type), intent(in) :: grid
  type(time_type), intent(in) :: day
  integer, dimension(:), allocatable :: axes

  ! Use the hor_grid and z_grid components of vardesc to determine the 
  ! desired axes to register the diagnostic field for.
  select case (var_desc%z_grid)

    case ("L")
      select case (var_desc%hor_grid)
        case ("q")
          allocate(axes(3)) ; axes(:) = grid%axesqL(:)
        case ("h")
          allocate(axes(3)) ; axes(:) = grid%axeshL(:)
        case ("u")
          allocate(axes(3)) ; axes(:) = grid%axesuL(:)
        case ("v")
          allocate(axes(3)) ; axes(:) = grid%axesvL(:)
        case ("z")
          allocate(axes(1)) ; axes(:) = grid%axeszL(:)
        case default
          call GOLD_error(FATAL, "ocean_register_diag: " // &
              "unknown hor_grid component "//trim(var_desc%hor_grid))
      end select

    case ("i")
      select case (var_desc%hor_grid)
        case ("q")
          allocate(axes(3)) ; axes(:) = grid%axesqi(:)
        case ("h")
          allocate(axes(3)) ; axes(:) = grid%axeshi(:)
        case ("u")
          allocate(axes(3)) ; axes(:) = grid%axesui(:)
        case ("v")
          allocate(axes(3)) ; axes(:) = grid%axesvi(:)
        case ("z")
          allocate(axes(1)) ; axes(:) = grid%axeszi(:)
        case default
          call GOLD_error(FATAL, "ocean_register_diag: " // &
            "unknown hor_grid component "//trim(var_desc%hor_grid))
      end select

    case ("1")
      allocate(axes(2))
      select case (var_desc%hor_grid)
        case ("q")
          axes(:) = grid%axesq1(:)
        case ("h")
          axes(:) = grid%axesh1(:)
        case ("u")
          axes(:) = grid%axesu1(:)
        case ("v")
          axes(:) = grid%axesv1(:)
        case default
          call GOLD_error(FATAL, "ocean_register_diag: " // &
            "unknown hor_grid component "//trim(var_desc%hor_grid))
      end select

    case default
      call GOLD_error(FATAL,&
        "ocean_register_diag: unknown z_grid component "//trim(var_desc%z_grid))
  end select

  ocean_register_diag = register_diag_field("ocean_model", trim(var_desc%name), axes, &
        day, trim(var_desc%longname), trim(var_desc%units), missing_value = -1.0e+34)

  if (allocated(axes)) deallocate(axes)

end function ocean_register_diag

subroutine diag_mediator_init(param_file, err_msg)
  type(param_file_type),      intent(in)  :: param_file
  character(len=*), optional, intent(out) :: err_msg

  integer :: ios, unit
  logical :: opened, new_file
  character(len=240) :: doc_file
  character(len=40)  :: mod  = "GOLD_diag_mediator" ! This module's name.

  call diag_manager_init(err_msg=err_msg)

  if (is_root_pe()) then
    doc_file = " "
    call read_param(param_file,"AVAILABLE_DIAGS_FILE", doc_file)

    call log_param(param_file, mod, "AVAILABLE_DIAGS_FILE", doc_file, &
                 "A file into which to write a list of all available \n"//&
                 "ocean diagnostics that can be included in a diag_table.", &
                 default=" ")
    if (len_trim(doc_file) > 0) then
      new_file = .true. ; if (doc_unit /= -1) new_file = .false.
    ! Find an unused unit number.
      do doc_unit=512,42,-1
        inquire( doc_unit, opened=opened)
        if (.not.opened) exit
      enddo

      if (opened) call GOLD_error(FATAL, &
          "diag_mediator_init failed to find an unused unit number.")

      if (new_file) then
        open(doc_unit, file=trim(doc_file), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='REPLACE', iostat=ios)
      else ! This file is being reopened, and should be appended.
        open(doc_unit, file=trim(doc_file), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='OLD', position='APPEND', iostat=ios)
      endif
      inquire(doc_unit, opened=opened)
      if ((.not.opened) .or. (ios /= 0)) then
        call GOLD_error(FATAL, "Failed to open available diags file "//trim(doc_file)//".")
      endif
    endif
  endif

end subroutine diag_mediator_init

subroutine diag_mediator_close_registration( )

  if (doc_unit > -1) then
    close(doc_unit) ; doc_unit = -2
  endif

end subroutine diag_mediator_close_registration

subroutine diag_mediator_end(time)
  type(time_type), intent(in) :: time

  call diag_manager_end(time)

  if (doc_unit > -1) then
    close(doc_unit) ; doc_unit = -3
  endif

end subroutine diag_mediator_end

end module GOLD_diag_mediator
