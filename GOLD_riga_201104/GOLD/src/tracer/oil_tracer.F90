module oil_tracer
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
!*  By Robert Hallberg, 2002                                           *
!*                                                                     *
!*    This file contains an example of the code that is needed to set  *
!*  up and use a set (in this case two) of dynamically passive tracers *
!*  for diagnostic purposes.  The tracers here are an ideal age tracer *
!*  that ages at a rate of 1/year once it is isolated from the surface,*
!*  and a vintage tracer, whose surface concentration grows exponen-   *
!*  with time with a 30-year timescale (similar to CFCs).              *
!*                                                                     *
!*    A single subroutine is called from within each file to register  *
!*  each of the tracers for reinitialization and advection and to      *
!*  register the subroutine that initializes the tracers and set up    *
!*  their output and the subroutine that does any tracer physics or    *
!*  chemistry along with diapycnal mixing (included here because some  *
!*  tracers may float or swim vertically or dye diapycnal processes).  *
!*                                                                     *
!*                                                                     *
!*  Macros written all in capital letters are defined in GOLD_memory.h.*
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, tr                                    *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use GOLD_diag_mediator, only : diag_ptrs
use GOLD_diag_to_Z, only : register_Z_tracer, diag_to_Z_CS
use GOLD_error_handler, only : GOLD_error, FATAL, WARNING
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_grid, only : ocean_grid_type
use GOLD_io, only : file_exists, read_data, slasher, vardesc
use GOLD_restart, only : register_restart_field, query_initialized, GOLD_restart_CS
use GOLD_sponge, only : set_up_sponge_field, sponge_CS
use GOLD_time_manager, only : time_type, get_time
use GOLD_tracer, only : register_tracer, advect_tracer_CS, tracer_vertdiff
use GOLD_tracer, only : add_tracer_diagnostics, add_tracer_OBC_values
use GOLD_tracer_Z_init, only : tracer_Z_init
use GOLD_variables, only : forcing, surface, ocean_OBC_type
use GOLD_variables, only : thermo_var_ptrs
use coupler_util, only : set_coupler_values, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <GOLD_memory.h>

public register_oil_tracer, initialize_oil_tracer
public oil_tracer_column_physics, oil_tracer_surface_state
public oil_stock, oil_tracer_end

! NTR_MAX is the maximum number of tracers in this module.
integer, parameter :: NTR_MAX = 20

type p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d

type, public :: oil_tracer_CS ; private
  integer :: ntr    ! The number of tracers that are actually used.
  logical :: coupled_tracers = .false.  ! These tracers are not offered to the
                                        ! coupler.
  character(len = 200) :: IC_file ! The file in which the age-tracer initial values
                    ! can be found, or an empty string for internal initialization.
  logical :: Z_IC_file ! If true, the IC_file is in Z-space.  The default is false.
  real :: oil_source_longitude, oil_source_latitude ! Lat,lon of source location (geographic)
  integer :: oil_source_i=-999, oil_source_j=-999 ! Local i,j of source location (computational)
  real :: oil_source_rate     ! Rate of oil injection (kg/s)
  real :: oil_start_year      ! The year in which tracers start aging, or at which the
                              ! surface value equals young_val, in years.
  real :: oil_end_year        ! The year in which tracers start aging, or at which the
                              ! surface value equals young_val, in years.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(advect_tracer_CS), pointer :: tr_adv_CSp => NULL()
  real, pointer :: tr(:,:,:,:) => NULL()   ! The array of tracers used in this
                                           ! subroutine, in g m-3?
  type(p3d), dimension(NTR_MAX) :: &
    tr_adx, &! Tracer zonal advective fluxes in g m-3 m3 s-1.
    tr_ady, &! Tracer meridional advective fluxes in g m-3 m3 s-1.
    tr_dfx, &! Tracer zonal diffusive fluxes in g m-3 m3 s-1.
    tr_dfy   ! Tracer meridional diffusive fluxes in g m-3 m3 s-1.
  real, dimension(NTR_MAX) :: &
    IC_val = 0.0, &    ! The (uniform) initial condition value.
    young_val = 0.0, & ! The value assigned to tr at the surface.
    land_val = -1.0, & ! The value of tr used where land is masked out.
    sfc_growth_rate    ! The exponential growth rate for the surface value,
                       ! in units of year-1.
  real, dimension(NTR_MAX) ::  oil_decay_days, & ! Decay time scale of oil (in days)
                               oil_decay_rate    ! Decay rate of oil (in s^-1) calculated from oil_decay_days
  integer, dimension(NTR_MAX) ::  oil_source_k ! Layer of source
  logical :: mask_tracers  ! If true, oils are masked out in massless layers.
  logical :: oil_may_reinit  ! If true, oil may go through the
                           ! initialization code if they are not found in the
                           ! restart files.
  integer, dimension(NTR_MAX) :: &
    ind_tr, &  ! Indices returned by aof_set_coupler_flux if it is used and the
               ! surface tracer concentrations are to be provided to the coupler.
    id_tracer = -1, id_tr_adx = -1, id_tr_ady = -1, &
    id_tr_dfx = -1, id_tr_dfy = -1

  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
                             ! ocean diagnostic fields and control variables.
  type(GOLD_restart_CS), pointer :: restart_CSp => NULL()

  type(vardesc) :: tr_desc(NTR_MAX)
end type oil_tracer_CS

contains

function register_oil_tracer(G, param_file, CS, diag, tr_adv_CSp, &
                                   restart_CS)
  type(ocean_grid_type),     intent(in) :: G
  type(param_file_type),     intent(in) :: param_file
  type(oil_tracer_CS),       pointer    :: CS
  type(diag_ptrs), target,   intent(in) :: diag
  type(advect_tracer_CS),    pointer    :: tr_adv_CSp
  type(GOLD_restart_CS),     pointer    :: restart_CS
! This subroutine is used to register tracer fields and subroutines
! to be used with GOLD.
! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  tr_adv_CSp - A pointer that is set to point to the control structure
!                  for the tracer advection and diffusion module.
!  (in)      restart_CS - A pointer to the restart control structure.

  character(len=128) :: version = '$Id: oil_tracer.F90,v 1.1.2.7 2010/09/15 21:13:31 rwh Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "oil_tracer" ! This module's name.
  character(len=200) :: inputdir ! The directory where the input files are.
  character(len=200) :: IC_file  ! The initial condition file without the path.
  character(len=3)   :: name_tag ! String for creating identifying oils
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_oil_tracer
  integer :: isd, ied, jsd, jed, nz, m, i, j
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke

  if (associated(CS)) then
    call GOLD_error(WARNING, "register_oil_tracer called with an "// &
                             "associated control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag

  CS%IC_file = ""  ; call read_param(param_file,"OIL_IC_FILE",CS%IC_file)
  IC_file = CS%IC_file
  CS%Z_IC_file = .false.
  call read_param(param_file,"OIL_IC_FILE_IS_Z",CS%Z_IC_file)
  if ((len_trim(CS%IC_file) > 0) .and. (scan(CS%IC_file,'/') == 0)) then
    ! Add the directory if CS%IC_file is not already a complete path.
    inputdir = "." ;  call read_param(param_file,"INPUTDIR",inputdir)
    CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
  else
    CS%IC_file = ""
  endif
  CS%mask_tracers = .false. ; CS%oil_may_reinit = .false.
  CS%oil_start_year = 0.0 ; CS%oil_end_year = 1.0E99
  CS%oil_source_k(:) = 0 ; CS%oil_decay_days(:) = 0. ; CS%oil_decay_rate(:) = 0.
  CS%oil_source_rate = 1.0
  call read_param(param_file,"MASK_MASSLESS_TRACERS",CS%mask_tracers)
  call read_param(param_file,"OIL_MAY_REINIT", CS%oil_may_reinit)
  call read_param(param_file, "OIL_SOURCE_LONGITUDE", CS%oil_source_longitude, .true.)
  call read_param(param_file, "OIL_SOURCE_LATITUDE", CS%oil_source_latitude, .true.)
  call read_param(param_file, "OIL_SOURCE_LAYER", CS%oil_source_k, .true.)
  call read_param(param_file, "OIL_DECAY_DAYS", CS%oil_decay_days, .true.)
  call read_param(param_file, "OIL_SOURCE_RATE", CS%oil_source_rate)
  call read_param(param_file, "OIL_DATED_START_YEAR", CS%oil_start_year)
  call read_param(param_file, "OIL_DATED_END_YEAR", CS%oil_end_year)

  CS%ntr = 0
  do m=1,NTR_MAX
    if (CS%oil_source_k(m)/=0) then
      write(name_tag(1:3),'("_",I2.2)') m
      CS%ntr = CS%ntr + 1
      CS%tr_desc(m) = vardesc("oil"//trim(name_tag),"Oil Tracer",'h','L','s',"kg/m3", 'f')
      CS%IC_val(m) = 0.0
      if (CS%oil_decay_days(m)>0.) then
        CS%oil_decay_rate(m)=1./(86400.0*CS%oil_decay_days(m))
      elseif (CS%oil_decay_days(m)<0.) then
        CS%oil_decay_rate(m)=-1.
      endif
    endif
  enddo

  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr)) ; CS%tr(:,:,:,:) = 0.0

  do m=1,CS%ntr
    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for the restart file.
    call register_restart_field(tr_ptr, tr_ptr, CS%tr_desc(m), &
                                .not.CS%oil_may_reinit,restart_CS)
    ! Register the tracer for horizontal advection & diffusion.
    call register_tracer(tr_ptr, CS%tr_desc(m)%name, param_file, tr_adv_CSp)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = aof_set_coupler_flux(trim(CS%tr_desc(m)%name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_oil_tracer")
  enddo

  CS%tr_adv_CSp => tr_adv_CSp
  CS%restart_CSp => restart_CS
  register_oil_tracer = .true.

  ! Write all relevant parameters to the model log.
  call log_version(param_file, mod, version, tagname, "")
  call log_param(param_file, mod, "INPUTDIR/CFC_IC_FILE", CS%IC_file)
  call log_param(param_file, mod, "OIL_IC_FILE", IC_file, &
                 "The file in which the oil tracer initial values can be \n"//&
                 "found, or an empty string for internal initialization.", &
                 default=" ")
  call log_param(param_file, mod, "OIL_IC_FILE_IS_Z", CS%Z_IC_file, &
                 "If true, OIL_IC_FILE is in depth space, not layer space", &
                 default=.false.)
  call log_param(param_file, mod, "MASK_MASSLESS_TRACERS", CS%mask_tracers, &
                 "If true, the tracers are masked out in massless layer. \n"//&
                 "This can be a problem with time-averages.", default=.false.)
  call log_param(param_file, mod, "OIL_MAY_REINIT", CS%oil_may_reinit, &
                 "If true, oil tracers may go through the initialization \n"//&
                 "code if they are not found in the restart files. \n"//&
                 "Otherwise it is a fatal error if the oil tracers are not \n"//&
                 "found in the restart files of a restared run.", &
                 default=.false.)
  call log_param(param_file, mod, "OIL_SOURCE_LONGITUDE", CS%oil_source_longitude, &
                 "The geographic longitude of the oil source.", units="degrees E")
  call log_param(param_file, mod, "OIL_SOURCE_LATITUDE", CS%oil_source_latitude, &
                 "The geographic latitude of the oil source.", units="degrees N")
  call log_param(param_file, mod, "OIL_SOURCE_LAYER", CS%oil_source_k(1:CS%ntr), &
                 "The layer into which the oil is introduced, or a \n"//&
                 "negative number for a vertically uniform source, \n"//&
                 "or 0 not to use this tracer.", units="Layer", default=0)
  call log_param(param_file, mod, "OIL_SOURCE_RATE", CS%oil_source_rate, &
                 "The rate of oil injection.", units="kg s-1", default=1.0)
  call log_param(param_file, mod, "OIL_DECAY_DAYS", CS%oil_decay_days(1:CS%ntr), &
                 "The decay timescale in days (if positive), or no decay \n"//&
                 "if 0, or use the temperature dependent decay rate of \n"//&
                 "Adcroft et al. (GRL, 2010) if negative.", units="days", &
                 default=0.0)
  call log_param(param_file, mod, "OIL_DECAY_RATE", CS%oil_decay_rate(1:CS%ntr))
  call log_param(param_file, mod, "OIL_DATED_START_YEAR", CS%oil_start_year, &
                 "The time at which the oil source starts", units="years", &
                 default=0.0)
  call log_param(param_file, mod, "OIL_DATED_END_YEAR", CS%oil_end_year, &
                 "The time at which the oil source ends", units="years", &
                 default=1.0e99)

end function register_oil_tracer

subroutine initialize_oil_tracer(restart, day, G, h, OBC, CS, sponge_CSp, &
                                       diag_to_Z_CSp)
  logical,                            intent(in) :: restart
  type(time_type), target,            intent(in) :: day
  type(ocean_grid_type),              intent(in) :: G
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in) :: h
  type(ocean_OBC_type),               pointer    :: OBC
  type(oil_tracer_CS),          pointer    :: CS
  type(sponge_CS),                    pointer    :: sponge_CSp
  type(diag_to_Z_CS),                 pointer    :: diag_to_Z_CSp
!   This subroutine initializes the CS%ntr tracer fields in tr(:,:,:,:)
! and it sets up the tracer output.

! Arguments: restart - .true. if the fields have already been read from
!                     a restart file.
!  (in)      day - Time of the start of the run.
!  (in)      G - The ocean's grid structure.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in/out)  CS - The control structure returned by a previous call to
!                 register_oil_tracer.
!  (in/out)  sponge_CSp - A pointer to the control structure for the sponges, if
!                         they are in use.  Otherwise this may be unassociated.
!  (in/out)  diag_to_Z_Csp - A pointer to the control structure for diagnostics
!                            in depth space.
  character(len=16) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for age tracer fluxes, either
                                ! years m3 s-1 or years kg s-1.
  logical :: OK
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: Isdq, Iedq, Jsdq, Jedq

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jed = G%Jedq

  ! Establish location of source
  do j=G%jsdq+1,G%jed ; do i=G%isdq+1,G%ied
    ! This test for i,j index is specific to a lat/lon (non-rotated grid).
    ! and needs to be generalized to work properly on the tri-polar grid.
    if (CS%oil_source_longitude<G%geolonq(i,j) .and. &
        CS%oil_source_longitude>=G%geolonq(i-1,j) .and. &
        CS%oil_source_latitude<G%geolatq(i,j) .and. &
        CS%oil_source_latitude>=G%geolatq(i,j-1) ) then
      CS%oil_source_i=i
      CS%oil_source_j=j
    endif
  enddo; enddo

  CS%Time => day

  do m=1,CS%ntr
    if ((.not.restart) .or. (CS%oil_may_reinit .and. .not. &
        query_initialized(CS%tr(:,:,:,m), CS%tr_desc(m)%name, CS%restart_CSp))) then

      if (len_trim(CS%IC_file) > 0) then
  !  Read the tracer concentrations from a netcdf file.
        if (.not.file_exists(CS%IC_file, G%Domain)) &
          call GOLD_error(FATAL, "initialize_oil_tracer: "// &
                                 "Unable to open "//CS%IC_file)

        if (CS%Z_IC_file) then
          OK = tracer_Z_init(CS%tr(:,:,:,m), h, CS%IC_file, CS%tr_desc(m)%name,&
                             G, -1e34, 0.0) ! CS%land_val(m))
          if (.not.OK) then
            OK = tracer_Z_init(CS%tr(:,:,:,m), h, CS%IC_file, &
                     trim(CS%tr_desc(m)%name)//"_z", G, -1e34, 0.0) ! CS%land_val(m))
            if (.not.OK) call GOLD_error(FATAL,"initialize_oil_tracer: "//&
                    "Unable to read "//trim(CS%tr_desc(m)%name)//" from "//&
                    trim(CS%IC_file)//".")
          endif
        else
          call read_data(CS%IC_file, trim(CS%tr_desc(m)%name), CS%tr(:,:,:,m), &
                         domain=G%Domain%mpp_domain)
        endif
      else
        do k=1,nz ; do j=js,je ; do i=is,ie
          if (G%hmask(i,j) < 0.5) then
            CS%tr(i,j,k,m) = CS%land_val(m)
          else
            CS%tr(i,j,k,m) = CS%IC_val(m)
          endif
        enddo ; enddo ; enddo
      endif

    endif ! restart
  enddo ! Tracer loop

  if (associated(OBC)) then
  ! All tracers but the first have 0 concentration in their inflows. As this
  ! is the default value, the following calls are unnecessary.
  ! do m=1,CS%ntr
  !  call add_tracer_OBC_values(trim(CS%tr_desc(m)%name), CS%advect_tracer_CSp, 0.0)
  ! enddo
  endif

  ! This needs to be changed if the units of tracer are changed above.
  if (G%Boussinesq) then ; flux_units = "years m3 s-1"
  else ; flux_units = "years kg s-1" ; endif

  do m=1,CS%ntr
    ! Register the tracer for the restart file.
    name = CS%tr_desc(m)%name ; longname = CS%tr_desc(m)%longname
    units = CS%tr_desc(m)%units
    CS%id_tracer(m) = register_diag_field("ocean_model", trim(name), G%axeshl, &
        day, trim(longname) , trim(units))
    CS%id_tr_adx(m) = register_diag_field("ocean_model", trim(name)//"_adx", &
        G%axesul, day, trim(longname)//" advective zonal flux" , &
        trim(flux_units))
    CS%id_tr_ady(m) = register_diag_field("ocean_model", trim(name)//"_ady", &
        G%axesvl, day, trim(longname)//" advective meridional flux" , &
        trim(flux_units))
    CS%id_tr_dfx(m) = register_diag_field("ocean_model", trim(name)//"_dfx", &
        G%axesul, day, trim(longname)//" diffusive zonal flux" , &
        trim(flux_units))
    CS%id_tr_dfy(m) = register_diag_field("ocean_model", trim(name)//"_dfy", &
        G%axesvl, day, trim(longname)//" diffusive zonal flux" , &
        trim(flux_units))
    if (CS%id_tr_adx(m) > 0) call safe_alloc_ptr(CS%tr_adx(m)%p,Isdq,Iedq,jsd,jed,nz)
    if (CS%id_tr_ady(m) > 0) call safe_alloc_ptr(CS%tr_ady(m)%p,isd,ied,Jsdq,Jedq,nz)
    if (CS%id_tr_dfx(m) > 0) call safe_alloc_ptr(CS%tr_dfx(m)%p,Isdq,Iedq,jsd,jed,nz)
    if (CS%id_tr_dfy(m) > 0) call safe_alloc_ptr(CS%tr_dfy(m)%p,isd,ied,Jsdq,Jedq,nz)

!    Register the tracer for horizontal advection & diffusion.
    if ((CS%id_tr_adx(m) > 0) .or. (CS%id_tr_ady(m) > 0) .or. &
        (CS%id_tr_dfx(m) > 0) .or. (CS%id_tr_dfy(m) > 0)) &
      call add_tracer_diagnostics(name, CS%tr_adv_CSp, CS%tr_adx(m)%p, &
                                  CS%tr_ady(m)%p,CS%tr_dfx(m)%p,CS%tr_dfy(m)%p)

    call register_Z_tracer(CS%tr(:,:,:,m), trim(name)//"_z", longname, units, &
                           day, G, diag_to_Z_CSp)
  enddo

end subroutine initialize_oil_tracer

subroutine oil_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, CS, tv)
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in) :: h_old, h_new, ea, eb
  type(forcing),                      intent(in) :: fluxes
  real,                               intent(in) :: dt
  type(ocean_grid_type),              intent(in) :: G
  type(oil_tracer_CS),                pointer    :: CS
  type(thermo_var_ptrs),              intent(in) :: tv
!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! This is a simple example of a set of advected passive tracers.

! Arguments: h_old -  Layer thickness before entrainment, in m or kg m-2.
!  (in)      h_new -  Layer thickness after entrainment, in m or kg m-2.
!  (in)      ea - an array to which the amount of fluid entrained
!                 from the layer above during this call will be
!                 added, in m or kg m-2.
!  (in)      eb - an array to which the amount of fluid entrained
!                 from the layer below during this call will be
!                 added, in m or kg m-2.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_oil_tracer.
!
! The arguments to this subroutine are redundant in that
!     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]

  real :: Isecs_per_year = 1.0 / (365.0*86400.0)
  real :: year, h_total, ldecay
  integer :: secs, days
  integer :: i, j, k, is, ie, js, je, nz, m, k_max
  real, allocatable :: local_tr(:,:,:)
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  do m=1,CS%ntr
    call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G)
  enddo

  !   Set the surface value of tracer 1 to increase exponentially
  ! with a 30 year time scale.
  call get_time(CS%Time, secs, days)
  year = (86400.0*days + real(secs)) * Isecs_per_year

  ! Decay tracer (limit decay rate to 1/dt - just in case)
  do m=2,CS%ntr
    do k=1,nz ; do j=js,je ; do i=is,ie
      !CS%tr(i,j,k,m) = CS%tr(i,j,k,m) - dt*CS%oil_decay_rate(m)*CS%tr(i,j,k,m) ! Simple
      !CS%tr(i,j,k,m) = CS%tr(i,j,k,m) - min(dt*CS%oil_decay_rate(m),1.)*CS%tr(i,j,k,m) ! Safer
      if (CS%oil_decay_rate(m)>0.) then
        CS%tr(i,j,k,m) = G%hmask(i,j)*max(1.-dt*CS%oil_decay_rate(m),0.)*CS%tr(i,j,k,m) ! Safest
      elseif (CS%oil_decay_rate(m)<0.) then
        ldecay = 12.*(3.0**(-(tv%T(i,j,k)-20.)/10.)) ! Timescale in days
        ldecay = 1./(86400.*ldecay) ! Rate in s^-1
        CS%tr(i,j,k,m) = G%hmask(i,j)*max(1.-dt*ldecay,0.)*CS%tr(i,j,k,m)
      endif
    enddo ; enddo ; enddo
  enddo

  ! Add oil at the source location
  if (year>=CS%oil_start_year .and. year<=CS%oil_end_year .and. &
      CS%oil_source_i>-999 .and. CS%oil_source_j>-999) then
    i=CS%oil_source_i ; j=CS%oil_source_j
    k_max=nz ; h_total=0.
    do k=nz, 2, -1
      h_total = h_total + h_new(i,j,k)
      if (h_total<10.) k_max=k-1 ! Find bottom most interface that is 10 m above bottom
    enddo
    do m=1,CS%ntr
      k=CS%oil_source_k(m)
      if (k>0) then
        k=min(k,k_max) ! Only insert k or first layer with interface 10 m above bottom
        CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + CS%oil_source_rate*dt/((h_new(i,j,k)+G%h_subroundoff) &
                                           * G%DXDYh(i,j) )
      elseif (k<0) then
        h_total=G%h_subroundoff
        do k=1, nz
          h_total = h_total + h_new(i,j,k)
        enddo
        do k=1, nz
          CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + CS%oil_source_rate*dt/(h_total &
                                           * G%DXDYh(i,j) )
        enddo
      endif
    enddo
  endif

  allocate(local_tr(G%isd:G%ied,G%jsd:G%jed,nz))
  do m=1,CS%ntr
    if (CS%id_tracer(m)>0) then
      if (CS%mask_tracers) then
        do k=1,nz ; do j=js,je ; do i=is,ie
          if (h_new(i,j,k) < 1.1*G%Angstrom) then
            local_tr(i,j,k) = CS%land_val(m)
          else
            local_tr(i,j,k) = CS%tr(i,j,k,m)
          endif
        enddo ; enddo ; enddo
      else
        do k=1,nz ; do j=js,je ; do i=is,ie
          local_tr(i,j,k) = CS%tr(i,j,k,m)
        enddo ; enddo ; enddo
      endif ! CS%mask_tracers
      call post_data(CS%id_tracer(m),local_tr,CS%diag)
    endif ! CS%id_tracer(m)>0
    if (CS%id_tr_adx(m)>0) &
      call post_data(CS%id_tr_adx(m),CS%tr_adx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_ady(m)>0) &
      call post_data(CS%id_tr_ady(m),CS%tr_ady(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfx(m)>0) &
      call post_data(CS%id_tr_dfx(m),CS%tr_dfx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfy(m)>0) &
      call post_data(CS%id_tr_dfy(m),CS%tr_dfy(m)%p(:,:,:),CS%diag)
  enddo
  deallocate(local_tr)

end subroutine oil_tracer_column_physics

function oil_stock(h, stocks, G, CS, names, units, stock_index)
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in)    :: h
  real, dimension(:),                 intent(out)   :: stocks
  type(ocean_grid_type),              intent(in)    :: G
  type(oil_tracer_CS),          pointer       :: CS
  character(len=*), dimension(:),     intent(out)   :: names
  character(len=*), dimension(:),     intent(out)   :: units
  integer, optional,                  intent(in)    :: stock_index
  integer                                           :: oil_stock
! This function calculates the mass-weighted integral of all tracer stocks,
! returning the number of stocks it has calculated.  If the stock_index
! is present, only the stock corresponding to that coded index is returned.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (out)     stocks - the mass-weighted integrated amount of each tracer,
!                     in kg times concentration units.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_oil_tracer.
!  (out)     names - the names of the stocks calculated.
!  (out)     units - the units of the stocks calculated.
!  (in,opt)  stock_index - the coded index of a specific stock being sought.
! Return value: the number of stocks calculated here.

  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  oil_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,CS%ntr
    names(m) = CS%tr_desc(m)%name ; units(m) = trim(CS%tr_desc(m)%units)//" kg"
    stocks(m) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      stocks(m) = stocks(m) + CS%tr(i,j,k,m) * &
                             (G%hmask(i,j) * G%DXDYh(i,j) * h(i,j,k))
    enddo ; enddo ; enddo
    stocks(m) = G%H_to_kg_m2 * stocks(m)
  enddo
  oil_stock = CS%ntr

end function oil_stock

subroutine oil_tracer_surface_state(state, h, G, CS)
  type(surface),                      intent(inout) :: state
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in)    :: h
  type(ocean_grid_type),              intent(in)    :: G
  type(oil_tracer_CS),          pointer       :: CS
!   This particular tracer package does not report anything back to the coupler.
! The code that is here is just a rough guide for packages that would.
! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_oil_tracer.
  integer :: m, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  
  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,CS%ntr
      !   This call loads the surface vlues into the appropriate array in the
      ! coupler-type structure.
      call set_coupler_values(CS%tr(:,:,1,m), state%tr_fields, CS%ind_tr(m), &
                              ind_csurf, is, ie, js, je)
    enddo
  endif

end subroutine oil_tracer_surface_state

subroutine oil_tracer_end(CS)
  type(oil_tracer_CS), pointer :: CS
  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    do m=1,CS%ntr
      if (associated(CS%tr_adx(m)%p)) deallocate(CS%tr_adx(m)%p)
      if (associated(CS%tr_ady(m)%p)) deallocate(CS%tr_ady(m)%p)
      if (associated(CS%tr_dfx(m)%p)) deallocate(CS%tr_dfx(m)%p)
      if (associated(CS%tr_dfy(m)%p)) deallocate(CS%tr_dfy(m)%p)
    enddo

    deallocate(CS)
  endif
end subroutine oil_tracer_end

end module oil_tracer
