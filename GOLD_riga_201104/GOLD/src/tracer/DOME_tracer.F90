module DOME_tracer
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
!*  up and use a set (in this case eleven) of dynamically passive      *
!*  tracers.  These tracers dye the inflowing water or water initially *
!*  within a range of latitudes or water initially in a range of       *
!*  depths.                                                            *
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
use GOLD_restart, only : register_restart_field, GOLD_restart_CS
use GOLD_sponge, only : set_up_sponge_field, sponge_CS
use GOLD_time_manager, only : time_type, get_time
use GOLD_tracer, only : register_tracer, advect_tracer_CS, tracer_vertdiff
use GOLD_tracer, only : add_tracer_diagnostics, add_tracer_OBC_values
use GOLD_variables, only : forcing, surface, ocean_OBC_type

use coupler_util, only : set_coupler_values, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <GOLD_memory.h>

public register_DOME_tracer, initialize_DOME_tracer
public DOME_tracer_column_physics, DOME_tracer_surface_state, DOME_tracer_end

! ntr is the number of tracers in this module.
integer, parameter :: ntr = 11

type p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d

type, public :: DOME_tracer_CS ; private
  logical :: coupled_tracers = .false.  ! These tracers are not offered to the
                                        ! coupler.
  character(len = 200) :: tracer_IC_file ! The full path to the IC file, or " "
                                   ! to initialize internally.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(advect_tracer_CS), pointer :: tr_adv_CSp => NULL()
  real, pointer :: tr(:,:,:,:) => NULL()   ! The array of tracers used in this
                                           ! subroutine, in g m-3?
  real, pointer :: tr_aux(:,:,:,:) => NULL() ! The masked tracer concentration
                                             ! for output, in g m-3.
  type(p3d), dimension(NTR) :: &
    tr_adx, &! Tracer zonal advective fluxes in g m-3 m3 s-1.
    tr_ady, &! Tracer meridional advective fluxes in g m-3 m3 s-1.
    tr_dfx, &! Tracer zonal diffusive fluxes in g m-3 m3 s-1.
    tr_dfy   ! Tracer meridional diffusive fluxes in g m-3 m3 s-1.
  real :: land_val(NTR) = -1.0 ! The value of tr used where land is masked out.
  logical :: mask_tracers  ! If true, tracers are masked out in massless layers.
  logical :: use_sponge

  integer, dimension(NTR) :: ind_tr ! Indices returned by aof_set_coupler_flux
             ! if it is used and the surface tracer concentrations are to be
             ! provided to the coupler.

  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
                             ! ocean diagnostic fields and control variables.
  integer, dimension(NTR) :: id_tracer = -1, id_tr_adx = -1, id_tr_ady = -1
  integer, dimension(NTR) :: id_tr_dfx = -1, id_tr_dfy = -1

  type(vardesc) :: tr_desc(NTR)
end type DOME_tracer_CS

contains

function register_DOME_tracer(G, param_file, CS, diag, tr_adv_CSp, &
                                      restart_CS)
  type(ocean_grid_type),   intent(in) :: G
  type(param_file_type),   intent(in) :: param_file
  type(DOME_tracer_CS),    pointer    :: CS
  type(diag_ptrs), target, intent(in) :: diag
  type(advect_tracer_CS),  pointer    :: tr_adv_CSp
  type(GOLD_restart_CS),   pointer    :: restart_CS
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
  character(len=80)  :: name, longname
  character(len=128) :: version = '$Id: DOME_tracer.F90,v 1.1.2.3.2.11 2010/10/01 16:11:09 rwh Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "DOME_tracer" ! This module's name.
  character(len=200) :: inputdir, tracer_IC_file
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_DOME_tracer
  integer :: isd, ied, jsd, jed, nz, m
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke

  if (associated(CS)) then
    call GOLD_error(WARNING, "DOME_register_tracer called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag
  CS%tracer_IC_file = " " ; tracer_IC_file = " "
  call read_param(param_file,"DOME_TRACER_IC_FILE",tracer_IC_file)
  if (len_trim(tracer_IC_file) >= 1) then
    inputdir = "." ;  call read_param(param_file,"INPUTDIR",inputdir)
    CS%tracer_IC_file = trim(slasher(inputdir))//trim(tracer_IC_file)
  endif
  CS%use_sponge = .false. ;  call read_param(param_file,"SPONGE",CS%use_sponge)

  allocate(CS%tr(isd:ied,jsd:jed,nz,NTR)) ; CS%tr(:,:,:,:) = 0.0
  if (CS%mask_tracers) then
    allocate(CS%tr_aux(isd:ied,jsd:jed,nz,NTR)) ; CS%tr_aux(:,:,:,:) = 0.0
  endif

  do m=1,NTR
    CS%tr_desc(m) = vardesc("tr","Tracer",'h','L','s',"kg kg-1", 'd')
    if (m < 10) then ; write(name,'("tr_D",I1.1)') m
    else ; write(name,'("tr_D",I2.2)') m ; endif
    write(longname,'("Concentration of DOME Tracer ",I2.2)') m
    CS%tr_desc(m)%name = name
    CS%tr_desc(m)%longname = longname
    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for the restart file.
    call register_restart_field(tr_ptr, tr_ptr, CS%tr_desc(m),.true.,restart_CS)
    ! Register the tracer for horizontal advection & diffusion.
    call register_tracer(tr_ptr, CS%tr_desc(m)%name, param_file, tr_adv_CSp)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = aof_set_coupler_flux(trim(CS%tr_desc(m)%name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_DOME_tracer")
  enddo

  ! Write all relevant parameters to the model log.
  call log_version(param_file, mod, version, tagname, "")
  call log_param(param_file, mod, "DOME_TRACER_IC_FILE", tracer_IC_file, &
                 "The name of a file from which to read the initial \n"//&
                 "conditions for the DOME tracers, or blank to initialize \n"//&
                 "them internally.", default=" ")
  if (len_trim(tracer_IC_file) >= 1) &
    call log_param(param_file, mod, "INPUTDIR/DOME_TRACER_IC_FILE", &
                   CS%tracer_IC_file)
  call log_param(param_file, mod, "SPONGE", CS%use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. \n"//&
                 "The exact location and properties of those sponges are \n"//&
                 "specified from GOLD_initialization.F90.", default=.false.)
  CS%tr_adv_CSp => tr_adv_CSp
  register_DOME_tracer = .true.
end function register_DOME_tracer

subroutine initialize_DOME_tracer(restart, day, G, h, OBC, CS, sponge_CSp, &
                                  diag_to_Z_CSp)
  logical,                            intent(in) :: restart
  type(time_type), target,            intent(in) :: day
  type(ocean_grid_type),              intent(in) :: G
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in) :: h
  type(ocean_OBC_type),               pointer    :: OBC
  type(DOME_tracer_CS),               pointer    :: CS
  type(sponge_CS),                    pointer    :: sponge_CSp
  type(diag_to_Z_CS),                 pointer    :: diag_to_Z_CSp
!   This subroutine initializes the NTR tracer fields in tr(:,:,:,:)
! and it sets up the tracer output.

! Arguments: restart - .true. if the fields have already been read from
!                     a restart file.
!  (in)      day - Time of the start of the run.
!  (in)      G - The ocean's grid structure.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in/out)  CS - The control structure returned by a previous call to
!                 DOME_register_tracer.
!  (in/out)  sponge_CSp - A pointer to the control structure for the sponges, if
!                         they are in use.  Otherwise this may be unassociated.
!  (in/out)  diag_to_Z_Csp - A pointer to the control structure for diagnostics
!                            in depth space.
  real, allocatable :: temp(:,:,:)
  real, pointer, dimension(:,:,:) :: &
    OBC_tr1_u => NULL(), & ! These arrays should be allocated and set to
    OBC_tr1_v => NULL()    ! specify the values of tracer 1 that should come
                           ! in through u- and v- points through the open
                           ! boundary conditions, in the same units as tr.
  character(len=16) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for tracer fluxes, usually
                            ! kg(tracer) kg(water)-1 m3 s-1 or kg(tracer) s-1.
  real, pointer :: tr_ptr(:,:,:) => NULL()
  real :: PI     ! 3.1415926... calculated as 4*atan(1)
  real :: tr_y   ! Initial zonally uniform tracer concentrations.
  real :: dist2  ! The distance squared from a line, in m2.
  real :: h_neglect         ! A thickness that is so small it is usually lost
                            ! in roundoff and can be neglected, in m.
  real :: e(SZK_(G)+1), e_top, e_bot, d_tr
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: Isdq, Iedq, Jsdq, Jedq

  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jed = G%Jedq
  h_neglect = G%H_subroundoff

  CS%Time => day

  if (.not.restart) then
    if (len_trim(CS%tracer_IC_file) >= 1) then
      !  Read the tracer concentrations from a netcdf file.
      if (.not.file_exists(CS%tracer_IC_file, G%Domain)) &
        call GOLD_error(FATAL, "DOME_initialize_tracer: Unable to open "// &
                        CS%tracer_IC_file)
      do m=1,NTR
        call read_data(CS%tracer_IC_file, trim(CS%tr_desc(m)%name), &
                       CS%tr(:,:,:,m), domain=G%Domain%mpp_domain)
      enddo
    else
      do m=1,NTR
        do k=1,nz ; do j=js,je ; do i=is,ie
          CS%tr(i,j,k,m) = 1.0e-20 ! This could just as well be 0.
        enddo ; enddo ; enddo
      enddo

!    This sets a stripe of tracer across the basin.
      do m=2,NTR ; do j=js,je ; do i=is,ie
        tr_y = 0.0
        if ((m <= 6) .and. (G%geolath(i,j) > (300.0+50.0*real(m-1))) .and. &
            (G%geolath(i,j) < (350.0+50.0*real(m-1)))) tr_y = 1.0
        do k=1,nz
!      This adds the stripes of tracer to every layer.
            CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + tr_y
        enddo
      enddo; enddo; enddo

      if (NTR > 7) then
        do j=js,je ; do i=is,ie
          e(nz+1) = -G%D(i,j)
          do k=nz,1,-1
            e(k) = e(k+1) + h(i,j,k)*G%H_to_m
            do m=7,NTR
              e_top = -600.0*real(m-1) + 3000.0
              e_bot = -600.0*real(m-1) + 2700.0
              if (e_top < e(k)) then
                if (e_top < e(k+1)) then ; d_tr = 0.0
                elseif (e_bot < e(k+1)) then
                  d_tr = (e_top-e(k+1)) / ((h(i,j,k)+h_neglect)*G%H_to_m)
                else ; d_tr = (e_top-e_bot) / ((h(i,j,k)+h_neglect)*G%H_to_m)
                endif
              elseif (e_bot < e(k)) then
                if (e_bot < e(k+1)) then ; d_tr = 1.0
                else ; d_tr = (e(k)-e_bot) / ((h(i,j,k)+h_neglect)*G%H_to_m)
                endif
              else
                d_tr = 0.0
              endif
              if (h(i,j,k) < 2.0*G%Angstrom) d_tr=0.0
              CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + d_tr
            enddo
          enddo
        enddo ; enddo
      endif

    endif
  endif ! restart

  if ( CS%use_sponge ) then
!   If sponges are used, this example damps tracers in sponges in the
! northern half of the domain to 1 and tracers in the southern half
! to 0.  For any tracers that are not damped in the sponge, the call
! to set_up_sponge_field can simply be omitted.
    if (.not.associated(sponge_CSp)) &
      call GOLD_error(FATAL, "DOME_initialize_tracer: "// &
        "The pointer to sponge_CSp must be associated if SPONGE is defined.")

    allocate(temp(G%isd:G%ied,G%jsd:G%jed,nz))
    do k=1,nz ; do j=js,je ; do i=is,ie
      if (G%geolath(i,j) > 700.0 .and. (k > nz/2)) then
        temp(i,j,k) = 1.0
      else
        temp(i,j,k) = 0.0
      endif
    enddo ; enddo ; enddo

!   do m=1,NTR
    do m=1,1
      ! This is needed to force the compiler not to do a copy in the sponge
      ! calls.  Curses on the designers and implementers of Fortran90.
      tr_ptr => CS%tr(:,:,:,m)
      call set_up_sponge_field(temp,tr_ptr,nz,sponge_CSp)
    enddo
    deallocate(temp)
  endif

  if (associated(OBC)) then
    if (OBC%apply_OBC_v) then
      allocate(OBC_tr1_v(G%isd:G%ied,G%jsd:G%jed,nz))
      do k=1,nz ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
        if (k < nz/2) then ; OBC_tr1_v(i,j,k) = 0.0
        else ; OBC_tr1_v(i,j,k) = 1.0 ; endif
      enddo ; enddo ; enddo
      call add_tracer_OBC_values(trim(CS%tr_desc(1)%name), CS%tr_adv_CSp, &
                                 0.0, OBC_in_v=OBC_tr1_v)
    else
      ! This is not expected in the DOME example.
      call add_tracer_OBC_values(trim(CS%tr_desc(1)%name), CS%tr_adv_CSp, 0.0)
    endif
    ! All tracers but the first have 0 concentration in their inflows. As this
    ! is the default value, the following calls are unnecessary.
    do m=2,NTR
      call add_tracer_OBC_values(trim(CS%tr_desc(m)%name), CS%tr_adv_CSp, 0.0)
    enddo
  endif

  ! This needs to be changed if the units of tracer are changed above.
  if (G%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
  else ; flux_units = "kg s-1" ; endif

  do m=1,NTR
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

end subroutine initialize_DOME_tracer

subroutine DOME_tracer_column_physics(h_old, h_new,  ea,  eb, fluxes, dt, G, CS)
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in) :: h_old, h_new, ea, eb
  type(forcing),                      intent(in) :: fluxes
  real,                               intent(in) :: dt
  type(ocean_grid_type),              intent(in) :: G
  type(DOME_tracer_CS),               pointer    :: CS
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
!                 DOME_register_tracer.
!
! The arguments to this subroutine are redundant in that
!     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]

  real :: b1(SZI_(G))          ! b1 and c1 are variables used by the
  real :: c1(SZI_(G),SZK_(G))  ! tridiagonal solver.
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(CS)) return

  do m=1,NTR
    call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G)
  enddo

  if (CS%mask_tracers) then
    do m = 1,NTR ; if (CS%id_tracer(m) > 0) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        if (h_new(i,j,k) < 1.1*G%Angstrom) then
          CS%tr_aux(i,j,k,m) = CS%land_val(m)
        else
          CS%tr_aux(i,j,k,m) = CS%tr(i,j,k,m)
        endif
      enddo ; enddo ; enddo
    endif ; enddo
  endif

  do m=1,NTR
    if (CS%mask_tracers) then
      if (CS%id_tracer(m)>0) &
        call post_data(CS%id_tracer(m),CS%tr_aux(:,:,:,m),CS%diag)
    else
      if (CS%id_tracer(m)>0) &
        call post_data(CS%id_tracer(m),CS%tr(:,:,:,m),CS%diag)
    endif
    if (CS%id_tr_adx(m)>0) &
      call post_data(CS%id_tr_adx(m),CS%tr_adx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_ady(m)>0) &
      call post_data(CS%id_tr_ady(m),CS%tr_ady(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfx(m)>0) &
      call post_data(CS%id_tr_dfx(m),CS%tr_dfx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_dfy(m)>0) &
      call post_data(CS%id_tr_dfy(m),CS%tr_dfy(m)%p(:,:,:),CS%diag)
  enddo

end subroutine DOME_tracer_column_physics

subroutine DOME_tracer_surface_state(state, h, G, CS)
  type(surface),                      intent(inout) :: state
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in)    :: h
  type(ocean_grid_type),              intent(in)    :: G
  type(DOME_tracer_CS),               pointer       :: CS
!   This particular tracer package does not report anything back to the coupler.
! The code that is here is just a rough guide for packages that would.
! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 DOME_register_tracer.
  integer :: i, j, m, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  
  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,ntr
      !   This call loads the surface vlues into the appropriate array in the
      ! coupler-type structure.
      call set_coupler_values(CS%tr(:,:,1,1), state%tr_fields, CS%ind_tr(m), &
                              ind_csurf, is, ie, js, je)
    enddo
  endif

end subroutine DOME_tracer_surface_state

subroutine DOME_tracer_end(CS)
  type(DOME_tracer_CS), pointer :: CS
  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    if (associated(CS%tr_aux)) deallocate(CS%tr_aux)
    do m=1,NTR
      if (associated(CS%tr_adx(m)%p)) deallocate(CS%tr_adx(m)%p)
      if (associated(CS%tr_ady(m)%p)) deallocate(CS%tr_ady(m)%p)
      if (associated(CS%tr_dfx(m)%p)) deallocate(CS%tr_dfx(m)%p)
      if (associated(CS%tr_dfy(m)%p)) deallocate(CS%tr_dfy(m)%p)
    enddo

    deallocate(CS)
  endif
end subroutine DOME_tracer_end

end module DOME_tracer
