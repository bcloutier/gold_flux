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
!----------------------------------------------------------------
! <CONTACT EMAIL="Niki.Zadeh@noaa.gov"> Niki Zadeh 
! </CONTACT>
! 
! <REVIEWER EMAIL="William.Cooke@noaa.gov"> William Cooke
! </REVIEWER>
!
! <OVERVIEW>
!  This module drives the generic version of tracers TOPAZ and CFC
! </OVERVIEW>
!----------------------------------------------------------------

#include <GOLD_memory.h>

module GOLD_generic_tracer

#ifdef _USE_GENERIC_TRACER
#include <fms_platform.h>

  use mpp_mod,        only: stdout, mpp_error, FATAL,WARNING,NOTE 
  use field_manager_mod, only: fm_get_index,fm_string_len

  use generic_tracer, only: generic_tracer_register, generic_tracer_get_diag_list
  use generic_tracer, only: generic_tracer_init, generic_tracer_source, generic_tracer_register_diag
  use generic_tracer, only: generic_tracer_coupler_get, generic_tracer_coupler_set
  use generic_tracer, only: generic_tracer_end, generic_tracer_get_list, do_generic_tracer
  use generic_tracer, only: generic_tracer_update_from_bottom,generic_tracer_vertdiff_G

  use g_tracer_utils,   only: g_tracer_get_name,g_tracer_set_values,g_tracer_set_common,g_tracer_get_common
  use g_tracer_utils,   only: g_tracer_get_next,g_tracer_type,g_tracer_is_prog,g_tracer_flux_init
  use g_tracer_utils,   only: g_tracer_send_diag,g_tracer_get_values
  use g_tracer_utils,   only: g_tracer_get_pointer,g_tracer_get_alias,g_diag_type

  use GOLD_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
  use GOLD_diag_mediator, only : diag_ptrs
  use GOLD_diag_to_Z, only : register_Z_tracer, diag_to_Z_CS
  use GOLD_error_handler, only : GOLD_error, FATAL, WARNING, NOTE, is_root_pe
  use GOLD_file_parser, only : read_param, param_file_type, log_param, log_version
  use GOLD_grid, only : ocean_grid_type
  use GOLD_io, only : file_exists, read_data, slasher, vardesc
  use GOLD_restart, only : register_restart_field, query_initialized, GOLD_restart_CS
  use GOLD_sponge, only : set_up_sponge_field, sponge_CS
  use GOLD_time_manager, only : time_type, get_time
  use GOLD_tracer, only : register_tracer, advect_tracer_CS, tracer_vertdiff
  use GOLD_tracer, only : add_tracer_diagnostics, add_tracer_OBC_values
  use GOLD_tracer_Z_init, only : tracer_Z_init
  use GOLD_variables, only : forcing, surface, ocean_OBC_type, thermo_var_ptrs
  use GOLD_variables, only : optics_type


  implicit none ; private
  logical :: g_registered = .false.

  public register_GOLD_generic_tracer, initialize_GOLD_generic_tracer
  public GOLD_generic_tracer_column_physics, GOLD_generic_tracer_surface_state
  public end_GOLD_generic_tracer, GOLD_generic_tracer_get
  public GOLD_generic_tracer_stock
  public GOLD_generic_flux_init
  public GOLD_generic_tracer_min_max

  type, public :: GOLD_generic_tracer_CS ; private
     character(len = 200) :: IC_file ! The file in which the generic tracer initial values can
                       ! be found, or an empty string for internal initialization.
     logical :: Z_IC_file ! If true, the generic_tracer IC_file is in Z-space.  The default is false..
     real :: tracer_IC_val = 0.0    ! The initial value assigned to tracers.
     real :: tracer_land_val = -1.0 ! The values of tracers used where  land is masked out.
     logical :: tracers_may_reinit  ! If true, tracers may go through the
                              ! initialization code if they are not found in the
                              ! restart files.

     type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
                                ! ocean diagnostic fields and control variables.
     type(GOLD_restart_CS), pointer :: restart_CSp => NULL()

     real    :: Rho_0
     integer :: nkml
     !   The following pointer will be directed to the first element of the
     ! linked list of generic tracers.
     type(g_tracer_type), pointer :: g_tracer_list => NULL()
     !   The following pointer will be directed to the first element of the
     ! linked list of generic diagnostics fields that must be Z registered by GOLD.
     type(g_diag_type), pointer :: g_diag_list => NULL()

  end type GOLD_generic_tracer_CS

  character(len=128) :: version = '$Id: GOLD_generic_tracer.F90,v 1.1.4.16 2010/09/24 19:11:44 aja Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'

contains

  ! <SUBROUTINE NAME="register_GOLD_generic_tracer">
  !  <OVERVIEW>
  !   Initialize phase I: Add the generic tracers
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine:
  !     Initializes the generic tracer packages and adds their tracers to the list
  !     Adds the tracers in the list of generic tracers to the set of GOLD tracers (i.e., GOLD-register them)
  !     Register these tracers for restart
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call register_GOLD_generic_tracer(G, param_file, CS, diag, tr_adv_CSp, restart_CS)
  !  </TEMPLATE>
  ! </SUBROUTINE>

  function register_GOLD_generic_tracer(G, param_file, CS, diag, tr_adv_CSp, restart_CS)
    type(ocean_grid_type), intent(in)   :: G
    type(param_file_type), intent(in)   :: param_file
    type(GOLD_generic_tracer_CS),   pointer      :: CS
    type(diag_ptrs), target, intent(in) :: diag
    type(advect_tracer_CS), pointer     :: tr_adv_CSp
    type(GOLD_restart_CS),   pointer     :: restart_CS
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
    logical :: register_GOLD_generic_tracer


    character(len=fm_string_len), parameter :: sub_name = 'register_GOLD_generic_tracer'
    character(len=200) :: inputdir ! The directory where NetCDF input files are.
    ! These can be overridden later in via the field manager?

    integer :: ntau, k,i,j,axes(3)
    type(g_tracer_type), pointer      :: g_tracer,g_tracer_next
    character(len=fm_string_len)      :: g_tracer_name,longname,units
    real, dimension(:,:,:,:), pointer   :: tr_field 
    real, dimension(:,:,:), pointer     :: tr_ptr 
    real, dimension(G%isd:G%ied, G%jsd:G%jed,G%ke)         :: grid_tmask
    character(len=200) :: IC_file  ! The initial condition file without the path.
    integer, dimension(G%isd:G%ied, G%jsd:G%jed)           :: grid_kmt
    type(vardesc) :: var_desc

    register_GOLD_generic_tracer = .false.
    if (associated(CS)) then
       call mpp_error(WARNING, "register_GOLD_generic_tracer called with an "// &
            "associated control structure.")
       return
    endif
    allocate(CS)


    !Register all the generic tracers used and create the list of them. 
    !This can be called by ALL PE's. No array fields allocated.
    if (.not. g_registered) then 
       call generic_tracer_register
       g_registered = .true.
    endif

    !nnz: Trick to get the time. I don't know how it works! time=diag%time_end
    CS%diag => diag  

    CS%IC_file = ""  ; call read_param(param_file,"GENERIC_TRACER_IC_FILE",CS%IC_file)
    CS%Z_IC_file = .false. ; call read_param(param_file,"GENERIC_TRACER_IC_FILE_IS_Z",CS%Z_IC_file)
    IC_file = CS%IC_file
    if ((len_trim(CS%IC_file) > 0) .and. (scan(CS%IC_file,'/') == 0)) then
      ! Add the directory if CS%IC_file is not already a complete path.
      inputdir = "." ;  call read_param(param_file,"INPUTDIR",inputdir)
      CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
    endif
    CS%tracers_may_reinit = .false.
    call read_param(param_file,"TRACERS_MAY_REINIT", CS%tracers_may_reinit)

    call log_version(param_file, sub_name, version, tagname, "")
    call log_param(param_file, sub_name, "INPUTDIR/GENERIC_TRACER_IC_FILE", CS%IC_file)
    call log_param(param_file, sub_name, "GENERIC_TRACER_IC_FILE", IC_file, &
                 "The file in which the generic trcer initial values can \n"//&
                 "be found, or an empty string for internal initialization.", &
                 default=" ")
    call log_param(param_file, sub_name, "GENERIC_TRACER_IC_FILE_IS_Z", CS%Z_IC_file, &
                 "If true, GENERIC_TRACER_IC_FILE is in depth space, not \n"//&
                 "layer space.",default=.false.)
    call log_param(param_file, sub_name, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code \n"//&
                 "if they are not found in the restart files.  Otherwise \n"//&
                 "it is a fatal error if tracers are not found in the \n"//&
                 "restart files of a restared run.", default=.false.)

    CS%restart_CSp => restart_CS


    ntau=1 ! GOLD needs the fields at only one time step 


    !At this point G%hmask is not allocated. 
    !postpone diag_registeration to initialize_GOLD_generic_tracer

    !!nnz: G%axeshl is not defined either. 
    !Fields cannot be diag registered as they are allocated and have to registered later.
    grid_tmask(:,:,:) = 0.0
    grid_kmt(:,:) = 0.0
    axes(:) = -1

    !
    ! Initialize all generic tracers
    !
    call generic_tracer_init(G%isc,G%iec,G%jsc,G%jec,G%isd,G%ied,G%jsd,G%jed,&
         G%ke,ntau,axes,grid_tmask,grid_kmt,CS%diag%time_end)


    !
    ! GOLD-register the generic tracers
    !

    !Get the tracer list
    call generic_tracer_get_list(CS%g_tracer_list)
    if(.NOT. associated(CS%g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")
    !For each tracer name get its T_prog index and get its fields
    
    g_tracer=>CS%g_tracer_list
    do
       call g_tracer_get_alias(g_tracer,g_tracer_name)

       call g_tracer_get_pointer(g_tracer,g_tracer_name,'field',tr_field)
       call g_tracer_get_values(g_tracer,g_tracer_name,'longname', longname)
       call g_tracer_get_values(g_tracer,g_tracer_name,'units',units )

       !nnz: Hard coded stuff. Need get/set routines
       var_desc = vardesc(g_tracer_name,longname,'h','L','s',units,'f')
       !!nnz: GOLD field is 3D. Does this affect performance? Need it be override field?
       tr_ptr => tr_field(:,:,:,1)
       ! Register tracer for restart file. 
       ! mandatory field in restart file is set to .false. 
       ! 2008/12/08 jgj: change default to true, so all fields must be present in restart.
       ! 2010/02/04 jgj: if tracers_may_reinit is true, tracers may go through 
       ! initialization code if not found in restart
       call register_restart_field(tr_ptr, tr_ptr, var_desc, .not.CS%tracers_may_reinit, restart_CS)

       ! Register prognastic tracer for horizontal advection & diffusion.
       if(g_tracer_is_prog(g_tracer)) call register_tracer(tr_ptr, g_tracer_name, param_file, tr_adv_CSp)

       !traverse the linked list till hit NULL
       call g_tracer_get_next(g_tracer, g_tracer_next)
       if(.NOT. associated(g_tracer_next)) exit
       g_tracer=>g_tracer_next  

    enddo

    !Rho_0 is a property of GOLD and not a property of any tracer
    !Read it directly from GOLD_input
    CS%Rho_0=1.0
    call read_param(param_file,"RHO_0",CS%Rho_0,.true.)
    !So is NKML
    CS%nkml=1
    call read_param(param_file,"NKML",CS%nkml,.true.)

    register_GOLD_generic_tracer = .true.
  end function register_GOLD_generic_tracer

  ! <SUBROUTINE NAME="initialize_GOLD_generic_tracer">
  !  <OVERVIEW>
  !   Initialize phase II:  Initialize required variables for generic tracers
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   There are some steps of initialization that cannot be done in register_GOLD_generic_tracer
  !   This is the place and time to do them:
  !       Set the grid mask and initial time for all generic tracers.
  !       Diag_register them.
  !       Z_diag_register them.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call initialize_GOLD_generic_tracer(restart, day, G, h, OBC, CS, sponge_CSp, diag_to_Z_CSp) 
 ! </SUBROUTINE>
  subroutine initialize_GOLD_generic_tracer(restart, day, G, h, OBC, CS, sponge_CSp, &
       diag_to_Z_CSp)
    logical,                            intent(in) :: restart
    type(time_type), target,            intent(in) :: day
    type(ocean_grid_type),              intent(in) :: G
    real, dimension(NXMEM_,NYMEM_,NZ_), intent(in) :: h
    type(ocean_OBC_type),               pointer    :: OBC
    type(GOLD_generic_tracer_CS),                pointer    :: CS
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
    !                 register_GOLD_generic_tracer.
    !  (in/out)  sponge_CSp - A pointer to the control structure for the sponges, if
    !                         they are in use.  Otherwise this may be unassociated.
    !  (in/out)  diag_to_Z_Csp - A pointer to the control structure for diagnostics
    !                            in depth space.

    character(len=fm_string_len), parameter :: sub_name = 'initialize_GOLD_generic_tracer'
    logical :: OK
    integer :: i, j, k, isc, iec, jsc, jec, nk
    type(g_tracer_type), pointer    :: g_tracer,g_tracer_next
    type(g_diag_type)  , pointer    :: g_diag,g_diag_next
    character(len=fm_string_len)      :: g_tracer_name, longname, units
    real, dimension(:,:,:,:), pointer   :: tr_field 
    real, dimension(:,:,:), pointer     :: tr_ptr 
    real,    dimension(G%isd:G%ied, G%jsd:G%jed,1:G%ke) :: grid_tmask
    integer, dimension(G%isd:G%ied, G%jsd:G%jed)        :: grid_kmt

    !! 2010/02/04  Add code to re-initialize Generic Tracers if needed during a model simulation
    !! By default, restart cpio should not contain a Generic Tracer IC file and step below will be skipped.
    !! Ideally, the generic tracer IC file should have the tracers on Z levels.

    isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nk = G%ke

    !Get the tracer list
    if(.NOT. associated(CS%g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")
    !For each tracer name get its  fields
    g_tracer=>CS%g_tracer_list  

    do
      call g_tracer_get_alias(g_tracer,g_tracer_name)
      call g_tracer_get_pointer(g_tracer,g_tracer_name,'field',tr_field)
      tr_ptr => tr_field(:,:,:,1)

      if (.not.restart .or. (CS%tracers_may_reinit .and. &
          .not.query_initialized(tr_ptr, g_tracer_name, CS%restart_CSp))) then

        if (len_trim(CS%IC_file) > 0) then
          !  Read the tracer concentrations from a netcdf file.
          if (.not.file_exists(CS%IC_file)) call GOLD_error(FATAL, &
                  "initialize_GOLD_Generic_tracer: Unable to open "//CS%IC_file)
          if (CS%Z_IC_file) then
            OK = tracer_Z_init(tr_ptr, h, CS%IC_file, g_tracer_name, G)
            if (.not.OK) then
              OK = tracer_Z_init(tr_ptr, h, CS%IC_file, trim(g_tracer_name)//"_z", G)
              if (.not.OK) call GOLD_error(FATAL,"initialize_GOLD_Generic_tracer: "//&
                      "Unable to read "//trim(g_tracer_name)//" from "//&
                      trim(CS%IC_file)//".")
            endif
            call GOLD_error(NOTE,"initialize_GOLD_generic_tracer: "//&
                            "initialized generic tracer "//trim(g_tracer_name)//&
                            " using Generic Tracer File on Z: "//CS%IC_file)
          else
            ! native grid  
            call GOLD_error(NOTE,"initialize_GOLD_generic_tracer: "//&
                  "Using Generic Tracer IC file on native grid "//trim(CS%IC_file)//".")
            call read_data(CS%IC_file, trim(g_tracer_name), tr_ptr, domain=G%Domain%mpp_domain)
          endif
        else
          call GOLD_error(FATAL,"initialize_GOLD_generic_tracer: "//&
                  "check Generic Tracer IC filename "//trim(CS%IC_file)//".")
        endif

      endif

      !traverse the linked list till hit NULL
      call g_tracer_get_next(g_tracer, g_tracer_next)
      if(.NOT. associated(g_tracer_next)) exit
      g_tracer=>g_tracer_next  
    enddo
    !! end section to re-initialize generic tracers


    !Now we can reset the grid mask, axes and time to their true values
    !Note that grid_tmask must be set correctly on the data domain boundary 
    !so that coast mask can be deduced from it.
    grid_tmask(:,:,:) = 0.0
    grid_kmt(:,:) = 0
    do j = G%jsd, G%jed ; do i = G%isd, G%ied  
       if (G%hmask(i,j) .gt. 0) then
          grid_tmask(i,j,:) = 1.0
          grid_kmt(i,j) = G%ke ! Tell the code that a layer thicker than 1m is the bottom layer.
       endif
    enddo ; enddo 
    
    call g_tracer_set_common(G%isc,G%iec,G%jsc,G%jec,G%isd,G%ied,G%jsd,G%jed,G%ke,1,G%axeshl,grid_tmask,grid_kmt,day) 

    ! Register generic tracer modules diagnostics

    call generic_tracer_register_diag()

    ! Register Z diagnostic output.
    !Get the tracer list
    if(.NOT. associated(CS%g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")
    !For each tracer name get its  fields
    g_tracer=>CS%g_tracer_list  
    do
       call g_tracer_get_alias(g_tracer,g_tracer_name)

       call g_tracer_get_pointer(g_tracer,g_tracer_name,'field',tr_field)
       tr_ptr => tr_field(:,:,:,1)
       call g_tracer_get_values(g_tracer,g_tracer_name,'longname', longname)
       call g_tracer_get_values(g_tracer,g_tracer_name,'units',units )

       call register_Z_tracer(tr_ptr, trim(g_tracer_name)//"_z",longname , units, &
            day, G, diag_to_Z_CSp)

       !traverse the linked list till hit NULL
       call g_tracer_get_next(g_tracer, g_tracer_next)
       if(.NOT. associated(g_tracer_next)) exit
       g_tracer=>g_tracer_next  

    enddo

    !For each special diagnostics name get its  fields
    !Get the diag list
    call generic_tracer_get_diag_list(CS%g_diag_list)
    if(associated(CS%g_diag_list)) then
       g_diag=>CS%g_diag_list  
       do    
          if(g_diag%Z_diag .ne. 0) &
               call register_Z_tracer(g_diag%field_ptr, trim(g_diag%name)//"_z",g_diag%longname , g_diag%units, &
               day, G, diag_to_Z_CSp)
          
          !traverse the linked list till hit NULL
          g_diag=>g_diag%next  
          if(.NOT. associated(g_diag)) exit
          
       enddo
    endif

  end subroutine initialize_GOLD_generic_tracer

  ! <SUBROUTINE NAME="GOLD_generic_tracer_column_physics">
  !  <OVERVIEW>
  !   Column physics for generic tracers.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine does:
  !       Get the coupler values for generic tracers that exchange with atmosphere
  !       Update generic tracer concentration fields from sources and sinks.
  !       Vertically diffuse generic tracer concentration fields.
  !       Update generic tracers from bottom and their bottom reservoir.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call GOLD_generic_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, CS, tv, optics)
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine GOLD_generic_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, CS, tv, optics)
    real, dimension(NXMEM_,NYMEM_,NZ_), intent(in) :: h_old, h_new, ea, eb
    type(forcing),                      intent(in) :: fluxes
    real,                               intent(in) :: dt
    type(ocean_grid_type),              intent(in) :: G
    type(GOLD_generic_tracer_CS),       pointer    :: CS
    type(thermo_var_ptrs),              intent(in) :: tv
    type(optics_type),                  intent(in) :: optics
    !   This subroutine applies diapycnal diffusion and any other column
    ! tracer physics or chemistry to the tracers from this file.
    ! CFCs are relatively simple, as they are passive tracers. with only a surface
    ! flux as a source.

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
    !                 register_GOLD_generic_tracer.
    !
    ! The arguments to this subroutine are redundant in that
    !     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]
    character(len=fm_string_len), parameter :: sub_name = 'GOLD_generic_tracer_column_physics'

    type(g_tracer_type), pointer  :: g_tracer, g_tracer_next
    character(len=fm_string_len)  :: g_tracer_name
    real, dimension(:,:), pointer :: stf_array,trunoff_array,runoff_tracer_flux_array

    real, dimension(G%isd:G%ied,G%jsd:G%jed,G%ke) :: rho_dzt, dzt
    real, dimension(G%isd:G%ied,G%jsd:G%jed)      :: hblt_depth
    integer :: i, j, k, isc, iec, jsc, jec, nk

    isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nk = G%ke

    !Get the tracer list
    if(.NOT. associated(CS%g_tracer_list)) call mpp_error(FATAL,&
         trim(sub_name)//": No tracer in the list.")
    !
    !Update the fields of the generic tracers from Tr(n)
    !
    !nnz: This step is not necessary in GOLD code because the tracers
    !are "registered" for restart and advection. I.e., the ocean model
    !has the pointers to %field arrays and use them directly.
    !I think this is a much better design than T_prog(n) of MOM which
    !get updated by ocean model and have to be copied back and for
    !to the tracers fields. 
    !The problem with MOM is T_prog%field is not a pointer but is an allocatable array. 


    !
    !Extract the tracer surface fields from coupler and update tracer fields from sources 
    !
    call generic_tracer_coupler_get(fluxes%tr_fluxes)

    !
    !Add contribution of river to surface flux
    !
    g_tracer=>CS%g_tracer_list  
    do
       if(_ALLOCATED(g_tracer%trunoff)) then
          call g_tracer_get_alias(g_tracer,g_tracer_name)
          call g_tracer_get_pointer(g_tracer,g_tracer_name,'stf',   stf_array)
          call g_tracer_get_pointer(g_tracer,g_tracer_name,'trunoff',trunoff_array)
          call g_tracer_get_pointer(g_tracer,g_tracer_name,'runoff_tracer_flux',runoff_tracer_flux_array)
          !nnz: Why is fluxes%river = 0?
          runoff_tracer_flux_array = trunoff_array * fluxes%liq_runoff
          stf_array = stf_array + runoff_tracer_flux_array
       endif

       !traverse the linked list till hit NULL
       call g_tracer_get_next(g_tracer, g_tracer_next)
       if(.NOT. associated(g_tracer_next)) exit
       g_tracer=>g_tracer_next  

    enddo

    !
    !Prepare input arrays for source update
    !

    rho_dzt(:,:,:) = G%H_to_kg_m2 * G%Angstrom
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       rho_dzt(i,j,k) = G%H_to_kg_m2 * h_old(i,j,k)
    enddo; enddo ; enddo !}

    dzt(:,:,:) = 1.0
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       dzt(i,j,k) = G%H_to_m * h_old(i,j,k)
    enddo; enddo ; enddo !}


    ! Boussinesq model 
    hblt_depth(:,:) = G%H_to_m * G%Angstrom
    do j=jsc,jec ; do i=isc,iec ; 
       hblt_depth(i,j) = G%H_to_m * h_old(i,j,1)
    enddo; enddo
    do k=2,CS%nkml ; do j=jsc,jec ; do i=isc,iec
       hblt_depth(i,j) = hblt_depth(i,j) + G%H_to_m * h_old(i,j,k)
    enddo; enddo ; enddo


    !
    !Calculate tendencies (i.e., field changes at dt) from the sources / sinks
    !

    call generic_tracer_source(tv%T,tv%S,rho_dzt,dzt,hblt_depth,G%isd,G%jsd,1,dt,&
         G%dxdyh,CS%diag%time_end,&
         optics%nbands, optics%max_wavelength_band, optics%sw_pen_band, optics%opacity_band)

    !
    !Update Tr(n)%field from explicit vertical diffusion
    !

    ! Use a tridiagonal solver to determine the concentrations after the
    ! surface source is applied and diapycnal advection and diffusion occurs.

    call generic_tracer_vertdiff_G(h_old, ea, eb, dt, G%kg_m2_to_H, G%m_to_H, 1) !Last arg is tau which is always 1 for GOLD

    ! Update bottom fields after vertical processes

    call generic_tracer_update_from_bottom(dt, 1, CS%diag%time_end) !Second arg is tau which is always 1 for GOLD

  end subroutine GOLD_generic_tracer_column_physics

  ! <SUBROUTINE NAME="GOLD_generic_tracer_stock">
  !  <OVERVIEW>
  !   Calculate the mass-weighted integral of tracer concentrations.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !     This subroutine calculates mass-weighted integral on the PE either
  !   of all available tracer concentrations, or of a tracer that is 
  !   being requested specifically, returning the number of stocks it has
  !   calculated.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   ns = GOLD_generic_tracer_stock(h, stocks, G, CS, names, units, stock_index)
  !  </TEMPLATE>
  ! </SUBROUTINE>

  function GOLD_generic_tracer_stock(h, stocks, G, CS, names, units, stock_index)
    real, dimension(NXMEM_,NYMEM_,NZ_), intent(in)    :: h
    real, dimension(:),                 intent(out)   :: stocks
    type(ocean_grid_type),              intent(in)    :: G
    type(GOLD_generic_tracer_CS),       pointer       :: CS
    character(len=*), dimension(:),     intent(out)   :: names
    character(len=*), dimension(:),     intent(out)   :: units
    integer, optional,                  intent(in)    :: stock_index
    integer                                           :: GOLD_generic_tracer_stock
  ! This function calculates the mass-weighted integral of all tracer stocks,
  ! returning the number of stocks it has calculated.  If the stock_index
  ! is present, only the stock corresponding to that coded index is returned.

  ! Arguments: h - Layer thickness, in m or kg m-2.
  !  (out)     stocks - the mass-weighted integrated amount of each tracer,
  !                     in kg times concentration units.
  !  (in)      G - The ocean's grid structure.
  !  (in)      CS - The control structure returned by a previous call to
  !                 register_GOLD_generic_tracer.
  !  (out)     names - the names of the stocks calculated.
  !  (out)     units - the units of the stocks calculated.
  !  (in,opt)  stock_index - the coded index of a specific stock being sought.
  ! Return value: the number of stocks calculated here.
    type(g_tracer_type), pointer  :: g_tracer, g_tracer_next
    real, dimension(:,:,:,:), pointer   :: tr_field 
    real, dimension(:,:,:), pointer     :: tr_ptr 
    character(len=fm_string_len), parameter :: sub_name = 'GOLD_generic_tracer_stock'
 
    integer :: i, j, k, is, ie, js, je, nz, m
    is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

    GOLD_generic_tracer_stock = 0
    if (.not.associated(CS)) return

    if (present(stock_index)) then ; if (stock_index > 0) then
      ! Check whether this stock is available from this routine.

      ! No stocks from this routine are being checked yet.  Return 0.
      return
    endif ; endif

    if(.NOT. associated(CS%g_tracer_list)) return ! No stocks.

    m=1 ; g_tracer=>CS%g_tracer_list
    do
      call g_tracer_get_alias(g_tracer,names(m))
      call g_tracer_get_values(g_tracer,names(m),'units',units(m))
      units(m) = trim(units(m))//" kg"
      call g_tracer_get_pointer(g_tracer,names(m),'field',tr_field)

      stocks(m) = 0.0
      tr_ptr => tr_field(:,:,:,1)
      do k=1,nz ; do j=js,je ; do i=is,ie
        stocks(m) = stocks(m) + tr_ptr(i,j,k) * &
                               (G%hmask(i,j) * G%DXDYh(i,j) * h(i,j,k))
      enddo ; enddo ; enddo
      stocks(m) = G%H_to_kg_m2 * stocks(m)

      !traverse the linked list till hit NULL
      call g_tracer_get_next(g_tracer, g_tracer_next)
      if(.NOT. associated(g_tracer_next)) exit
      g_tracer=>g_tracer_next
      m = m+1
    enddo

    GOLD_generic_tracer_stock = m

  end function GOLD_generic_tracer_stock

  ! <SUBROUTINE NAME="GOLD_generic_tracer_min_max">
  !  <OVERVIEW>
  !   Find the min and max of tracer concentrations.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !     This subroutine find the global min and max of either
  !   of all available tracer concentrations, or of a tracer that is 
  !   being requested specifically, returning the number of tracers it has gone through.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   ns = GOLD_generic_tracer_min_max(do_minmax, gmin, gmax, igmin, jgmin, kgmin, igmax, jgmax, kgmax , G, CS, names, units)
  !  </TEMPLATE>
  ! </SUBROUTINE>

  function GOLD_generic_tracer_min_max(ind_start, got_minmax, gmin, gmax, xgmin, ygmin, zgmin, xgmax, ygmax, zgmax , G, CS, names, units)
    use mpp_utilities_mod, only: mpp_array_global_min_max
    integer,                            intent(in)    :: ind_start
    logical, dimension(:),              intent(out)   :: got_minmax
    real, dimension(:),                 intent(out)   :: gmin,gmax
    real, dimension(:),                 intent(out)   :: xgmin, ygmin, zgmin, xgmax, ygmax, zgmax
    type(ocean_grid_type),              intent(in)    :: G
    type(GOLD_generic_tracer_CS),       pointer       :: CS
    character(len=*), dimension(:),     intent(out)   :: names
    character(len=*), dimension(:),     intent(out)   :: units
    integer                                           :: GOLD_generic_tracer_min_max
  ! This function calculates the mass-weighted integral of all tracer stocks,
  ! returning the number of stocks it has calculated.  If the stock_index
  ! is present, only the stock corresponding to that coded index is returned.

  ! Arguments: h - Layer thickness, in m or kg m-2.
  !  (out)     gmin , gmax - the global minimum and maximum of each tracer,
  !                     in kg times concentration units.
  !  (in)      G - The ocean's grid structure.
  !  (in)      CS - The control structure returned by a previous call to
  !                 register_GOLD_generic_tracer.
  !  (out)     names - the names of the stocks calculated.
  !  (out)     units - the units of the stocks calculated.
  !  (in,opt)  trace_index - the coded index of a specific tracer being sought.
  ! Return value: the number of tracers done here.

    type(g_tracer_type), pointer  :: g_tracer, g_tracer_next
    real, dimension(:,:,:,:), pointer   :: tr_field 
    real, dimension(:,:,:), pointer     :: tr_ptr 
    character(len=fm_string_len), parameter :: sub_name = 'GOLD_generic_tracer_min_max'

    real, dimension(:,:,:),pointer :: grid_tmask
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau

    integer :: i, j, k, is, ie, js, je, nz, m
    real, allocatable, dimension(:) :: geo_z
    
    is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

    GOLD_generic_tracer_min_max = 0
    if (.not.associated(CS)) return

    if(.NOT. associated(CS%g_tracer_list)) return ! No stocks.


    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask) 
    
    !Don't know how to get the depth corresponding to k
    !Fake it for now.
    allocate(geo_z(nk))
    do k=1,nk ;  geo_z(k) = G%Rlay(k) ; enddo

       
    m=ind_start ; g_tracer=>CS%g_tracer_list
    do
      call g_tracer_get_alias(g_tracer,names(m))
      call g_tracer_get_values(g_tracer,names(m),'units',units(m))
      units(m) = trim(units(m))//" kg"
      call g_tracer_get_pointer(g_tracer,names(m),'field',tr_field)

      gmin(m) = -1.0
      gmax(m) = -1.0

      tr_ptr => tr_field(:,:,:,1)


      call mpp_array_global_min_max(tr_ptr, grid_tmask,isd,jsd,isc,iec,jsc,jec,nk , gmin(m), gmax(m), &
                                    G%geolonh,G%geolath,geo_z,xgmin(m), ygmin(m), zgmin(m), xgmax(m), ygmax(m), zgmax(m))

      got_minmax(m) = .true.


      !traverse the linked list till hit NULL
      call g_tracer_get_next(g_tracer, g_tracer_next)
      if(.NOT. associated(g_tracer_next)) exit
      g_tracer=>g_tracer_next
      m = m+1
    enddo

    GOLD_generic_tracer_min_max = m

  end function GOLD_generic_tracer_min_max


  ! <SUBROUTINE NAME="GOLD_generic_tracer_surface_state">
  !  <OVERVIEW>
  !   Calculate the surface state and set coupler values
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine calculates the surface state and set coupler values for 
  !   those generic tracers that havd flux exchange with atmosphere.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call GOLD_generic_tracer_surface_state(state, h, G, CS) 
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine GOLD_generic_tracer_surface_state(state, h, G, CS)
    type(surface),                      intent(inout) :: state
    real, dimension(NXMEM_,NYMEM_,NZ_), intent(in) :: h
    type(ocean_grid_type),              intent(in) :: G
    type(GOLD_generic_tracer_CS),                pointer    :: CS
    !   This subroutine sets up the fields that the coupler needs to calculate the
    ! CFC fluxes between the ocean and atmosphere.
    ! Arguments: state - A structure containing fields that describe the
    !                    surface state of the ocean.
    !  (in)      h - Layer thickness, in m.
    !  (in)      G - The ocean's grid structure.
    !  (in)      CS - The control structure returned by a previous call to
    !                 register_GOLD_generic_tracer.

    character(len=fm_string_len), parameter :: sub_name = 'GOLD_generic_tracer_surface_state'
    real, dimension(G%isd:G%ied,G%jsd:G%jed,1:G%ke,1) :: rho0
    type(g_tracer_type), pointer :: g_tracer 

    !Set coupler values
    !nnz: fake rho0
    rho0=1.0

    call generic_tracer_coupler_set(state%tr_fields,&
         ST=state%SST,&
         SS=state%SSS,&
         rho=rho0,& !nnz: required for MOM
         ilb=G%isd, jlb=G%jsd,&
         tau=1)

    !Output diagnostics via diag_manager for all tracers in this module
    if(.NOT. associated(CS%g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         "No tracer in the list.")
    call g_tracer_send_diag(CS%g_tracer_list, CS%diag%time_end, tau=1)

  end subroutine GOLD_generic_tracer_surface_state

!ALL PE subroutine on Ocean!  Due to otpm design the fluxes should be initialized like this on ALL PE's!
  subroutine GOLD_generic_flux_init

    integer :: ind
    character(len=fm_string_len)   :: g_tracer_name,longname, package,units,old_package,file_in,file_out
    real :: const_init_value
    character(len=fm_string_len), parameter :: sub_name = 'GOLD_generic_flux_init'
    type(g_tracer_type), pointer :: g_tracer_list,g_tracer,g_tracer_next

    if (.not. g_registered) then 
       call generic_tracer_register
       g_registered = .true.
    endif

    call generic_tracer_get_list(g_tracer_list)
    if(.NOT. associated(g_tracer_list)) then
       call mpp_error(WARNING, trim(sub_name)// ": No generic tracer in the list.")
       return
    endif

    g_tracer=>g_tracer_list  
    do

       call g_tracer_flux_init(g_tracer)
    
       !traverse the linked list till hit NULL
       call g_tracer_get_next(g_tracer, g_tracer_next)
       if(.NOT. associated(g_tracer_next)) exit
       g_tracer=>g_tracer_next  

    enddo
    
  end subroutine GOLD_generic_flux_init



  subroutine GOLD_generic_tracer_get(name,member,array, CS)
    character(len=*),         intent(in)  :: name
    character(len=*),         intent(in)  :: member
    real, dimension(:,:,:),   intent(out) :: array
    type(GOLD_generic_tracer_CS), pointer :: CS
    !  (in)      CS - The control structure returned by a previous call to
    !                 register_GOLD_generic_tracer.
   
    real, dimension(:,:,:),   pointer :: array_ptr
    character(len=fm_string_len), parameter :: sub_name = 'GOLD_generic_tracer_get'

    call g_tracer_get_pointer(CS%g_tracer_list,name,member,array_ptr)
    array(:,:,:) = array_ptr(:,:,:)

  end subroutine GOLD_generic_tracer_get 

  ! <SUBROUTINE NAME="end_GOLD_generic_tracer">
  !  <OVERVIEW>
  !   Ends the generic tracer module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Call the end for generic tracer module and deallocate all temp arrays
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call end_GOLD_generic_tracer(CS) 
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine end_GOLD_generic_tracer(CS)
    type(GOLD_generic_tracer_CS), pointer :: CS
    !   This subroutine deallocates the memory owned by this module.
    ! Argument: CS - The control structure returned by a previous call to
    !                register_GOLD_generic_tracer.
    integer :: m
    
    call generic_tracer_end

    if (associated(CS)) then
       deallocate(CS)
    endif
  end subroutine end_GOLD_generic_tracer

#endif /* _USE_GENERIC_TRACER */
end module GOLD_generic_tracer
