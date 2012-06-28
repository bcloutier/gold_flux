module GOLD_sponge
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
!*  By Robert Hallberg, March 1999-June 2000                           *
!*                                                                     *
!*    This program contains the subroutines that implement sponge      *
!*  regions, in which the stratification and water mass properties     *
!*  are damped toward some profiles.  There are three externally       *
!*  callable subroutines in this file.                                 *
!*                                                                     *
!*    initialize_sponge determines the mapping from the model          *
!*  variables into the arrays of damped columns.  This remapping is    *
!*  done for efficiency and to conserve memory.  Only columns which    *
!*  have positive inverse damping times and which are deeper than a    *
!*  supplied depth are placed in sponges.  The inverse damping         *
!*  time is also stored in this subroutine, and memory is allocated    *
!*  for all of the reference profiles which will subsequently be       *
!*  provided through calls to set_up_sponge_field.  The three          *
!*  arguments are a two-dimensional array containing the damping       *
!*  rates, the number of reference profiles which will be provided,    *
!*  and the minimum bottom depth over which to apply the damping.      *
!*                                                                     *
!*    set_up_sponge_field is called to provide a reference profile     *
!*  and the location of the field that will be damped back toward      *
!*  that reference profile.  A third argument, the number of layers    *
!*  in the field is also provided, but currently only Rml is set to    *
!*  be applied at fewer than nz layers.                                *
!*    The first call to set_up_sponge_field must be to register the    *
!*  interface heights, and the location of the field for the first     *
!*  call is ignored, since the location of the thickness can vary.     *
!*  If BULKMIXEDLAYER is defined, the second call to set_up_sponge-    *
!*  _field must be for Rml.  The order of any other calls to           *
!*  set_up_sponge_field does not matter.                               *
!*                                                                     *
!*    Apply_sponge damps all of the fields that have been registered   *
!*  with set_up_sponge_field toward their reference profiles.  The     *
!*  four arguments are the thickness to be damped, the amount of time  *
!*  over which the damping occurs, and arrays to which the movement    *
!*  of fluid into a layer from above and below will be added. The      *
!*  effect on momentum of the sponge may be accounted for later using  *
!*  the movement of water recorded in these later arrays.              *
!*                                                                     *
!*    All of the variables operated upon in this file are defined at   *
!*  the thickness points.                                              *
!*                                                                     *
!*  Macros written all in capital letters are defined in GOLD_memory.h.*
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, T, S, Rml, Iresttime, ea, eb          *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_coms, only : sum_across_PEs
use GOLD_error_handler, only : GOLD_error, FATAL, NOTE, WARNING, is_root_pe
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_grid, only : ocean_grid_type

! Planned extension:  Support for time varying sponge targets.

implicit none ; private

#include <GOLD_memory.h>

public set_up_sponge_field, initialize_sponge, apply_sponge

type :: p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d
type :: p2d
  real, dimension(:,:), pointer :: p => NULL()
end type p2d

type, public :: sponge_CS ; private
  logical :: bulkmixedlayer  ! If true, a refined bulk mixed layer is used with
                             ! nkml sublayers and nkbl buffer layer.
  integer :: nkml            ! The number of mixed layer sublayers.
  integer :: nkbl            ! The number of buffer layer sublayers.
  integer :: nz              ! The total number of layers.
  integer :: num_col         ! The number of sponge points within the
                             ! computational domain.
  integer :: fldno = 0       ! The number of fields which have already been
                             ! registered by calls to set_up_sponge_field
  integer, pointer :: col_i(:) => NULL()  ! Arrays containing the i- and j- indicies
  integer, pointer :: col_j(:) => NULL()  ! of each of the columns being damped.
  real, pointer :: Iresttime_col(:) => NULL()  ! The inverse restoring time of
                             ! each column.
  type(p3d) :: var(MAX_FIELDS)  ! Pointers to the fields that are being damped.
  type(p2d) :: Ref_val(MAX_FIELDS)  ! The values to which the fields are damped.
end type sponge_CS

contains

subroutine initialize_sponge(Iresttime, int_height, G, param_file, CS)
  real,                  intent(in) :: Iresttime(NXMEM_,NYMEM_)
  real,                  intent(in) :: int_height(NXMEM_,NYMEM_,NZp1_)
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  type(sponge_CS),       pointer    :: CS
! This subroutine determines the number of points which are within
! sponges in this computational domain.  Only points that have
! positive values of Iresttime and which hmask indicates are ocean
! points are included in the sponges.

! Arguments: Iresttime - The inverse of the restoring time, in s-1.
!  (in)      int_height - The interface heights to damp back toward, in m.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
  character(len=128) :: version = '$Id: GOLD_sponge.F90,v 13.0.2.1.2.10 2011/10/07 22:29:44 Robert.Hallberg Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "GOLD_sponge"  ! This module's name.
  logical :: use_sponge
  integer :: i, j, k, col, total_sponge_cols

  if (associated(CS)) then
    call GOLD_error(WARNING, "initialize_sponge called with an associated "// &
                            "control structure.")
    return
  endif

  use_sponge = .false. ; call read_param(param_file,"SPONGE",use_sponge)

  if (.not.use_sponge) return
  allocate(CS)

  CS%nz = G%ke
  CS%bulkmixedlayer = .false.
  call read_param(param_file,"BULKMIXEDLAYER",CS%bulkmixedlayer)
  if (CS%bulkmixedlayer) then
    CS%nkml = 1 ; call read_param(param_file,"NKML",CS%nkml)
    CS%nkbl = 1 ; call read_param(param_file,"NKBL",CS%nkbl)
  endif

  CS%num_col = 0 ; CS%fldno = 1
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((Iresttime(i,j)>0.0) .and. (G%hmask(i,j)>0)) &
      CS%num_col = CS%num_col + 1
  enddo ; enddo

  if (CS%num_col > 0) then

    allocate(CS%Iresttime_col(CS%num_col)) ; CS%Iresttime_col = 0.0
    allocate(CS%col_i(CS%num_col))         ; CS%col_i = 0
    allocate(CS%col_j(CS%num_col))         ; CS%col_j = 0

    col = 1
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if ((Iresttime(i,j)>0.0) .and. (G%hmask(i,j)>0)) then
        CS%col_i(col) = i ; CS%col_j(col) = j
        CS%Iresttime_col(col) = Iresttime(i,j)
        col = col +1
      endif
    enddo ; enddo

    allocate(CS%Ref_val(1)%p(CS%nz+1,CS%num_col))
    do col=1,CS%num_col ; do k=1,CS%nz+1
      CS%Ref_val(1)%p(k,col) = int_height(CS%col_i(col),CS%col_j(col),k)
    enddo ; enddo

  endif

  total_sponge_cols = CS%num_col
  call sum_across_PEs(total_sponge_cols)

  ! Write all relevant parameters to the model log.
  call log_version(param_file, mod, version, tagname)
  call log_param(param_file, mod, "SPONGE", use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. \n"//&
                 "The exact location and properties of those sponges are \n"//&
                 "specified from GOLD_initialization.F90.", default=.false.)
  call log_param(param_file, mod, "BULKMIXEDLAYER", CS%bulkmixedlayer, &
                 "If defined, use a refined Kraus-Turner-like bulk mixed \n"//&
                 "layer with transitional buffer layers.", default=.false.)
  call log_param(param_file, mod, "Total sponge columns", total_sponge_cols, &
                 "The total number of columns where sponges are applied.")
  if (CS%bulkmixedlayer) then
    call log_param(param_file, mod, "NKML", CS%nkml, &
                 "The number of sublayers within the mixed layer if \n"//&
                 "BULKMIXEDLAYER is true.", units="nondim", default=1)
    call log_param(param_file, mod, "NKBL", CS%nkbl, &
                 "The number of layers that are used as variable density \n"//&
                 "buffer layers if BULKMIXEDLAYER is true.", units="nondim", &
                 default=1)
  endif

end subroutine initialize_sponge

subroutine set_up_sponge_field(sp_val, f_ptr, nlay, CS)
  real,         intent(in) :: sp_val(NXMEM_,NYMEM_,NZ_)
  real, target, intent(in) :: f_ptr(NXMEM_,NYMEM_,NZ_)
  integer,      intent(in) :: nlay
  type(sponge_CS), pointer :: CS
! This subroutine stores the reference profile for the variable
! whose address is given by f_ptr. nlay is the number of layers in
! this variable.  The first call to this subroutine must be to
! register the interface depth profiles, and the second must be to
! register the mixed layer buoyancy (Rml) profiles if BULKMIXEDLAYER
! is defined.  Subsequent calls can be made in any order.

! Arguments: sp_val - The reference profiles of the quantity being
!                     registered.
!  (in)      f_ptr - a pointer to the field which will be damped.
!  (in)      nlay - the number of layers in this quantity.
!  (in/out)  CS - A pointer to the control structure for this module that is
!                 set by a previous call to initialize_sponge.

  integer :: k, col
  character(len=256) :: mesg ! String for error messages

  if (.not.associated(CS)) return

  CS%fldno = CS%fldno + 1

  if (CS%fldno > MAX_FIELDS) then
    write(mesg,'("Increase MAX_FIELDS to at least ",I3," in GOLD_memory.h or decrease &
           &the number of fields to be damped in the call to &
           &initialize_sponge." )') CS%fldno
    call GOLD_error(FATAL,"set_up_sponge_field: "//mesg)
  endif

  allocate(CS%Ref_val(CS%fldno)%p(CS%nz,CS%num_col))
  CS%Ref_val(CS%fldno)%p(:,:) = 0.0
  do col=1,CS%num_col
    do k=1,nlay
      CS%Ref_val(CS%fldno)%p(k,col) = sp_val(CS%col_i(col),CS%col_j(col),k)
    enddo
    do k=nlay+1,CS%nz
      CS%Ref_val(CS%fldno)%p(k,col) = 0.0
    enddo
  enddo

  CS%var(CS%fldno)%p => f_ptr

  if ( CS%bulkmixedlayer ) then
    if ((CS%fldno==2) .and. (nlay>1)) then
      if (is_root_pe()) call GOLD_error(NOTE, "set_up_sponge_field: &
           &Caution: Only 1 layer should be used for the reference value of the&
           & mixed layer density.  The second call to set_up_sponge_field &
           & should set the reference mixed layer density.")
    elseif ((CS%fldno>2) .and. (nlay/=CS%nz)) then
      write(mesg,'("Danger:  Sponge reference fields require nz (",I3,") layers,&
          & except for the second call to set_up_sponge_field, which initializes&
          & the reference mixed layer density and requires 1 layer.")' ) &
            CS%nz
      if (is_root_pe()) call GOLD_error(WARNING,"set_up_sponge_field: "//mesg)
    endif
  elseif (nlay/=CS%nz) then
    write(mesg,'("Danger: Sponge reference fields require nz (",I3,") layers.&
        & A field with ",I3," layers was passed to set_up_sponge_field.")') &
          CS%nz, nlay
    if (is_root_pe()) call GOLD_error(WARNING,"set_up_sponge_field: "//mesg)
  endif

end subroutine set_up_sponge_field

subroutine apply_sponge(h, dt, G, ea, eb, CS)
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(inout) :: h
  real,                               intent(in)    :: dt
  type(ocean_grid_type),              intent(in)    :: G
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(out)   :: ea
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(out)   :: eb
  type(sponge_CS),                    pointer       :: CS

! This subroutine applies damping to the layers thicknesses, mixed
! layer buoyancy, and a variety of tracers for every column where
! there is damping.

! Arguments: h -  Layer thickness, in m.
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.
!  (out)     ea - an array to which the amount of fluid entrained
!                 from the layer above during this call will be
!                 added, in m.
!  (out)     eb - an array to which the amount of fluid entrained
!                 from the layer below during this call will be
!                 added, in m.
!  (in)      CS - A pointer to the control structure for this module that is
!                 set by a previous call to initialize_sponge.

  real :: damp     ! The timestep times the local damping
                   ! coefficient.  Nondimensional.
  real :: e(SZ3_(h)+1)  ! The interface heights, in m or kg m-2, usually negative.
  real :: e0       ! The height of the free surface in m.
  real :: e_str    ! A nondimensional amount by which the reference
                   ! profile must be stretched for the free surfaces
                   ! heights in the two profiles to agree.
  real :: w        ! The thickness of water moving upward through an
                   ! interface within 1 timestep, in m.
  real :: wm       ! wm is w if w is negative and 0 otherwise, in m.
  real :: wb       ! w at the interface below a layer, in m.
  real :: wpb      ! wpb is wb if wb is positive and 0 otherwise, m.
  real :: I1pdamp  ! I1pdamp is 1/(1 + damp).  Nondimensional.
  integer :: c, i, j, k, nz, m, nkmb

  if (.not.associated(CS)) return
  nz = G%ke
  if (CS%bulkmixedlayer) nkmb = CS%nkml + CS%nkbl

  do c=1,CS%num_col
! c is an index for the next 3 lines but a multiplier for the rest of the loop
! Therefore we use c as per C code and increment the index where necessary.
    i = CS%col_i(c) ; j = CS%col_j(c)
    damp = dt*CS%Iresttime_col(c)

    e(1) = 0.0 ; e0 = 0.0
    do k=1,nz
      e(k+1) = e(k) - h(i,j,k)
    enddo
    e_str = e(nz+1) / CS%Ref_val(1)%p(nz+1,c)

    if ( CS%bulkmixedlayer ) then
      I1pdamp = 1.0 / (1.0 + damp)
      do k=1,nkmb
        CS%var(2)%p(i,j,k) = I1pdamp * &
            (CS%var(2)%p(i,j,k) + CS%Ref_val(2)%p(1,c)*damp)
      enddo
      do k=1,nkmb
        do m=3,CS%fldno
          CS%var(m)%p(i,j,k) = I1pdamp * &
              (CS%var(m)%p(i,j,k) + CS%Ref_val(m)%p(k,c)*damp)
        enddo
      enddo

      wpb = 0.0; wb = 0.0
      do k=nz,nkmb+1,-1
        if (G%Rlay(k) > CS%var(2)%p(i,j,1)) then
          w = MIN((((e(k)-e0) - e_str*CS%Ref_val(1)%p(k,c)) * damp), &
                    (wb+h(i,j,k) - G%Angstrom))
          wm = 0.5*(w-ABS(w))
          do m=3,CS%fldno
            CS%var(m)%p(i,j,k) = (h(i,j,k)*CS%var(m)%p(i,j,k) + &
                     CS%Ref_val(m)%p(k,c)*(damp*h(i,j,k) + wpb - wm)) / &
                     (h(i,j,k)*(1.0 + damp) + wpb - wm)
          enddo
          eb(i,j,k) = eb(i,j,k) + wpb
          ea(i,j,k) = ea(i,j,k) - wm
          h(i,j,k) = h(i,j,k) + (wb - w)
          wb = w
          wpb = w - wm
        else
          do m=3,CS%fldno
            CS%var(m)%p(i,j,k) = I1pdamp * &
              (CS%var(m)%p(i,j,k) + CS%Ref_val(m)%p(k,c)*damp)
          enddo
          w = wb + (h(i,j,k) - G%Angstrom)
          wm = 0.5*(w-ABS(w))
          eb(i,j,k) = eb(i,j,k) + wpb
          ea(i,j,k) = ea(i,j,k) - wm
          h(i,j,k)  = h(i,j,k)  + (wb - w)
          wb = w
          wpb = w - wm
        endif
      enddo

      if (wb < 0) then
        do k=nkmb,1,-1
          w = MIN((wb + (h(i,j,k) - G%Angstrom)),0.0)
          h(i,j,k)  = h(i,j,k)  + (wb - w)
          ea(i,j,k) = ea(i,j,k) - w
          wb = w
        enddo
      else
        w = wb
        do k=CS%nkml,nkmb
          eb(i,j,k) = eb(i,j,k) + w
        enddo

        k = CS%nkml
        h(i,j,k) = h(i,j,k) + w
        do m=3,CS%fldno
          CS%var(m)%p(i,j,k) = (CS%var(m)%p(i,j,k)*h(i,j,k) + &
                                CS%Ref_val(m)%p(k,c)*w) / (h(i,j,k) + w)
        enddo
      endif

      do k=1,nkmb
        CS%var(2)%p(i,j,k) = I1pdamp * &
            (CS%var(2)%p(i,j,k) + CS%Ref_val(2)%p(1,c)*damp)
        do m=3,CS%fldno
          CS%var(m)%p(i,j,k) = I1pdamp * &
              (CS%var(m)%p(i,j,k) + CS%Ref_val(m)%p(CS%nkml,c)*damp)
        enddo
      enddo

    else                                          ! not BULKMIXEDLAYER

      wpb = 0.0
      wb = 0.0
      do k=nz,1,-1
        w = MIN((((e(k)-e0) - e_str*CS%Ref_val(1)%p(k,c)) * damp), &
                  (wb+h(i,j,k) - G%Angstrom))
        wm = 0.5*(w-ABS(w))
        do m=2,CS%fldno
          CS%var(m)%p(i,j,k) = (h(i,j,k)*CS%var(m)%p(i,j,k) + &
              CS%Ref_val(m)%p(k,c) * (damp*h(i,j,k) + wpb - wm)) / &
                     (h(i,j,k)*(1.0 + damp) + wpb - wm)
        enddo
        eb(i,j,k) = eb(i,j,k) + wpb
        ea(i,j,k) = ea(i,j,k) - wm
        h(i,j,k)  = h(i,j,k)  + (wb - w)
        wb = w
        wpb = w - wm
      enddo

    endif                                         ! end BULKMIXEDLAYER
  enddo ! end of c loop

end subroutine apply_sponge

end module GOLD_sponge
