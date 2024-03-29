module GOLD_interface_heights

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
!*  By Robert Hallberg, February 2001                                  *
!*                                                                     *
!*    The subroutines here calculate the heights of interfaces in a    *
!*  way that is consistent with the calculations in the pressure       *
!*  gradient acceleration calculation.  In a Boussinseq model this is  *
!*  pretty simple vertical sum, but in a non-Boussinesq model it uses  *
!*  integrals of the equation of state.                                *
!*                                                                     *
!*  Macros written all in capital letters are defined in GOLD_memory.h.*
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, f                                     *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, D                                     *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_error_handler, only : GOLD_error, FATAL
use GOLD_file_parser, only : log_version
use GOLD_grid, only : ocean_grid_type
use GOLD_variables, only : thermo_var_ptrs
use GOLD_EOS, only : int_specific_vol_dp

implicit none ; private

#include <GOLD_memory.h>

public find_eta

interface find_eta
  module procedure find_eta_2d, find_eta_3d
end interface find_eta

contains

subroutine find_eta_3d(h, tv, G_Earth, G, eta, eta_bt, halo_size)
  real, dimension(NXMEM_,NYMEM_,NZ_),       intent(in)  :: h
  type(thermo_var_ptrs),                    intent(in)  :: tv
  real,                                     intent(in)  :: G_Earth
  type(ocean_grid_type),                    intent(in)  :: G
  real, dimension(NXMEM_,NYMEM_,NZp1_),     intent(out) :: eta
  real, dimension(NXMEM_,NYMEM_), optional, intent(in)  :: eta_bt
  integer,                        optional, intent(in)  :: halo_size
!   This subroutine determines the heights of all interfaces between layers,
! using the appropriate form for consistency with the calculation of the
! pressure gradient forces.  Additionally, these height may be dilated for
! consistency with the corresponding time-average quantity from the barotropic
! calculation.

! Arguments: h - Layer thickness, in m.  Intent in.
!  (in)      tv - a structure pointing to various thermodynamic variables.
!  (in)      G_Earth - The Earth's gravitational acceleration, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (out)     eta - The free surface height realtive to mean sea level, in m.
!  (in,opt)  eta_bt - An optional barotropic variable that gives the "correct"
!                     free surface height (Boussinesq) or total water column
!                     mass per unit aread (non-Boussinesq).  This is used to
!                     dilate the layer thicknesses when calculating interface
!                     heights, in m or kg m-2.
!  (in,opt)  halo_size - The width of halo points on which to calculate eta.

  real :: p(SZ1_(h),SZ2_(h),SZ3_(h)+1)
  real :: dz_geo(SZ1_(h),SZ2_(h))      ! The change in geopotential height
                                       ! across a layer, in m2 s-2.
  real :: dilate(SZ1_(h)), htot(SZ1_(h))
  real :: I_gEarth
  integer i, j, k, isv, iev, jsv, jev, nz, halo

  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  
  isv = G%isc-halo ; iev = G%iec+halo ; jsv = G%jsc-halo ; jev = G%jec+halo
  nz = G%ke

  if ((isv<G%isd) .or. (iev>G%ied) .or. (jsv<G%jsd) .or. (jev>G%jed)) &
    call GOLD_error(FATAL,"find_eta called with an overly large halo_size.")

  I_gEarth = 1.0 / G_Earth

  do j=jsv,jev ; do i=isv,iev ; eta(i,j,nz+1) = -G%D(i,j) ; enddo ; enddo

  if (G%Boussinesq) then
    do k=nz,1,-1 ; do j=jsv,jev ; do i=isv,iev
      eta(i,j,k) = eta(i,j,k+1) + h(i,j,k)
    enddo ; enddo ; enddo
    if (present(eta_bt)) then
      ! Dilate the water column to agree with the free surface height
      ! that is used for the dynamics.
      do j=jsv,jev
        do i=isv,iev
          dilate(i) = (eta_bt(i,j) + G%D(i,j)) / (eta(i,j,1) + G%D(i,j))
        enddo
        do k=1,nz ; do i=isv,iev
          eta(i,j,k) = dilate(i) * (eta(i,j,k) + G%D(i,j)) - G%D(i,j)
        enddo ; enddo
      enddo
    endif
  else
    if (associated(tv%eqn_of_state)) then
      ! ### THIS SHOULD BE P_SURF, IF AVAILABLE.
      do j=jsv,jev ; do i=isv,iev ; p(i,j,1) = 0.0 ; enddo ; enddo
      do k=2,nz+1 ; do j=jsv,jev ; do i=isv,iev
        p(i,j,k) = p(i,j,k-1) + G_Earth*G%H_to_kg_m2*h(i,j,k-1)
      enddo ; enddo ; enddo

      do k=nz,1,-1
        call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p(:,:,k), p(:,:,k+1), &
                                 0.0, G, tv%eqn_of_state, dz_geo, halo_size=halo)
        do j=jsv,jev ; do i=isv,iev
          eta(i,j,k) = eta(i,j,k+1) + I_gEarth * dz_geo(i,j)
        enddo ; enddo
      enddo
    elseif (associated(tv%Rml)) then
      do k=nz,tv%nk_Rml+1,-1 ; do j=jsv,jev ; do i=isv,iev
        eta(i,j,k) = eta(i,j,k+1) + G%H_to_kg_m2*h(i,j,k)/G%Rlay(k)
      enddo ; enddo ; enddo
      do k=tv%nk_Rml,1,-1 ; do j=jsv,jev ; do i=isv,iev
        eta(i,j,k) = eta(i,j,k+1) + G%H_to_kg_m2*h(i,j,k)/tv%Rml(i,j,k)
      enddo ; enddo ; enddo
    else
      do k=nz,1,-1 ; do j=jsv,jev ; do i=isv,iev
        eta(i,j,k) = eta(i,j,k+1) + G%H_to_kg_m2*h(i,j,k)/G%Rlay(k)
      enddo ; enddo ; enddo
    endif
    if (present(eta_bt)) then
      ! Dilate the water column to agree with the free surface height
      ! from the time-averaged barotropic solution.
      do j=jsv,jev
        do i=isv,iev ; htot(i) = G%H_subroundoff ; enddo
        do k=1,nz ; do i=isv,iev ; htot(i) = htot(i) + h(i,j,k) ; enddo ; enddo
        do i=isv,iev ; dilate(i) = eta_bt(i,j) / htot(i) ; enddo
        do k=1,nz ; do i=isv,iev
          eta(i,j,k) = dilate(i) * (eta(i,j,k) + G%D(i,j)) - G%D(i,j)
        enddo ; enddo
      enddo
    endif
  endif

end subroutine find_eta_3d

subroutine find_eta_2d(h, tv, G_Earth, G, eta, eta_bt, halo_size)
  real, dimension(NXMEM_,NYMEM_,NZ_),       intent(in)  :: h
  type(thermo_var_ptrs),                    intent(in)  :: tv
  real,                                     intent(in)  :: G_Earth
  type(ocean_grid_type),                    intent(in)  :: G
  real, dimension(NXMEM_,NYMEM_),           intent(out) :: eta
  real, dimension(NXMEM_,NYMEM_), optional, intent(in)  :: eta_bt
  integer,                        optional, intent(in)  :: halo_size
!   This subroutine determines the free surface height, using the appropriate
! form for consistency with the calculation of the pressure gradient forces.
! Additionally, the sea surface height may be adjusted for consistency with the
! corresponding time-average quantity from the barotropic calculation.

! Arguments: h - Layer thickness, in m.  Intent in.
!  (in)      tv - a structure pointing to various thermodynamic variables.
!  (in)      G_Earth - The Earth's gravitational acceleration, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (out)     eta - The free surface height relative to mean sea level, in m.
!  (in,opt)  eta_bt - An optional barotropic variable that gives the "correct"
!                     free surface height (Boussinesq) or total water column
!                     mass per unit aread (non-Boussinesq), in m or kg m-2.
!  (in,opt)  halo_size - The width of halo points on which to calculate eta.

  real, dimension(SZ1_(h),SZ2_(h)) :: &
    p_top, &   ! The pressure at the interface above a layer, in Pa.
    p_bot, &   ! The pressure at the interface below a layer, in Pa.
    dz_geo     ! The change in geopotential height across a layer, in m2 s-2.
  real :: htot(SZ1_(h))  ! The sum of all layers' thicknesses, in kg m-2 or m.
  real :: I_gEarth
  integer i, j, k, is, ie, js, je, nz, halo

  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  nz = G%ke

  I_gEarth = 1.0 / G_Earth

  do j=js,je ; do i=is,ie ; eta(i,j) = -G%D(i,j) ; enddo ; enddo

  if (G%Boussinesq) then
    if (present(eta_bt)) then
      do j=js,je ; do i=is,ie
        eta(i,j) = eta_bt(i,j)
      enddo ; enddo
    else
      do k=1,nz ; do j=js,je ; do i=is,ie
        eta(i,j) = eta(i,j) + h(i,j,k)
      enddo ; enddo ; enddo
    endif
  else
    if (associated(tv%eqn_of_state)) then
      do j=js,je ; do i=is,ie
        p_top(i,j) = 0.0 ; p_bot(i,j) = 0.0
      enddo ; enddo
      do k=1,nz
        do j=js,je ; do i=is,ie
          p_top(i,j) = p_bot(i,j)
          p_bot(i,j) = p_top(i,j) + G_Earth*G%H_to_kg_m2*h(i,j,k)
        enddo ; enddo
        call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p_top, p_bot, 0.0, &
                                 G, tv%eqn_of_state, dz_geo, halo_size=halo)
        do j=js,je ; do i=is,ie
          eta(i,j) = eta(i,j) + I_gEarth * dz_geo(i,j)
        enddo ; enddo
      enddo
    elseif (associated(tv%Rml)) then
      do k=1,tv%nk_Rml ; do j=js,je ; do i=is,ie
        eta(i,j) = eta(i,j) + G%H_to_kg_m2*h(i,j,k)/tv%Rml(i,j,k)
      enddo ; enddo ; enddo
      do k=tv%nk_Rml+1,nz ; do j=js,je ; do i=is,ie
        eta(i,j) = eta(i,j) + G%H_to_kg_m2*h(i,j,k)/G%Rlay(k)
      enddo ; enddo ; enddo
    else
      do k=1,nz ; do j=js,je ; do i=is,ie
        eta(i,j) = eta(i,j) + G%H_to_kg_m2*h(i,j,k)/G%Rlay(k)
      enddo ; enddo ; enddo
    endif
    if (present(eta_bt)) then
      !   Dilate the water column to agree with the the time-averaged column
      ! mass from the barotropic solution.
      do j=js,je
        do i=is,ie ; htot(i) = G%H_subroundoff ; enddo
        do k=1,nz ; do i=is,ie ; htot(i) = htot(i) + h(i,j,k) ; enddo ; enddo
        do i=is,ie
          eta(i,j) = (eta_bt(i,j) / htot(i)) * (eta(i,j) + G%D(i,j)) - G%D(i,j)
        enddo
      enddo
    endif
  endif
  
end subroutine find_eta_2d

end module GOLD_interface_heights
