module GOLD_continuity_Hallberg95
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
!*  By Robert Hallberg, April 1994                                     *
!*                                                                     *
!*    This program contains the subroutine that advects layer          *
!*  thickness.  The scheme is described in Hsu and Arakawa, Mon.       *
!*  Wea. Rev. 1990, and is based in large measure on an earlier        *
!*  advection scheme by Takacs, Mon. Wea. Rev. 1985.  This scheme      *
!*  has been modified by R. Hallberg to be positive definite and is    *
!*  second order accurate in space and time.  In the limit of          *
!*  uniform velocities, the scheme is third order accurate in space    *
!*  and time.  This scheme does a predictor-corrector time step,       *
!*  first zonally, then meridionally.                                  *
!*                                                                     *
!*    continuity takes 7 arguments, described below.  Of these h and   *
!*  and hin may refer to the same location.                            *
!*                                                                     *
!*    The advection scheme is stable for a CFL number (based on the    *
!*  sum of velocities and internal gravity wave group velocities) of   *
!*  less than 1.  Positive definiteness is guaranteed only for a CFL   *
!*  number based on velocities (not wave group velocities) of less     *
!*  than 0.37609, but is not typically a limitation.                   *
!*                                                                     *
!*  Macros written all in capital letters are defined in GOLD_memory.h.*
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, vh                                    *
!*    j    x ^ x ^ x   At >:  u, uh                                    *
!*    j    > o > o >   At o:  h, hin                                   *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use GOLD_diag_mediator, only : time_type, diag_ptrs
use GOLD_domains, only : pass_var, pass_vector
use GOLD_error_handler, only : GOLD_error, FATAL, WARNING
use GOLD_file_parser, only : log_version, param_file_type
use GOLD_grid, only : ocean_grid_type
use GOLD_variables, only : ocean_OBC_type, OBC_SIMPLE

implicit none ; private

#include <GOLD_memory.h>

public continuity_orig, continuity_orig_init, continuity_orig_end

integer :: id_clock_update, id_clock_correct, id_clock_pass

type, public :: continuity_orig_CS ; private
  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
                             ! ocean diagnostic fields.
end type continuity_orig_CS

contains

subroutine continuity_orig(u, v, hin, h, uh, vh, dt, G, CS, uhbt, vhbt, OBC)
  real, intent(in),  dimension(NXMEMQ_,NYMEM_,NZ_) :: u
  real, intent(in),  dimension(NXMEM_,NYMEMQ_,NZ_) :: v
  real, intent(in),  dimension(NXMEM_,NYMEM_,NZ_)  :: hin
  real, intent(out), dimension(NXMEM_,NYMEM_,NZ_)  :: h
  real, intent(out), dimension(NXMEMQ_,NYMEM_,NZ_) :: uh
  real, intent(out), dimension(NXMEM_,NYMEMQ_,NZ_) :: vh
  real, intent(in)                                 :: dt
  type(ocean_grid_type), intent(inout)             :: G
  type(continuity_orig_CS), pointer                :: CS
  real, intent(in), optional, dimension(NXMEMQ_,NYMEM_) :: uhbt
  real, intent(in), optional, dimension(NXMEM_,NYMEMQ_) :: vhbt
  type(ocean_OBC_type), pointer, optional         :: OBC
!    This subroutine time steps the layer thicknesses.
!  A positive definite scheme based on Hsu & Arakawa (1990) is used.

! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      hin - Initial layer thickness, in m.
!  (out)     h - Final layer thickness, in m.
!  (out)     uh - Volume flux through zonal faces = u*h*dy, m3 s-1.
!  (out)     vh - Volume flux through meridional faces = v*h*dx,
!                  in m3 s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 continuity_orig_init.
!  (in)      uhbt - The summed volume flux through zonal faces, m3 s-1.
!  (in)      vhbt - The summed volume flux through meridional faces, m3 s-1.
!  (in)      OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!
!    This subroutine time steps the layer thicknesses.
!  A positive definite scheme based on Hsu & Arakawa (1990) is used.

  real :: hpi(SZI_(G))  ! The predicted layer thickness, in m.
  real :: up(SZIQ_(G))  ! Positive u velocity, in m s-1.
  real :: um(SZIQ_(G))  ! Negative u velocity, in m s-1.
  real :: hpj(SZJ_(G))  ! The predicted layer thickness, in m.
  real :: vp(SZJ_(G))   ! Positive v velocity, in m s-1.
  real :: vm(SZJ_(G))   ! Negative v velocity, in m s-1.
  real, dimension(SZI_(G),SZJ_(G),SZ3_(u)) :: &
    Area_u, Area_v
  real, dimension(SZI_(G),SZJ_(G)) :: &
    uhbt_resid, vhbt_resid
  real, dimension(SZI_(G)) :: &
    Area_sum, u_cor, v_cor, uhbt_res, vhbt_res

  logical, dimension(SZJ_(G)) :: do_uh, do_vh
  logical :: apply_OBC_u = .false., apply_OBC_v = .false.
  real, dimension(SZI_(G),SZJ_(G)) :: h_sum
  real :: vol_out
  real, dimension(SZI_(G)) :: &
    udy_cor, vdx_cor

  real :: a, gp, gm,g0, g1  ! a, gp, gm, g0, & g1 are temporary
                            ! variables in the advection scheme.
  real :: gam               ! 0 for upwind, 1 for third-order.
                            ! Note that this is 1-gamma of Hsu & Arakawa (1990).
  real :: vgeom_v           ! The ratio of the geometric mean velocity to
                            ! the velocity at a point, nondim.
  real :: max_vg_v          ! A limit that is applied to vgeom_v to avoid
                            ! either a very large or a negative effective
                            ! thickness, nondim.
  real :: h_neglect2        ! The square of a thickness that is so small it is
                            ! usually lost in roundoff and can be neglected, in m.
  real, parameter :: C1_6 = 1.0 / 6.0
  real :: Cdt_6             ! Cdt_6 is dt/6.
  real :: Idt               ! Idt is 1/dt.
  integer :: is, ie, js, je, nz
  integer :: i, j, k, it
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(CS)) call GOLD_error(FATAL, &
         "GOLD_continuity_Hallberg95: Module must be initialized before it is used.")
  Cdt_6 = dt*C1_6
  h_neglect2 = G%H_subroundoff**2

  if (present(OBC)) then ; if (associated(OBC)) then
    apply_OBC_u = OBC%apply_OBC_u ; apply_OBC_v = OBC%apply_OBC_v
  endif ; endif

!    First, advect zonally.
  call cpu_clock_begin(id_clock_update)
  do k=1,nz
    do j=js-2,je+2

!    Use upwind differencing as a predictor.
      do i=is-2,ie+1
        up(i) = (u(i,j,k) + ABS(u(i,j,k))) * 0.5
        um(i) = u(i,j,k) - up(i)
        uh(i,j,k) = G%dy_u(i,j) * (up(i)*hin(i,j,k) + um(i)*hin(i+1,j,k))
      enddo
      if (apply_OBC_u) then ; do i=is-2,ie+1
        if (OBC%OBC_mask_u(i,j) .and. (OBC%OBC_kind_u(I,j) == OBC_SIMPLE)) &
          uh(i,j,k) = OBC%uh(i,j,k)
      enddo ; endif

      do i=is-1,ie+1
        hpi(i) = hin(i,j,k) - dt * G%IDXDYh(i,j) * (uh(i,j,k) - uh(i-1,j,k))
      enddo

      do i=is-1,ie
!    Calculate the corrected thickness fluxes using the Hsu/Arakawa Corrector.
        if (um(i) < 0.0) then
          g1 = (hpi(i) - hin(i+1,j,k))
          g0 = (hpi(i) + hin(i+2,j,k)) - (hin(i+1,j,k) + hpi(i+1))
          gam = 1.0 - ((g0*g0 + g1*g1) / &
                 (((g0*g0 + g1*g1) + hin(i,j,k)*hin(i+1,j,k)) + h_neglect2) )**2
          a = C1_6 - um(i) * Cdt_6 * G%IDXu(i,j)

          if (um(i+1) == 0.0) then
            vgeom_v = 0.0
          else
            if (((hin(i+2,j,k) - hpi(i+1)) * gam) <= 0.5*hin(i+1,j,k)) then
              max_vg_v = 6.0 ! This is a very weak limit. 0.5 = 3.0/6.0.
            else ! Limit the velocity ratio to avoid negative effective thickness.
              max_vg_v = 3.0*hin(i+1,j,k) / ((hin(i+2,j,k) - hpi(i+1)) * gam)
            endif
            vgeom_v = max_vg_v
            if (um(i+1) > um(i)*max_vg_v**2) vgeom_v = sqrt(um(i+1)/um(i))
          endif
          Area_u(i,j,k) = G%dy_u(i,j) * (hin(i+1,j,k) + gam * ((0.5-a) * &
              (hpi(i)-hin(i+1,j,k)) + (a*vgeom_v) * (hpi(i+1)-hin(i+2,j,k))))
          uh(i,j,k) = um(i) * Area_u(i,j,k)
        else
          g1 = (hin(i,j,k) - hpi(i+1))
          g0 = ((hin(i-1,j,k) + hpi(i+1)) - (hin(i,j,k) + hpi(i)) )
          gam = 1.0 - ((g0*g0 + g1*g1) / &
                 (((g0*g0 + g1*g1) + hin(i,j,k)*hin(i+1,j,k)) + h_neglect2) )**2
          a = C1_6 + up(i) * Cdt_6 * G%IDXu(i,j)
          if (up(i-1) == 0.0) then ; vgeom_v = 0.0
          else
            if (((hin(i-1,j,k) - hpi(i)) * gam) <= 0.5*hin(i,j,k)) then
              max_vg_v = 6.0 ! This is a very weak limit. 0.5 = 3.0/6.0.
            else ! Limit the velocity ratio to avoid negative effective thickness.
              max_vg_v = 3.0*hin(i,j,k) / ((hin(i-1,j,k) - hpi(i)) * gam)
            endif
            vgeom_v = max_vg_v
            if (up(i-1) < up(i)*max_vg_v**2) vgeom_v = sqrt(up(i-1)/up(i))
          endif
          Area_u(i,j,k) = G%dy_u(i,j) * (hin(i,j,k) + gam * ((0.5-a) * &
              (hpi(i+1)-hin(i,j,k)) + (a * vgeom_v) * (hpi(i) - hin(i-1,j,k))))
          uh(i,j,k) = up(i) * Area_u(i,j,k)
        endif
      enddo
      if (apply_OBC_u) then ; do i=is-1,ie
        if (OBC%OBC_mask_u(i,j) .and. (OBC%OBC_kind_u(I,j) == OBC_SIMPLE)) &
          uh(i,j,k) = OBC%uh(i,j,k)
      enddo ; endif

    enddo ! j loop

!   Calculate the corrected thickness change.
    do j=js-2,je+2 ; do i=is,ie
      h(i,j,k) = hin(i,j,k) - dt* G%IDXDYh(i,j) * (uh(i,j,k) - uh(i-1,j,k))
    enddo ; enddo

!   Now, advect meridionally.
    do i=is,ie

!   Use upwind differencing as a predictor.
      do j=js-2,je+1
        vp(j) = (v(i,j,k) + ABS(v(i,j,k))) * 0.5
        vm(j) = v(i,j,k) - vp(j)
        vh(i,j,k) = G%dx_v(i,j) * (vp(j)*h(i,j,k) + vm(j)*h(i,j+1,k))
      enddo
      if (apply_OBC_v) then ; do j=js-2,je+1
        if (OBC%OBC_mask_v(i,j) .and. (OBC%OBC_kind_v(i,J) == OBC_SIMPLE)) &
          vh(i,j,k) = OBC%vh(i,j,k)
      enddo ; endif

      do j=js-1,je+1
        hpj(j) = h(i,j,k) - dt * G%IDXDYh(i,j) * (vh(i,j,k) - vh(i,j-1,k))
      enddo

!    Calculate fluxes and use the Hsu/Arakawa Corrector.

      do j=js-1,je
        if (vm(j) < 0.0) then
          g1 = (hpj(j) - h(i,j+1,k))
          g0 = (hpj(j) + h(i,j+2,k)) - (h(i,j+1,k) + hpj(j+1))
          gam = 1.0 - ((g0*g0 + g1*g1) / &
                 (((g0*g0 + g1*g1) + h(i,j,k)*h(i,j+1,k)) + h_neglect2) )**2
          a  = c1_6 - vm(j) * Cdt_6 * G%IDYv(i,j)
          if (vm(j+1) == 0.0) then
            vgeom_v = 0.0
          else
            if (((hin(i,j+2,k) - hpj(j+1)) * gam) <= 0.5*hin(i,j+1,k)) then
              max_vg_v = 6.0 ! This is a very weak limit. 0.5 = 3.0/6.0.
            else ! Limit the velocity ratio to avoid negative effective thickness.
              max_vg_v = 3.0*hin(i,j+1,k) / ((hin(i,j+2,k) - hpj(j+1)) * gam)
            endif
            vgeom_v = max_vg_v
            if (vm(j+1) > vm(j)*max_vg_v**2) vgeom_v = sqrt(vm(j+1)/vm(j))
          endif
          Area_v(i,j,k) = G%dx_v(i,j) * (h(i,j+1,k) + gam * ((0.5 - a) * &
              (hpj(j) - h(i,j+1,k)) + (a * vgeom_v) * (hpj(j+1) - h(i,j+2,k))))
          vh(i,j,k) = vm(j) * Area_v(i,j,k)
        else
          g1 = (h(i,j,k) - hpj(j+1))
          g0 = (h(i,j-1,k) + hpj(j+1)) - (h(i,j,k) + hpj(j))
          gam = 1.0 - ((g0*g0 + g1*g1) / &
                 (((g0*g0 + g1*g1) + h(i,j,k)*h(i,j+1,k)) + h_neglect2) )**2
          a  = c1_6 + vp(j) * Cdt_6 * G%IDYv(i,j)
          if (vp(j-1) == 0.0) then
            vgeom_v = 0.0
          else
            if (((hin(i,j-1,k) - hpj(j)) * gam) <= 0.5*hin(i,j,k)) then
              max_vg_v = 6.0 ! This is a very weak limit. 0.5 = 3.0/6.0.
            else ! Limit the velocity ratio to avoid negative effective thickness for a < 1/4.
              max_vg_v = 3.0*hin(i,j,k) / ((hin(i,j-1,k) - hpj(j)) * gam)
            endif
            vgeom_v = max_vg_v
            if (vp(j-1) < vp(j)*max_vg_v**2) vgeom_v = sqrt(vp(j-1)/vp(j))
          endif
          Area_v(i,j,k) = G%dx_v(i,j) * (h(i,j,k) + gam * ((0.5 - a) * &
              (hpj(j+1) - h(i,j,k)) + (a * vgeom_v) * (hpj(j) - h(i,j-1,k))))
          vh(i,j,k) = vp(j) * Area_v(i,j,k)
        endif
      enddo
      if (apply_OBC_v) then ; do j=js-1,je
        if (OBC%OBC_mask_v(i,j) .and. (OBC%OBC_kind_v(i,J) == OBC_SIMPLE)) &
          vh(i,j,k) = OBC%vh(i,j,k)
      enddo ; endif
    enddo ! i loop

!    Calculate the corrected thickness change.
    do j=js,je ; do i=is,ie
      h(i,j,k) = h(i,j,k) - dt * G%IDXDYh(i,j) * (vh(i,j,k) - vh(i,j-1,k))
      if (h(i,j,k) < G%Angstrom) h(i,j,k) = G%Angstrom
    enddo ; enddo
  enddo ! k loop

  if (ASSOCIATED(CS%diag%uh_lay)) then
    do k=1,nz ; do j=js,je ; do i=is-1,ie
      CS%diag%uh_lay(i,j,k) = uh(i,j,k)
    enddo ; enddo ; enddo
  endif
  if (ASSOCIATED(CS%diag%vh_lay)) then
    do k=1,nz ; do j=js-1,je ; do i=is,ie
      CS%diag%vh_lay(i,j,k) = vh(i,j,k)
    enddo ; enddo ; enddo
  endif
  call cpu_clock_end(id_clock_update)

  if (present(uhbt) .and. present(vhbt)) then
    ! Correct the fluxes so that they agree with the
    ! barotropic prescription, as suggested by P. Schopf.
    Idt = 1.0/dt

    call cpu_clock_begin(id_clock_pass)
    call pass_var(h, G%Domain)
    call cpu_clock_end(id_clock_pass)

    call cpu_clock_begin(id_clock_correct)
    do it=1,2
      do j=js,je
        do i=is-1,ie
          uhbt_res(i) = uhbt(i,j) - uh(i,j,1)
          Area_sum(i) = Area_u(i,j,1)
        enddo
        do k=2,nz ; do i=is-1,ie
          uhbt_res(i) = uhbt_res(i) - uh(i,j,k)
          Area_sum(i) = Area_sum(i) + Area_u(i,j,k)
        enddo ; enddo
        do i=is-1,ie
          u_cor(i) = 0.0
          if (Area_sum(i) > 0.0) u_cor(i) = uhbt_res(i) / Area_sum(i)
          if (abs(u_cor(i)) >  0.1*G%DXu(i,j)*Idt) then
            if (u_cor(i) > 0.0) then ; u_cor(i) = 0.1*G%DXu(i,j)*Idt
            else ; u_cor(i) = -0.1*G%DXu(i,j)*Idt ; endif
          endif
        enddo

        ! Again, following Paul's advice, if the corrective flow is in the
        ! opposite direction of the original flow, use the smaller of the
        ! original thickness and the upwind thickness.
        do k=1,nz ; do i=is-1,ie
          if ((it == 1) .and. (u_cor(i) * u(i,j,k) <= 0.0)) then
            if (u_cor(i) < 0.0) then
              if (Area_u(i,j,k) > G%dy_u(i,j)*h(i+1,j,k)) &
                Area_u(i,j,k) = G%dy_u(i,j)*h(i+1,j,k)
            else
              if (Area_u(i,j,k) > G%dy_u(i,j)*h(i,j,k)) &
                Area_u(i,j,k) = G%dy_u(i,j)*h(i,j,k)
            endif
          endif
          uh(i,j,k) = uh(i,j,k) + u_cor(i) * Area_u(i,j,k)
        enddo ; enddo
      enddo

      do j=js-1,je
        do i=is,ie
          vhbt_res(i) = vhbt(i,j) - vh(i,j,1)
          Area_sum(i) = Area_v(i,j,1)
        enddo
        do k=2,nz ; do i=is,ie
          vhbt_res(i) = vhbt_res(i) - vh(i,j,k)
          Area_sum(i) = Area_sum(i) + Area_v(i,j,k)
        enddo ; enddo
        do i=is,ie
          v_cor(i) = 0.0
          if (Area_sum(i) > 0.0) v_cor(i) = vhbt_res(i) / Area_sum(i)
          if (abs(v_cor(i)) >  0.1*G%DYv(i,j)*Idt) then
            if (v_cor(i) > 0.0) then ; v_cor(i) = 0.1*G%DYv(i,j)*Idt
            else ; v_cor(i) = -0.1*G%DYv(i,j)*Idt ; endif
          endif
        enddo
        do k=1,nz ; do i=is,ie
          if ((it == 1) .and. (v_cor(i) * v(i,j,k) <= 0.0)) then
            if (v_cor(i) < 0.0) then
              if (Area_v(i,j,k) > G%dx_v(i,j)*h(i,j+1,k)) &
                Area_v(i,j,k) = G%dx_v(i,j)*h(i,j+1,k)
            else
              if (Area_v(i,j,k) > G%dx_v(i,j)*h(i,j,k)) &
                Area_v(i,j,k) = G%dx_v(i,j)*h(i,j,k)
            endif
          endif
          vh(i,j,k) = vh(i,j,k) + v_cor(i) * Area_v(i,j,k)
        enddo ; enddo
      enddo

  !   Calculate the corrected thickness change.
      do k=1,nz ; do j=js,je ; do i=is,ie
        h(i,j,k) = hin(i,j,k) - dt*G%IDXDYh(i,j) * &
            ((uh(i,j,k) - uh(i-1,j,k)) + (vh(i,j,k) - vh(i,j-1,k)))
        if (h(i,j,k) < G%Angstrom) h(i,j,k) = G%Angstrom
      enddo ; enddo ; enddo

    enddo ! it Iterations...
    call cpu_clock_end(id_clock_correct)


    do it=3,4
      ! These are the desperation iterations - move as much water as
      ! possible using upwind fluxes.  Hopefully this is used rarely!

      call cpu_clock_begin(id_clock_pass)
      call pass_var(h, G%Domain)
      call pass_vector(uh, vh, G%Domain)
      call cpu_clock_end(id_clock_pass)

      call cpu_clock_begin(id_clock_correct)
      do j=js-1,je+1
        do i=is-1,ie+1 ; h_sum(i,j) = h(i,j,1) ; enddo
        do k=2,nz ; do i=is-1,ie+1
          h_sum(i,j) = h_sum(i,j) + h(i,j,k)
        enddo ; enddo
      enddo
      do j=js-1,je+1
        do_uh(j) = .false.
        do i=is-2,ie+1
          uhbt_resid(i,j) = uhbt(i,j) - uh(i,j,1)
        enddo
        do k=2,nz ; do i=is-2,ie+1
          uhbt_resid(i,j) = uhbt_resid(i,j) - uh(i,j,k)
        enddo ; enddo
        do i=is-1,ie
          if (1.0e-10*(G%DXDYh(i,j)*h_sum(i,j) + G%DXDYh(i+1,j)*h_sum(i+1,j)) < &
              abs(dt*uhbt_resid(i,j))) then
            do_uh(j) = .true.
          else
            uhbt_resid(i,j) = 0.0
          endif
        enddo
      enddo
      do j=js-2,je+1
        do_vh(j) = .false.
        do i=is-1,ie+1
          vhbt_resid(i,j) = vhbt(i,j) - vh(i,j,1)
        enddo
        do k=2,nz ; do i=is-1,ie+1
          vhbt_resid(i,j) = vhbt_resid(i,j) - vh(i,j,k)
        enddo ; enddo
        if ((j>=js-1) .and. (j<=je)) then ; do i=is-1,ie
          if (1.0e-10*(G%DXDYh(i,j)*h_sum(i,j) + G%DXDYh(i,j+1)*h_sum(i,j+1)) < &
              abs(dt*vhbt_resid(i,j))) then
            do_vh(j) = .true.
          else
            vhbt_resid(i,j) = 0.0
          endif
        enddo ; endif
      enddo

      do j=js,je ; if (do_uh(j)) then
        do i=is-1,ie
          if (uhbt_resid(i,j) > 0.0) then
            vol_out = dt*uhbt_resid(i,j)
            if (uhbt_resid(i-1,j) < 0.0) vol_out = vol_out - dt*uhbt_resid(i-1,j)
            if (vhbt_resid(i,j-1) < 0.0) vol_out = vol_out - dt*vhbt_resid(i,j-1)
            if (vhbt_resid(i,j) > 0.0) vol_out = vol_out + dt*vhbt_resid(i,j)
            if (vol_out <= 0.99*G%DXDYh(i,j)*h_sum(i,j)) then
              udy_cor(i) = uhbt_resid(i,j) / h_sum(i,j)
            else
              udy_cor(i) = uhbt_resid(i,j)*0.99*G%DXDYh(i,j) / vol_out
            endif
          elseif (uhbt_resid(i,j) < 0.0) then
            vol_out = -dt*uhbt_resid(i,j)
            if (uhbt_resid(i+1,j) > 0.0) vol_out = vol_out + dt*uhbt_resid(i+1,j)
            if (vhbt_resid(i+1,j-1) < 0.0) vol_out = vol_out - dt*vhbt_resid(i+1,j-1)
            if (vhbt_resid(i+1,j) > 0.0) vol_out = vol_out + dt*vhbt_resid(i+1,j)
            if (vol_out <= 0.99*G%DXDYh(i,j)*h_sum(i+1,j)) then
              udy_cor(i) = uhbt_resid(i,j) / h_sum(i+1,j)
            else
              udy_cor(i) = uhbt_resid(i,j)*0.99*G%DXDYh(i+1,j) / vol_out
            endif
          else
            udy_cor(i) = 0.0
          endif
        enddo
        do k=1,nz ; do i=is-1,ie
          if (udy_cor(i) < 0.0) then
            uh(i,j,k) = uh(i,j,k) + udy_cor(i) * h(i+1,j,k)
          else if (udy_cor(i) > 0.0) then
            uh(i,j,k) = uh(i,j,k) + udy_cor(i) * h(i,j,k)
          endif
        enddo ; enddo
      endif ; enddo
      do j=js-1,je ; if (do_vh(j)) then
        do i=is,ie
          if (vhbt_resid(i,j) > 0.0) then
            vol_out = dt*vhbt_resid(i,j)
            if (uhbt_resid(i-1,j) < 0.0) vol_out = vol_out - dt*uhbt_resid(i-1,j)
            if (uhbt_resid(i,j) > 0.0) vol_out = vol_out + dt*uhbt_resid(i,j)
            if (vhbt_resid(i,j-1) < 0.0) vol_out = vol_out - dt*vhbt_resid(i,j-1)
            if (vol_out <= 0.99*G%DXDYh(i,j)*h_sum(i,j)) then
              vdx_cor(i) = vhbt_resid(i,j) / h_sum(i,j)
            else
              vdx_cor(i) = vhbt_resid(i,j)*0.99*G%DXDYh(i,j) / vol_out
            endif
          elseif (vhbt_resid(i,j) < 0.0) then
            vol_out = -dt*vhbt_resid(i,j)
            if (uhbt_resid(i-1,j+1) < 0.0) vol_out = vol_out - dt*uhbt_resid(i-1,j+1)
            if (uhbt_resid(i,j+1) > 0.0) vol_out = vol_out + dt*uhbt_resid(i,j+1)
            if (vhbt_resid(i,j+1) > 0.0) vol_out = vol_out + dt*vhbt_resid(i,j+1)
            if (vol_out <= 0.99*G%DXDYh(i,j)*h_sum(i,j+1)) then
              vdx_cor(i) = vhbt_resid(i,j) / h_sum(i,j+1)
            else
              vdx_cor(i) = vhbt_resid(i,j)*0.99*G%DXDYh(i,j+1) / vol_out
            endif
          else
            vdx_cor(i) = 0.0
          endif
        enddo
        do k=1,nz ; do i=is,ie
          if (vdx_cor(i) < 0.0) then
            vh(i,j,k) = vh(i,j,k) + vdx_cor(i) * h(i,j+1,k)
          else if (vdx_cor(i) > 0.0) then
            vh(i,j,k) = vh(i,j,k) + vdx_cor(i) * h(i,j,k)
          endif
        enddo ; enddo
      endif ; enddo

      if (apply_OBC_v) then
        do k=1,nz ; do j=js-1,je ; do i=is,ie
          if (OBC%OBC_mask_v(i,j) .and. (OBC%OBC_kind_v(i,J) == OBC_SIMPLE)) &
            vh(i,j,k) = OBC%vh(i,j,k)
        enddo ; enddo ; enddo
      endif
      if (apply_OBC_u) then
        do k=1,nz ; do j=js,je ; do i=is-1,ie
          if (OBC%OBC_mask_u(i,j) .and. (OBC%OBC_kind_u(I,j) == OBC_SIMPLE)) &
            uh(i,j,k) = OBC%uh(i,j,k)
        enddo ; enddo ; enddo
      endif

      do k=1,nz ; do j=js,je ; do i=is,ie
        h(i,j,k) = hin(i,j,k) - dt*G%IDXDYh(i,j) * &
            ((uh(i,j,k) - uh(i-1,j,k)) + (vh(i,j,k) - vh(i,j-1,k)))
        if (h(i,j,k) < G%Angstrom) h(i,j,k) = G%Angstrom
      enddo ; enddo ; enddo
      call cpu_clock_end(id_clock_correct)
    enddo ! iterations 3 and 4.

  endif ! Correction of the fluxes.

end subroutine continuity_orig


subroutine continuity_orig_init(Time, G, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ptrs), target, intent(inout) :: diag
  type(continuity_orig_CS), pointer      :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
  character(len=128) :: version = '$Id: GOLD_continuity_Hallberg95.F90,v 1.1.2.2.2.7 2010/05/24 19:11:44 rwh Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "GOLD_continuity_Hallberg95" ! This module's name.

  if (associated(CS)) then
    call GOLD_error(WARNING, "continuity_orig_init called with associated control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag

  id_clock_update = cpu_clock_id('(Ocean continuity update)', grain=CLOCK_ROUTINE)
  id_clock_correct = cpu_clock_id('(Ocean continuity correction)', grain=CLOCK_ROUTINE)
  id_clock_pass = cpu_clock_id('(Ocean continuity halo updates)', grain=CLOCK_ROUTINE)

  call log_version(param_file, mod, version, tagname)

end subroutine continuity_orig_init

subroutine continuity_orig_end(CS)
  type(continuity_orig_CS), pointer :: CS
  deallocate(CS)
end subroutine continuity_orig_end

end module GOLD_continuity_Hallberg95
