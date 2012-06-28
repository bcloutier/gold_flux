module GOLD_PressureForce_Mont
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
!*  By Robert Hallberg, April 1994 - June 2000                         *
!*                                                                     *
!*    This file contains the subroutine that calculates the hori-      *
!*  zontal accelerations due to pressure gradients.  A second-order    *
!*  accurate, centered scheme is used.  If a split time stepping       *
!*  scheme is used, the vertical decomposition into barotropic and     *
!*  baroclinic contributions described by Hallberg (J Comp Phys 1997)  *
!*  is used.  With a nonlinear equation of state, compressibility is   *
!*  added along the lines proposed by Sun et al. (JPO 1999), but with  *
!*  compressibility coefficients based on a fit to a user-provided     *
!*  reference profile.                                                 *
!*                                                                     *
!*    PressureForce takes 8 arguments, which are described below. If   *
!*  a non-split time stepping scheme is used, the last three arguments *
!*  are ignored.                                                       *
!*                                                                     *
!*  Macros written all in capital letters are defined in GOLD_memory.h.*
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, f                                     *
!*    j+1  > o > o >   At ^:  v, PFv                                   *
!*    j    x ^ x ^ x   At >:  u, PFu                                   *
!*    j    > o > o >   At o:  h, D, M, e, p, pbce, gtot, T, S, rho_st  *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1                                                  *
!*           i  i+1                                                    *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_CompressComp, only : uncompress_e_rho, Grad_z_estar
use GOLD_CompressComp, only : uncompress_p_alpha, Grad_p_pstar, Compress_CS
use GOLD_diag_mediator, only : post_data, register_diag_field
use GOLD_diag_mediator, only : safe_alloc_ptr, diag_ptrs, time_type
use GOLD_error_handler, only : GOLD_error, GOLD_mesg, FATAL, WARNING, is_root_pe
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_grid, only : ocean_grid_type
use GOLD_tidal_forcing, only : calc_tidal_forcing, tidal_forcing_CS
use GOLD_variables, only : thermo_var_ptrs
use GOLD_EOS, only : calculate_density, calculate_density_derivs
use GOLD_EOS, only : int_specific_vol_dp
implicit none ; private

#include <GOLD_memory.h>

public PressureForce_Mont_Bouss, PressureForce_Mont_nonBouss, Set_pbce_Bouss
public Set_pbce_nonBouss, PressureForce_Mont_init, PressureForce_Mont_end

type, public :: PressureForce_Mont_CS ; private
  logical :: bulkmixedlayer ! If true, a refined bulk mixed layer is used with
                            ! tv%nk_Rml variable density mixed & buffer layers.
  logical :: tides          ! If true, apply tidal momentum forcing.
  real    :: Rho0           !   The density used in the Boussinesq
                            ! approximation, in kg m-3.
  real    :: Rho_atm        !   The assumed atmospheric density, in kg m-3.
                            ! By default, Rho_atm is 0.
  real    :: GFS_scale
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
                             ! ocean diagnostic fields.
  integer :: id_PFu_bc = -1, id_PFv_bc = -1, id_eta = -1, id_bott_press = -1
  integer :: id_e_tidal = -1
  type(Compress_CS), pointer :: compress_CSp => NULL()
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL()
end type PressureForce_Mont_CS

contains

subroutine PressureForce_Mont_nonBouss(h, tv, PFu, PFv, G, CS, p_atm, pbce, eta)
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in)   :: h
  type(thermo_var_ptrs), intent(inout)             :: tv
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(out) :: PFu
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(out) :: PFv
  type(ocean_grid_type),               intent(in)  :: G
  type(PressureForce_Mont_CS),         pointer     :: CS
  real, dimension(:,:),               optional, pointer     :: p_atm
  real, dimension(NXMEM_,NYMEM_,NZ_), optional, intent(out) :: pbce
  real, dimension(NXMEM_,NYMEM_),     optional, intent(out) :: eta

!    This subroutine determines the acceleration due to pressure forces in a
!  non-Boussinesq fluid using the compressibility compensated (if appropriate)
!  Montgomery-potential form described in Hallberg (Ocean Mod., 2005).
!    To work, the following fields must be set outside of the usual
!  ie to ie, je to je range before this subroutine is called:
!  h[ie+1] and h[je+1] and (if BULKMIXEDLAYER is set) Rml[ie+1] and
!  Rml[je+1] and (if tv%form_of_EOS is set) T[ie+1], S[ie+1],
!  T[je+1], and S[je+1].
! Arguments: h - Layer thickness, in m.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (out)     PFu - Zonal acceleration due to pressure gradients
!                  (equal to -dM/dx) in m s-2.
!  (out)     PFv - Meridional acceleration due to pressure
!                  gradients (equal to -dM/dy) in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 PressureForce_init.
!  (in)      p_atm - The pressure at the ice-ocean or atmosphere-ocean
!                    interface in Pa.
!  (out)     pbce - The baroclinic pressure anomaly in each layer
!                   due to free surface height anomalies, in m s-2.
!                   pbce points to a space with nz layers or NULL.
!  (out)     eta - The free surface height used to calculate PFu and PFv, in m,
!                  with any tidal contributions or compressibility compensation.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    M, &          ! The Montgomery potential, M = (p/rho + gz) , in m2 s-2.
    alpha_star, & ! Compression adjusted specific volume, in m3 kg-1.
    dalpha_star   ! The change in alpha_star across a layer, in m3 kg-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: p ! Interface pressure in Pa.
                ! p may be adjusted (with a nonlinear equation of state) so that
                ! its derivative compensates for the adiabatic compressibility
                ! in seawater, but p will still be close to the pressure.
  real, dimension(SZIQ_(G),SZJ_(G),SZK_(G)+1) :: gx_p
  real, dimension(SZI_(G),SZJQ_(G),SZK_(G)+1) :: gy_p
       ! gx_p and gy_p are the zonal and meridional gradients of the
       ! compressibility compensated pressure along isobars.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), target :: &
    T_tmp, &    ! Temporary array of temperatures where layers that are lighter
                ! than the mixed layer have the mixed layer's properties, in C.
    S_tmp       ! Temporary array of salinities where layers that are lighter
                ! than the mixed layer have the mixed layer's properties, in psu.

  real, dimension(SZI_(G),SZJ_(G)) :: &
    dM, &         !   A barotropic correction to the Montgomery potentials to
                  ! enable the use of a reduced gravity form of the equations,
                  ! in m2 s-2.
    dp_star, &    ! Layer thickness after compensation for compressibility, in Pa.
    e_tidal, &    !   Bottom geopotential anomaly due to tidal forces from
                  ! astronomical sources and self-attraction and loading, in m.
    geopot_bot, & !   Bottom geopotential relative to time-mean sea level,
                  ! including any tidal contributions, in units of m2 s-2.
    dz_geo, &     !   The change in geopotential across a layer, in m2 s-2.
    Rho_cv_BL, &  !   The coordinate potential density in the deepest variable
                  ! density near-surface layer, in kg m-3.
    SSH           !   Sea surface height anomalies, in m.
  real :: p_ref(SZI_(G))     !   The pressure used to calculate the coordinate
                             ! density, in Pa (usually 2e7 Pa = 2000 dbar).
  real :: PFu_bc, PFv_bc     ! The pressure gradient force due to along-layer
                             ! compensated density gradients, in m s-2.
  real :: PFu_3, PFv_3       ! The pressure gradient forces due to along-isobar
                             ! compensated pressure gradients, in m s-2.
  real :: dp_neglect         ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected, in Pa.
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: use_EOS         ! If true, density is calculated from T & S using
                             ! an equation of state.
  logical :: is_split        ! A flag indicating whether the pressure
                             ! gradient terms are to be split into
                             ! barotropic and baroclinic pieces.
  type(thermo_var_ptrs) :: tv_tmp! A structure of temporary T & S.

  real :: I_gEarth
  real :: dalpha
  real :: Pa_to_H   ! A factor to convert from Pa to the thicknesss units (H).
  real :: alpha_Lay(SZK_(G)) ! The specific volume of each layer, in kg m-3.
  real :: dalpha_int(SZK_(G)+1) ! The change in specific volume across each
                             ! interface, in kg m-3.
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  integer :: i, j, k
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke ; nkmb=tv%nk_Rml
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq

  use_p_atm = .false.
  if (present(p_atm)) then ; if (associated(p_atm)) use_p_atm = .true. ; endif
  is_split = .false. ; if (present(pbce)) is_split = .true.
  use_EOS = associated(tv%eqn_of_state)

  if (.not.associated(CS)) call GOLD_error(FATAL, &
      "GOLD_PressureForce: Module must be initialized before it is used.")

  I_gEarth = 1.0 / G%g_Earth
  dp_neglect = G%H_to_Pa * G%H_subroundoff
  do k=1,nz ; alpha_Lay(k) = 1.0 / G%Rlay(k) ; enddo
  do k=2,nz ; dalpha_int(k) = alpha_Lay(k-1) - alpha_Lay(k) ; enddo

  if (use_p_atm) then
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1 ; p(i,j,1) = p_atm(i,j) ; enddo ; enddo
  else
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1 ; p(i,j,1) = 0.0 ; enddo ; enddo
  endif
  do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    p(i,j,k+1) = p(i,j,k) + G%H_to_Pa * h(i,j,k)
  enddo ; enddo ; enddo

  if (present(eta)) then
    Pa_to_H = 1.0 / G%H_to_Pa
    if (use_p_atm) then ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta(i,j) = (p(i,j,nz+1) - p_atm(i,j))*Pa_to_H ! eta has the same units as h.
    enddo ; enddo ; else ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta(i,j) = p(i,j,nz+1)*Pa_to_H ! eta has the same units as h.
    enddo ; enddo ; endif
  endif

  if (CS%tides) then
    !   Determine the sea surface height anomalies, to enable the calculation
    ! of self-attraction and loading.
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      SSH(i,j) = -G%D(i,j)
    enddo ; enddo
    if (use_EOS) then
      do k=1,nz
        call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p(:,:,k), p(:,:,k+1), &
                                 0.0, G, tv%eqn_of_state, dz_geo, halo_size=1)
        do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
          SSH(i,j) = SSH(i,j) + I_gEarth * dz_geo(i,j)
        enddo ; enddo
      enddo
    elseif (CS%bulkmixedlayer) then
      do k=1,nkmb ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        SSH(i,j) = SSH(i,j) + G%H_to_kg_m2*h(i,j,k)/tv%Rml(i,j,k)
      enddo ; enddo ; enddo
      do k=nkmb+1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        SSH(i,j) = SSH(i,j) + G%H_to_kg_m2*h(i,j,k)*alpha_Lay(k)
      enddo ; enddo ; enddo
    else
      do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        SSH(i,j) = SSH(i,j) + G%H_to_kg_m2*h(i,j,k)*alpha_Lay(k)
      enddo ; enddo ; enddo
    endif

    call calc_tidal_forcing(CS%Time, SSH, e_tidal, G, CS%tides_CSp)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      geopot_bot(i,j) = -G%g_Earth*(e_tidal(i,j) + G%D(i,j))
    enddo ; enddo
  else
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      geopot_bot(i,j) = -G%g_Earth*G%D(i,j)
    enddo ; enddo
  endif

  if (use_EOS) then
    call Grad_p_pstar(p, gx_p, gy_p, nz, G, CS%compress_CSp)

    !   Calculate compression compensated in-situ specific volumes (alpha_star),
    ! and make compensating adjustments to p.

    !   With a bulk mixed layer, replace the T & S of any layers that are
    ! lighter than the the buffer layer with the properties of the buffer
    ! layer.  These layers will be massless anyway, and it avoids any
    ! formal calculations with hydrostatically unstable profiles.
    if (CS%bulkmixedlayer) then
      tv_tmp%T => T_tmp ; tv_tmp%S => S_tmp
      tv_tmp%eqn_of_state => tv%eqn_of_state
      do k=1,nkmb ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
      enddo ; enddo ; enddo
      do i=Isq,Ieq+1 ; p_ref(i) = tv%P_Ref ; enddo
      do j=Jsq,Jeq+1
        call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_ref, &
                        Rho_cv_BL(:,j), Isq, Ieq-Isq+2, tv%eqn_of_state)
      enddo
      do k=nkmb+1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        if (G%Rlay(k) < Rho_cv_BL(i,j)) then
          tv_tmp%T(i,j,k) = tv%T(i,j,nkmb) ; tv_tmp%S(i,j,k) = tv%S(i,j,nkmb)
        else
          tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
        endif
      enddo ; enddo ; enddo
    else
      tv_tmp%T => tv%T ; tv_tmp%S => tv%S
      tv_tmp%eqn_of_state => tv%eqn_of_state
    endif

    call uncompress_p_alpha(p,tv_tmp,alpha_star,G,CS%compress_CSp, dalpha_star)

    do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      ! Here dalpha = alpha_bot - alpha_top
      dalpha = dalpha_star(i,j,k) - alpha_star(i,j,k)
      alpha_star(i,j,k) = 0.5*(alpha_star(i,j,k) + dalpha_star(i,j,k))
      dalpha_star(i,j,k) = dalpha
    enddo ; enddo ; enddo
  else if (CS%bulkmixedlayer) then
    do k=1,nkmb ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      alpha_star(i,j,k) = 1.0 / tv%Rml(i,j,k)
    enddo ; enddo ; enddo
    do k=nkmb+1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      if (G%Rlay(k) < tv%Rml(i,j,nkmb)) then
        alpha_star(i,j,k) = alpha_star(i,j,nkmb)
      else
        alpha_star(i,j,k) = alpha_Lay(k)
      endif
    enddo ; enddo ; enddo
  endif                                               ! use_EOS

  if (use_EOS .or. CS%bulkmixedlayer) then
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      M(i,j,nz) = geopot_bot(i,j) + p(i,j,nz+1) * alpha_star(i,j,nz)
    enddo ; enddo
    do k=nz-1,1,-1 ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      M(i,j,k) = M(i,j,k+1) + p(i,j,k+1) * (alpha_star(i,j,k) - alpha_star(i,j,k+1))
    enddo ; enddo ; enddo   
  else ! not use_EOS .or. bulkmixedlayer
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      M(i,j,nz) = geopot_bot(i,j) + p(i,j,nz+1) * alpha_Lay(nz)
    enddo ; enddo
    do k=nz-1,1,-1 ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      M(i,j,k) = M(i,j,k+1) + p(i,j,k+1) * dalpha_int(k+1)
    enddo ; enddo ; enddo   
  endif ! use_EOS .or. bulkmixedlayer

  if (CS%GFS_scale < 1.0) then
    ! Adjust the Montgomery potential to make this a reduced gravity model.
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      dM(i,j) = (CS%GFS_scale - 1.0) * M(i,j,1)
    enddo ; enddo
    do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      M(i,j,k) = M(i,j,k) + dM(i,j)
    enddo ; enddo ; enddo

    !   Could instead do the following, to avoid taking small differences
    ! of large numbers...
!   do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
!     M(i,j,1) = CS%GFS_scale * M(i,j,1)
!   enddo ; enddo
!   if (use_EOS .or. CS%bulkmixedlayer) then
!     do k=2,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
!       M(i,j,k) = M(i,j,k-1) - p(i,j,k) * (alpha_star(i,j,k-1) - alpha_star(i,j,k))
!     enddo ; enddo ; enddo   
!   else ! not use_EOS .or. bulkmixedlayer
!     do k=2,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
!        M(i,j,k) = M(i,j,k-1) - p(i,j,k) * dalpha_int(k)
!     enddo ; enddo ; enddo   
!   endif ! use_EOS .or. bulkmixedlayer

  endif

  if (use_EOS) then
    do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      M(i,j,k) = M(i,j,k) + 0.125 * dalpha_star(i,j,k) * (p(i,j,k+1) - p(i,j,k))
    enddo ; enddo ; enddo
  endif

  ! Note that ddM/dPb = alpha_star(i,j,1)
  if (present(pbce)) then
    call Set_pbce_nonBouss(p, tv_tmp, G, G%g_Earth, CS%GFS_scale, pbce, &
                           alpha_star, dalpha_star)
  endif

!    Calculate the pressure force. On a Cartesian grid,
!      PFu = - dM/dx   and  PFv = - dM/dy.
  if (use_EOS .or. CS%bulkmixedlayer) then
    PFu_3 = 0.0 ; PFv_3 = 0.0
    do k=1,nz
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        dp_star(i,j) = (p(i,j,k+1) - p(i,j,k)) + dp_neglect
      enddo ; enddo
      do j=js,je ; do I=Isq,Ieq
        ! PFu_bc = p* grad alpha*
        PFu_bc = (alpha_star(i+1,j,k) - alpha_star(i,j,k)) * (G%IDXu(I,j) * &
          ((dp_star(i,j) * dp_star(i+1,j) + (p(i,j,k) * dp_star(i+1,j) + &
           p(i+1,j,k) * dp_star(i,j))) / (dp_star(i,j) + dp_star(i+1,j))))
        if (use_EOS) then
        ! PFu_3 = alpha* grad p*
          PFu_3 = 0.25 * ((alpha_star(i+1,j,k) + alpha_star(i,j,k)) * &
                  (gx_p(i,j,k) + gx_p(i,j,k+1)))
        endif
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * G%IDXu(I,j) + &
                     (PFu_bc + PFu_3)
        if (ASSOCIATED(CS%diag%PFu_bc)) CS%diag%PFu_bc(i,j,k) = PFu_bc
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv_bc = (alpha_star(i,j+1,k) - alpha_star(i,j,k)) * (G%IDYv(i,J) * &
          ((dp_star(i,j) * dp_star(i,j+1) + (p(i,j,k) * dp_star(i,j+1) + &
          p(i,j+1,k) * dp_star(i,j))) / (dp_star(i,j) + dp_star(i,j+1))))
        if (use_EOS) then
          PFv_3 = 0.25 * ((alpha_star(i,j+1,k) + alpha_star(i,j,k)) * &
                  (gy_p(i,j,k) + gy_p(i,j,k+1)))
        endif
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * G%IDYv(i,J) + &
                     (PFv_bc + PFv_3)
        if (ASSOCIATED(CS%diag%PFv_bc)) CS%diag%PFv_bc(i,j,k) = PFv_bc
      enddo ; enddo
    enddo ! k-loop
  else ! .not. use_EOS
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * G%IDXu(I,j)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * G%IDYv(i,J)
      enddo ; enddo
    enddo
  endif ! use_EOS

  if (CS%id_eta>0) call post_data(CS%id_eta, CS%diag%eta, CS%diag)
  if (CS%id_bott_press>0) call post_data(CS%id_bott_press, &
                                         p(:,:,nz+1), CS%diag)
  if (CS%id_e_tidal>0) call post_data(CS%id_e_tidal, e_tidal, CS%diag)   

end subroutine PressureForce_Mont_nonBouss

subroutine PressureForce_Mont_Bouss(h, tv, PFu, PFv, G, CS, p_atm, pbce, eta)
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in)   :: h
  type(thermo_var_ptrs), intent(inout)             :: tv
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(out) :: PFu
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(out) :: PFv
  type(ocean_grid_type),               intent(in)  :: G
  type(PressureForce_Mont_CS),         pointer     :: CS
  real, dimension(:,:),               optional, pointer     :: p_atm
  real, dimension(NXMEM_,NYMEM_,NZ_), optional, intent(out) :: pbce
  real, dimension(NXMEM_,NYMEM_),     optional, intent(out) :: eta

!    This subroutine determines the acceleration due to pressure
!  forces.
!    To work, the following fields must be set outside of the usual
!  ie to ie, je to je range before this subroutine is called:
!   h[ie+1] and h[je+1] and (if BULKMIXEDLAYER is set) Rml[ie+1] and
!   Rml[je+1] and (if tv%form_of_EOS is set) T[ie+1], S[ie+1],
!   T[je+1], and S[je+1].
! Arguments: h - Layer thickness, in m.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (out)     PFu - Zonal acceleration due to pressure gradients
!                  (equal to -dM/dx) in m s-2.
!  (out)     PFv - Meridional acceleration due to pressure
!                  gradients (equal to -dM/dy) in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 PressureForce_init.
!  (in)      p_atm - the pressure at the ice-ocean or atmosphere-ocean
!                    interface in Pa.
!  (out)     pbce - the baroclinic pressure anomaly in each layer
!                   due to free surface height anomalies, in m s-2.
!  (out)     eta - the free surface height used to calculate PFu and PFv, in m,
!                  with any tidal contributions or compressibility compensation.

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    M, &        ! The Montgomery potential, M = (p/rho + gz) , in m2 s-2.
    rho_star, & ! In-situ density divided by the derivative with depth of the
                ! corrected e times (G_Earth/Rho0).  In units of m s-2.
    drho_star   ! Difference over a layer of the in-situ density divided by the
                ! derivative with depth of the corrected e times
                ! (G_Earth/Rho0), in m s-2.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: e ! Interface height in m.
                ! e may be adjusted (with a nonlinearequation of state) so that
                ! its derivative compensates for the adiabatic compressibility
                ! in seawater, but e will still be close to the interface depth.
  real, dimension(SZIQ_(G),SZJ_(G),SZK_(G)+1) :: gx_e
  real, dimension(SZI_(G),SZJQ_(G),SZK_(G)+1) :: gy_e
       ! gx_e and gy_e are the zonal and meridional gradients of the
       ! compressibility compensated interface height relative to geopotentials.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), target :: &
    T_tmp, &    ! Temporary array of temperatures where layers that are lighter
                ! than the mixed layer have the mixed layer's properties, in C.
    S_tmp       ! Temporary array of salinities where layers that are lighter
                ! than the mixed layer have the mixed layer's properties, in psu.

  real :: Rho_cv_BL(SZI_(G),SZJ_(G))  !   The coordinate potential density in
                ! the deepest variable density near-surface layer, in kg m-3.
  real :: h_star(SZI_(G),SZJ_(G)) ! Layer thickness after compensation
                             ! for compressibility, in m.
  real :: e_tidal(SZI_(G),SZJ_(G)) ! Bottom geopotential anomaly due to tidal
                             ! forces from astronomical sources and self-
                             ! attraction and loading, in m.
  real :: p_ref(SZI_(G))     !   The pressure used to calculate the coordinate
                             ! density, in Pa (usually 2e7 Pa = 2000 dbar).
  real :: I_Rho0             ! 1/Rho0.
  real :: G_Rho0             ! G_Earth / Rho0 in m4 s-2 kg-1.
  real :: PFu_bc, PFv_bc     ! The pressure gradient force due to along-layer
                             ! compensated density gradients, in m s-2.
  real :: PFu_3, PFv_3       ! The third pressure gradient term due to horizontal
                             ! gradients of compensated geopotential, introduced by
                             ! Hallberg (Ocean Mod., 2005).
  real :: dr                 ! Temporary variables.
  real :: h_neglect          ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected, in m.
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: use_EOS         ! If true, density is calculated from T & S using
                             ! an equation of state.
  logical :: is_split        ! A flag indicating whether the pressure
                             ! gradient terms are to be split into
                             ! barotropic and baroclinic pieces.
  type(thermo_var_ptrs) :: tv_tmp! A structure of temporary T & S.
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  integer :: i, j, k

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke ; nkmb=tv%nk_Rml
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq

  use_p_atm = .false.
  if (present(p_atm)) then ; if (associated(p_atm)) use_p_atm = .true. ; endif
  is_split = .false. ; if (present(pbce)) is_split = .true.
  use_EOS = associated(tv%eqn_of_state)

  if (.not.associated(CS)) call GOLD_error(FATAL, &
       "GOLD_PressureForce: Module must be initialized before it is used.")

  h_neglect = G%H_subroundoff * G%H_to_m
  I_Rho0 = 1.0/CS%Rho0
  G_Rho0 = G%g_Earth/G%Rho0

  if (CS%tides) then
    !   Determine the surface height anomaly for calculating self attraction
    ! and loading.  This should really be based on bottom pressure anomalies,
    ! but that is not yet implemented, and the current form is correct for
    ! barotropic tides.
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e(i,j,1) = -1.0*G%D(i,j)
    enddo ; enddo
    do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e(i,j,1) = e(i,j,1) + h(i,j,k)
    enddo ; enddo ; enddo
    call calc_tidal_forcing(CS%Time, e(:,:,1), e_tidal, G, CS%tides_CSp)
  endif

!    Here layer interface heights, e, are calculated.
  if (CS%tides) then
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e(i,j,nz+1) = -1.0*G%D(i,j) - e_tidal(i,j)
    enddo ; enddo
  else
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e(i,j,nz+1) = -1.0*G%D(i,j)
    enddo ; enddo
  endif
  do k=nz,1,-1 ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    e(i,j,k) = e(i,j,k+1) + h(i,j,k)
  enddo ; enddo ; enddo
  if (use_EOS) then
    call Grad_z_estar(e, gx_e, gy_e, nz, G, CS%compress_CSp)

!   Calculate compression compensated in-situ densities (rho_star),
! and make compensating adjustments to e.

! With a bulk mixed layer, replace the T & S of any layers that are
! lighter than the the buffer layer with the properties of the buffer
! layer.  These layers will be massless anyway, and it avoids any
! formal calculations with hydrostatically unstable profiles.

    if (CS%bulkmixedlayer) then
      tv_tmp%T => T_tmp ; tv_tmp%S => S_tmp
      tv_tmp%eqn_of_state => tv%eqn_of_state
      do k=1,nkmb ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
      enddo ; enddo ; enddo
      do i=Isq,Ieq+1 ; p_ref(i) = tv%P_Ref ; enddo
      do j=Jsq,Jeq+1
        call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_ref, &
                        Rho_cv_BL(:,j), Isq, Ieq-Isq+2, tv%eqn_of_state)
      enddo
      do k=nkmb+1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        if (G%Rlay(k) < Rho_cv_BL(i,j)) then
          tv_tmp%T(i,j,k) = tv%T(i,j,nkmb) ; tv_tmp%S(i,j,k) = tv%S(i,j,nkmb)
        else
          tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
        endif
      enddo ; enddo ; enddo
    else
      tv_tmp%T => tv%T ; tv_tmp%S => tv%S
      tv_tmp%eqn_of_state => tv%eqn_of_state
    endif

    call uncompress_e_rho(e,tv_tmp,rho_star,G,CS%compress_CSp, drho_star)

    do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      dr = drho_star(i,j,k) - rho_star(i,j,k)
      rho_star(i,j,k) = 0.5*(rho_star(i,j,k) + drho_star(i,j,k))
      drho_star(i,j,k) = dr
    enddo ; enddo ; enddo
  else if (CS%bulkmixedlayer) then
    do k=1,nkmb ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      rho_star(i,j,k) = G_Rho0*tv%Rml(i,j,k)
    enddo ; enddo ; enddo
    do k=nkmb+1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      if (G%Rlay(k) < tv%Rml(i,j,nkmb)) then
        rho_star(i,j,k) = rho_star(i,j,nkmb)
      else
        rho_star(i,j,k) = G_Rho0*G%Rlay(k)
      endif
    enddo ; enddo ; enddo
  endif                                               ! use_EOS

!    Here the layer Montgomery potentials, M, are calculated.
  if (use_EOS .or. CS%bulkmixedlayer) then
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      M(i,j,1) = CS%GFS_scale * (rho_star(i,j,1) * e(i,j,1))
      if (use_p_atm) M(i,j,1) = M(i,j,1) + p_atm(i,j) * I_Rho0
    enddo ; enddo
    do k=2,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      M(i,j,k) = M(i,j,k-1) + (rho_star(i,j,k) - rho_star(i,j,k-1)) * e(i,j,k)
    enddo ; enddo ; enddo
    ! Note that this cannot be combined with the previous loop.
    if (use_EOS) then
      do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        M(i,j,k) = M(i,j,k) + 0.125 * drho_star(i,j,k) * (e(i,j,k+1) - e(i,j,k))
      enddo ; enddo ; enddo
    endif
  else ! not use_EOS .or. bulkmixedlayer
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      M(i,j,1) = G%g_prime(1) * e(i,j,1)
      if (use_p_atm) M(i,j,1) = M(i,j,1) + p_atm(i,j) * I_Rho0
    enddo ; enddo
    do k=2,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      M(i,j,k) = M(i,j,k-1) + G%g_prime(k) * e(i,j,k)
    enddo ; enddo ; enddo
  endif ! use_EOS .or. bulkmixedlayer

  if (present(pbce)) then
    call Set_pbce_Bouss(e, tv_tmp, G, G%g_Earth, CS%Rho0, CS%GFS_scale, pbce, &
                        rho_star, drho_star)
  endif

!    Calculate the pressure force. On a Cartesian grid,
!      PFu = - dM/dx   and  PFv = - dM/dy.
  if (use_EOS .or. CS%bulkmixedlayer) then
    PFu_3 = 0.0 ; PFv_3 = 0.0
    do k=1,nz
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        h_star(i,j) = (e(i,j,k) - e(i,j,k+1)) + h_neglect
      enddo ; enddo
      do j=js,je ; do I=Isq,Ieq
        PFu_bc = -1.0*(rho_star(i+1,j,k) - rho_star(i,j,k)) * (G%IDXu(I,j) * &
          ((h_star(i,j) * h_star(i+1,j) - (e(i,j,k) * h_star(i+1,j) + &
          e(i+1,j,k) * h_star(i,j))) / (h_star(i,j) + h_star(i+1,j))))
        if (use_EOS) then
          PFu_3 = -0.25 * ((rho_star(i+1,j,k) + rho_star(i,j,k)) * &
                   (gx_e(i,j,k) + gx_e(i,j,k+1)))
        endif
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * G%IDXu(I,j) + &
                     (PFu_bc - PFu_3)
        if (ASSOCIATED(CS%diag%PFu_bc)) CS%diag%PFu_bc(i,j,k) = PFu_bc
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv_bc = -1.0*(rho_star(i,j+1,k) - rho_star(i,j,k)) * (G%IDYv(i,J) * &
          ((h_star(i,j) * h_star(i,j+1) - (e(i,j,k) * h_star(i,j+1) + &
          e(i,j+1,k) * h_star(i,j))) / (h_star(i,j) + h_star(i,j+1))))
        if (use_EOS) then
          PFv_3 = -0.25 * ((rho_star(i,j+1,k) + rho_star(i,j,k)) * &
                   (gy_e(i,j,k) + gy_e(i,j,k+1)))
        endif
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * G%IDYv(i,J) + &
                     (PFv_bc - PFv_3)
        if (ASSOCIATED(CS%diag%PFv_bc)) CS%diag%PFv_bc(i,j,k) = PFv_bc
      enddo ; enddo
    enddo ! k-loop
  else ! .not. use_EOS
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * G%IDXu(I,j)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * G%IDYv(i,J)
      enddo ; enddo
    enddo
  endif ! use_EOS

  if (present(eta)) then
    if (CS%tides) then
    ! eta is the sea surface height relative to a time-invariant geoid, for
    ! comparison with what is used for eta in btstep.  See how e was calculated
    ! about 200 lines above.
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = e(i,j,1) + e_tidal(i,j)
      enddo ; enddo
    else
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = e(i,j,1)
      enddo ; enddo
    endif
  endif

! Here the pressure gradient accelerations are offered for averaging.
  if (CS%id_PFu_bc>0) call post_data(CS%id_PFu_bc, CS%diag%PFu_bc, CS%diag)
  if (CS%id_PFv_bc>0) call post_data(CS%id_PFv_bc, CS%diag%PFv_bc, CS%diag)
  if (CS%id_e_tidal>0) call post_data(CS%id_e_tidal, e_tidal, CS%diag)

end subroutine PressureForce_Mont_Bouss

subroutine Set_pbce_Bouss(e, tv, G, g_Earth, Rho0, GFS_scale, pbce, rho_star, drho_star)
  real, dimension(NXMEM_,NYMEM_,NZp1_), intent(in)  :: e
  type(thermo_var_ptrs),                intent(in)  :: tv
  type(ocean_grid_type),                intent(in)  :: G
  real,                                 intent(in)  :: g_Earth
  real,                                 intent(in)  :: Rho0
  real,                                 intent(in)  :: GFS_scale
!  type(PressureForce_Mont_CS),          pointer     :: CS
  real, dimension(NXMEM_,NYMEM_,NZ_),   intent(out) :: pbce
  real, dimension(NXMEM_,NYMEM_,NZ_), optional, intent(in) :: rho_star
  real, dimension(NXMEM_,NYMEM_,NZ_), optional, intent(in) :: drho_star
!    This subroutine determines the partial derivative of the acceleration due 
!  to pressure forces with the free surface height.
! Arguments: e - Interface height, in m.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (in)      G - The ocean's grid structure.
!  (in)      g_Earth - The gravitational acceleration, in m s-2.
!  (in)      Rho0 - The "Boussinesq" ocean density, in kg m-3.
!  (in)      Rho_atmos - The atmospheric density, in kg m-3.  Typically this
!                        should be 0, but it could be set to a significant
!                        fraction of the ocean density to use the model in a
!                        reduced gravity mode.
!  (in)      CS - The control structure returned by a previous call to
!                 PressureForce_init.
!  (out)     pbce - the baroclinic pressure anomaly in each layer
!                   due to free surface height anomalies, in m s-2.
!  (in,opt)  rho_star - The layer densities (maybe compressibility compensated),
!                       times g/rho_0, in m s-2.
!  (in,opt)  drho_star - The difference in rho_star beteen the top and bottom
!                        of a layer from uncompensated compressibility, in m s-2.
   
  real :: Ihtot(SZI_(G))     ! The inverse of the sum of the layer
                             ! thicknesses, in m-1.
  real :: press(SZI_(G))     ! Interface pressure, in Pa.
  real :: T_int(SZI_(G))     ! Interface temperature in C.
  real :: S_int(SZI_(G))     ! Interface salinity in PSU.
  real :: dR_dT(SZI_(G))     ! Partial derivatives of density with temperature
  real :: dR_dS(SZI_(G))     ! and salinity in kg m-3 K-1 and kg m-3 PSU-1.
  real :: rho_in_situ(SZI_(G)) !In-situ density at the top of a layer.
  real :: G_Rho0             ! g_Earth / Rho0 in m4 s-2 kg-1.
  real :: Rho0xG             ! g_Earth * Rho0 in kg s-2 m-2.
  logical :: use_EOS         ! If true, density is calculated from T & S using
                             ! an equation of state.
  logical :: Bulkmixedlayer  ! True if the mixed layer density is a spatially
                             ! varying array.
  real :: h_neglect          ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected, in m.
  integer :: Isq, Ieq, Jsq, Jeq, nz, i, j, k
    
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq ; nz = G%ke
 
  Rho0xG = Rho0*g_Earth
  G_Rho0 = g_Earth/Rho0
  use_EOS = associated(tv%eqn_of_state)
  bulkmixedlayer = associated(tv%Rml)
  h_neglect = G%H_subroundoff * G%H_to_m

  if (use_EOS .or. bulkmixedlayer) then
    if (present(rho_star)) then
      do j=Jsq,Jeq+1
        do i=Isq,Ieq+1
          Ihtot(i) = 1.0 / ((e(i,j,1)-e(i,j,nz+1)) + h_neglect)
          pbce(i,j,1) = GFS_scale * rho_star(i,j,1)
        enddo
        do k=2,nz ; do i=Isq,Ieq+1
          pbce(i,j,k) = pbce(i,j,k-1) + (rho_star(i,j,k)-rho_star(i,j,k-1)) * &
                        ((e(i,j,k) - e(i,j,nz+1)) * Ihtot(i))
        enddo ; enddo
        if (use_EOS .and. present(drho_star)) then
          do k=1,nz ; do i=Isq,Ieq+1
            pbce(i,j,k) = pbce(i,j,k) + 0.125 * drho_star(i,j,k) * &
                                        (e(i,j,k+1) - e(i,j,k))*Ihtot(i)
          enddo ; enddo
        endif
      enddo ! end of j loop
    else
      do j=Jsq,Jeq+1
        do i=Isq,Ieq+1
          Ihtot(i) = 1.0 / ((e(i,j,1)-e(i,j,nz+1)) + h_neglect)
          press(i) = -Rho0xG*e(i,j,1)
        enddo
        call calculate_density(tv%T(:,j,1), tv%S(:,j,1), press, rho_in_situ, &
                               Isq, Ieq-Isq+2, tv%eqn_of_state)
        do i=Isq,Ieq+1
          pbce(i,j,1) = G_Rho0*(GFS_scale * rho_in_situ(i))
        enddo
        do k=2,nz
          do i=Isq,Ieq+1
            press(i) = -Rho0xG*e(i,j,k)
            T_int(i) = 0.5*(tv%T(i,j,k-1)+tv%T(i,j,k))
            S_int(i) = 0.5*(tv%S(i,j,k-1)+tv%S(i,j,k))
          enddo                        
          call calculate_density_derivs(T_int, S_int, press, dR_dT, dR_dS, &
                                        Isq, Ieq-Isq+2, tv%eqn_of_state)
          do i=Isq,Ieq+1
            pbce(i,j,k) = pbce(i,j,k-1) + G_Rho0 * &
               ((e(i,j,k) - e(i,j,nz+1)) * Ihtot(i)) * &
               (dR_dT(i)*(tv%T(i,j,k)-tv%T(i,j,k-1)) + &
                dR_dS(i)*(tv%S(i,j,k)-tv%S(i,j,k-1)))
          enddo                        
        enddo
      enddo ! end of j loop
    endif
  else ! not use_EOS .or. bulkmixedlayer
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1
        Ihtot(i) = 1.0 / ((e(i,j,1)-e(i,j,nz+1)) + h_neglect)
        pbce(i,j,1) = G%g_prime(1)
      enddo
      do k=2,nz ; do i=Isq,Ieq+1
        pbce(i,j,k) = pbce(i,j,k-1) + &
                      G%g_prime(k) * ((e(i,j,k) - e(i,j,nz+1)) * Ihtot(i))
     enddo ; enddo
    enddo ! end of j loop
  endif ! use_EOS .or. bulkmixedlayer
                      
end subroutine Set_pbce_Bouss

subroutine Set_pbce_nonBouss(p, tv, G, g_Earth, GFS_scale, pbce, alpha_star, dalpha_star)
  real, dimension(NXMEM_,NYMEM_,NZp1_), intent(in)  :: p
  type(thermo_var_ptrs),                intent(in)  :: tv
  type(ocean_grid_type),                intent(in)  :: G
  real,                                 intent(in)  :: g_Earth
  real,                                 intent(in)  :: GFS_scale
!  type(PressureForce_Mont_CS),          pointer     :: CS
  real, dimension(NXMEM_,NYMEM_,NZ_),   intent(out) :: pbce
  real, dimension(NXMEM_,NYMEM_,NZ_), optional, intent(in) :: alpha_star
  real, dimension(NXMEM_,NYMEM_,NZ_), optional, intent(in) :: dalpha_star
!    This subroutine determines the partial derivative of the acceleration due 
!  to pressure forces with the column mass.
! Arguments: p - Interface pressures, in Pa.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (in)      G - The ocean's grid structure.
!  (in)      g_Earth - The gravitational acceleration, in m s-2.
!  (in)      Rho_atmos - The atmospheric density, in kg m-3.  Typically this
!                        should be 0, but it could be set to a significant
!                        fraction of the ocean density to use the model in a
!                        reduced gravity mode.
!  (in)      CS - The control structure returned by a previous call to
!                 PressureForce_init.
!  (out)     pbce - the baroclinic pressure anomaly in each layer
!                   due to total thickness anomalies, in m4 s-2 kg-1.
!  (in,opt)  alpha_star - The layer specific volumes (maybe compressibility
!                         compensated), in m3 kg-1.
!  (in,opt)  dalpha_star - The difference in alpha_star between the top and
!                          bottom of a layer from uncompensated compressibility, in m3 kg-1.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    dpbce, &      !   A barotropic correction to the pbce to enable the use of
                  ! a reduced gravity form of the equations, in m4 s-2 kg-1.
    C_htot        ! dP_dH divided by the total ocean pressure, m2 kg-1.
  real :: T_int(SZI_(G))     ! Interface temperature in C.
  real :: S_int(SZI_(G))     ! Interface salinity in PSU.
  real :: dR_dT(SZI_(G))     ! Partial derivatives of density with temperature
  real :: dR_dS(SZI_(G))     ! and salinity in kg m-3 K-1 and kg m-3 PSU-1.
  real :: rho_in_situ(SZI_(G)) !In-situ density at an interface, in kg m-3.
  real :: alpha_Lay(SZK_(G)) ! The specific volume of each layer, in kg m-3.
  real :: dalpha_int(SZK_(G)+1) ! The change in specific volume across each
                             ! interface, in kg m-3.
  real :: dP_dH   !   A factor that converts from thickness to pressure,
                  ! usually in Pa m2 kg-1.
  real :: dp_neglect         ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected, in Pa.
  logical :: use_EOS         ! If true, density is calculated from T & S using
                             ! an equation of state.
  logical :: Bulkmixedlayer  ! True if the mixed layer density is a spatially
                             ! varying array.
  integer :: Isq, Ieq, Jsq, Jeq, nz, i, j, k
    
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq ; nz = G%ke

  use_EOS = associated(tv%eqn_of_state)
  bulkmixedlayer = associated(tv%Rml)

  dP_dH = g_Earth * G%H_to_kg_m2
  dp_neglect = dP_dH * G%H_subroundoff

  do k=1,nz ; alpha_Lay(k) = 1.0 / G%Rlay(k) ; enddo
  do k=2,nz ; dalpha_int(k) = alpha_Lay(k-1) - alpha_Lay(k) ; enddo

  if (use_EOS .or. bulkmixedlayer) then
    if (present(alpha_star)) then
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        C_htot(i,j) = dP_dH / ((p(i,j,nz+1)-p(i,j,1)) + dp_neglect)
        pbce(i,j,nz) = dP_dH * alpha_star(i,j,nz)
      enddo ; enddo
      do k=nz-1,1,-1 ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        pbce(i,j,k) = pbce(i,j,k+1) + ((p(i,j,k+1)-p(i,j,1)) * C_htot(i,j)) * &
            (alpha_star(i,j,k) - alpha_star(i,j,k+1))
      enddo ; enddo ; enddo
    else
      do j=Jsq,Jeq+1
        call calculate_density(tv%T(:,j,nz), tv%S(:,j,nz), p(:,j,nz+1), &
                               rho_in_situ, Isq, Ieq-Isq+2, tv%eqn_of_state)
        do i=Isq,Ieq+1
          C_htot(i,j) = dP_dH / ((p(i,j,nz+1)-p(i,j,1)) + dp_neglect)
          pbce(i,j,nz) = dP_dH / rho_in_situ(i)
        enddo
      enddo
      do k=nz-1,1,-1 ; do j=Jsq,Jeq+1
        do i=Isq,Ieq+1
          T_int(i) = 0.5*(tv%T(i,j,k)+tv%T(i,j,k+1))
          S_int(i) = 0.5*(tv%S(i,j,k)+tv%S(i,j,k+1))
        enddo                        
        call calculate_density(T_int, S_int, p(:,j,k+1), rho_in_situ, &
                                    Isq, Ieq-Isq+2, tv%eqn_of_state)
        call calculate_density_derivs(T_int, S_int, p(:,j,k+1), dR_dT, dR_dS, &
                                    Isq, Ieq-Isq+2, tv%eqn_of_state)
        do i=Isq,Ieq+1
          pbce(i,j,k) = pbce(i,j,k+1) + ((p(i,j,k+1)-p(i,j,1))*C_htot(i,j)) * &
              ((dR_dT(i)*(tv%T(i,j,k+1)-tv%T(i,j,k)) + &
                dR_dS(i)*(tv%S(i,j,k+1)-tv%S(i,j,k))) / rho_in_situ(i)**2)
        enddo
      enddo ; enddo
    endif
  else ! not use_EOS .or. bulkmixedlayer
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      C_htot(i,j) = dP_dH / ((p(i,j,nz+1)-p(i,j,1)) + dp_neglect)
      pbce(i,j,nz) = dP_dH * alpha_Lay(nz)
    enddo ; enddo
    do k=nz-1,1,-1 ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      pbce(i,j,k) = pbce(i,j,k+1) + ((p(i,j,k+1)-p(i,j,1))*C_htot(i,j)) * &
          dalpha_int(k+1)
    enddo ; enddo ; enddo
  endif ! use_EOS .or. bulkmixedlayer

  if (GFS_scale < 1.0) then
    ! Adjust the Montgomery potential to make this a reduced gravity model.
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      dpbce(i,j) = (GFS_scale - 1.0) * pbce(i,j,1)
    enddo ; enddo
    do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      pbce(i,j,k) = pbce(i,j,k) + dpbce(i,j)
    enddo ; enddo ; enddo
  endif
  
  if (use_EOS .and. present(dalpha_star)) then
    do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      pbce(i,j,k) = pbce(i,j,k) + 0.125 * dalpha_star(i,j,k) * &
                           ((p(i,j,k+1)-p(i,j,k))*C_htot(i,j))
    enddo ; enddo ; enddo
  endif

end subroutine Set_pbce_nonBouss


subroutine PressureForce_Mont_init(Time, G, param_file, diag, CS, compress_CSp, &
                              tides_CSp)
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ptrs), target, intent(inout) :: diag
  type(PressureForce_Mont_CS),  pointer  :: CS
  type(Compress_CS),      optional, pointer :: compress_CSp
  type(tidal_forcing_CS), optional, pointer :: tides_CSp
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module.
!  (in)      compress_CSp - a pointer to the control structure of the
!                           compressibility compensation module.
!  (in)      tides_CSp - a pointer to the control structure of the tide module.
  logical :: use_temperature, use_EOS
  character(len=128) :: version = '$Id: GOLD_PressureForce_Montgomery.F90,v 1.1.2.1.2.15 2011/05/12 22:08:05 Robert.Hallberg Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod   ! This module's name.

  if (associated(CS)) then
    call GOLD_error(WARNING, "PressureForce_init called with an associated "// &
                            "control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag ; CS%Time => Time
  if (present(compress_CSp)) then
    if (associated(compress_CSp)) CS%compress_CSp => compress_CSp
  endif
  if (present(tides_CSp)) then
    if (associated(tides_CSp)) CS%tides_CSp => tides_CSp
  endif


  call read_param(param_file,"RHO_0",CS%Rho0,.true.)
  ! use_EOS is only an option if TEMPERATURE is defined.
  call read_param(param_file,"TEMPERATURE",use_temperature,.true.)
  use_EOS = use_temperature
  if (use_temperature) then
    ! This checks both the old and new names for the option of using an EOS.
    call read_param(param_file,"NONLINEAR_EOS",use_EOS)
    call read_param(param_file,"USE_EOS",use_EOS)
  endif

  CS%bulkmixedlayer = .false. ; CS%tides = .false.
  call read_param(param_file,"TIDES",CS%tides)
  call read_param(param_file,"BULKMIXEDLAYER",CS%bulkmixedlayer)

  if (CS%bulkmixedlayer .or. use_EOS) then
    CS%id_PFu_bc = register_diag_field('ocean_model', 'PFu_bc', G%axesul, Time, &
         'Density Gradient Zonal Pressure Force Accel.', "meter second-2")
    CS%id_PFv_bc = register_diag_field('ocean_model', 'PFv_bc', G%axesvl, Time, &
         'Density Gradient Meridional Pressure Force Accel.', "meter second-2")
    if (CS%id_PFu_bc > 0) then
      call safe_alloc_ptr(diag%PFu_bc,G%Isdq,G%Iedq,G%jsd,G%jed,G%ke)
      diag%PFu_bc(:,:,:) = 0.0
    endif
    if (CS%id_PFv_bc > 0) then
      call safe_alloc_ptr(diag%PFv_bc,G%isd,G%ied,G%Jsdq,G%Jedq,G%ke)
      diag%PFv_bc(:,:,:) = 0.0
    endif
  endif

  if (CS%tides) then
    CS%id_e_tidal = register_diag_field('ocean_model', 'e_tidal', G%axesh1, &
        Time, 'Tidal Forcing Astronomical and SAL Height Anomaly', 'meter')
  endif

  CS%GFS_scale = 1.0
  if (G%g_prime(1) /= G%g_Earth) CS%GFS_scale = G%g_prime(1) / G%g_Earth

  mod = "GOLD_PressureForce_Mont"
  call log_version(param_file, mod, version, tagname)
  call log_param(param_file, mod, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3")
  call log_param(param_file, mod, "TIDES", CS%tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)
  call log_param(param_file, mod, "BULKMIXEDLAYER", CS%bulkmixedlayer, &
                 "If true, use a Kraus-Turner-like bulk mixed layer \n"//&
                 "with transitional buffer layers.  Layers 1 through  \n"//&
                 "NKML+NKBL have variable densities. There must be at \n"//&
                 "least NKML+NKBL+1 layers if BULKMIXEDLAYER is true.", &
                 default=.false.)
  call log_param(param_file, mod, "GFS / G_EARTH", CS%GFS_scale)

end subroutine PressureForce_Mont_init


subroutine PressureForce_Mont_end(CS)
  type(PressureForce_Mont_CS), pointer :: CS
  if (associated(CS)) deallocate(CS)
end subroutine PressureForce_Mont_end

end module GOLD_PressureForce_Mont
