module GOLD_thickness_diffuse
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
!*  By Robert Hallberg, July 1999                                      *
!*                                                                     *
!*    The subroutine in this file implements horizontal interface      *
!*  depth diffusion.  In a Z-coordinate model, this would be described *
!*  as Gent-McWilliams diffusion, with an identical interpretation     *
!*  of the diffusion coefficient.  The coefficient is locally limited  *
!*  to guarantee numerical stability, and are vertically uniform.      *
!*                                                                     *
!*  Macros written all in capital letters are defined in GOLD_memory.h.*
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, vh, vav                               *
!*    j    x ^ x ^ x   At >:  u, uh, uav                               *
!*    j    > o > o >   At o:  h                                        *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_diag_mediator, only : post_data, query_averaging_enabled, diag_ptrs
use GOLD_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use GOLD_error_handler, only : GOLD_error, FATAL, WARNING
use GOLD_EOS, only : calculate_density, calculate_density_derivs
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_grid, only : ocean_grid_type
use GOLD_interface_heights, only : find_eta
use GOLD_lateral_mixing_coeffs, only : VarMix_CS
use GOLD_variables, only : thermo_var_ptrs
use GOLD_EOS, only : calculate_density
use GOLD_MEKE_types, only : MEKE_type

implicit none ; private

public thickness_diffuse, thickness_diffuse_init, thickness_diffuse_end
public vert_fill_TS

#include <GOLD_memory.h>

type, public :: thickness_diffuse_CS ; private
  real    :: Khth           ! The background interface depth diffusivity in m2 s-1.
  real    :: Khth_Slope_Cff ! Slope dependence coefficient of Khth in m2 s-1.
  real    :: Khth_Min       ! Minimum value of Khth in m2 s-1.
  real    :: Khth_Max       ! Maximum value of Khth in m2 s-1, or 0 for no max.
  real    :: slope_max      ! Slopes steeper than this are limited in some way.
  real    :: kappa_smooth   ! A diffusivity that is used to interpolate more
                            ! sensible values of T & S into thin layers.
  logical :: thickness_diffuse ! If true, interfaces heights are diffused
                             ! with a coefficient of Khth.
  logical :: diffuse_isopycnals ! If true, it is isopycnal surfaces, not model
                            ! interfaces, that are smoothed.  This currently
                            ! only changes things if bulkmixedlayer is defined.
  logical :: full_thickness_diffuse  ! If true, use the new full-column code
                            ! to calculate an overturning residual circulation.
  logical :: bulkmixedlayer ! If true, a refined bulk mixed layer is used.
  integer :: nkmb           ! The combined number of mixed- and buffer-layers.
  integer :: nkml           ! The number of layers within the mixed layer.
  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
                             ! ocean diagnostic fields.
  real, pointer :: GMwork(:,:) => NULL()  ! Work by thick. diff. in W m-2.
  integer :: id_uhGM = -1, id_vhGM = -1, id_GMwork = -1, id_KH_u = -1, id_KH_v = -1
 ! integer :: id_sfn_slope_x = -1, id_sfn_slope_y = -1, id_sfn_x = -1, id_sfn_y = -1
end type thickness_diffuse_CS

contains

subroutine thickness_diffuse(h, uhtr, vhtr, tv, dt, G, MEKE, VarMix, CS)
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(inout) :: h
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(inout) :: uhtr
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(inout) :: vhtr
  type(thermo_var_ptrs),               intent(in)    :: tv
  real,                                intent(in)    :: dt
  type(ocean_grid_type),               intent(in)    :: G
  type(MEKE_type),                     pointer       :: MEKE
  type(VarMix_CS),                     pointer       :: VarMix
  type(thickness_diffuse_CS),          pointer       :: CS
!    This subroutine does interface depth diffusion.  The fluxes are
!  limited to give positive definiteness, and the diffusivities are
!  limited to guarantee stability.

! Arguments: h - Layer thickness, in m.
!  (in/out)  uhtr - Accumulated zonal mass fluxes in m3.
!  (in/out)  vhtr - Accumulated meridional mass fluxes in m3.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      VarMix - A structure containing a number of fields related to
!                     variable lateral mixing.
!  (in)      MEKE - A structure containing information about the Mesoscale Eddy
!                   Kinetic Energy parameterization; this might be unassociated.
!  (in)      CS - The control structure returned by a previous call to
!                 thickness_diffuse_init.


  real :: uhD(SZIQ_(G), SZJ_(G), SZK_(G)) ! uhD & vhD are the diffusive u*h &
  real :: vhD(SZI_(G), SZJQ_(G), SZK_(G)) ! v*h fluxes, in m3 s-1.

  real :: KH_u(SZIQ_(G), SZJ_(G))   ! KH_u & KH_v are the horizontal
  real :: KH_v(SZI_(G), SZJQ_(G))   ! thickness diffusivities at u & v (m2/s)
                                    ! grid points, in m2 s-1.
  real :: Khth_Loc      ! Locally calculated thickness mixing coefficient (m2/s)
  logical :: use_VarMix, Resoln_scaled
  integer :: i, j, k, is, ie, js, je, nz


  if (.not. ASSOCIATED(CS)) call GOLD_error(FATAL, "GOLD_thickness_diffuse:"// &
         "Module must be initialized before it is used.")
  if ((.not.CS%thickness_diffuse) .or. &
       .not.( CS%Khth > 0.0 .or. associated(VarMix) .or. associated(MEKE) ) ) return

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (ASSOCIATED(MEKE)) then
    if (ASSOCIATED(MEKE%GM_src)) then
      do j=js,je ; do i=is,ie ; MEKE%GM_src(i,j) = 0. ; enddo ; enddo
    endif
  endif

  use_VarMix = .false. ; Resoln_scaled = .false.
  if (Associated(VarMix)) then
    use_VarMix = VarMix%use_variable_mixing
    Resoln_scaled = VarMix%Resoln_scaled_KhTh
  endif

  ! Set the diffusivities.
  do j=js,je ; do I=is-1,ie
    Khth_Loc = CS%Khth
    if (use_VarMix) &
      Khth_Loc = Khth_Loc + CS%KHTH_Slope_Cff*VarMix%L2u(I,j)*VarMix%SN_u(I,j)
    if (associated(MEKE%Kh)) &
      Khth_Loc = Khth_Loc + MEKE%KhTh_fac*sqrt(MEKE%Kh(i,j)*MEKE%Kh(i+1,j))
    if (Resoln_scaled) &
      Khth_Loc = Khth_Loc * 0.5*(VarMix%Res_fn_h(i,j) + VarMix%Res_fn_h(i+1,j))
    if (CS%Khth_Max > 0) then
      Khth_Loc = max(CS%Khth_min, min(Khth_Loc,CS%Khth_Max))
    else
      Khth_Loc = max(CS%Khth_min, Khth_Loc)
    endif
    KH_u(I,j) = 0.2/(dt*(G%IDXu(I,j)*G%IDXu(I,j) + G%IDYu(I,j)*G%IDYu(I,j)))
    if (KH_u(I,j) > Khth_Loc) KH_u(I,j) = Khth_Loc
  enddo ; enddo

  do J=js-1,je ; do i=is,ie
    Khth_Loc = CS%Khth
    if (use_VarMix) &
      Khth_Loc = Khth_Loc + CS%KHTH_Slope_Cff*VarMix%L2v(i,J)*VarMix%SN_v(i,J)
    if (associated(MEKE%Kh)) &
      Khth_Loc = Khth_Loc + MEKE%KhTh_fac*sqrt(MEKE%Kh(i,j)*MEKE%Kh(i,j+1))
    if (Resoln_scaled) &
      Khth_Loc = Khth_Loc * 0.5*(VarMix%Res_fn_h(i,j) + VarMix%Res_fn_h(i,j+1))
    if (CS%Khth_Max > 0) then
      Khth_Loc = max(CS%Khth_min, min(Khth_Loc,CS%Khth_Max))
    else
      Khth_Loc = max(CS%Khth_min, Khth_Loc)
    endif
    KH_v(i,J) = 0.2/(dt*(G%IDXv(i,J)*G%IDXv(i,J) + G%IDYv(i,J)*G%IDYv(i,J)))
    if (KH_v(i,J) > Khth_Loc) KH_v(i,J) = Khth_Loc
  enddo ; enddo


  if (CS%full_thickness_diffuse) then
    call thickness_diffuse_full(h, Kh_u, Kh_v, tv, uhD, vhD, dt, G, MEKE, CS)
  else
    call thickness_diffuse_orig(h, Kh_u, Kh_v, tv, uhD, vhD, dt, G, MEKE, CS)
  endif


  do k=1,nz
    do j=js,je ; do I=is-1,ie
      uhtr(I,j,k) = uhtr(I,j,k) + uhD(I,j,k)*dt
      if (ASSOCIATED(CS%diag%uhGM)) CS%diag%uhGM(I,j,k) = uhD(I,j,k)
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      vhtr(i,J,k) = vhtr(i,J,k) + vhD(i,J,k)*dt
      if (ASSOCIATED(CS%diag%vhGM)) CS%diag%vhGM(i,J,k) = vhD(i,J,k)
    enddo ; enddo
    do j=js,je ; do i=is,ie
      h(i,j,k) = h(i,j,k) - dt * G%IDXDYh(i,j) * &
          ((uhD(I,j,k) - uhD(I-1,j,k)) + (vhD(i,J,k) - vhD(i,J-1,k)))
      if (h(i,j,k) < G%Angstrom) h(i,j,k) = G%Angstrom
    enddo ; enddo
  enddo

! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_uhGM > 0) call post_data(CS%id_uhGM, CS%diag%uhGM, CS%diag)
    if (CS%id_vhGM > 0) call post_data(CS%id_vhGM, CS%diag%vhGM, CS%diag)
    if (CS%id_GMwork > 0) call post_data(CS%id_GMwork, CS%GMwork, CS%diag)
    if (CS%id_KH_u > 0) call post_data(CS%id_KH_u, KH_u, CS%diag)
    if (CS%id_KH_v > 0) call post_data(CS%id_KH_v, KH_v, CS%diag)
  endif

end subroutine thickness_diffuse

subroutine thickness_diffuse_orig(h, Kh_u, Kh_v, tv, uhD, vhD, dt, G, MEKE, CS)
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(in)  :: h
  real, dimension(NXMEMQ_,NYMEM_),     intent(in)  :: Kh_u
  real, dimension(NXMEM_,NYMEMQ_),     intent(in)  :: Kh_v
  type(thermo_var_ptrs),               intent(in)  :: tv
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(out) :: uhD
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(out) :: vhD
  real,                                intent(in)  :: dt
  type(ocean_grid_type),               intent(in)  :: G
  type(MEKE_type),                     pointer     :: MEKE
  type(thickness_diffuse_CS),          pointer     :: CS
!    This subroutine does interface depth diffusion.  The fluxes are
!  limited to give positive definiteness, and the diffusivities are
!  limited to guarantee stability.

! Arguments: h - Layer thickness, in m.
!  (in)      Kh_u - Thickness diffusivity at u points, in m2 s-1.
!  (in)      Kh_v - Thickness diffusivity at v points, in m2 s-1.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (out)     uhD - Zonal mass fluxes in m3 s-1.
!  (out)     vhD - Meridional mass fluxes in m3 s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      MEKE - A structure containing information about the Mesoscale Eddy
!                   Kinetic Energy parameterization; this might be unassociated.
!  (in)      CS - The control structure returned by a previous call to
!                 thickness_diffuse_init.
  real :: e(SZI_(G), SZJ_(G), SZK_(G)+1) ! The heights of the interfaces, relative
                                    ! to mean sea level, in m, pos. upward.
  real :: h_avail(SZI_(G), SZJ_(G), SZK_(G)) ! The mass available for diffusion
                                    ! out of each face, divided by dt, in m3 s-1.
  real :: uhtot(SZIQ_(G), SZJ_(G))  ! The vertical sum of uh, in m3 s-1.
  real :: vhtot(SZI_(G), SZJQ_(G))  ! The vertical sum of vh, in m3 s-1.
  real :: E_x(SZIQ_(G), SZJ_(G))  ! X-slope of interface at u points (for diagnostics)
  real :: E_y(SZI_(G), SZJQ_(G))  ! Y-slope of interface at u points (for diagnostics)
  real :: GMwork(SZI_(G), SZJ_(G))  ! Energetic work removed by thickness mixing (W m-2)
  real :: uhDtent, vhDtent          ! The tentative (unlimited) diffusive
                                    ! fluxes, in m3 s-1.
  real :: I4dt                      ! 1 / 4 dt
  ! The remaining arrays are only used with the "diffuse_isopycnals" option.
  real, dimension(SZI_(G), SZJ_(G)) :: &
    Rho_MLb, &          !   The potential density at the bottom of the mixed
                        ! layer in kg m-3.
    Rho_lay             !   The potential density of a layer, in kg m-3
  real :: p_ref(SZI_(G))!   An array full of the reference pressure, in Pa.

  real, dimension(SZIQ_(G), SZJ_(G)) :: h_mlu_avail_w, h_mlu_avail_e
  real, dimension(SZI_(G), SZJQ_(G)) :: h_mlv_avail_s, h_mlv_avail_n
                        !   The available volume within the mixed layer to
                        ! the west and east (_w and _e) of u points or to the
                        ! south and north (_s and _n) of v points, divided by
                        ! the timestep, in m3 s-1.
  real :: fr_u(SZIQ_(G), SZJ_(G)), fr_v(SZI_(G), SZJQ_(G))
                        !   The fraction of the available volume in the upwind
                        ! column of the mixed layer that is moved across a
                        ! cell face, nondimensional.
  integer :: kbu(SZIQ_(G), SZJ_(G)), kbv(SZI_(G), SZJQ_(G))
                        !   Layers deeper than these values have top interfaces
                        ! that approximate isopycnals.  In these layers and
                        ! above, the bolus velocity is constant with depth.
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  logical :: find_work  ! If true, find the change in energy due to the fluxes.
  logical :: domore_k
  integer :: is, ie, js, je, nz, nkmb
  integer :: i, j, k, kb_max

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  nkmb = tv%nk_Rml

  I4dt = 0.25 / dt
  use_EOS = associated(tv%eqn_of_state)

  !   Determine the thickness above which the interface slopes are not reliable
  ! estimates of isopycnal slopes.  The coordinate variable (often sigma2) is
  ! used, and while this may not be the most accurate, it is probably not as
  ! bad an approximation to the physics of the ocean as is diffusing interface
  ! height in the first place.
  do j=js,je ; do I=is-1,ie ; kbu(I,j) = -1 ; enddo ; enddo
  do J=js-1,je ; do i=is,ie ; kbv(i,J) = -1 ; enddo ; enddo

  if (CS%bulkmixedlayer .and. CS%diffuse_isopycnals) then
    do j=js,je ; do I=is-1,ie ; if (G%umask(I,j)<0.5) kbu(I,j) = 0 ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; if (G%vmask(i,J)<0.5) kbv(i,J) = 0 ; enddo ; enddo

    ! Find the maximum density in the mixed and buffer layers.
    k = nkmb
    if (use_EOS) then
      p_ref(:) = tv%P_Ref
      do j=js-1,je+1
        call calculate_density(tv%T(:,j,1),tv%S(:,j,1),p_ref,Rho_MLb(:,j), &
                               is-1,ie-is+3, tv%eqn_of_state)
      enddo
      do k=2,nkmb ; do j=js-1,je+1
        call calculate_density(tv%T(:,j,k),tv%S(:,j,k),p_ref,Rho_lay(:,j), &
                               is-1,ie-is+3, tv%eqn_of_state)
        do i=is-1,ie+1
          if ((Rho_MLb(i,j) < Rho_lay(i,j)) .and. (h(i,j,k) > 2.0*G%Angstrom)) &
              Rho_MLb(i,j) = Rho_lay(i,j)
        enddo
      enddo ; enddo
    else
      do j=js-1,je+1 ; do i=is-1,ie+1 ; Rho_MLb(i,j) = tv%Rml(i,j,1) ; enddo ; enddo
      do k=2,nkmb ; do j=js-1,je+1 ; do i=is-1,ie+1
        if ((Rho_MLb(i,j) < tv%Rml(i,j,k)) .and. (h(i,j,k) > 2.0*G%Angstrom)) &
            Rho_MLb(i,j) = tv%Rml(i,j,k)
      enddo ; enddo ; enddo
    endif

    !   Find the uppermost interior layer that is denser than the mixed layers.
    ! The interfaces below this layer are considered to approximate isopycnals
    ! well enough for this parameterization.
    do k=nkmb+1,nz-1
      domore_k = .false.
      if (use_EOS) then
        do j=js-1,je+1
          call calculate_density(tv%T(:,j,k),tv%S(:,j,k),p_ref,Rho_lay(:,j), &
                                 is-1, ie-is+3, tv%eqn_of_state)
        enddo
        do j=js,je ; do I=is-1,ie ; if (kbu(I,j) < 0) then
          if (min(Rho_lay(i,j),Rho_lay(i+1,j)) > max(Rho_MLb(i,j), Rho_MLb(i+1,j))) then
            kbu(I,j) = k
          else ; domore_k = .true. ; endif
        endif ; enddo ; enddo
        do J=js-1,je ; do i=is,ie ; if (kbv(i,J) < 0) then
          if (min(Rho_lay(i,j),Rho_lay(i,j+1)) > max(Rho_MLb(i,j), Rho_MLb(i,j+1))) then
            kbv(i,J) = k
          else ; domore_k = .true. ; endif
        endif ; enddo ; enddo
      else ! Linear EOS, so G%Rlay is appropriate for all layers.
        do j=js,je ; do I=is-1,ie ; if (kbu(I,j) < 0) then
          if (G%Rlay(k) > max(Rho_MLb(i,j), Rho_MLb(i+1,j))) then ; kbu(I,j) = k
          else ; domore_k = .true. ; endif
        endif ; enddo ; enddo
        do J=js-1,je ; do i=is,ie ; if (kbv(i,J) < 0) then
          if (G%Rlay(k) > max(Rho_MLb(i,j), Rho_MLb(i,j+1))) then ; kbv(i,J) = k
          else ; domore_k = .true. ; endif
        endif ; enddo ; enddo
      endif
      if (.not.domore_k) exit ! At all points credible interfaces have been found.
    enddo

    ! Find the shallowest layer that contains any interior points.
    kb_max = 0
    do j=js,je ; do I=is-1,ie
      if (kbu(I,j) < 0) kbu(I,j) = nz
      if (kbu(I,j) > kb_max) kb_max = kbu(I,j)
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      if (kbv(i,J) < 0) kbv(i,J) = nz
      if (kbv(i,J) > kb_max) kb_max = kbv(i,J)
    enddo ; enddo
  endif ! Bulk mixed layer and diffuse_isopycnals

  call find_eta(h, tv, G%g_Earth, G, e, halo_size=1)

  do j=js,je ; do I=is-1,ie
    uhtot(I,j) = 0.0
    h_mlu_avail_w(I,j) = 0.0 ; h_mlu_avail_e(I,j) = 0.0
  enddo ; enddo

  do J=js-1,je ; do i=is,ie
    vhtot(i,J) = 0.0
    h_mlv_avail_s(i,J) = 0.0 ; h_mlv_avail_n(i,J) = 0.0
  enddo ; enddo

  do k=nz,1,-1
    do j=js-1,je+1 ; do i=is-1,ie+1
      h_avail(i,j,k) =  max(I4dt*G%DXDYh(i,j)*(h(i,j,k)-G%Angstrom),0.0)
    enddo ; enddo
    do j=js,je ; do I=is-1,ie
      if (k>kbu(I,j)) then
        uhDtent = (KH_u(I,j)*(e(i,j,k)-e(i+1,j,k))*G%dy_u(I,j)*G%IDXu(I,j)) * &
                   G%m_to_H - uhtot(I,j)
        if (uhDtent >= 0.0) then
          uhD(I,j,k) = min(uhDtent,h_avail(i,j,k))
        else
          uhD(I,j,k) = max(uhDtent,-h_avail(i+1,j,k))
        endif

        uhtot(I,j) = uhtot(I,j)  + uhD(I,j,k)
      else
        h_mlu_avail_w(I,j) = h_mlu_avail_w(I,j) + h_avail(i,j,k)
        h_mlu_avail_e(I,j) = h_mlu_avail_e(I,j) + h_avail(i+1,j,k)
      endif
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      if (k>kbv(i,J)) then
        vhDtent = (KH_v(i,J)*(e(i,j,k)-e(i,j+1,k))*G%dx_v(i,J)*G%IDYv(i,J)) * &
                   G%m_to_H - vhtot(i,J)
        if (vhDtent >= 0.0) then
          vhD(i,J,k) = min(vhDtent,h_avail(i,j,k))
        else
          vhD(i,J,k) = max(vhDtent,-h_avail(i,j+1,k))
        endif

        vhtot(i,J) = vhtot(i,J)  + vhD(i,J,k)
      else
        h_mlv_avail_s(i,J) = h_mlv_avail_s(i,J) + h_avail(i,j,k)
        h_mlv_avail_n(i,J) = h_mlv_avail_n(i,J) + h_avail(i,j+1,k)
      endif
    enddo ; enddo
  enddo

  if (CS%bulkmixedlayer .and. CS%diffuse_isopycnals) then
    do j=js,je ; do I=is-1,ie
      !   The first term is proportional to the sea surface height gradient,
      ! and should be quite small.
      uhDtent = (KH_u(I,j)*(e(i,j,1)-e(i+1,j,1))*G%dy_u(I,j)*G%IDXu(I,j)) * &
                 G%m_to_H - uhtot(I,j)
      if (uhDtent >= 0) then
        if (uhDtent >= h_mlu_avail_w(I,j)) then ; fr_u(I,j) = 1.0
        else ; fr_u(I,j) = uhDtent / h_mlu_avail_w(I,j) ; endif
      else
        if (-uhDtent >= h_mlu_avail_e(I,j)) then ; fr_u(I,j) = -1.0
        else ; fr_u(I,j) = uhDtent / h_mlu_avail_e(I,j) ; endif
      endif
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      vhDtent = (KH_v(i,J)*(e(i,j,1)-e(i,j+1,1))*G%dx_v(i,J)*G%IDYv(i,J)) * &
                 G%m_to_H - vhtot(i,J)
      if (vhDtent >= 0) then
        if (vhDtent >= h_mlv_avail_s(i,J)) then ; fr_v(i,J) = 1.0
        else ; fr_v(i,J) = vhDtent / h_mlv_avail_s(i,J) ; endif
      else
        if (-vhDtent >= h_mlv_avail_n(i,J)) then ; fr_v(i,J) = -1.0
        else ; fr_v(i,J) = vhDtent / h_mlv_avail_n(i,J) ; endif
      endif
    enddo ; enddo
    !   It is not expected that fr_u or fr_v will be 1 or -1 when uhDtent is
    ! nonzero.  This should not happen with a vertically constant interface
    ! height diffusivity.
    do k=1,kb_max
      do j=js,je ; do I=is-1,ie ; if (k<=kbu(I,j)) then
          if (fr_u(I,j) >= 0.0) then
            uhD(I,j,k) = fr_u(I,j) * h_avail(i,j,k)
          else
            uhD(I,j,k) = fr_u(I,j) * h_avail(i+1,j,k)
          endif
        endif ; enddo ; enddo
      do J=js-1,je ; do i=is,ie ; if (k<=kbv(i,J)) then
          if (fr_v(i,J) >= 0.0) then
            vhD(i,J,k) = fr_v(i,J) * h_avail(i,j,k)
          else
            vhD(i,J,k) = fr_v(i,J) * h_avail(i,j+1,k)
          endif
        endif ; enddo ; enddo
    enddo
  endif


! Diagnostic of work done by thickness mixing flux
  find_work = .false.
  if (ASSOCIATED(MEKE)) find_work = ASSOCIATED(MEKE%GM_src)
  find_work = (ASSOCIATED(CS%GMwork) .or. find_work)

  if (find_work) then
    do j=js-1,je+1 ; do i=is-1,ie+1
      GMwork(i,j) = 0.0
    enddo ; enddo
    do j=js,je ; do I=is-1,ie
      uhtot(I,j) = 0.0
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      vhtot(i,J) = 0.0
    enddo ; enddo
    do k=nz,1,-1
      do j=js,je ; do I=is-1,ie
        E_x(I,j) = (e(i+1,j,k)-e(i,j,k))*G%IDXu(I,j)
        uhtot(I,j) = uhtot(I,j)  + uhD(I,j,k)*G%IDYu(I,j)
      enddo ; enddo
      do J=js-1,je ; do i=is,ie
        E_y(i,J) = (e(i,j+1,k)-e(i,j,k))*G%IDYv(i,J)
        vhtot(i,J) = vhtot(i,J)  + vhD(i,J,k)*G%IDXv(i,J)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        GMwork(i,j) = GMwork(i,j) + 0.5 * G%g_prime(k) * G%H_to_kg_m2 * &
               (( uhtot(i,j)*E_x(i,j) + uhtot(i-1,j  )*E_x(i-1, j ) ) + &
                ( vhtot(i,j)*E_y(i,j) + vhtot(i  ,j-1)*E_y(i  ,j-1) ))
      enddo ; enddo
    enddo
    if (ASSOCIATED(MEKE)) then ; if (ASSOCIATED(MEKE%GM_src)) then
      do j=js,je ; do i=is,ie ; MEKE%GM_src(i,j) = GMwork(i,j) ; enddo ; enddo
    endif ; endif
    if (ASSOCIATED(CS%GMwork)) then ; do j=js,je ; do i=is,ie
      CS%GMwork(i,j) = GMwork(i,j)
    enddo ; enddo ; endif
  endif

end subroutine thickness_diffuse_orig

subroutine thickness_diffuse_full(h, Kh_u, Kh_v, tv, uhD, vhD, dt, G, MEKE, CS)
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(in)  :: h
  real, dimension(NXMEMQ_,NYMEM_),     intent(in)  :: Kh_u
  real, dimension(NXMEM_,NYMEMQ_),     intent(in)  :: Kh_v
  type(thermo_var_ptrs),               intent(in)  :: tv
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(out) :: uhD
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(out) :: vhD
  real,                                intent(in)  :: dt
  type(ocean_grid_type),               intent(in)  :: G
  type(MEKE_type),                     pointer     :: MEKE
  type(thickness_diffuse_CS),          pointer     :: CS
!    This subroutine does interface depth diffusion.  The fluxes are
!  limited to give positive definiteness, and the diffusivities are
!  limited to guarantee stability.

! Arguments: h - Layer thickness, in m.
!  (in)      Kh_u - Thickness diffusivity at u points, in m2 s-1.
!  (in)      Kh_v - Thickness diffusivity at v points, in m2 s-1.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (out)     uhD - Zonal mass fluxes in m3 s-1.
!  (out)     vhD - Meridional mass fluxes in m3 s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      MEKE - A structure containing information about the Mesoscale Eddy
!                   Kinetic Energy parameterization; this might be unassociated.
!  (in)      CS - The control structure returned by a previous call to
!                 thickness_diffuse_init.


  real, dimension(SZI_(G), SZJ_(G), SZK_(G)) :: &
    T, &          ! The temperature (or density) in C, with the values in
                  ! in massless layers filled vertically by diffusion.
    S, &          ! The filled salinity, in PSU, with the values in
                  ! in massless layers filled vertically by diffusion.
    Rho, &        ! Density itself, when a nonlinear equation of state is
                  ! not in use.
    h_avail, &    ! The mass available for diffusion out of each face, divided
                  ! by dt, in m3 s-1.
    h_frac        ! The fraction of the mass in the column above the bottom
                  ! interface of a layer that is within a layer, ND. 0<h_frac<=1
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)+1) :: &
    e, &          ! The interfaces heights relative to mean sea level, in m.
    pres, &       ! The pressure at an interface, in Pa.
    h_avail_rsum  ! The running sum of h_avail above an interface, in m3 s-1.
  real, dimension(SZIQ_(G), SZJ_(G)) :: &
    drho_dT_u, &  ! The derivatives of density with temperature and
    drho_dS_u     ! salinity at u points, in kg m-3 K-1 and kg m-3 psu-1.
  real, dimension(SZI_(G), SZJQ_(G)) :: &
    drho_dT_v, &  ! The derivatives of density with temperature and
    drho_dS_v     ! salinity at v points, in kg m-3 K-1 and kg m-3 psu-1.
  real :: uhtot(SZIQ_(G), SZJ_(G))  ! The vertical sum of uhD, in m3 s-1.
  real :: vhtot(SZI_(G), SZJQ_(G))  ! The vertical sum of vhD, in m3 s-1.
  real, dimension(SZIQ_(G)) :: &
    T_u, S_u, &   ! Temperature, salinity, and pressure on the interface at
    pres_u        ! the u-point in the horizontal.
  real, dimension(SZI_(G)) :: &
    T_v, S_v, &   ! Temperature, salinity, and pressure on the interface at
    pres_v        ! the v-point in the horizontal.
  real :: Work_u(SZIQ_(G), SZJ_(G)) ! The work being done by the thickness
  real :: Work_v(SZI_(G), SZJQ_(G)) ! diffusion integrated over a cell, in W.
  real :: Work_h                    ! The work averaged over an h-cell in W m-2.
  real :: I4dt                      ! 1 / 4 dt
  real :: drdiA, drdiB  ! Along layer zonal- and meridional- potential density
  real :: drdjA, drdjB  ! gradients in the layers above (A) and below(B) the
                        ! interface times the grid spacing, in kg m-3.
  real :: drdkL, drdkR  ! Vertical density differences across an interface,
                        ! in kg m-3.
  real :: hg2A, hg2B, hg2L, hg2R
  real :: haA, haB, haL, haR
  real :: dzaL, dzaR
  real :: wtA, wtB, wtL, wtR
  real :: drdx, drdy, drdz  ! Zonal, meridional, and vertical density gradients,
                            ! in units of kg m-4.
  real :: Sfn_est       ! Two preliminary estimates (before limiting) of the
  real :: Sfn_unlim     ! overturning streamfunction, both in m3 s-1.
  real :: Sfn           ! The overturning streamfunction, in m3 s-1.
  real :: Sfn_safe      ! The streamfunction that goes linearly back to 0 at the
                        ! top.  This is a good thing to use when the slope is
                        ! so large as to be meaningless.
  real :: Slope         ! The slope of density surfaces, calculated in a way
                        ! that is always between -1 and 1.
  real :: mag_grad2     ! The squared magnitude of the 3-d density gradient, in kg2 m-8.
  real :: slope2_Ratio  ! The ratio of the slope squared to slope_max squared.
  real :: I_slope_max2  ! The inverse of slope_max squared, nondimensional.
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected, in H.
  real :: h_neglect2    ! h_neglect^2, in H2.
  real :: dz_neglect    ! A thickness in m that is so small it is usually lost
                        ! in roundoff and can be neglected, in m.
  real :: G_scale       ! The gravitational accerlation times the conversion
                        ! factor from thickness to m, in m s-2 or m4 s-2 kg-1.
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  logical :: find_work  ! If true, find the change in energy due to the fluxes.
  integer :: nk_linear  ! The number of layers over which the streamfunction
                        ! goes to 0.

! Diagnostics that should be eliminated altogether later...
 ! real, dimension(SZIQ_(G), SZJ_(G), SZK_(G)+1) :: sfn_x, sfn_slope_x
 ! real, dimension(SZI_(G), SZJQ_(G), SZK_(G)+1) :: sfn_y, sfn_slope_y

  integer :: is, ie, js, je, nz, nkmb, Isdq
  integer :: i, j, k
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke ; Isdq = G%Isdq
  nkmb = tv%nk_Rml

  I4dt = 0.25 / dt
  I_slope_max2 = 1.0 / (CS%slope_max**2)
  G_scale = G%g_Earth * G%H_to_m
  h_neglect = G%H_subroundoff ; h_neglect2 = h_neglect**2
  dz_neglect = G%H_subroundoff*G%H_to_m

  use_EOS = associated(tv%eqn_of_state)

  nk_linear = 1 ; if (CS%bulkmixedlayer) nk_linear = CS%nkml

  find_work = .false.
  if (ASSOCIATED(MEKE)) find_work = ASSOCIATED(MEKE%GM_src)
  find_work = (ASSOCIATED(CS%GMwork) .or. find_work)

  if (use_EOS) then
    call vert_fill_TS(h, tv%T, tv%S, CS%kappa_smooth, dt, T, S, G, 1)
  elseif (CS%bulkmixedlayer) then
    ! Use T as a temporary array and S as a scratch array.
    do k=1,nkmb ; do j=js-1,je+1 ; do i=is-1,ie+1
      T(i,j,k) = tv%Rml(i,j,k)
    enddo ; enddo ; enddo
    do k=nkmb+1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
      T(i,j,k) = G%Rlay(k)
    enddo ; enddo ; enddo
    call vert_fill_TS(h, T, T, CS%kappa_smooth, dt, Rho, S, G, 1)
  endif

  ! Find the maximum and minimum permitted streamfunction.
  do j=js-1,je+1 ; do i=is-1,ie+1
    h_avail_rsum(i,j,1) = 0.0
    pres(i,j,1) = 0.0  ! ### This should be atmospheric pressure.

    h_avail(i,j,1) = max(I4dt*G%DXDYh(i,j)*(h(i,j,1)-G%Angstrom),0.0)
    h_avail_rsum(i,j,2) = h_avail(i,j,1)
    h_frac(i,j,1) = 1.0
    pres(i,j,2) = pres(i,j,1) + G%H_to_Pa*h(i,j,1)
  enddo ; enddo
  do k=2,nz
    do j=js-1,je+1 ; do i=is-1,ie+1
      h_avail(i,j,k) = max(I4dt*G%DXDYh(i,j)*(h(i,j,k)-G%Angstrom),0.0)
      h_avail_rsum(i,j,k+1) = h_avail_rsum(i,j,k) + h_avail(i,j,k)
      h_frac(i,j,k) = 0.0 ; if (h_avail(i,j,k) > 0.0) &
        h_frac(i,j,k) = h_avail(i,j,k) / h_avail_rsum(i,j,k+1)
      pres(i,j,k+1) = pres(i,j,k) + G%H_to_Pa*h(i,j,k)
    enddo ; enddo
  enddo

  call find_eta(h, tv, G%g_Earth, G, e, halo_size=1)

  do j=js,je ; do I=is-1,ie
    uhtot(I,j) = 0.0 ; Work_u(I,j) = 0.0
    ! sfn_x(I,j,1) = 0.0 ; sfn_slope_x(I,j,1) = 0.0
    ! sfn_x(I,j,nz+1) = 0.0 ; sfn_slope_x(I,j,nz+1) = 0.0
  enddo ; enddo
  do J=js-1,je ; do i=is,ie
    vhtot(i,J) = 0.0 ; Work_v(i,J) = 0.0
    ! sfn_y(i,J,1) = 0.0 ; sfn_slope_y(i,J,1) = 0.0
    ! sfn_y(i,J,nz+1) = 0.0 ; sfn_slope_y(i,J,nz+1) = 0.0
  enddo ; enddo

  do k=nz,2,-1
    if (find_work .and. .not.(use_EOS .or. CS%bulkmixedlayer)) then
      drdiA = 0.0 ; drdiB = 0.0 ; drdjA = 0.0 ; drdjB = 0.0
      drdkL = G%g_prime(k) ; drdkR = G%g_prime(k)
    endif

    ! Calculate the zonal fluxes and gradients.
    if (use_EOS .and. ((k > nk_linear) .or. find_work)) then
      do j=js,je
        do I=is-1,ie
          pres_u(I) = 0.5*(pres(i,j,k) + pres(i+1,j,k))
          T_u(I) = 0.25*((T(i,j,k) + T(i+1,j,k)) + (T(i,j,k-1) + T(i+1,j,k-1)))
          S_u(I) = 0.25*((S(i,j,k) + S(i+1,j,k)) + (S(i,j,k-1) + S(i+1,j,k-1)))
        enddo
        call calculate_density_derivs(T_u, S_u, pres_u, drho_dT_u(:,j), &
                     drho_dS_u(:,j), (is-Isdq+1)-1, ie-is+2, tv%eqn_of_state)
      enddo
    endif

    do j=js,je ; do I=is-1,ie
      if ((k > nk_linear) .or. find_work) then
        if (use_EOS) then
          ! Estimate the horizontal density gradients along layers.
          drdiA = drho_dT_u(I,j) * (T(i+1,j,k-1)-T(i,j,k-1)) + &
                  drho_dS_u(I,j) * (S(i+1,j,k-1)-S(i,j,k-1))
          drdiB = drho_dT_u(I,j) * (T(i+1,j,k)-T(i,j,k)) + &
                  drho_dS_u(I,j) * (S(i+1,j,k)-S(i,j,k))

          ! Estimate the vertical density gradients times the grid spacing.
          drdkL = (drho_dT_u(I,j) * (T(i,j,k)-T(i,j,k-1)) + &
                   drho_dS_u(I,j) * (S(i,j,k)-S(i,j,k-1)))
          drdkR = (drho_dT_u(I,j) * (T(i+1,j,k)-T(i+1,j,k-1)) + &
                   drho_dS_u(I,j) * (S(i+1,j,k)-S(i+1,j,k-1)))
        elseif (CS%bulkmixedlayer) then
          drdiA = Rho(i+1,j,k-1) - Rho(i,j,k-1)
          drdiB = Rho(i+1,j,k) - Rho(i,j,k)
          drdkL = Rho(i,j,k) - Rho(i,j,k-1)
          drdkR = Rho(i+1,j,k) - Rho(i+1,j,k-1)
        endif
      endif


      if (k > nk_linear) then
        if (use_EOS .or. CS%bulkmixedlayer) then
          hg2A = h(i,j,k-1)*h(i+1,j,k-1) + h_neglect2
          hg2B = h(i,j,k)*h(i+1,j,k) + h_neglect2
          hg2L = h(i,j,k-1)*h(i,j,k) + h_neglect2
          hg2R = h(i+1,j,k-1)*h(i+1,j,k) + h_neglect2
          haA = 0.5*(h(i,j,k-1) + h(i+1,j,k-1))
          haB = 0.5*(h(i,j,k) + h(i+1,j,k)) + h_neglect
          haL = 0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect
          haR = 0.5*(h(i+1,j,k-1) + h(i+1,j,k)) + h_neglect
          if (G%Boussinesq) then
            dzaL = haL * G%H_to_m ; dzaR = haR * G%H_to_m
          else
            dzaL = 0.5*(e(i,j,k-1) - e(i,j,k+1)) + dz_neglect
            dzaR = 0.5*(e(i+1,j,k-1) - e(i+1,j,k+1)) + dz_neglect
          endif
          ! Use the harmonic mean thicknesses to weight the horizontal gradients.
          ! These unnormalized weights have been rearranged to minimize divisions.
          wtA = hg2A*haB ; wtB = hg2B*haA
          wtL = hg2L*(haR*dzaR) ; wtR = hg2R*(haL*dzaL)

          drdz = (wtL * drdkL + wtR * drdkR) / (dzaL*wtL + dzaR*wtR)
          ! The expression for drdz above is mathematically equivalent to:
          !   drdz = ((hg2L/haL) * drdkL/dzaL + (hg2R/haR) * drdkR/dzaR) / &
          !          ((hg2L/haL) + (hg2R/haR))
          ! This is the gradient of density along geopotentials.
          drdx = ((wtA * drdiA + wtB * drdiB) / (wtA + wtB) - &
                  drdz * (e(i,j,k)-e(i+1,j,k))) * G%IDXu(I,j)

          ! This estimate of slope is accurate for small slopes, but bounded
          ! to be between -1 and 1.
          mag_grad2 = drdx**2 + drdz**2
          if (mag_grad2 > 0.0) then
            Slope = drdx / sqrt(mag_grad2)
            slope2_Ratio = Slope**2 * I_slope_max2
          else ! Just in case mag_grad2 = 0 ever.
            Slope = 0.0
            slope2_Ratio = 1.0e20  ! Force the use of the safe streamfunction.
          endif

          ! Estimate the streamfunction at each interface.
          Sfn_unlim = -((KH_u(I,j)*G%dy_u(I,j))*Slope) * G%m_to_H
          if (uhtot(I,j) <= 0.0) then
            ! The transport that must balance the transport below is positive.
            Sfn_safe = uhtot(I,j) * (1.0 - h_frac(i,j,k))
          else !  (uhtot(I,j) > 0.0)
            Sfn_safe = uhtot(I,j) * (1.0 - h_frac(i+1,j,k))
          endif

          ! Avoid moving dense water upslope from below the level of
          ! the bottom on the receiving side.
          if (Sfn_unlim > 0.0) then ! The flow below this interface is positive.
            if (e(i,j,k) < e(i+1,j,nz+1)) then
              Sfn_unlim = 0.0 ! This is not uhtot, because it may compensate for
                              ! deeper flow in very unusual cases.
            elseif (e(i+1,j,nz+1) > e(i,j,k+1)) then
              ! Scale the transport with the fraction of the donor layer above
              ! the bottom on the receiving side.
              Sfn_unlim = Sfn_unlim * ((e(i,j,k) - e(i+1,j,nz+1)) / &
                                       ((e(i,j,k) - e(i,j,k+1)) + dz_neglect))
            endif
          else
            if (e(i+1,j,k) < e(i,j,nz+1)) then ; Sfn_unlim = 0.0
            elseif (e(i,j,nz+1) > e(i+1,j,k+1)) then
              Sfn_unlim = Sfn_unlim * ((e(i+1,j,k) - e(i,j,nz+1)) / &
                                     ((e(i+1,j,k) - e(i+1,j,k+1)) + dz_neglect))
            endif
          endif

          Sfn_est = (Sfn_unlim + slope2_Ratio*Sfn_safe) / (1.0 + slope2_Ratio)
        else  ! The layers are constant density.
          Sfn_est = ((KH_u(I,j)*G%dy_u(I,j)) * &
                     ((e(i,j,k)-e(i+1,j,k))*G%IDXu(I,j))) * G%m_to_H
        endif

        ! Make sure that there is enough mass above to allow the streamfunction
        ! to satisfy the boundary condition of 0 at the surface.
        Sfn = min(max(Sfn_est, -h_avail_rsum(i,j,k)), h_avail_rsum(i+1,j,k))

        ! The actual transport is limited by the mass available in the two
        ! neighboring grid cells.
        uhD(I,j,k) = max(min((Sfn - uhtot(I,j)), h_avail(i,j,k)), &
                         -h_avail(i+1,j,k))

 !       sfn_x(I,j,k) = max(min(Sfn, uhtot(I,j)+h_avail(i,j,k)), &
 !                          uhtot(I,j)-h_avail(i+1,j,k))
 !       sfn_slope_x(I,j,k) = max(uhtot(I,j)-h_avail(i+1,j,k), &
 !                                min(uhtot(I,j)+h_avail(i,j,k), &
 !             min(h_avail_rsum(i+1,j,k), max(-h_avail_rsum(i,j,k), &
 !             (KH_u(I,j)*G%dy_u(I,j)) * ((e(i,j,k)-e(i+1,j,k))*G%IDXu(I,j)) )) ))
      else ! k <= nk_linear
        ! Balance the deeper flow with a return flow uniformly distributed
        ! though the remaining near-surface layers.  This is the same as
        ! using Sfn_safe above.  There is no need to apply the limiters in
        ! this case.
        if (uhtot(I,j) <= 0.0) then
          uhD(I,j,k) = -uhtot(I,j) * h_frac(i,j,k)
        else !  (uhtot(I,j) > 0.0)
          uhD(I,j,k) = -uhtot(I,j) * h_frac(i+1,j,k)
        endif

 !       sfn_x(I,j,k) = sfn_x(I,j,k+1) + uhD(I,j,k)
 !       if (sfn_slope_x(I,j,k+1) <= 0.0) then
 !         sfn_slope_x(I,j,k) = sfn_slope_x(I,j,k+1) * (1.0 - h_frac(i,j,k))
 !       else
 !         sfn_slope_x(I,j,k) = sfn_slope_x(I,j,k+1) * (1.0 - h_frac(i+1,j,k))
 !       endif
      endif

      uhtot(I,j) = uhtot(I,j) + uhD(I,j,k)

      if (find_work) then
        !   This is the energy tendency based on the original profiles, and does
        ! not include any nonlinear terms due to a finite time step (which would
        ! involve interactions between the fluxes through the different faces.
        !   A second order centered estimate is used for the density transfered
        ! between water columns.

        Work_u(I,j) = Work_u(I,j) + G_scale * &
          ( uhtot(I,j) * (drdkR * e(i+1,j,k) - drdkL * e(i,j,k)) - &
            (uhD(I,j,k) * drdiB) * 0.25 * &
            ((e(i,j,k) + e(i,j,k+1)) + (e(i+1,j,k) + e(i+1,j,k+1))) )
      endif

    enddo ; enddo

    ! Calculate the meridional fluxes and gradients.
    if (use_EOS .and. ((k > nk_linear) .or. find_work)) then
      do J=js-1,je
        do i=is,ie
          pres_v(i) = 0.5*(pres(i,j,k) + pres(i,j+1,k))
          T_v(i) = 0.25*((T(i,j,k) + T(i,j+1,k)) + (T(i,j,k-1) + T(i,j+1,k-1)))
          S_v(i) = 0.25*((S(i,j,k) + S(i,j+1,k)) + (S(i,j,k-1) + S(i,j+1,k-1)))
        enddo
        call calculate_density_derivs(T_v, S_v, pres_v, drho_dT_v(:,J), &
                     drho_dS_v(:,J), is, ie-is+1, tv%eqn_of_state)
      enddo
    endif
    do J=js-1,je ; do i=is,ie
      if ((k > nk_linear) .or. find_work) then
        if (use_EOS) then
          ! Estimate the horizontal density gradients along layers.
          drdjA = drho_dT_v(i,J) * (T(i,j+1,k-1)-T(i,j,k-1)) + &
                  drho_dS_v(i,J) * (S(i,j+1,k-1)-S(i,j,k-1))
          drdjB = drho_dT_v(i,J) * (T(i,j+1,k)-T(i,j,k)) + &
                  drho_dS_v(i,J) * (S(i,j+1,k)-S(i,j,k))

          ! Estimate the vertical density gradients times the grid spacing.
          drdkL = (drho_dT_v(i,J) * (T(i,j,k)-T(i,j,k-1)) + &
                   drho_dS_v(i,J) * (S(i,j,k)-S(i,j,k-1)))
          drdkR = (drho_dT_v(i,J) * (T(i,j+1,k)-T(i,j+1,k-1)) + &
                   drho_dS_v(i,J) * (S(i,j+1,k)-S(i,j+1,k-1)))
        elseif (CS%bulkmixedlayer) then
          drdjA = Rho(i,j+1,k-1) - Rho(i,j,k-1)
          drdjB = Rho(i,j+1,k) - Rho(i,j,k)
          drdkL = Rho(i,j,k) - Rho(i,j,k-1)
          drdkR = Rho(i,j+1,k) - Rho(i,j+1,k-1)
        endif
      endif

      if (k > nk_linear) then
        if (use_EOS .or. CS%bulkmixedlayer) then
          hg2A = h(i,j,k-1)*h(i,j+1,k-1) + h_neglect2
          hg2B = h(i,j,k)*h(i,j+1,k) + h_neglect2
          hg2L = h(i,j,k-1)*h(i,j,k) + h_neglect2
          hg2R = h(i,j+1,k-1)*h(i,j+1,k) + h_neglect2
          haA = 0.5*(h(i,j,k-1) + h(i,j+1,k-1)) + h_neglect
          haB = 0.5*(h(i,j,k) + h(i,j+1,k)) + h_neglect
          haL = 0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect
          haR = 0.5*(h(i,j+1,k-1) + h(i,j+1,k)) + h_neglect
          if (G%Boussinesq) then
            dzaL = haL * G%H_to_m ; dzaR = haR * G%H_to_m
          else
            dzaL = 0.5*(e(i,j,k-1) - e(i,j,k+1)) + dz_neglect
            dzaR = 0.5*(e(i,j+1,k-1) - e(i,j+1,k+1)) + dz_neglect
          endif
          ! Use the harmonic mean thicknesses to weight the horizontal gradients.
          ! These unnormalized weights have been rearranged to minimize divisions.
          wtA = hg2A*haB ; wtB = hg2B*haA
          wtL = hg2L*(haR*dzaR) ; wtR = hg2R*(haL*dzaL)

          drdz = (wtL * drdkL + wtR * drdkR) / (dzaL*wtL + dzaR*wtR)
          ! The expression for drdz above is mathematically equivalent to:
          !   drdz = ((hg2L/haL) * drdkL/dzaL + (hg2R/haR) * drdkR/dzaR) / &
          !          ((hg2L/haL) + (hg2R/haR))
          ! This is the gradient of density along geopotentials.
          drdy = ((wtA * drdjA + wtB * drdjB) / (wtA + wtB) - &
                  drdz * (e(i,j,k)-e(i,j+1,k))) * G%IDYv(i,J)

          ! This estimate of slope is accurate for small slopes, but bounded
          ! to be between -1 and 1.
          mag_grad2 = drdy**2 + drdz**2
          if (mag_grad2 > 0.0) then
            Slope = drdy / sqrt(mag_grad2)
            slope2_Ratio = Slope**2 * I_slope_max2
          else ! Just in case mag_grad2 = 0 ever.
            Slope = 0.0
            slope2_Ratio = 1.0e20  ! Force the use of the safe streamfunction.
          endif

          ! Estimate the streamfunction at each interface.
          Sfn_unlim = -((KH_v(i,J)*G%dx_v(i,J))*Slope) * G%m_to_H
          if (vhtot(i,J) <= 0.0) then
            Sfn_safe = vhtot(i,J) * (1.0 - h_frac(i,j,k))
          else !  (vhtot(I,j) > 0.0)
            Sfn_safe = vhtot(i,J) * (1.0 - h_frac(i,j+1,k))
          endif

          ! Avoid moving dense water upslope from below the level of
          ! the bottom on the receiving side.
          if (Sfn_unlim > 0.0) then ! The flow below this interface is positive.
            if (e(i,j,k) < e(i,j+1,nz+1)) then
              Sfn_unlim = 0.0 ! This is not vhtot, because it may compensate for
                              ! deeper flow in very unusual cases.
            elseif (e(i,j+1,nz+1) > e(i,j,k+1)) then
              ! Scale the transport with the fraction of the donor layer above
              ! the bottom on the receiving side.
              Sfn_unlim = Sfn_unlim * ((e(i,j,k) - e(i,j+1,nz+1)) / &
                                       ((e(i,j,k) - e(i,j,k+1)) + dz_neglect))
            endif
          else
            if (e(i,j+1,k) < e(i,j,nz+1)) then ; Sfn_unlim = 0.0
            elseif (e(i,j,nz+1) > e(i,j+1,k+1)) then
              Sfn_unlim = Sfn_unlim * ((e(i,j+1,k) - e(i,j,nz+1)) / &
                                     ((e(i,j+1,k) - e(i,j+1,k+1)) + dz_neglect))
            endif
          endif

          ! Estimate the streamfunction at each interface.
          Sfn_est = (Sfn_unlim + slope2_Ratio*Sfn_safe) / (1.0 + slope2_Ratio)
        else  ! The layers are constant density.
          Sfn_est = ((KH_v(i,J)*G%dx_v(i,J)) * &
                     ((e(i,j,k)-e(i,j+1,k))*G%IDYv(i,J))) * G%m_to_H
        endif

        ! Make sure that there is enough mass above to allow the streamfunction
        ! to satisfy the boundary condition of 0 at the surface.
        Sfn = min(max(Sfn_est, -h_avail_rsum(i,j,k)), h_avail_rsum(i,j+1,k))

        ! The actual transport is limited by the mass available in the two
        ! neighboring grid cells.
        vhD(i,J,k) = max(min((Sfn - vhtot(i,J)), h_avail(i,j,k)), &
                         -h_avail(i,j+1,k))

  !      sfn_y(i,J,k) = max(min(Sfn, vhtot(i,J)+h_avail(i,j,k)), &
  !                         vhtot(i,J)-h_avail(i,j+1,k))
  !      sfn_slope_y(i,J,k) = max(vhtot(i,J)-h_avail(i,j+1,k), &
  !                               min(vhtot(i,J)+h_avail(i,j,k), &
  !            min(h_avail_rsum(i,j+1,k), max(-h_avail_rsum(i,j,k), &
  !            (KH_v(i,J)*G%dx_v(i,J)) * ((e(i,j,k)-e(i,j+1,k))*G%IDYv(i,J)) )) ))
      else  ! k <= nk_linear
        ! Balance the deeper flow with a return flow uniformly distributed
        ! though the remaining near-surface layers.
        if (vhtot(i,J) <= 0.0) then
          vhD(i,J,k) = -vhtot(i,J) * h_frac(i,j,k)
        else !  (vhtot(i,J) > 0.0)
          vhD(i,J,k) = -vhtot(i,J) * h_frac(i,j+1,k)
        endif

  !      sfn_y(i,J,k) = sfn_y(i,J,k+1) + vhD(i,J,k)
  !      if (sfn_slope_y(i,J,k+1) <= 0.0) then
  !        sfn_slope_y(i,J,k) = sfn_slope_y(i,J,k+1) * (1.0 - h_frac(i,j,k))
  !      else
  !        sfn_slope_y(i,J,k) = sfn_slope_y(i,J,k+1) * (1.0 - h_frac(i,j+1,k))
  !      endif
      endif

      vhtot(i,J) = vhtot(i,J)  + vhD(i,J,k)

      if (find_work) then
        !   This is the energy tendency based on the original profiles, and does
        ! not include any nonlinear terms due to a finite time step (which would
        ! involve interactions between the fluxes through the different faces.
        !   A second order centered estimate is used for the density transfered
        ! between water columns.

        Work_v(i,J) = Work_v(i,J) + G_scale * &
          ( vhtot(i,J) * (drdkR * e(i,j+1,k) - drdkL * e(i,j,k)) - &
           (vhD(i,J,k) * drdjB) * 0.25 * &
           ((e(i,j,k) + e(i,j,k+1)) + (e(i,j+1,k) + e(i+1,j,k+1))) )
      endif
    enddo ; enddo
  enddo ! k-loop

  ! In layer 1, enforce the boundary conditions that Sfn(z=0) = 0.0
  if (.not.find_work .or. .not.(use_EOS .or. CS%bulkmixedlayer)) then
    do j=js,je ; do I=is-1,ie ; uhD(I,j,1) = -uhtot(I,j) ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; vhD(i,J,1) = -vhtot(i,J) ; enddo ; enddo
  else
    if (use_EOS) then ; do j=js,je
      do I=is-1,ie
        pres_u(I) = 0.5*(pres(i,j,1) + pres(i+1,j,1))
        T_u(I) = 0.5*(T(i,j,1) + T(i+1,j,1))
        S_u(I) = 0.5*(S(i,j,1) + S(i+1,j,1))
      enddo
      call calculate_density_derivs(T_u, S_u, pres_u, drho_dT_u(:,j), &
                   drho_dS_u(:,j), (is-Isdq+1)-1, ie-is+2, tv%eqn_of_state)
    enddo ; endif
    do j=js,je ; do I=is-1,ie
      uhD(I,j,1) = -uhtot(I,j)

      if (use_EOS) then
        drdiB = drho_dT_u(I,j) * (T(i+1,j,1)-T(i,j,1)) + &
                drho_dS_u(I,j) * (S(i+1,j,1)-S(i,j,1))
      elseif (CS%bulkmixedlayer) then
        drdiB = Rho(i+1,j,1) - Rho(i,j,1)
      endif
      Work_u(I,j) = Work_u(I,j) + G_scale * ( (uhD(I,j,1) * drdiB) * 0.25 * &
          ((e(i,j,1) + e(i,j,2)) + (e(i+1,j,1) + e(i+1,j,2))) )

    enddo ; enddo

    if (use_EOS) then ; do J=js-1,je
      do i=is,ie
        pres_v(i) = 0.5*(pres(i,j,1) + pres(i,j+1,1))
        T_v(i) = 0.5*(T(i,j,1) + T(i,j+1,1))
        S_v(i) = 0.5*(S(i,j,1) + S(i,j+1,1))
      enddo
      call calculate_density_derivs(T_v, S_v, pres_v, drho_dT_v(:,J), &
                   drho_dS_v(:,J), is, ie-is+1, tv%eqn_of_state)
    enddo ; endif
    do J=js-1,je ; do i=is,ie
      vhD(i,J,1) = -vhtot(i,J)

      if (use_EOS) then
        drdjB = drho_dT_v(i,J) * (T(i,j+1,1)-T(i,j,1)) + &
                drho_dS_v(i,J) * (S(i,j+1,1)-S(i,j,1))
      elseif (CS%bulkmixedlayer) then
        drdjB = Rho(i,j+1,1) - Rho(i,j,1)
      endif
      Work_v(i,J) = Work_v(i,J) - G_scale * ( (vhD(i,J,1) * drdjB) * 0.25 * &
          ((e(i,j,1) + e(i,j,2)) + (e(i,j+1,1) + e(i,j+1,2))) )
    enddo ; enddo
  endif

  if (find_work) then ; do j=js,je ; do i=is,ie
    ! Note that the units of Work_v and Work_u are W, while Work_h is W m-2.
    Work_h = 0.25 * G%IDXDYh(i,j) * &
      ((Work_u(I-1,j) + Work_u(I,j)) + (Work_v(i,J-1) + Work_v(i,J)))
    if (ASSOCIATED(CS%GMwork)) CS%GMwork(i,j) = Work_h
    if (ASSOCIATED(MEKE)) then ; if (ASSOCIATED(MEKE%GM_src)) then
      MEKE%GM_src(i,j) = MEKE%GM_src(i,j) + Work_h
    endif ; endif
  enddo ; enddo ; endif

 ! if (CS%id_sfn_x > 0) call post_data(CS%id_sfn_x, sfn_x, CS%diag)
 ! if (CS%id_sfn_y > 0) call post_data(CS%id_sfn_y, sfn_y, CS%diag)
 ! if (CS%id_sfn_slope_x > 0) call post_data(CS%id_sfn_slope_x, sfn_slope_x, CS%diag)
 ! if (CS%id_sfn_slope_y > 0) call post_data(CS%id_sfn_slope_y, sfn_slope_y, CS%diag)

end subroutine thickness_diffuse_full

subroutine vert_fill_TS(h, T_in, S_in, kappa, dt, T_f, S_f, G, halo_here)
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(in)    :: h
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(in)    :: T_in
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(in)    :: S_in
  real,                                intent(in)    :: kappa
  real,                                intent(in)    :: dt
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(out)   :: T_f
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(out)   :: S_f
  type(ocean_grid_type),               intent(in)    :: G
  integer,                   optional, intent(in)    :: halo_here
!    This subroutine fills massless layers with sensible values of two
!*  tracer arrays (nominally temperature and salinity) by diffusing
!*  vertically with a (small?) constant diffusivity.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (in)      T_in - The input temperature, in K.
!  (in)      S_in - The input salinity, in psu.
!  (in)      kappa - The diapycnal diffusivity, in m2 s-1.
!  (in)      dt - Time increment in s.
!  (out)     T_f - The filled temperature, in K.
!  (out)     S_f - The filled salinity, in psu.
!  (in)      G - The ocean's grid structure.
!  (in,opt)  halo_here - the number of halo points to work on, 0 by default.

  real :: ent(SZI_(G),SZK_(G))     ! The diffusive entrainment (kappa*dt)/dz
                                   ! between layers in a timestep in m or kg m-2.
  real :: b1(SZI_(G)), d1(SZI_(G)) ! b1, c1, and d1 are variables used by the
  real :: c1(SZI_(G),SZK_(G))      ! tridiagonal solver.
  real :: kap_dt_x2                ! The product of 2*kappa*dt, converted to
                                   ! the same units as h, in m2 or kg2 m-4.
  real :: h0                       ! A negligible thickness, in m or kg m-2, to
                                   ! allow for zero thicknesses.
  integer :: i, j, k, is, ie, js, je, nz, halo

  halo=0 ; if (present(halo_here)) halo = max(halo_here,0)

  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  nz = G%ke

  kap_dt_x2 = (2.0*kappa*dt)*G%m_to_H**2
  h0 = 1.0e-16*sqrt(kappa*dt)*G%m_to_H

  if (kap_dt_x2 <= 0.0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      T_f(i,j,k) = T_in(i,j,k) ; S_f(i,j,k) = S_in(i,j,k)
    enddo ; enddo ; enddo
  else
    do j=js,je
      do i=is,ie
        ent(i,2) = kap_dt_x2 / ((h(i,j,1)+h(i,j,2)) + h0)
        b1(i) = 1.0 / (h(i,j,1)+ent(i,2))
        d1(i) = b1(i) * h(i,j,1)
        T_f(i,j,1) = (b1(i)*h(i,j,1))*T_in(i,j,1)
        S_f(i,j,1) = (b1(i)*h(i,j,1))*S_in(i,j,1)
      enddo
      do k=2,nz-1 ; do i=is,ie
        ent(i,k+1) = kap_dt_x2 / ((h(i,j,k)+h(i,j,k+1)) + h0)
        c1(i,k) = ent(i,k) * b1(i)
        b1(i) = 1.0 / ((h(i,j,k) + d1(i)*ent(i,k)) + ent(i,k+1))
        d1(i) = b1(i) * (h(i,j,k) + d1(i)*ent(i,k))
        T_f(i,j,k) = b1(i) * (h(i,j,k)*T_in(i,j,k) + ent(i,k)*T_f(i,j,k-1))
        S_f(i,j,k) = b1(i) * (h(i,j,k)*S_in(i,j,k) + ent(i,k)*S_f(i,j,k-1))
      enddo ; enddo
      do i=is,ie
        c1(i,nz) = ent(i,nz) * b1(i)
        b1(i) = 1.0 / (h(i,j,nz) + d1(i)*ent(i,nz))
        T_f(i,j,nz) = b1(i) * (h(i,j,nz)*T_in(i,j,nz) + ent(i,nz)*T_f(i,j,nz-1))
        S_f(i,j,nz) = b1(i) * (h(i,j,nz)*S_in(i,j,nz) + ent(i,nz)*S_f(i,j,nz-1))
      enddo
      do k=nz-1,1,-1 ; do i=is,ie
        T_f(i,j,k) = T_f(i,j,k) + c1(i,k+1)*T_f(i,j,k+1)
        S_f(i,j,k) = S_f(i,j,k) + c1(i,k+1)*S_f(i,j,k+1)
      enddo ; enddo
    enddo
  endif

end subroutine vert_fill_TS


subroutine thickness_diffuse_init(Time, G, param_file, diag, CS)
  type(time_type),       intent(in) :: Time
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  type(diag_ptrs), target, intent(inout) :: diag
  type(thickness_diffuse_CS),     pointer    :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
  integer :: nkbl

  character(len=128) :: version = '$Id: GOLD_thickness_diffuse.F90,v 13.0.2.4.2.41 2011/04/18 12:01:06 Robert.Hallberg Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "GOLD_thickness_diffuse" ! This module's name.
  character(len=48)  :: flux_units

  if (associated(CS)) then
    call GOLD_error(WARNING, &
      "Thickness_diffuse_init called with an associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag
  CS%thickness_diffuse = .false. ; CS%Khth = 0.0
  CS%full_thickness_diffuse = .false.
  call read_param(param_file,"THICKNESSDIFFUSE",CS%thickness_diffuse)
  if (CS%thickness_diffuse) call read_param(param_file,"KHTH",CS%Khth,.true.)
  CS%KHTH_Slope_Cff = 0.; call read_param(param_file,"KHTH_SLOPE_CFF",CS%KHTH_Slope_Cff)
  CS%KHTH_Min = 0.0 ; call read_param(param_file,"KHTH_MIN",CS%KHTH_Min)
  CS%KHTH_Max = 0.0 ; call read_param(param_file,"KHTH_MAX",CS%KHTH_Max)
  call read_param(param_file,"FULL_THICKNESSDIFFUSE",CS%full_thickness_diffuse)

  CS%diffuse_isopycnals = .true.
  call read_param(param_file,"DIFFUSE_ISOPYCNALS",CS%diffuse_isopycnals)
  CS%bulkmixedlayer = .false.
  call read_param(param_file,"BULKMIXEDLAYER",CS%bulkmixedlayer)
  if (CS%bulkmixedlayer) then
    CS%nkml = 1 ; call read_param(param_file,"NKML",CS%nkml)
  endif
  CS%slope_max = 1.0e-2 ; CS%kappa_smooth = 1.0e-6
  call read_param(param_file,"KHTH_SLOPE_MAX",CS%slope_max)
  call read_param(param_file,"KD_SMOOTH",CS%kappa_smooth)

  if (G%Boussinesq) then ; flux_units = "meter3 second-1"
  else ; flux_units = "kilogram second-1" ; endif

  CS%id_uhGM = register_diag_field('ocean_model', 'uhGM', G%axesul, Time, &
           'Time Mean Diffusive Zonal Thickness Flux', flux_units)
  if (CS%id_uhGM > 0) call safe_alloc_ptr(diag%uhGM,G%Isdq,G%Iedq,G%jsd,G%jed,G%ke)
  CS%id_vhGM = register_diag_field('ocean_model', 'vhGM', G%axesvl, Time, &
           'Time Mean Diffusive Meridional Thickness Flux', flux_units)
  if (CS%id_vhGM > 0) call safe_alloc_ptr(diag%vhGM,G%isd,G%ied,G%Jsdq,G%Jedq,G%ke)
  CS%id_GMwork = register_diag_field('ocean_model', 'GMwork', G%axesh1, Time, &
           'Time Mean Integral Work done by Diffusive Thickness Flux', 'Watt meter-2')
  if (CS%id_GMwork > 0) call safe_alloc_ptr(CS%GMwork,G%isd,G%ied,G%jsd,G%jed)
  CS%id_KH_u = register_diag_field('ocean_model', 'KHTH_u', G%axesu1, Time, &
           'Thickness Diffusivity at U-point', 'meter second-2')
  CS%id_KH_v = register_diag_field('ocean_model', 'KHTH_v', G%axesv1, Time, &
           'Thickness Diffusivity at V-point', 'meter second-2')

 ! CS%id_sfn_x =  register_diag_field('ocean_model', 'sfn_x', G%axesui, Time, &
 !          'Parameterized Zonal Overturning Streamfunction', 'meter3 second-1')
 ! CS%id_sfn_y =  register_diag_field('ocean_model', 'sfn_y', G%axesvi, Time, &
 !          'Parameterized Meridional Overturning Streamfunction', 'meter3 second-1')
 ! CS%id_sfn_slope_x =  register_diag_field('ocean_model', 'sfn_sl_x', G%axesui, Time, &
 !          'Parameterized Zonal Overturning Streamfunction from Interface Slopes', 'meter3 second-1')
 ! CS%id_sfn_slope_y =  register_diag_field('ocean_model', 'sfn_sl_y', G%axesvi, Time, &
 !          'Parameterized Meridional Overturning Streamfunction from Interface Slopes', 'meter3 second-1')

  ! Write all relevant parameters to the model log.
  call log_version(param_file, mod, version, tagname, "")
  call log_param(param_file, mod, "THICKNESSDIFFUSE", CS%thickness_diffuse, &
                 "If true, interfaces heights are diffused with a \n"//&
                 "coefficient of KHTH.", default=.false.)
  call log_param(param_file, mod, "FULL_THICKNESSDIFFUSE", CS%full_thickness_diffuse, &
                 "If true, use the new full-column code to calculate the \n"//&
                 "overturning residual circulation.  ###THIS DEFAULT SHOULD BE CHANGED", &
                 default=.false.)
  if (CS%thickness_diffuse) call log_param(param_file, mod, "KHTH", CS%Khth, &
                 "The background horizontal thickness diffusivity.", &
                 units = "m2 s-1", default=0.0)
  call log_param(param_file, mod, "KHTH_SLOPE_CFF", CS%KHTH_Slope_Cff, &
                 "The nondimensional coefficient in the Visbeck formula \n"//&
                 "for the interface depth diffusivity", units="nondim", &
                 default=0.0)
  call log_param(param_file, mod, "KHTH_MIN", CS%KHTH_Min, &
                 "The minimum horizontal thickness diffusivity.", &
                 units = "m2 s-1", default=0.0)
  call log_param(param_file, mod, "KHTH_MAX", CS%KHTH_Max, &
                 "The maximum horizontal thickness diffusivity.", &
                 units = "m2 s-1", default=0.0)
  call log_param(param_file, mod, "DIFFUSE_ISOPYCNALS", CS%diffuse_isopycnals, &
                 "If true, , within the mixed layer, the overturning \n"//&
                 "velocities are constant.  Otherwise the mass fluxes \n"//&
                 "flatten all interfaces between layers.", default=.true.)
  call log_param(param_file, mod, "BULKMIXEDLAYER", CS%bulkmixedlayer, &
                 "If defined, use a refined Kraus-Turner-like bulk mixed \n"//&
                 "layer with transitional buffer layers.")
  call log_param(param_file, mod, "KHTH_SLOPE_MAX", CS%slope_max, &
                 "A slope beyond which the calculated isopycnal slope is \n"//&
                 "not reliable and is scaled away. This is used with \n"//&
                 "FULL_THICKNESSDIFFUSE.", units="nondim", default=0.01)
  call log_param(param_file, mod, "KD_SMOOTH", CS%kappa_smooth, &
                 "A diapycnal diffusivity that is used to interpolate \n"//&
                 "more sensible values of T & S into thin layers.", &
                 default=1.0e-6)
  if (CS%bulkmixedlayer) then
    call log_param(param_file, mod, "NKML",CS%nkml, &
                 "The number of sublayers within the mixed layer if \n"//&
                 "BULKMIXEDLAYER is true.", units="nondim", default=1)
  endif

end subroutine thickness_diffuse_init

subroutine thickness_diffuse_end(CS)
  type(thickness_diffuse_CS), pointer :: CS
  if(associated(CS)) deallocate(CS)
end subroutine thickness_diffuse_end

end module GOLD_thickness_diffuse
