module GOLD_mixed_layer_restrat
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
!
!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, June 2002, 2006                                *
!*                                                                     *
!*    The subroutine in this file implements a parameterization of     *
!*  unresolved viscous mixed layer restratification of the mixed layer *
!*  as described in Hallberg (Aha Hulikoa, 2003).  This is based on the*
!*  subinertial mixed layer theory of Young (JPO, 1994).  There is no  *
!*  net horizontal volume transport due to this parameterization, and  *
!*  no direct effect below the mixed layer. If OLD_RESTRAT_PARAM is    *
!*  defined, the strength of the parameterization is determined by the *
!*  coefficient ML_RESTRAT_COEF from GOLD_input, with units of s, that *
!*  essentially replaces f^-1 in describing the strength of the        *
!*  unresolved overturning circulation. As this circulation should go  *
!*  as the square of the horizontal density gradient, this may         *
!*  justifiably be larger than f^-1.  There should probably be a       *
!*  dependence on both latitude and the ratio of the grid-spacing to   *
!*  the deformation radius, but this older theory does not include     *
!*  these factors.                                                     *
!*                                                                     *
!*    A newer option, developed in detail by Fox-Kemper, et al., is    *
!*  also available with a run-time switch, to set the restratification *
!*  timescale to agree with his high-resolution studies of mixed layer *
!*  restratification.  In this case FOX_KEMPER_ML_RESTRAT_COEF is a    *
!*  nondimensional number of order a few tens, proportional to the     *
!*  ratio of the defromation radius to the dominant horizontal         *
!*  lengthscale of the submesoscale mixed layer instabilities.         *
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
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_grid, only : ocean_grid_type
use GOLD_variables, only : thermo_var_ptrs, forcing
use GOLD_EOS, only : calculate_density

implicit none ; private

#include <GOLD_memory.h>

public mixedlayer_restrat, mixedlayer_restrat_init

type, public :: mixedlayer_restrat_CS ; private
  logical :: bulkmixedlayer  ! If true, a refined bulk mixed layer is used.
  real    :: ml_restrat_coef ! (New param) A nondimensional factor by which the 
                             ! instability is enhanced over what would be
                             ! predicted based on the resolved  gradients.  This
                             ! increases with grid spacing^2, up to something
                             ! of order 500.
                             ! (Old param) A coefficient in s (perhaps OMEGA-1)
                             ! relating the mixed layer restratification to the
                             ! horizontal density gradients.
  integer :: nkml            ! The number of sublayers in the mixed layer.
  logical :: old_restrat     ! If true, the restratification uses the older
                             ! form, predating the Fox-Kemper et al. form.
  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
                             ! ocean diagnostic fields and control variables.
  integer :: id_urestrat_time , id_vrestrat_time 
  integer :: id_uhml = -1, id_vhml = -1
end type mixedlayer_restrat_CS

contains

subroutine mixedlayer_restrat(h, uhtr, vhtr, tv, fluxes, dt, G, CS)
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(inout) :: h
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(inout) :: uhtr
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(inout) :: vhtr
  type(thermo_var_ptrs),               intent(in)    :: tv
  type(forcing),                       intent(in)    :: fluxes
  real,                                intent(in)    :: dt
  type(ocean_grid_type),               intent(in)    :: G
  type(mixedlayer_restrat_CS),         pointer       :: CS
!    This subroutine does interface depth diffusion.  The fluxes are
!  limited to give positive definiteness, and the diffusivities are
!  limited to guarantee stability.

! Arguments: h - Layer thickness, in m or kg m-2.  (Intent in/out.)
!                The units of h are referred to as H below.
!  (in/out)  uhtr - Accumulated zonal mass fluxes in m3 or kg.
!  (in/out)  vhtr - Accumulated meridional mass fluxes in m3 or kg.
!  (in)      tv - A structure containing the thermobaric variables.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 mixedlayer_restrat_init.

  real :: uhml(SZIQ_(G),SZJ_(G),SZK_(G)) ! The zonal and meridional mixed layer
  real :: vhml(SZI_(G),SZJQ_(G),SZK_(G)) ! fluxes, in m3 s-1 or kg s-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    h_avail       ! The volume available for diffusion out of each face of each
                  ! sublayer of the mixed layer, divided by dt, in units
                  ! of H m2 s-1 (i.e., m3 s-1 or kg s-1).
  real, dimension(SZI_(G),SZJ_(G)) :: &
    htot, &       ! The sum of the thicknesses of layers in the mixed layer, H.
    Rml_av        ! g_Rho0 times the average mixed layer density, in m s-2.
  real :: g_Rho0  ! G_Earth/Rho0 in m4 s-2 kg-1.
  real :: Rho0(SZI_(G)) ! Potential density relative to the surface, in kg m-3.
  real :: p0(SZI_(G))   ! A pressure of 0, in Pa.

  real :: h_vel         ! htot interpolated onto velocity points, in m (not H).
  real :: absf          ! The absolute value of f, interpolated to velocity
                        ! points, in s-1.
  real :: u_star        ! The surface friction velocity, interpolated to velocity
                        ! points, in m s-1.
  real :: mom_mixrate   ! The rate at which momentum is homogenized within the
                        ! mixed layer in s-1.
  real :: timescale     ! The mixing growth timescale in s.
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected, in H.
  real :: dz_neglect    ! A thickness in m that is so small it is usually lost
                        ! in roundoff and can be neglected, in m.
  real :: I4dt          ! 1 / 4 dt
  real :: I2htot        ! Twice the total mixed layer thickness,
  real :: z_topx2       ! depth of the top of a layer, and
  real :: hx2           ! layer thickness, all at velocity points and in H.
  real :: a(SZK_(G))    ! A nondimensional value relating the overall flux
                        ! magnitudes (uDml & vDml) to the realized flux in a
                        ! layer.  The vertical sum of a() through the pieces of
                        ! the mixed layer must be 0.
  real :: uDml(SZIQ_(G))  ! The zonal and meridional volume fluxes in the upper
  real :: vDml(SZI_(G))   ! half of the mixed layer in H m2 s-1 (m3 s-1 or kg s-1).
  real :: utimescale_diag(SZIQ_(G),SZJ_(G)) ! The restratification timescales
  real :: vtimescale_diag(SZI_(G),SZJQ_(G)) ! in the zonal and meridional
                                            ! directions, in s, stored in 2-D
                                            ! arrays for diagnostic purposes.
  real, pointer, dimension(:,:,:) :: Rml
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq

  if (.not. associated(CS)) call GOLD_error(FATAL, "GOLD_mixedlayer_restrat: "// &
         "Module must be initialized before it is used.")
  if ((.not.CS%bulkmixedlayer) .or. (CS%nkml<2) .or. (CS%ml_restrat_coef<=0.0)) &
    return

  Rml => tv%Rml
  uDml(:) = 0.0 ; vDml(:) = 0.0
  I4dt = 0.25 / dt
  g_Rho0 = G%g_Earth/G%Rho0
  use_EOS = associated(tv%eqn_of_state)
  h_neglect = G%H_subroundoff
  dz_neglect = G%H_subroundoff*G%H_to_m

  ! Fix this later for CS%nkml >= 3.

  if (use_EOS) p0(:) = 0.0

  do j=js-1,je+1
    do i=is-1,ie+1
      htot(i,j) = 0.0 ; Rml_av(i,j) = 0.0
    enddo
    ! A bug was fixed here... nkml was double counted before, effectively
    ! increasing CONSTANT by a factor of 2.25.
    do k=1,CS%nkml
      if (use_EOS) then
        call calculate_density(tv%T(:,j,k),tv%S(:,j,k),p0,Rho0(:),is-1,ie-is+3,tv%eqn_of_state)
        do i=is-1,ie+1
          Rml_av(i,j) = Rml_av(i,j) + h(i,j,k)*Rho0(i)
          htot(i,j) = htot(i,j) + h(i,j,k)
          h_avail(i,j,k) = max(I4dt*G%DXDYh(i,j)*(h(i,j,k)-G%Angstrom),0.0)
        enddo
      else
        do i=is-1,ie+1
          Rml_av(i,j) = Rml_av(i,j) + h(i,j,k)*Rml(i,j,k)
          htot(i,j) = htot(i,j) + h(i,j,k)
          h_avail(i,j,k) = max(I4dt*G%DXDYh(i,j)*(h(i,j,k)-G%Angstrom),0.0)
        enddo
      endif
    enddo

    do i=is-1,ie+1
      Rml_av(i,j) = (g_Rho0*Rml_av(i,j)) / (htot(i,j) + h_neglect)
    enddo
  enddo

! TO DO:
!   1. Mixing extends below the mixing layer to the mixed layer.  Find it!
!   2. Add exponential tail to streamfunction?

  do j=js,je ; do i=is,ie ; utimescale_diag(i,j) = 0.0 ; enddo ; enddo
  do j=js,je ; do i=is,ie ; vtimescale_diag(i,j) = 0.0 ; enddo ; enddo

!   U - Component
  do j=js,je ; do I=is-1,ie
    h_vel = 0.5*(htot(i,j) + htot(i+1,j)) * G%H_to_m
    if (CS%old_restrat) then
      timescale = CS%ml_restrat_coef
    else
      u_star = 0.5*(fluxes%ustar(i,j) + fluxes%ustar(i+1,j))
      absf = 0.5*(abs(G%f(I,J-1)) + abs(G%f(I,J)))
      ! peak ML visc: u_star * 0.41 * (h_ml*u_star)/(absf*h_ml + 4.0*u_star)
      ! momentum mixing rate: pi^2*visc/h_ml^2 
      ! 0.41 is the von Karmen constant, 9.8696 = pi^2.
      mom_mixrate = (0.41*9.8696)*u_star**2 / &
                    (absf*h_vel**2 + 4.0*(h_vel+dz_neglect)*u_star)
      timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)

      timescale = timescale * CS%ml_restrat_coef
!        timescale = timescale*(2?)*(L_def/L_MLI)*min(EKE/MKE,1.0 + G%DYv(i,j)**2/L_def**2))
    endif
    utimescale_diag(I,j) = timescale

    uDml(I) = timescale * G%umask(I,j)*G%DYu(I,j)* &
        G%IDXu(I,j)*(Rml_av(i+1,j)-Rml_av(i,j)) * (h_vel**2 * G%m_to_H)

    if (uDml(i) == 0) then
      do k=1,CS%nkml ; uhml(I,j,k) = 0.0 ; enddo
    else
      I2htot = 1.0 / (htot(i,j) + htot(i+1,j) + h_neglect)
      z_topx2 = 0.0
      ! a(k) relates the sublayer transport to uDml with a linear profile.
      ! The sum of a through the mixed layers must be 0.
      do k=1,CS%nkml
        hx2 = (h(i,j,k) + h(i+1,j,k) + h_neglect)
        a(k) = (hx2 * I2htot) * (2.0 - 4.0*(z_topx2+0.5*hx2)*I2htot)
        z_topx2 = z_topx2 + hx2
        if (a(k)*uDml(I) > 0.0) then
          if (a(k)*uDml(I) > h_avail(i,j,k)) uDml(I) = h_avail(i,j,k) / a(k)
        else
          if (-a(k)*uDml(I) > h_avail(i+1,j,k)) uDml(I) = -h_avail(i+1,j,k)/a(k)
        endif
      enddo
      do k=1,CS%nkml
        uhml(I,j,k) = a(k)*uDml(I)
        uhtr(I,j,k) = uhtr(I,j,k) + uhml(I,j,k)*dt
      enddo
    endif
  enddo ; enddo

!  V- component
  do J=js-1,je ; do i=is,ie
    h_vel = 0.5*(htot(i,j) + htot(i,j+1)) * G%H_to_m
    if (CS%old_restrat) then
      timescale = CS%ml_restrat_coef
    else
      u_star = 0.5*(fluxes%ustar(i,j) + fluxes%ustar(i,j+1))
      absf = 0.5*(abs(G%f(I-1,J)) + abs(G%f(I,J)))
      ! peak ML visc: u_star * 0.41 * (h_ml*u_star)/(absf*h_ml + 4.0*u_star)
      ! momentum mixing rate: pi^2*visc/h_ml^2 
      ! 0.41 is the von Karmen constant, 9.8696 = pi^2.
      mom_mixrate = (0.41*9.8696)*u_star**2 / &
                    (absf*h_vel**2 + 4.0*(h_vel+dz_neglect)*u_star)
      timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)

      timescale = timescale * CS%ml_restrat_coef
!        timescale = timescale*(2?)*(L_def/L_MLI)*min(EKE/MKE,1.0 + G%DYv(i,j)**2/L_def**2))
    endif
    vtimescale_diag(i,J) = timescale

    vDml(i) = timescale * G%vmask(i,J)*G%DXv(i,J)* &
        G%IDYv(i,J)*(Rml_av(i,j+1)-Rml_av(i,j)) * (h_vel**2 * G%m_to_H)
    if (vDml(i) == 0) then
      do k=1,CS%nkml ; vhml(i,J,k) = 0.0 ; enddo
    else
      I2htot = 1.0 / (htot(i,j) + htot(i,j+1) + h_neglect)
      z_topx2 = 0.0
      ! a relates the sublayer transport to uDml with a linear profile.
      ! The sum of a through the mixed layers must be 0.
      do k=1,CS%nkml
        hx2 = (h(i,j,k) + h(i,j+1,k) + h_neglect)
        a(k) = (hx2 * I2htot) * (2.0 - 4.0*(z_topx2+0.5*hx2)*I2htot)
        z_topx2 = z_topx2 + hx2
        if (a(k)*vDml(i) > 0.0) then
          if (a(k)*vDml(i) > h_avail(i,j,k)) vDml(i) = h_avail(i,j,k) / a(k)
        else
          if (-a(k)*vDml(i) > h_avail(i,j+1,k)) vDml(i) = -h_avail(i,j+1,k)/a(k)
        endif
      enddo
      do k=1,CS%nkml
        vhml(i,J,k) = a(k)*vDml(i)
        vhtr(i,J,k) = vhtr(i,J,k) + vhml(i,J,k)*dt
      enddo
    endif
  enddo ; enddo

  do k=1,CS%nkml ; do j=js,je ; do i=is,ie
    h(i,j,k) = h(i,j,k) - dt*G%IDXDYh(i,j) * &
        ((uhml(I,j,k) - uhml(I-1,j,k)) + (vhml(i,J,k) - vhml(i,J-1,k)))
  enddo ; enddo ; enddo

! Offer fields for averaging.
  if (query_averaging_enabled(CS%diag) .and. &
      ((CS%id_urestrat_time > 0)  .or. (CS%id_vrestrat_time > 0))) then
    call post_data(CS%id_urestrat_time, utimescale_diag, CS%diag)
    call post_data(CS%id_vrestrat_time, vtimescale_diag, CS%diag)
  endif
  if (query_averaging_enabled(CS%diag) .and. &
      ((CS%id_uhml>0) .or. (CS%id_vhml>0))) then
    do k=CS%nkml+1,nz
      do j=js,je ; do I=Isq,Ieq ; uhml(I,j,k) = 0.0 ; enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie ; vhml(i,J,k) = 0.0 ; enddo ; enddo
    enddo
    if (CS%id_uhml > 0) call post_data(CS%id_uhml, uhml, CS%diag)
    if (CS%id_vhml > 0) call post_data(CS%id_vhml, vhml, CS%diag)
  endif

end subroutine mixedlayer_restrat

subroutine mixedlayer_restrat_init(Time, G, param_file, diag, CS)
  type(time_type),             intent(in)    :: Time
  type(ocean_grid_type),       intent(in)    :: G
  type(param_file_type),       intent(in)    :: param_file
  type(diag_ptrs), target,     intent(inout) :: diag
  type(mixedlayer_restrat_CS), pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                  for this module
  character(len=128) :: version = '$Id: GOLD_mixed_layer_restrat.F90,v 13.0.2.3.2.10 2010/08/31 21:25:28 rwh Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "GOLD_mixed_layer_restrat"  ! This module's name.
  character(len=48)  :: flux_units

  if (associated(CS)) then
    call GOLD_error(WARNING, "mixedlayer_restrat_init called with an "// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  if (G%Boussinesq) then ; flux_units = "meter3 second-1"
  else ; flux_units = "kilogram second-1" ; endif

  CS%diag => diag
  CS%bulkmixedlayer = .false.
  CS%nkml = 1 ; CS%ml_restrat_coef = 0.0 ; CS%old_restrat = .true.
  call read_param(param_file,"BULKMIXEDLAYER",CS%bulkmixedlayer)
  if (CS%bulkmixedlayer) then
    call read_param(param_file,"NKML",CS%nkml)
    call read_param(param_file,"OLD_RESTRAT_PARAM",CS%old_restrat)
    if (CS%old_restrat) then
      call read_param(param_file,"ML_RESTRAT_COEF",CS%ml_restrat_coef)
    else
      call read_param(param_file,"FOX_KEMPER_ML_RESTRAT_COEF",CS%ml_restrat_coef)
    endif  
  endif

  CS%id_uhml = register_diag_field('ocean_model', 'uhml', G%axesul, Time, &
      'Zonal Thickness Flux to Restratify Mixed Layer', flux_units)
  CS%id_vhml = register_diag_field('ocean_model', 'vhml', G%axesvl, Time, &
      'Meridional Thickness Flux to Restratify Mixed Layer', flux_units)
  CS%id_urestrat_time = register_diag_field('ocean_model', 'MLu_restrat_time', G%axesu1, Time, &
      'Mixed Layer Zonal Restratification Timescale', 'second')
  CS%id_vrestrat_time = register_diag_field('ocean_model', 'MLv_restrat_time', G%axesu1, Time, &
      'Mixed Layer Meridional Restratification Timescale', 'second')

  ! Write all relevant parameters to the model log.
  call log_version(param_file, mod, version, tagname, "")
  call log_param(param_file, mod, "BULKMIXEDLAYER", CS%bulkmixedlayer, &
                 "If defined, use a refined Kraus-Turner-like bulk mixed \n"//&
                 "layer with transitional buffer layers.")
  if (CS%bulkmixedlayer) then
    call log_param(param_file, mod, "NKML", CS%nkml, &
                 "The number of sublayers within the mixed layer if \n"//&
                 "BULKMIXEDLAYER is true.", units="nondim", default=1)
    call log_param(param_file, mod, "OLD_RESTRAT_PARAM", CS%old_restrat, &
                 "If true, use the older version of the mixed layer \n"//&
                 "restratification code that predates the Fox-Kemper \n"//&
                 "et al. parameterization. ###THIS DEFAULT SHOULD BE CHANGED.###", &
                 default=.true.)
    if (CS%old_restrat) then
      call log_param(param_file, mod, "ML_RESTRAT_COEF", CS%ml_restrat_coef, &
                 "A coefficient in s (perhaps OMEGA^-1) relating the \n"//&
                 "mixed layer restratification to the horizontal density \n"//&
                 "gradients. This is only used with OLD_RESTRAT_PARAM defined.", &
                 units="s", default=0.0)
    else
      call log_param(param_file, mod, "FOX_KEMPER_ML_RESTRAT_COEF", &
                     CS%ml_restrat_coef, &
                 "A nondimensional coefficient that is proportional to \n"//&
                 "the ratio of the deformation radius to the dominant \n"//&
                 "lengthscale of the submesoscale mixed layer \n"//&
                 "instabilities, times the minimum of the ratio of the \n"//&
                 "mesoscale eddy kinetic energy to the large-scale \n"//&
                 "geostrophic kinetic energy or 1 plus the square of the \n"//&
                 "grid spacing over the deformation radius, as detailed \n"//&
                 "by Fox-Kemper et al. (2010)", units="nondim", default=0.0)
    endif
  endif

end subroutine mixedlayer_restrat_init

end module GOLD_mixed_layer_restrat
