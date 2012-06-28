module GOLD_entrain_H2000
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
!*  By Robert Hallberg, September 1997 - July 2000                     *
!*                                                                     *
!*    This file contains the subroutines that implement diapycnal      *
!*  mixing and advection in isopycnal layers.  The main subroutine,    *
!*  calculate_entrainment, returns the entrainment by each layer       *
!*  across the interfaces above and below it.  These are calculated    *
!*  subject to the constraints that no layers can be driven to neg-    *
!*  ative thickness and that the each layer maintains its target       *
!*  density, using the scheme described in Hallberg (MWR 2000). There  *
!*  may or may not be a bulk mixed layer above the isopycnal layers.   *
!*                                                                     *
!*    The dual-stream entrainment scheme of MacDougall and Dewar       *
!*  (JPO 1997) is used for combined diapycnal advection and diffusion, *
!*  modified as described in Hallberg (MWR 2000) to be solved          *
!*  implicitly in time.  Any profile of diffusivities may be used.     *
!*  Diapycnal advection is fundamentally the residual of diapycnal     *
!*  diffusion, so the fully implicit upwind differencing scheme that   *
!*  is used is entirely appropriate.  The downward buoyancy flux in    *
!*  each layer is determined from an implicit calculation based on     *
!*  the previously calculated flux of the layer above and an estimated *
!*  flux in the layer below.  This flux is subject to the following    *
!*  conditions:  (1) the flux in the top and bottom layers are set by  *
!*  the boundary conditions, and (2) no layer may be driven below an   *
!*  Angstrom thickness.  If there is a bulk mixed layer, the buffer    *
!*  layer is treated as a fixed density layer with vanishingly small   *
!*  diffusivity.                                                       *
!*                                                                     *
!*    In addition, the model may adjust the fluxes to drive the layer  *
!*  densities (sigma 2?) back toward their target values.              *
!*                                                                     *
!*    The model now uses a dynamic determination of the number of      *
!*  iterations each column requires, both in the inner iterations      *
!*  (with the shears assumed fixed) and the outer iterations on the    *
!*  velocity shears.  In a brief 1-degree global test with a 10%       *
!*  tolerance for the discrepancies between the values of F used to    *
!*  set the shears and the value after the inner iteration, going from *
!*  a limit of 50 iterations to 20 saves about 3% of the total number  *
!*  of iterations, while going down to 10 saves a further 3%.  It is   *
!*  also possible to same roughly 5% of the total number of iterations *
!*  by using the value of F from the last time step as the starting    *
!*  guess, although this has not been implemented because it was felt  *
!*  to be too small a benefit for its impact on the code flow and      *
!*  memory use.                                                        *
!*                                                                     *
!*    Calculate_Rino_flux calculates a measure of the entrainment      *
!*  rate due to the Richardson number dependent mixing.  This is an    *
!*  implicit calculation for the layer, but ignores the effects of     *
!*  entrainment by neighboring layers.  These effects are accounted    *
!*  for in its calling subroutine, calculate_entrainment.              *
!*                                                                     *
!*    Determine_Kd is where the user should change the code to         *
!*  alter the specification of the non-Richardson number-dependent     *
!*  diapycnal diffusivity.                                             *
!*                                                                     *
!*  Macros written all in capital letters are defined in GOLD_memory.h.*
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, buoy, T, S, Rml, ea, eb, etc.         *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use GOLD_diag_mediator, only : diag_ptrs, time_type
use GOLD_error_handler, only : GOLD_error, is_root_pe, FATAL, WARNING, NOTE
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_grid, only : ocean_grid_type
use GOLD_variables, only : forcing, thermo_var_ptrs
use GOLD_EOS, only : calculate_density, calculate_density_derivs

implicit none ; private

#include <GOLD_memory.h>

public Entrainment_H2000, entrain_H2000_init, entrain_H2000_end

type, public :: entrain_H2000_CS ; private
  logical :: bulkmixedlayer  ! If true, a refined bulk mixed layer is used with
                             ! tv%nk_Rml variable density mixed & buffer layers.
  logical :: correct_density ! If true, the layer densities are restored toward
                             ! their target variables by the diapycnal mixing.
  logical :: RiNo_mix        ! If true, use Richardson number dependent mixing.
  real    :: Kd              ! The interior diapycnal diffusivity in m2 s-1.
  real    :: Kv              ! The interior vertical viscosity in m2 s-1.
  integer :: max_ent_it      ! The maximum number of iterations that may be
                             ! used to calculate the diapycnal entrainment.
  real    :: RiNo_crit       ! The critical shear Richardson number for
                             ! shear-entrainment. The theoretical value is 1.
  real    :: Shearmix_rate   ! A nondimensional rate scale for shear-driven
                             ! entrainment.  The original value is 0.1, but
                             ! maybe it should be lower.
  real    :: RiNo_crit_eq    ! The critical shear Richardson number for
                             ! shear-entrainment that is used near the equator.
  real    :: Shearmix_rate_eq  ! The value of Shearmix_rate that is used at
                             ! the equator.
  real    :: Shearmix_lat_eq ! The latitude for transitions from the equatorial
                             ! values of shear mixing to midlatitude values,
                             ! in degrees.
  integer :: max_RiNo_it     ! The maximum number of iterations that may be used
                             ! to estimate the Richardson number dependent mixing.
  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
                             ! ocean diagnostic fields.
  integer :: id_Kd = -1
end type entrain_H2000_CS

contains

subroutine Entrainment_H2000(u, v, h, tv, fluxes, dt, G, CS, ea, eb, kb, &
                             Kd_Lay, Kd_int)
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(in)  :: u
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(in)  :: v
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(in)  :: h
  type(thermo_var_ptrs),               intent(in)  :: tv
  type(forcing),                       intent(in)  :: fluxes
  real,                                intent(in)  :: dt
  type(ocean_grid_type),               intent(in)  :: G
  type(entrain_H2000_CS),              pointer     :: CS
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(out) :: ea, eb
  integer, dimension(NXMEM_,NYMEM_),   intent(out) :: kb
  real, dimension(NXMEM_,NYMEM_,NZ_),   optional, intent(in) :: Kd_Lay
  real, dimension(NXMEM_,NYMEM_,NZp1_), optional, intent(in) :: Kd_int
!   This subroutine calculates ea and eb, the rates at which a layer
! entrains from the layers above and below.  The entrainment rates
! are proportional to the buoyancy flux in a layer and inversely
! proportional to the density differences between layers.  The
! scheme that is used here is described in detail in Hallberg, Mon.
! Wea. Rev. 2000.

! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      fluxes - A structure of surface fluxes that may be used.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      dt - The time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 entrain_H2000_init.
!  (out)     ea - The amount of fluid entrained from the layer above within
!                 this time step, in the same units as h, m or kg m-2.
!  (out)     eb - The amount of fluid entrained from the layer below within
!                 this time step, in the same units as h, m or kg m-2.
!  (out)     kb - The index of the lightest layer denser than the buffer layer.
! At least one of the two arguments must be present.
!  (in,opt)  Kd_Lay - The diapycnal diffusivity of layers, in m2 s-1.
!  (in,opt)  Kd_int - The diapycnal diffusivity of interfaces, in m2 s-1.

! In the comments below, H is used for the units of h, m or kg m-2.
  real, dimension(SZ1_(h),SZ3_(h)) :: &
    Kd, &   ! The diapycnal diffusivity of layers, converted to units of
            ! H2 s-1, that is m2 s-1 or kg2 m-4 s-1.
    F, &    ! The density flux through a layer within a time step divided by the
            ! density difference across the interface below the layer, in H.
    maxF, & ! maxF is the maximum value of F that will not deplete all of the
            ! layers above or below a layer within a timestep, in H.
    minF, & ! minF is the minimum flux that should be expected in the absence of
            ! interactions between layers, in H.
    Fprev, &! The previous estimate of F, in H.
    dFdfm, &! The partial derivative of F with respect to changes in F of the
            ! neighboring layers.  Nondimensional.
    fm_ent,&! The value of fm that was used to calculate shears, in H.

    h_guess, & ! An estimate of the layer thicknesses after entrainment, but
            ! before the entrainments are adjusted to drive the layer
            ! densities toward their target values, in H.

    FRi0, & ! The density flux through a layer due to shear Richardson number
            ! dependent entrainment, in the absence of entrainment by the
            ! neighboring layers, divided by the density difference across the
            ! interface below the layer and I2p2gkp1_gk, in H.
    R0, &   ! Potential density relative to the surface, in kg m-3.
    F_Ri_est, & ! The value of F that was last used to calculate FRi0, in H.
    F_best, & ! The most-nearly self consistent estimate of F that has yet
            ! been found, in H.  Consistency is determined by the mismatch
            ! between the value of F used to calculate FRi0 and the value that
            ! is returned from the inner iteration.
    maxEnt  ! maxEnt is the maximum value of entrainment from below (with
            ! compensating entrainment from above to keep the layer density
            ! from changing) that will not deplete all of the layers above or
            ! below a layer within a timestep, in H.
  real :: FRi           ! The Richardson number dependent mixing, in H.
  real :: dFRi_dfm      ! The partial derivative of FRi with fm (divided by I2p2...

  real :: fm, fr, fk    ! Work variables with units of H, H, and H2.
  real :: F_est, F1     ! Various temporary estimates of F, in H.

  real :: b1(SZ1_(h))         ! b1 and c1 are variables used by the
  real :: c1(SZ1_(h),SZ3_(h)) ! tridiagonal solver.

  real :: htot(SZ1_(h)) ! The total thickness above or below a layer in H.
  real :: mFkb(SZ1_(h)) ! The total thickness in the mixed and buffer layers
                        ! times gml_gkbp1, in H.
  real :: Rcv(SZ1_(h))  ! Value of the coordinate variable (potential density)
                        ! based on the simulated T and S and P_Ref, kg m-3.
  real :: pres(SZ1_(h)) ! Reference pressure (tv%P_Ref) in Pa.
  real :: p_0(SZ1_(h))  ! An array of 0 pressures.

  real, dimension(SZ1_(h),SZ3_(h)) :: &
    I2p2gkbp1_gml, &    ! 1 / (2 + 2 * g_kb+1 / g_ml). Nondimensional.
    F_start_ml, &       ! The input value of F at which a buffer or mixed
                        ! layer starts to be entrained, in H.
    Fout_start          ! The final value of F when a buffer or mixed
                        ! layer starts to be entrained, in H.
  real, dimension(SZ3_(h)) :: &
    gkp1_gk, &    ! The reduced gravity of the interface below a layer divided
                  ! by the reduced gravity of the interface above it. Nondim.
    gk_gkp1       ! The reduced gravity of the interface above a layer divided
                  ! by the reduced gravity of the interface below it. Nondim.
  real, dimension(SZ1_(h)) :: &
    gkbp1_gkb, &  ! gkbp1_gkb is the inverse of the ratio of the reduced gravity
                  ! at the interface between the buffer layer and the next
                  ! denser layer and the reduced gravity of the next deeper
                  ! interface. Nondimensional.
    gkb_gkbp1     ! g_kb / g_kb+1. Nondimensional.
  real, dimension(SZ1_(h),SZ3_(h)) :: &
    gkbp1_gml, &  ! gkbp1_gml is the inverse of the ratio of the density
                  ! difference between the mixed layer and the layer below the
                  ! buffer layer to the density difference across the next
                  ! deeper interface. Nondimensional.
    gml_gkbp1     ! 1 / gkbp1_gml. Nondimensional.
  real :: gratsdt(SZ3_(h))       ! 2*dt*(2 + g_k+1 / g_k + g_k / g_k+1). In s.
  real :: I2p2gkp1_gk(SZ3_(h))   ! 1 / (2 + 2 * g_k+1 / g_k). Nondimensional.
  real :: gratskbdt(SZ1_(h))     ! 2*dt*(2 + g_kb+1/g_kb + g_kb/g_kb+1), in s.
  real :: I2p2gkbp1_gkb(SZ1_(h)) ! 1 / (2 + 2 * g_kb+1 / g_kb). Nondimensional.
  real :: I2gprime(SZ1_(h),SZ3_(h)) ! A work-space for Calculate_Rino_flux that
                                    ! is retained between iterations.

  real :: I2p2gk_here         ! Either I2p2gkbp1_gml or I2p2gkbp1_gkb, nondimensional.
  logical :: flux_from_Ri     ! If true, the Richardson number flux is used to determine
                              ! the fluxes.
  real :: F_cor               ! A correction to the amount of F that is used to
                              ! entrain from the layer above, in H.
  real :: F_cor_max           ! The maximum permitted value of F_cor to avoid
                              ! giving the current layer negative thickness, in H.
  real :: ent_min_rem         ! The remaining entrainment from above that must be
                              ! found to keep the current layer above its minimum
                              ! thickness, in H.

  real :: glkdt_gkp1(SZ3_(h)) ! The arithmetic mean g' times the time step over
                          ! g' of the interface below a layer, in s.
  real :: g_R0            ! g_R0 is g/Rho in m4 kg-1 s-2.
  real :: Ig1             ! 1.0/g_prime(1), in s2 m-1.
  real :: z, Idenom       ! Temporary variables in units of H and H-1.
  real :: eps, tmp        ! Nondimensional temporary variables.
  real :: a(SZ3_(h)), a_0(SZ3_(h)) ! Nondimensional temporary variables.
  real :: I_Drho          ! Temporary variable in m3 kg-1.
  real :: Kd_here         ! The effective diapycnal diffusivity, in H2 s-1.
  real :: Angstrom        !   The minimum layer thickness, in H.

  ! The following group of variables are used when iterating to determine a
  ! self-consistent Richardson-number-dependent mixing.
  real, dimension(SZ1_(h)) :: &
    sd, &   ! sd is the summed difference between entrainment (F) that is
            ! calculated and the estimate of F (F_Ri_est) that was used to set
            ! values related to the Richardson number dependent mixing.  These
            ! must agree for the solution to be self-consistent, in which case
            ! sd is a sum of zeros.  sd has units of H.
    sad, &  ! sad is the summed absolute difference between F and F_Ri_est.
            ! sad is minimized to find self-consistent profiles.  sad is in H.
    sad_best, & ! sad_best is the best value of sad found so far, in H.
    sumF    ! sumF is the sum of F in H.  A good estimate has a small value of
            ! sad compared with sumF.
  real, dimension(3,SZ1_(h),SZ3_(h)) :: F_Ri_est_prev ! Three previous estimates
                              ! of entrainment that give estimates that are
                              ! systematically too large (1), too small (2), or
                              ! of ambiguous sign (3).
  integer, dimension(3,SZ1_(h)) :: prev_set ! The iteration number of a previous
                              ! entrainment estimate that is stored in
                              ! F_Ri_est_prev, or 0 if none has been stored.
  real, dimension(3,SZ1_(h)) :: sd_prev, sad_prev ! The values of sd and sad for
                              ! the entrainment profiles stored in F_Ri_est_prev.
  real, dimension(SZ1_(h)) :: w1, w3 ! w1, w3 are weights to be given to the
                                     ! most recent estimate of entrainment and
                                     ! the returned calculation of entrainment
                                     ! in finding the next estimate.
  real, dimension(3,SZ1_(h)) :: w2   ! w2 are weights of the 3 previous estimates
                                     ! of entrainment for the next estimate.
  integer, dimension(SZ1_(h)) :: m_set ! The index (1-3) of F_Ri_est_prev in which
                                       ! to store the current F_Ri_est profile.
  integer, dimension(SZ1_(h)) :: set_best ! The iteration at which the best-yet
                                          ! profile was found.

  logical :: do_i(SZ1_(h)), did_i(SZ1_(h)), reiterate, do_i_Ri(SZ1_(h)), inject_best(SZ1_(h))
  logical :: use_EOS      !   If true, density is calculated from temperature
                          ! and salnity using an equation of state.
  logical :: correct_density ! If true, the layer densities are restored toward
                             ! their target variables by the diapycnal mixing.
  integer :: it, it2, i, j, k, is, ie, js, je, nz, max_it2, K2, k3, kmb
  integer :: m1, m2
  real :: I_s0, I_s1, I_s2, I_s3
  integer :: kb_min       ! The minimum value of kb in the current j-row.
  integer :: kb_min_act   ! The minimum active value of kb in the current j-row.
  integer :: is1, ie1     ! The minimum and maximum active values of i in the current j-row.
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not. associated(CS)) call GOLD_error(FATAL, &
      "GOLD_entrain_H2000: Module must be initialized before it is used.")

  if (.not.(present(Kd_Lay) .or. present(Kd_int))) call GOLD_error(FATAL, &
      "GOLD_entrain_H2000: Either Kd_Lay or Kd_int must be present in call.")

  use_EOS = associated(tv%eqn_of_state)
  correct_density = (use_EOS .and. CS%correct_density)
  g_R0 = G%g_Earth/G%Rho0
  Angstrom = G%Angstrom

  if (CS%bulkmixedlayer) then
    kmb = tv%nk_Rml ; K2 = kmb + 1 ; kb_min = K2
  else
    kmb = 0 ; K2 = 2 ; kb_min = K2 ; kb(:,:) = 1
    F_start_ml(:,:) = 0.0  ! Initialized only to avoid comparison with a NaN.
  endif

  if (correct_density .or. CS%bulkmixedlayer) pres(:) = tv%P_Ref
  maxEnt(:,:) = 0.0
  FRi = 0.0

! gkp1_gk(:) = 0.0 ; gk_gkp1(:) = 0.0
! gkbp1_gkb(:) = 0.0 ; gkb_gkbp1(:) = 0.0
! gkbp1_gml(:,:) = 0.0 ; gml_gkbp1(:,:) = 0.0

  !   Set up a number of ratios of density differences between
  ! layers for convenience in later calculations.
  do k=2,nz-1
    gk_gkp1(k) = G%g_prime(k) / G%g_prime(k+1)
  enddo
  gk_gkp1(nz) = 2.0
  do k=2,nz
    gkp1_gk(k) = 1.0 / gk_gkp1(k)
    glkdt_gkp1(k) = 0.5*dt*(1.0+gk_gkp1(k))
    gratsdt(k) = 2.0*dt*(2.0+(gkp1_gk(k)+gk_gkp1(k)))
    I2p2gkp1_gk(k) = 0.5/(1.0+gkp1_gk(k))
  enddo
  Ig1 = 1.0 / G%g_prime(2)

  if (CS%bulkmixedlayer .and. use_EOS) then
    do i=is,ie ; p_0(i) = 0.0 ; enddo
  endif

  do j=js,je
    if (present(Kd_Lay)) then
      do k=1,nz ; do i=is,ie
        Kd(i,k) = G%m_to_H**2 * Kd_Lay(i,j,k)
      enddo ; enddo
    else ! Kd_int must be present, or there already would have been an error.
      do k=1,nz ; do i=is,i
        Kd(i,k) = 0.5 * G%m_to_H**2 * (Kd_int(i,j,k)+Kd_int(i,j,k+1))
      enddo ; enddo
    endif

    do i=is,ie ; do_i_Ri(i) = (G%hmask(i,j) > 0.5) ; enddo

    if (CS%bulkmixedlayer) then

      if (use_EOS) then ; do k=1,nz
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_0, R0(:,k), &
                               is, ie-is+1, tv%eqn_of_state)
      enddo ; endif

      ! Determine kb, the index of the layer below the deepest buffer layer.
      do i=is,ie
        !   Determine the next denser layer than the buffer layer in the
        ! coordinate density (sigma-2).
        do k=kmb+1,nz-1 ; if (tv%Rml(i,j,kmb) <= G%Rlay(k)) exit ; enddo
        kb(i,j) = k

        if (use_EOS) then
          !   Backtrack, in case there are massive layers above that are stable
          ! in sigma-0.
          do k=kb(i,j)-1,kmb+1,-1
            if (R0(i,kmb) > R0(i,k)) exit
            if (h(i,j,k)>2.0*Angstrom) kb(i,j) = k
          enddo
        endif
      enddo


      kb_min = nz
      eps = 0.1
      do i=is,ie
        if (kb(i,j) <= nz-1) then
!   Set up appropriately limited ratios of the reduced gravities of the
! interfaces above and below the buffer layer and the next denser layer.
          k = kb(i,j)

          I_Drho = g_R0 / G%g_prime(k+1)
          do k3=1,kmb
            a(k3) = (G%Rlay(k) - tv%Rml(i,j,k3)) * I_Drho
          enddo
          if (use_EOS .and. (a(kmb) < 2.0*eps*gk_gkp1(k))) then
!   If the buffer layer nearly matches the density of the layer below in the
! coordinate variable (sigma-2), use the sigma-0-based density ratio if it is
! greater (and stable).
            if ((R0(i,k) > R0(i,kmb)) .and. &
                (R0(i,k+1) > R0(i,k))) then
              I_Drho = 1.0 / (R0(i,k+1)-R0(i,k))
              a_0(kmb) = min((R0(i,k)-R0(i,kmb)) * I_Drho, gk_gkp1(k))
              if (a_0(kmb) > a(kmb)) then
                do k3=1,kmb-1
                  a_0(k3) = a_0(kmb) + (R0(i,kmb)-R0(i,k3)) * I_Drho
                enddo
                if (a(kmb) <= eps*gk_gkp1(k)) then
                  do k3=1,kmb ; a(k3) = a_0(k3) ; enddo
                else
! Alternative...  tmp = 0.5*(1.0 - cos(PI*(a(K2)/(eps*gk_gkp1(k)) - 1.0)) )
                  tmp = a(kmb)/(eps*gk_gkp1(k)) - 1.0
                  do k3=1,kmb ; a(k3) = tmp*a(k3) + (1.0-tmp)*a_0(k3) ; enddo
                endif
              endif
            endif
          endif

          gkb_gkbp1(i) = MAX(a(kmb),1e-5)
          gkbp1_gkb(i) = 1.0 / gkb_gkbp1(i)
          I2p2gkbp1_gkb(i) = 0.5/(1.0+gkbp1_gkb(i))
          gratskbdt(i) = 2.0*dt*(2.0+(gkbp1_gkb(i)+gkb_gkbp1(i)))

          do k3=kmb-1,1,-1
!           gml_gkbp1(i,k3) = MAX(a(k3),1e-5)
            ! Deliberately treat convective instabilities of the upper mixed 
            ! and buffer layers with respect to the deepest buffer layer as
            ! though they don't exist.  They will be eliminated by the upcoming
            ! call to the mixedlayer code anyway.
            gml_gkbp1(i,k3) = MAX(a(k3),gkb_gkbp1(i))
            gkbp1_gml(i,k3) = 1.0 / gml_gkbp1(i,k3)
            I2p2gkbp1_gml(i,k3) = 0.5/(1.0+gkbp1_gml(i,k3))
          enddo

          F_start_ml(i,kmb) = 0.0 ; Fout_start(i,kmb) = 0.0
          F_start_ml(i,kmb-1) = (h(i,j,kmb)-Angstrom)*gkb_gkbp1(i)
          Fout_start(i,kmb-1) = F_start_ml(i,kmb-1)
          do k3=kmb-1,2,-1
            F_start_ml(i,k3-1) = F_start_ml(i,k3) + &
               (I2p2gkbp1_gkb(i) * 2.0*(1.0+gkbp1_gml(i,k3))) * &
               (h(i,j,k3)-Angstrom)*gml_gkbp1(i,k3)
            Fout_start(i,k3-1) = Fout_start(i,k3) + &
               (h(i,j,k3)-Angstrom)*gml_gkbp1(i,k3)
          enddo

          if (kb_min>k) kb_min = k
        endif
        maxF(i,1) = 0.0
      enddo
    else
! Without a bulk mixed layer, surface fluxes are applied in this
! subroutine.  (Otherwise, they are handled in mixedlayer.)
!   Initially the maximum flux in layer zero is given by the surface
! buoyancy flux.  It will later be limited if the surface flux is
! too large.  Here buoy is the surface buoyancy flux.
      if (ASSOCIATED(fluxes%buoy)) then
        do i=is,ie
          maxF(i,1) = (dt*fluxes%buoy(i,j)*Ig1) * G%m_to_H
        enddo
      else
        if (ASSOCIATED(fluxes%liq_precip) .or. ASSOCIATED(fluxes%evap) .or. &
            ASSOCIATED(fluxes%sens) .or. ASSOCIATED(fluxes%sw)) then
          if (is_root_pe()) call GOLD_error(NOTE, "Calculate_Entrainment: &
                &The code to handle evaporation and precipitation without &
                &a bulk mixed layer has not been implemented.")
          if (is_root_pe()) call GOLD_error(FATAL, &
               "Either define BULKMIXEDLAYER in GOLD_input or use fluxes%buoy &
               &and a linear equation of state to drive the model.")
        else
          do i=is,ie
            maxF(i,1) = 0.0
          enddo
        endif
      endif
    endif

! The following code (1) calculates the maximum flux, maxF, and
! several related quantities as appropriate.
    if (CS%bulkmixedlayer) then
      do i=is,ie
        htot(i) = h(i,j,kmb)
        mFkb(i) = gkb_gkbp1(i)*(h(i,j,kmb) - Angstrom)
      enddo
      do k=1,kmb-1 ; do i=is,ie
        htot(i) = htot(i) + h(i,j,k)
        mFkb(i) = mFkb(i) + gml_gkbp1(i,k)*(h(i,j,k) - Angstrom)
      enddo ; enddo
      do i=is,ie ; if (kb(i,j) < nz) then
        maxF(i,kb(i,j)) = mFkb(i)
        maxEnt(i,kb(i,j)) = mFkb(i)
      endif ; enddo
    else
      do i=is,ie
        htot(i) = h(i,j,1) - Angstrom
      enddo
    endif
    do k=kb_min,nz-1 ; do i=is,ie
      if (k > kb(i,j)) then
        maxF(i,k) = gk_gkp1(k)*(maxF(i,k-1) + htot(i))
!        if (CS%RiNo_mix) maxEnt(i,k) = gk_gkp1(k)*htot(i)
!  Alternative?
        maxEnt(i,k) = gk_gkp1(k)*(maxEnt(i,k-1) + htot(i))
        htot(i) = htot(i) + (h(i,j,k) - Angstrom)
      endif
    enddo ; enddo
    do i=is,ie
      maxF(i,nz) = 0.0
      if (.not.CS%bulkmixedlayer) &
        maxF(i,nz) = MIN(0.0, gk_gkp1(nz)*(maxF(i,nz-1) + htot(i)))
      htot(i) = h(i,j,nz) - Angstrom
    enddo
    do k=nz-1,kb_min,-1 ; do i=is,ie
      if (k>=kb(i,j)) then
        maxF(i,k) = MIN(maxF(i,k),gkp1_gk(k+1)*maxF(i,k+1) + htot(i))
 !       if (CS%RiNo_mix) maxEnt(i,k) = MIN(maxEnt(i,k),htot(i))
 !  Alternative?
        maxEnt(i,k) = MIN(maxEnt(i,k),gkp1_gk(k+1)*maxEnt(i,k+1) + htot(i))
        htot(i) = htot(i) + (h(i,j,k) - Angstrom)
      endif
    enddo ; enddo
    if (.not.CS%bulkmixedlayer) then
      do i=is,ie
        maxF(i,1) = MIN(maxF(i,1),gkp1_gk(2)*maxF(i,2) + htot(i))
      enddo
    endif

!   Estimate the Richardson number dependent entrainment rate.
    if (CS%RiNo_mix) &
      call Calculate_Rino_flux(u,v,h,tv,kb,j,dt,G,CS, gk_gkp1, gkp1_gk, &
                               gkb_gkbp1, gkbp1_gkb, gml_gkbp1, gkbp1_gml,&
                               .true.,F,maxEnt,do_i_Ri,  gratsdt, gratskbdt, &
                               I2gprime, FRi0, fm_ent)


!   The following code provides an initial estimate of the flux in
! each layer, F.  The initial guess for the layer diffusive flux is
! the smaller of a forward discretization or the maximum diffusive
! flux starting from zero thickness in one time step without
! considering adjacent layers.
    do i=is,ie
      F(i,1) = maxF(i,1)
      F(i,nz) = maxF(i,nz)
    enddo
    do k=nz-1,K2,-1
      do i=is,ie
        if ((k>=kb(i,j)) .and. (do_i_Ri(i))) then
!   Here the layer flux is estimated, assuming no entrainment from
! the surrounding layers.  The estimate is a forward (steady) flux,
! limited by the maximum flux for a layer starting with zero
! thickness.  This is often a good guess and leads to few iterations.
          if (k>kb(i,j)) then
            F1 = MIN(sqrt((dt*gk_gkp1(k))*Kd(i,k)), &
                     glkdt_gkp1(k)*(Kd(i,k)/h(i,j,k)))
            if (CS%RiNo_mix) F1 = MAX(I2p2gkp1_gk(k)*FRi0(i,k), F1)
          else
            F1 = MIN(sqrt((dt*gkb_gkbp1(i))*Kd(i,k)), &
                     (0.5*dt*(gkb_gkbp1(i)+1.0))*(Kd(i,k)/h(i,j,k)))
            if (CS%RiNo_mix) F1 = MAX(I2p2gkbp1_gkb(i)*FRi0(i,k), F1)
          endif
          F(i,k) = MIN(maxF(i,k), F1)

!   Calculate the minimum flux that can be expected if there is no entrainment
! from the neighboring layers.  The 0.9 is used to give used to give a 10%
! known error tolerance.
          fm = h(i,j,k)
          if (k>kb(i,j)) then
            fk = Kd(i,k) * gratsdt(k)
            F_est = I2p2gkp1_gk(k) * fk / (fm + sqrt(fm*fm + fk))
          else
            fk = Kd(i,k) * gratskbdt(i)
            F_est = I2p2gkbp1_gkb(i) * fk / (fm + sqrt(fm*fm + fk))
          endif

          minF(i,k) = MIN(maxF(i,k), 0.9*F_est)
          if (k==kb(i,j)) minF(i,k) = 0.0 ! BACKWARD COMPATIBILITY - DELETE LATER?
        else
          F(i,k) = 0.0
        endif
      enddo ! end of i loop
    enddo ! end of k loop


    max_it2 = 1
    if (CS%RiNo_mix) max_it2 = CS%max_RiNo_it
    do it2=1,max_it2
! This is the start of the section that might be iterated.
      if (CS%bulkmixedlayer) then
        kb_min_act = nz
        do i=is,ie
          do_i(i) = do_i_Ri(i)
          if (do_i_Ri(i) .and. (kb(i,j) < kb_min_act)) kb_min_act = kb(i,j)
        enddo
      else
        kb_min_act = kb_min
        do i=is,ie ; do_i(i) = do_i_Ri(i) ; enddo
      endif
      is1 = ie+1 ; ie1 = is-1
      do i=is,ie ; if (do_i_Ri(i)) then
        is1 = i ; exit
      endif ; enddo
      do i=ie,is,-1 ; if (do_i_Ri(i)) then
        ie1 = i ; exit
      endif ; enddo

      do it=0,CS%max_ent_it-1
        do i=is1,ie1 ; if (do_i(i)) then
          if (.not.CS%bulkmixedlayer) F(i,1) = MIN(F(i,1),maxF(i,1))
          b1(i) = 1.0
        endif ; enddo ! end of i loop

        do k=kb_min_act,nz-1 ; do i=is1,ie1 ; if (do_i(i)) then
          if (k >= kb(i,j)) then
! Calculate the flux in layer k.
            if (k==kb(i,j)) then
! Calculate the flux in the layer beneath the buffer layer.
              I2p2gk_here = I2p2gkbp1_gkb(i)
              fm = -h(i,j,k) + gkp1_gk(k+1)*F(i,k+1)
              fk = gratskbdt(i)*Kd(i,k)
            else
              Fprev(i,k) = F(i,k)
              I2p2gk_here = I2p2gkp1_gk(k)
              fm = (F(i,k-1) - h(i,j,k)) + gkp1_gk(k+1)*F(i,k+1)
              fk = gratsdt(k)*Kd(i,k)
            endif
            fr = sqrt(fm*fm + fk)

            if (CS%RiNo_mix) then
              if (fm <= fm_ent(i,k)) then
                FRi = I2p2gk_here*(FRi0(i,k) + 2.0*(fm + h(i,j,k)))
                dFRi_dfm = 2.0
              else
!                num = (FRi0(i,k) + 2.0*h(i,j,k)) * (0.25*FRi0(i,k) + 1.5*(h(i,j,k)+Angstrom) + fm_ent(i,k))
!                Idenom = 1.0/(0.25*FRi0(i,k) + 1.5*(h(i,j,k)+Angstrom) + 0.5*(fm_ent(i,k) + fm))
!                FRi = I2p2gk_here * (2.0*fm + num * Idenom)
!                dFRi_dfm = 2.0 - 0.5 * num * Idenom**2
                ! This sets the final thickness to be (z*(z+2fm'))/(z+fm')^2
                ! smaller than the thickness when fm = fm_ent, where fm' = fm-fm_ent. It limits
                ! the energy used by entrainment in this layer when fm is large
                ! to between 2 and 3 (for h >> fm_ent+FRi0) times the value at fm_ent.
                z = (fm_ent(i,k) + FRi0(i,k) + 3.0*h(i,j,k) + Angstrom)
                Idenom = 1.0 / (z + fm - fm_ent(i,k))
                FRi = I2p2gk_here * (2.0*fm + (FRi0(i,k) + 2.0*h(i,j,k)) *&
                           ((z * (z + 2.0*(fm - fm_ent(i,k)))) * Idenom**2) )
                dFRi_dfm = 2.0 - 2.0*(FRi0(i,k) + 2.0*h(i,j,k)) * &
                           (z*(fm - fm_ent(i,k)))*Idenom**3
              endif
            endif ! Otherwise FRi is still 0 from initialization.

            if (fm>=0) then
              F1 = I2p2gk_here * (fm+fr)
            else
              F1 = I2p2gk_here * (fk / (-1.0*fm+fr))
            endif
            F(i,k) = MIN(maxF(i,k),(MAX(FRi,F1)))

            if (FRi >= F1) then ; F1 = FRi ; flux_from_Ri = .true.
            else ; flux_from_Ri = .false. ; endif

            if (k == kb(i,j)) then ; if (F1 >= F_start_ml(i,kmb-1)) then
            ! In this case, gratskbdt should probably have been altered above, but
            ! has not been.  This is equivalent to using a smaller diffusivity than
            ! indicated at the base of the mixed layer, but as the mixed layer will
            ! be entraining next (probably quite vigorously), this is an omission I
            ! am willing to make.  RWH

            ! In this case, the layer entirely entrains one or more of the buffer
            ! and mixed layer, and F1 is increased as the layer entrains mass
            ! more slowly from the lighter layers.  This is effectively repeatedly
            ! solving (4.4) of Hallberg (MWR, 2000) with changing values of
            ! gamma/delta_rho.
              if (F1 >= F_start_ml(i,1)) then
                I2p2gk_here = I2p2gkbp1_gml(i,1)
                F1 = I2p2gkbp1_gml(i,1) * 2.0*(1.0+gkbp1_gkb(i)) * &
                      (F1-F_start_ml(i,1)) + Fout_start(i,1)
              else
                do k3=kmb-1,2,-1 ; if (F1 < F_start_ml(i,k3-1)) then
                  I2p2gk_here = I2p2gkbp1_gml(i,k3)
                  F1 = I2p2gkbp1_gml(i,k3) * 2.0*(1.0+gkbp1_gkb(i)) * &
                        (F1-F_start_ml(i,k3)) + Fout_start(i,k3)
                  exit
                endif ; enddo
              endif
            endif ; endif
            F(i,k) = MIN(maxF(i,k),F1)

            if (F(i,k) >= maxF(i,k)) then
              dFdfm(i,k) = 0.0
            else if (.not.flux_from_Ri) then
              dFdfm(i,k) = I2p2gk_here * ((fr + fm) / fr)
            else
              dFdfm(i,k) = I2p2gk_here * dFRi_dfm
            endif

            if ((k > kb(i,j)) .and. (k > K2)) then
              ! This is part of a tridiagonal solver for the actual flux.
              c1(i,k) = dFdfm(i,k-1)*(gkp1_gk(k)*b1(i))
              b1(i) = 1.0 / (1.0 - c1(i,k)*dFdfm(i,k))
              F(i,k) = MIN(b1(i)*(F(i,k)-Fprev(i,k)) + Fprev(i,k), &
                              maxF(i,k))
              if (F(i,k) >= maxF(i,k)) dFdfm(i,k) = 0.0
            endif

          endif
        endif ; enddo ; enddo

        do k=nz-2,kb_min_act,-1 ; do i=is1,ie1
          if (do_i(i) .and. (k >= kb(i,j))) &
            F(i,k) = MIN((F(i,k)+c1(i,k+1)*(F(i,k+1)-Fprev(i,k+1))),maxF(i,k))
        enddo ; enddo

! Determine whether to do another iteration.
        if (it < CS%max_ent_it-1) then
          reiterate = .false.
          do i=is1,ie1
            did_i(i) = do_i(i) ; do_i(i) = .false.
          enddo
          do k=kb_min_act,nz-1 ; do i=is1,ie1
            if (did_i(i) .and. (k >= kb(i,j))) then
              if (F(i,k) < minF(i,k)) then
                F(i,k) = minF(i,k)
                do_i(i) = .true. ; reiterate = .true.
              elseif (k > kb(i,j)) then
                  if ((h(i,j,k) + ((1.0+gkp1_gk(k))*F(i,k) - &
                       (F(i,k-1) + gkp1_gk(k+1)*F(i,k+1)))) < 0.9*Angstrom) then
                    do_i(i) = .true. ; reiterate = .true.
                  endif
                else ! (k == kb(i,j))
! A more complicated test is required for the layer beneath the buffer layer,
! since its flux may be partially used to entrain directly from the mixed layer.
! Negative fluxes should not occur with the bulk mixed layer.
                  if (h(i,j,k) + ((F(i,k) + &
                      ea_kb(F(i,k),h,i,j,kmb,gkb_gkbp1,gkbp1_gkb,&
                            gml_gkbp1,gkbp1_gml,G)) - &
                      gkp1_gk(k+1)*F(i,k+1)) < 0.9*Angstrom) then
                    do_i(i) = .true. ; reiterate = .true.
                  endif
                endif
              endif
          enddo ; enddo
          if (.not.reiterate) exit
        endif ! end of if (it < CS%max_ent_it-1)
      enddo ! end of it loop
! This is the end of the section that might be iterated.


      if (it == (CS%max_ent_it)) then
!   Limit the flux so that the layer below is not depleted.
! This should only be applied to the last iteration.
        do i=is1,ie1
          if (F(i,nz-1) < 0.0) F(i,nz-1) = 0.0
        enddo
        do k=nz-2,kb_min_act,-1 ; do i=is1,ie1
          if (do_i(i) .and. (k>=kb(i,j))) then
            F(i,k) = MIN(MAX(minF(i,k),F(i,k)), (gkp1_gk(k+1)*F(i,k+1) + &
                 MAX(((F(i,k+1)-gkp1_gk(k+2)*F(i,k+2)) + &
                      (h(i,j,k+1) - Angstrom)), 0.5*(h(i,j,k+1)-Angstrom))))
          endif
        enddo ; enddo

!   Limit the flux so that the layer above is not depleted.
        do k=kb_min_act+1,nz-1 ; do i=is1,ie1 ; if (do_i(i)) then
          if ((.not.CS%bulkmixedlayer) .or. (k > kb(i,j)+1)) then
            F(i,k) = MIN(F(i,k), gk_gkp1(k)*( ((F(i,k-1) + &
                gkp1_gk(k-1)*F(i,k-1)) - F(i,k-2)) + (h(i,j,k-1) - Angstrom)))
          else if (k == kb(i,j)+1) then
            F(i,k) = MIN(F(i,k), gk_gkp1(k)*( (F(i,k-1) + &
                ea_kb(F(i,k-1),h,i,j,kmb,gkb_gkbp1,gkbp1_gkb,gml_gkbp1,&
                      gkbp1_gml,G)) + (h(i,j,k-1) - Angstrom)))
          endif
        endif ; enddo ; enddo
      endif

!   Determine which columns require further iteration on the Richardson number
! and the value of F to use to calculate the shears.
      if (it2 == 1) then
        reiterate = .true.
        do k=1,nz ; do i=is1,ie1
          F_Ri_est(i,k) = F(i,k) ; F_Ri_est_prev(1,i,k) = 0.0
          F_Ri_est_prev(2,i,k) = 0.0 ; F_Ri_est_prev(3,i,k) = 0.0
        enddo ; enddo
        do i=is1,ie1
          sd(i) = 0.0 ; sad(i) = 1.0*G%m_to_H ; sad_best(i) = 1.0e30*G%m_to_H
          set_best(i) = 0
          prev_set(1,i) = 0 ; prev_set(2,i) = 0 ; prev_set(3,i) = 0
          sad_prev(1,i) = 0.0 ; sad_prev(2,i) = 0.0 ; sad_prev(3,i) = 0.0
        enddo
      elseif (it2<max_it2) then
        reiterate = .false.
        do i=is1,ie1
          sd(i) = 0.0 ; sad(i) = 0.0 ;  sumF(i) = 0.0
          did_i(i) = do_i_Ri(i) ; do_i_Ri(i) = .false.
        enddo
        do k=kb_min_act,nz-1 ; do i=is1,ie1
          if ((k>=kb(i,j)) .and. did_i(i)) then
            if (.not.do_i_Ri(i)) then
              if (abs(F(i,k) - F_Ri_est(i,k)) > (0.01*h(i,j,k) + &
                  0.05*(max(F(i,k),0.0) + MAX(F_Ri_est(i,k),0.0)))) then
                do_i_Ri(i) = .true. ; reiterate = .true.
              endif
            endif
            sd(i) = sd(i) + (F(i,k) - F_Ri_est(i,k))
            sad(i) = sad(i) + abs(F(i,k) - F_Ri_est(i,k))
            sumF(i) = sumF(i) + max(F(i,k),0.0)
          endif
        enddo ; enddo

        do i=is1,ie1 ; if (do_i_Ri(i)) then
          if (sad(i) < sad_best(i)) then
            sad_best(i) = sad(i) ; set_best(i) = it2
          endif
          if (mod(it2-set_best(i),6)==5) then
            inject_best(i) = .true.
          else
            inject_best(i) = .false.
          endif
          m_set(i) = -1
          ! These are the default values - split the difference between the
          ! current input and output.
          w1(i)=0.5 ; w2(1,i)=0.0 ; w2(2,i)=0.0 ; w2(3,i)=0.0 ; w3(i)=0.5
          if (abs(sd(i)) > 0.2*sad(i)) then
            ! Replace one of the older estimates, and use the false position
            ! method (if appropriate) to estimate the next value.
            if (sd(i) < 0) then
              m1 = 1 ; m2 = 2
            else
              m1 = 2 ; m2 = 1
            endif
            w2(3,i) = 0.0
            if ((prev_set(m1,i) < 1) .or. (sad(i) < sad_prev(m1,i))) then
              ! This is better than a previous one - keep it.
              m_set(i) = m1 ; prev_set(m1,i) = it2
              sd_prev(m1,i) = sd(i) ; sad_prev(m1,i) = sad(i)
              if (prev_set(m2,i) > 0) then
                w1(i) = abs(sd(i) / (sd(i) - sd_prev(m2,i)))
                w2(m2,i) = 1.0 - w1(i) ; w2(m1,i) = 0.0 ; w3(i) = 0.0
              endif
            elseif (prev_set(m2,i) > 0) then
              ! With a worse estimate, try to improve the opposite extreme point
              ! rather than going for the true solution
              w1(i) = sad_prev(m2,i) / (sad(i) + sad_prev(m2,i))
              w2(m2,i) = 1.0 - w1(i) ; w2(m1,i) = 0.0 ; w3(i) = 0.0
            endif
          else
            if ((prev_set(3,i) < 1) .or. (sad(i) < sad_prev(3,i))) then
              m_set(i) = 3 ; prev_set(3,i) = it2
              sd_prev(3,i) = sd(i) ; sad_prev(3,i) = sad(i)
            endif
            if ((prev_set(1,i) > 0) .and. (prev_set(2,i) > 0) .and. &
                (prev_set(3,i) < it2) .and. (sad(i) > min(sad_prev(1,i),sad_prev(2,i)))) then
              ! This estimate is not so good as some others - try reverting to
              ! a weighted mean of the current and previous best estimates.
              I_s0 = 1.0/sad(i) ; I_s1 = 1.0/sad_prev(1,i)
              I_s2 = 1.0/sad_prev(2,i) ; I_s3 = 1.0/sad_prev(3,i)
              Idenom = 1.0 / (I_s0 + I_s1 + I_s2 + I_s3)
              w1(i) = I_s0*Idenom  ; w2(1,i) = I_s1*Idenom
              w2(2,i) = I_s2*Idenom ; w2(3,i) = I_s3*Idenom; w3(i) = 0.0
            endif
          endif
        endif ; enddo

        do i=is1,ie1 ; if (do_i_Ri(i) .and. ((it2 > set_best(i) + 20) .or. &
                         ((it2 == max_it2) .and. (it2 /= set_best(i))) )) then
          ! It is time to stop working on this column and use the best estimate
          ! found so far.  (And, yes, the i & k loops are out of order.)
          do_i_Ri(i) = .false.
          do k=max(kb(i,j),kb_min_act),nz-1 ; F(i,k) = F_best(i,k) ; enddo
        endif ; enddo

        do k=kb_min_act,nz-1 ; do i=is1,ie1
          if ((k>=kb(i,j)) .and. do_i_Ri(i)) then
            tmp = F_Ri_est(i,k)
            if (inject_best(i)) then
              F_Ri_est(i,k) = F_best(i,k)
            else
              F_Ri_est(i,k) = w3(i)*F(i,k) + w1(i)*F_Ri_est(i,k) + &
                 w2(1,i)*F_Ri_est_prev(1,i,k) + w2(2,i)*F_Ri_est_prev(2,i,k) + &
                 w2(3,i)*F_Ri_est_prev(3,i,k)
            endif
            if (m_set(i) > 0) F_Ri_est_prev(m_set(i),i,k) = tmp
            if (set_best(i) == it2) F_best(i,k) = F(i,k)
          endif
        enddo ; enddo
      endif ! End of the elseif it2<max_it2 branch.

!   Calculate a revised estimate of the Richardson number dependent
! entrainment rate.
      if ((it2<max_it2) .and. reiterate) &
        call Calculate_Rino_flux(u,v,h,tv,kb,j,dt,G,CS, gk_gkp1, gkp1_gk, &
                                 gkb_gkbp1, gkbp1_gkb, gml_gkbp1, gkbp1_gml, &
                                 .false., F_Ri_est,maxEnt, do_i_Ri, gratsdt, &
                                 gratskbdt, I2gprime, FRi0, fm_ent)
      if (.not.reiterate) exit

    enddo ! end of it2 loop

    call F_to_ent(F, h, kb, kmb, j, G, CS, gk_gkp1, gkp1_gk, gkb_gkbp1, gkbp1_gkb, &
                  gml_gkbp1, gkbp1_gml, ea(:,j,:), eb(:,j,:))

    ! Calculate the layer thicknesses after the entrainment to constrain the corrective fluxes.
    if (correct_density) then
      do i=is,ie
        h_guess(i,1) = h(i,j,1) - Angstrom + (eb(i,j,1) - ea(i,j,2))
        h_guess(i,nz) = h(i,j,nz) - Angstrom + (ea(i,j,nz) - eb(i,j,nz-1))
        if (h_guess(i,1) < 0.0) h_guess(i,1) = 0.0
        if (h_guess(i,nz) < 0.0) h_guess(i,nz) = 0.0
      enddo
      do k=2,nz-1 ; do i=is,ie
        h_guess(i,k) = h(i,j,k) - Angstrom + ((ea(i,j,k) - eb(i,j,k-1)) + &
                   (eb(i,j,k) - ea(i,j,k+1)))
        if (h_guess(i,k) < 0.0) h_guess(i,k) = 0.0
      enddo ; enddo
      if (CS%bulkmixedlayer) then
        do k=nz-1,kb_min,-1
          call calculate_density(tv%T(is:ie,j,k), tv%S(is:ie,j,k), pres(is:ie), &
                                 Rcv(is:ie), 1, ie-is+1, tv%eqn_of_state)
          do i=is,ie
            if ((k>kb(i,j)) .and. (F(i,k) > 0.0)) then
              ! Within a time step, a layer may entrain no more than
              ! its thickness for correction.  This limitation should
              ! apply extremely rarely, but precludes undesirable
              ! behavior.
              F_cor = h(i,j,k) * MIN(gkp1_gk(k) , MAX(-1.0, &
                         (G%Rlay(k) - Rcv(i)) / (G%Rlay(k+1)-G%Rlay(k))) )

              ! Ensure that (1) Entrainments are positive, (2) Corrections in
              ! a layer cannot deplete the layer itself (very generously), and
              ! (3) a layer can take no more than a quarter the mass of its neighbor.
              if (F_cor > 0.0) then
                F_cor = MIN(F_cor, 0.9*F(i,k), gk_gkp1(k)*0.5*h_guess(i,k), &
                            0.25*h_guess(i,k+1))
              else
                F_cor = -MIN(-F_cor, 0.9*F(i,k), 0.5*h_guess(i,k), &
                             0.25*gk_gkp1(k)*h_guess(i,k-1) )
              endif

              ea(i,j,k) = ea(i,j,k) - gkp1_gk(k)*F_cor
              eb(i,j,k) = eb(i,j,k) + F_cor
            else if ((k==kb(i,j)) .and. (F(i,k) > 0.0)) then
              !   The layer beneath the buffer layer should not be lightened to
              ! avoid unphysical convective mixed layer deepening.
              F_cor = h(i,j,k) * MIN(gkbp1_gkb(i) , MAX(0.0, &
                         (G%Rlay(k) - Rcv(i)) / (G%Rlay(k+1)-G%Rlay(k))) )

              if (F_cor > 0.0) then
                F_cor = MIN(F_cor, 0.9*F(i,k), 0.25*h_guess(i,k+1))
                ent_min_rem = ea(i,j,k) - 0.5*h_guess(i,k)
                if (ent_min_rem > 0.0) then
!                  F_cor = MIN(F_cor, 0.5*gkb_gkbp1(i)*h_guess(i,k))
                  if (ent_min_rem <= h(i,j,kmb)-Angstrom) then
                    F_cor = MIN(F_cor, 0.5*gkb_gkbp1(i)*h_guess(i,k))
                  else
                    ! If convective instabilities within the mixed and buffer layers
                    ! are not being avoided or overruled, the following code will be
                    ! necessary to ensure that ea is reduced by no more than h_guess/2.
                    F_cor_max = F(i,k) - (h(i,j,kmb)-Angstrom)*gkb_gkbp1(i)
                    ent_min_rem = ent_min_rem - (h(i,j,kmb)-Angstrom)
                    do k3=kmb-1,1,-1
                      if (ent_min_rem < (h(i,j,k3)-Angstrom)) then
                        F_cor_max = F_cor_max - ent_min_rem*gml_gkbp1(i,k3)
                        ent_min_rem = 0.0 ; exit
                      else
                        F_cor_max = F_cor_max - (h(i,j,k3)-Angstrom)*gml_gkbp1(i,k3)
                        ent_min_rem = ent_min_rem - (h(i,j,k3)-Angstrom)
                      endif
                    enddo
                    F_cor = MIN(F_cor, F_cor_max)
                  endif
                endif

                ea(i,j,k) = ea_kb(F(i,k)-F_cor,h,i,j,kmb,gkb_gkbp1,gkbp1_gkb, &
                                  gml_gkbp1,gkbp1_gml,G)
                eb(i,j,k) = eb(i,j,k) + F_cor
              endif
            else if (k < kb(i,j)) then
              ! Repetitive, unless ea(kb) has been corrected.
              ea(i,j,k) = ea(i,j,k+1)
            endif
          enddo
        enddo
        do k=kb_min-1,K2,-1 ; do i=is,ie
          ea(i,j,k) = ea(i,j,k+1)
        enddo ; enddo
        do k=kmb,2,-1 ; do i=is,ie
          ! Repetitive, unless ea(kb) has been corrected.
          ea(i,j,k) = MAX(0.0,ea(i,j,k+1)-h(i,j,k)+Angstrom)
        enddo ; enddo
      else ! not bulkmixedlayer
        do k=K2,nz-1
          call calculate_density(tv%T(is:ie,j,k), tv%S(is:ie,j,k), pres(is:ie), &
                                 Rcv(is:ie), 1, ie-is+1, tv%eqn_of_state)
          do i=is,ie ; if (F(i,k) > 0.0) then
            ! Within a time step, a layer may entrain no more than
            ! its thickness for correction.  This limitation should
            ! apply extremely rarely, but precludes undesirable
            ! behavior.
            F_cor = h(i,j,k) * MIN(gkp1_gk(k) , MAX(-1.0, &
                       (G%Rlay(k) - Rcv(i)) / (G%Rlay(k+1)-G%Rlay(k))) )

            ! Ensure that (1) Entrainments are positive, (2) Corrections in
            ! a layer cannot deplete the layer itself (very generously), and
            ! (3) a layer can take no more than a quarter the mass of its neighbor.
            if (F_cor >= 0.0) then
              F_cor = MIN(F_cor, 0.9*F(i,k), 0.5*gkp1_gk(k)*h_guess(i,k), &
                          0.25*h_guess(i,k+1))
            else
              F_cor = -MIN(-F_cor, 0.9*F(i,k), 0.5*h_guess(i,k), &
                           0.25*gk_gkp1(k)*h_guess(i,k-1) )
            endif

            ea(i,j,k) = ea(i,j,k) - gkp1_gk(k)*F_cor
            eb(i,j,k) = eb(i,j,k) + F_cor
          endif ; enddo
        enddo
      endif

    endif   ! correct_density

    if (ASSOCIATED(CS%diag%Kd)) then
      do k=2,nz-1 ; do i=is,ie
        if (k<kb(i,j)) then ; Kd_here = 0.0
        elseif (k==kb(i,j)) then
          Kd_here = F(i,k) * (h(i,j,k) + ea(i,j,k) + eb(i,j,k) - ea(i,j,k+1) - eb(i,j,k-1)) / &
            (I2p2gkbp1_gkb(i) * gratskbdt(i))
        else
          Kd_here = F(i,k) * (h(i,j,k) + ea(i,j,k) + eb(i,j,k) - ea(i,j,k+1) - eb(i,j,k-1)) / &
            (I2p2gkp1_gk(k) * gratsdt(k))
        endif

        CS%diag%Kd(i,j,k) = G%H_to_m**2 * MAX(Kd(i,k),Kd_here)
      enddo ; enddo
      do i=is,ie
        CS%diag%Kd(i,j,1) = G%H_to_m**2 * Kd(i,1)
        CS%diag%Kd(i,j,nz) = G%H_to_m**2 * Kd(i,nz)
      enddo
    endif

  enddo ! end of j loop

  ! Offer diagnostic fields for averaging.
  if (CS%id_Kd > 0) call post_data(CS%id_Kd, CS%diag%Kd, CS%diag)

end subroutine Entrainment_H2000


subroutine F_to_ent(F, h, kb, kmb, j, G, CS, gk_gkp1, gkp1_gk, gkb_gkbp1, &
                    gkbp1_gkb, gml_gkbp1, gkbp1_gml, ea, eb, do_i_in)
  real, dimension(NXMEM_,NZ_),          intent(in)  :: F
  real, dimension(NXMEM_,NYMEM_,NZ_),   intent(in)  :: h
  integer, dimension(NXMEM_,NYMEM_),    intent(in)  :: kb
  integer,                              intent(in)  :: kmb, j
  type(ocean_grid_type),                intent(in)  :: G
  type(entrain_H2000_CS),               intent(in)  :: CS
  real, dimension(NZ_),                 intent(in)  :: gk_gkp1, gkp1_gk
  real, dimension(NXMEM_),              intent(in)  :: gkb_gkbp1, gkbp1_gkb
  real, dimension(NXMEM_,NZ_),          intent(in)  :: gml_gkbp1, gkbp1_gml
  real, dimension(NXMEM_,NZ_),          intent(out) :: ea, eb
  logical, dimension(NXMEM_), optional, intent(in)  :: do_i_in
!   This subroutine calculates the actual entrainments (ea and eb) and the
! amount of surface forcing that is applied to each layer if there is no bulk
! mixed layer.  ea and eb are in the same units as h - m or kg m-2.

  logical :: do_i(SZ1_(h))
  integer :: i, k, is, ie, nz

  is = G%isc ; ie = G%iec ; nz = G%ke

  if (present(do_i_in)) then
    do i=is,ie ; do_i(i) = do_i_in(i) ; enddo
    do i=G%isc,G%iec ; if (do_i(i)) then
      is = i ; exit
    endif ; enddo
    do i=G%iec,G%isc,-1 ; if (do_i(i)) then
      ie = i ; exit
    endif ; enddo
  else
    do i=is,ie ; do_i(i) = .true. ; enddo
  endif

  do i=is,ie
    ea(i,nz) = 0.0 ; eb(i,nz) = 0.0
  enddo
  if (CS%bulkmixedlayer) then
    do i=is,ie
      ea(i,1) = 0.0 ; eb(i,1) = 0.0
    enddo
    do k=nz-1,kmb+1,-1 ; do i=is,ie ; if (do_i(i)) then
      if (k > kb(i,j)) then
        ! With a bulk mixed layer, surface buoyancy fluxes are applied
        ! elsewhere, so F should always be nonnegative.
        ea(i,k) = gkp1_gk(k)*F(i,k)
        eb(i,k) = F(i,k)
      else if (k == kb(i,j)) then
        ea(i,k) = ea_kb(F(i,k),h,i,j,kmb,gkb_gkbp1,gkbp1_gkb,gml_gkbp1,gkbp1_gml,G)
        eb(i,k) = F(i,k)
      else
        ea(i,k) = ea(i,k+1)
        eb(i,k) = 0.0
      endif
    endif ; enddo ; enddo
    do k=kmb,2,-1 ; do i=is,ie ; if (do_i(i)) then
      ea(i,k) = MAX(0.0,ea(i,k+1)-h(i,j,k)+G%Angstrom)
      eb(i,k) = 0.0
    endif ; enddo ; enddo
  else                                          ! not BULKMIXEDLAYER
    ! Calculate the entrainment by each layer from above and below.
    ! Entrainment is always positive, but F may be negative due to
    ! surface buoyancy fluxes.
    do i=is,ie
      ea(i,1) = 0.0 ; eb(i,1) = MAX(F(i,1),0.0)
      ea(i,2) = gkp1_gk(2)*F(i,2) - MIN(F(i,1),0.0)
    enddo

    do k=2,nz-1 ; do i=is,ie
      eb(i,k) = MAX(F(i,k),0.0)
      ea(i,k+1) = gkp1_gk(k+1)*F(i,k+1) - (F(i,k)-eb(i,k))
      if (ea(i,k+1) < 0.0) then
        eb(i,k) = eb(i,k) - ea(i,k+1)
        ea(i,k+1) = 0.0
      endif
    enddo ; enddo
  endif                                         ! end BULKMIXEDLAYER
end subroutine F_to_ent

function ea_kb(F_kb,h,i,j,kmb,gkb_gkbp1,gkbp1_gkb,gml_gkbp1,gkbp1_gml,G)
  real :: ea_kb
  real, intent(in) :: F_kb
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in)  :: h
  integer, intent(in) :: i, j, kmb
  real, dimension(NXMEM_), intent(in) :: gkb_gkbp1, gkbp1_gkb
  real, dimension(NXMEM_,NZ_), intent(in)  :: gml_gkbp1, gkbp1_gml
  type(ocean_grid_type), intent(in)    :: G
! Arguments: F_kb - The value of F in the buffer layer, in m.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      i,j,k - indices of the uppermost massive interior layer.
!  (in)      kmb
!  (in)      gkbp1_gkb - the inverse of the ratio of the reduced gravity
!                         at the interface between the buffer layer and the
!                         next denser layer and the reduced gravity of the
!                         next deeper interface.  Nondimensional.
!  (in)      gkb_gkbp1 - g_kb / g_kb+1. Nondimensional.
!  (in)      gkbp1_gml - the inverse of the ratio of the density difference
!                        between the mixed layer and the layer below the
!                        buffer layer to the density difference across the
!                        next deeper interface. Nondimensional.
!  (in)      gml_gkbp1 -  1 / gkbp1_gml. Nondimensional.
!  (in)      G - The ocean's grid structure.

! This function returns the entrainment from above by the layer below
! the buffer layer.
  real :: F1, h_ent
  integer :: k3

! Here F1 is the flux that remains to be entrained from one of the
! buffer or mixed layers.  h_ent is the fluid entrained from each layer.
  h_ent = MIN((h(i,j,kmb)-G%Angstrom),gkbp1_gkb(i)*F_kb)
  F1 = F_kb - h_ent*gkb_gkbp1(i)
  ea_kb = h_ent
  do k3=kmb-1,1,-1
    if (F1 <= 0.0) exit
    h_ent = MIN((h(i,j,k3)-G%Angstrom),gkbp1_gml(i,k3)*F1)
    F1 = F1 - h_ent*gml_gkbp1(i,k3)
    ea_kb = ea_kb + h_ent
  enddo

end function ea_kb

subroutine Calculate_Rino_flux(u, v, h, tv, kb, j, dt, G, CS, gk_gkp1, gkp1_gk, &
                               gkb_gkbp1, gkbp1_gkb, gml_gkbp1, gkbp1_gml, internal_estimate, &
                               F, maxEnt, do_i, gratsdt, gratskbdt, I2gprime, FRi0, fm_ent)
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(in)    :: u
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(in)    :: v
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(in)    :: h
  type(thermo_var_ptrs),               intent(in)    :: tv
  integer, dimension(NXMEM_,NYMEM_),   intent(in)    :: kb
  integer,                             intent(in)    :: j
  real,                                intent(in)    :: dt
  type(ocean_grid_type),               intent(in)    :: G
  type(entrain_H2000_CS),              intent(in)    :: CS
  real, dimension(NZ_),                intent(in)    :: gk_gkp1, gkp1_gk
  real, dimension(NXMEM_),             intent(in)    :: gkb_gkbp1, gkbp1_gkb
  real, dimension(NXMEM_,NZ_),         intent(in)    :: gml_gkbp1, gkbp1_gml
  logical,                             intent(in)    :: internal_estimate
  real,    dimension(NXMEM_,NZ_),      intent(in)    :: F, maxEnt
  logical, dimension(NXMEM_),          intent(in)    :: do_i
  real,    dimension(NZ_),             intent(in)    :: gratsdt
  real,    dimension(NXMEM_),          intent(in)    :: gratskbdt
  real,    dimension(NXMEM_,NZ_),      intent(inout) :: I2gprime
  real,    dimension(NXMEM_,NZ_),      intent(out)   :: FRi0
  real,    dimension(NXMEM_,NZ_),      intent(out)   :: fm_ent
!   This subroutine calculates FRi0, a measure of the rate of Rich-
! ardson number dependent mixing.  Specifically, FRi0 is the density
! flux through a layer due to shear Richardson number dependent en-
! trainment, in the absence of entrainment by the neighboring layers,
! divided by the density difference across the interface below the
! layer and I2p2gkp1_gk (i.e. roughly multiplied by 4).

! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      kb - The index of the lightest layer denser than the
!                 buffer layer.  kb has NYMEM rows or is a NULL ptr.
!  (in)      j - The meridional index upon which to work.
!  (in)      dt - The time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - This module's control structure.
!  (in)      gkp1_gk - The reduced gravity of the interface below a layer
!                      divided by the reduced gravity of the interface above
!                      it. Nondimensional.
!  (in)      gk_gkp1 - The reduced gravity of the interface above a layer
!                      divided by the reduced gravity of the interface below
!                      it. Nondimensional.
!  (in)      gkbp1_gkb - the inverse of the ratio of the reduced gravity
!                         at the interface between the buffer layer and the
!                         next denser layer and the reduced gravity of the
!                         next deeper interface.  Nondimensional.
!  (in)      gkb_gkbp1 - g_kb / g_kb+1. Nondimensional.
!  (in)      gkbp1_gml - the inverse of the ratio of the density difference
!                        between the mixed layer and the layer below the
!                        buffer layer to the density difference across the
!                        next deeper interface. Nondimensional.
!  (in)      gml_gkbp1 -  1 / gkbp1_gml. Nondimensional.
!  (in)      internal_estimate - A flag that indicates that the
!                 entrainment must be estimated internally, rather
!                 than using ea and eb.
!  (in)      F - The layer buoyancy flux, in units of m or kg m-2.
!  (in)      maxEnt - maximum value of entrainment from below (with
!                     compensating entrainment from above to keep the
!                     layer density from changing) that will not
!                     deplete all of the layers above or below a
!                     layer within a timestep, in m or kg m-2.
!  (in)      do_i - A mask indicating which columns to work on.
!  (in)      gratsdt - 2*dt*(2 + g_k+1/g_k + g_k/g_k+1).
!  (in)      gratskbdt - 2*dt*(2 + g_kb+1/g_kb + g_kb/g_kb+1).
!  (inout)   I2gprime - a work-space for this subroutine, containing
!                       rho_0 / 2g divided by the in-situ potential density
!                       difference across a layer, in s2 m-1.
!  (out)     FRi0 - The density flux through a layer due to shear
!                   Richardson number dependent entrainment in the
!                   absence of entrainment by the neighboring layers,
!                   divided by the density difference across the
!                   interface below the layer and I2p2gkp1_gk, in m or kg m-2.

  real, dimension(SZ1_(h),SZ3_(h)) :: &
    duvert2, &  ! The squared velocity difference between adjacent layers,
                ! in m2 s-2.
    u_h, &      ! The zonal and meridional velocities at layer thickness points,
    v_h, &      ! used to calculate h_Ri, in m s-1.

    F_tmp, &    ! An internal estimate of F, in m or kg m-2.
    Ent_a, &    ! The layer entrainment from above in m or kg m-2.
    Ent_b       ! The layer entrainment from below in m or kg m-2.
  real, dimension(SZ1_(h)) :: &
    pressure, & ! The in-situ pressure at the layer center in Pa.
    drho_dT, &  ! The partial derivative of density with temperature and
    drho_dS     ! salinity at ambient pressure, temperature and salinity,
                ! in kg m-3 K-1 and kg m-3 psu-1.

  real, dimension(SZ1_(h)) :: &
    RiNo_crit, &  ! The critical Richardson number, using the
                  ! same definition of the critical Richardson
                  ! number as used in Hallberg (MWR 2000).
    shear_rate    ! A nondimensional scaling for the entrainment rate,
                  ! divided by 80 - see Hallberg (MWR 2000) for details.

  real :: gR0_2         ! 0.5 G times the conversion from thickness to layer
                        ! mass, in kg m-2 s-2 or m s-2.
  real :: g_R0          ! g_R0 is g/Rho in m4 kg-1 s-2.

  real :: dRa, dRb      ! The density differences across the interfaces above
                        ! and below a layer, in kg m-3.
  real :: h_Ri          ! The layer thickness divided by the layer
                        ! shear Richardson number, in m or kg m-2.
  real :: sh2           ! The summed squared velocity differences between
                        ! layers, in m2 s-2.
  real :: fm2, fr2, fk2 ! Work variables with units of H, H, H2, where H is
                        ! the units of h, m or kg m-2.
  real :: entmag        ! A scaling value for the Ri#-dependent layer
                        ! entrainment within a timestep, in m or kg m-2.
  real :: h_min         ! A thickness over which to average u_h and v_h, in
                        ! the same units as h, m or kg m-2.
  real :: z             ! A nondimensional temporary variable.
  real :: PI            ! 3.1415926... calculated as 4*atan(1)
  logical :: use_EOS    !   If true, density is calculated from temperature
                        ! and salnity using an equation of state.
  integer :: i, k, is, ie, nz, K2, kmb
  integer, save :: j_prev = -1000 ! The value of j from the previous call.

  is = G%isc ; ie = G%iec ; nz = G%ke
  use_EOS = associated(tv%eqn_of_state)
  if (CS%bulkmixedlayer) then ; kmb = tv%nk_Rml ; K2 = kmb + 1
  else ; kmb = 0 ; K2 = 2 ; endif
  gR0_2 = 0.5*G%H_to_Pa ; g_R0 = G%g_Earth/G%Rho0

  PI = 4.0*atan(1.0)
  do i=is,ie
    if (abs(G%geolath(i,j)) >= 2.0*CS%Shearmix_lat_eq) then
      RiNo_crit(i) = CS%RiNo_crit
      shear_rate(i) = 0.0125*CS%Shearmix_rate
    elseif (abs(G%geolath(i,j)) <= CS%Shearmix_lat_eq) then
      RiNo_crit(i) = CS%RiNo_crit_eq
      shear_rate(i) = 0.0125*CS%Shearmix_rate_eq
    else
      z = cos(PI*(abs(G%geolath(i,j))-CS%Shearmix_lat_eq) / CS%Shearmix_lat_eq)
      RiNo_crit(i) = z*CS%RiNo_crit_eq + (1.0-z)*CS%RiNo_crit
      shear_rate(i) = 0.0125*(z*CS%Shearmix_rate_eq + (1.0-z)*CS%Shearmix_rate)
    endif
  enddo

  if (internal_estimate) then
! If internal_estimate is true, thin layers entrain until they have
! a thickness h_min.  This is a crude estimate, but avoids viscous
! homogenization of velocities.
    h_min = G%m_to_H * sqrt(dt*(MAX(CS%Kv,CS%Kd)))
    do i=is,ie
      F_tmp(i,1) = 0.0 ; F_tmp(i,nz) = 0.0
    enddo
    do k=2,nz-1 ; do i=is,ie
      if (k<kb(i,j)) then
        F_tmp(i,k) = 0.0
      elseif (k==kb(i,j)) then
        F_tmp(i,k) = MIN(MAX(h_min-h(i,j,k),0.0)/(1.0+gkbp1_gkb(i)),maxEnt(i,k))
      else
        F_tmp(i,k) = MIN(MAX(h_min-h(i,j,k),0.0)/(1.0+gkp1_gk(k)),maxEnt(i,k))
      endif
    enddo ; enddo
    call F_to_ent(F_tmp, h, kb, kmb, j, G, CS, gk_gkp1, gkp1_gk, gkb_gkbp1, &
                  gkbp1_gkb, gml_gkbp1, gkbp1_gml, Ent_a, Ent_b, do_i)
    do k=K2,nz-1 ; do i=is,ie
      fm_ent(i,k) = -h(i,j,k) + gkp1_gk(k+1)*F_tmp(i,k+1)
      if (k > kb(i,j)) fm_ent(i,k) = fm_ent(i,k) + F_tmp(i,k-1)
    enddo ; enddo
  else
    call F_to_ent(F, h, kb, kmb, j, G, CS, gk_gkp1, gkp1_gk, gkb_gkbp1, &
                  gkbp1_gkb, gml_gkbp1, gkbp1_gml, Ent_a, Ent_b, do_i)
    do k=K2,nz-1 ; do i=is,ie ; if (do_i(i)) then
      fm_ent(i,k) = -h(i,j,k) + gkp1_gk(k+1)*F(i,k+1)
      if (k > kb(i,j)) fm_ent(i,k) = fm_ent(i,k) + F(i,k-1)
    endif ; enddo ; enddo
  endif ! end of internal_estimate branches


! Calculate u_h and v_h.
  call Estimate_u_h(u, v, h, j, Ent_a, Ent_b, do_i, u_h, v_h, G)

! Calculate the gprime across each layer.  This does not change with iterations.
  if ((use_EOS) .and. (j /= j_prev)) then
    do i=is,ie
      pressure(i) = gR0_2*h(i,j,1)
    enddo
    do k=2,K2-1 ; do i=is,ie
      pressure(i) = pressure(i) + gR0_2 * (h(i,j,k-1)+h(i,j,k))
    enddo ; enddo
    do k=K2,nz-1
      do i=is,ie
        pressure(i) = pressure(i) + gR0_2 * (h(i,j,k-1)+h(i,j,k))
      enddo
      call calculate_density_derivs(tv%T(:,j,k),tv%S(:,j,k),pressure,drho_dT, &
                               drho_dS,is,ie-is+1, tv%eqn_of_state)
      do i=is,ie ; if ((k >= kb(i,j)) .and. (do_i(i))) then
        if (k > kb(i,j)) then
          dRa = drho_dT(i) * (tv%T(i,j,k)-tv%T(i,j,k-1)) + &
                drho_dS(i) * (tv%S(i,j,k)-tv%S(i,j,k-1))
        else
          dRa = drho_dT(i) * (tv%T(i,j,k)-tv%T(i,j,kmb)) + &
                drho_dS(i) * (tv%S(i,j,k)-tv%S(i,j,kmb))
        endif
        if (dRa <= 0.0) dRa = 0.0
!   The following threshold, 1e-12 kg m-3, should guarantee no
! floating exceptions, provided that duvert < 1e4 m2 s-2, and other
! values are appropriate to the Earth's oceans.
        dRb = drho_dT(i) * (tv%T(i,j,k+1)-tv%T(i,j,k)) + &
              drho_dS(i) * (tv%S(i,j,k+1)-tv%S(i,j,k))
        if (dRb <= 1.0e-12) dRb = 1.0e-12

! It may be that we should be calculating these fluxes in such a way that
! h_Ri is calculating using the properties of whichever mixed or buffer layer
! would be entrained, both here and later where sh2 and duvert2 are calculated.
        if (k > kb(i,j)) then
          I2gprime(i,k) = 1.0 / (g_R0*(gkp1_gk(k)*dRa + dRb))
        else
          I2gprime(i,k) = 1.0 / (g_R0*(gkbp1_gkb(i)*dRa + dRb))
        endif
      endif ; enddo
    enddo
    j_prev = j
  endif                                             ! use_EOS.

! Calculate h_Ri.

  if (.not.CS%bulkmixedlayer) then
    do i=is,ie ; if (do_i(i)) then
      duvert2(i,2) = (u_h(i,2)-u_h(i,1)) * (u_h(i,2)-u_h(i,1)) + &
                     (v_h(i,2)-v_h(i,1)) * (v_h(i,2)-v_h(i,1))
    endif ; enddo
  endif                                             ! BULKMIXEDLAYER

  do k=K2,nz-1 ; do i=is,ie ; if ((k >= kb(i,j)) .and. (do_i(i))) then
    duvert2(i,k+1) = (u_h(i,k+1)-u_h(i,k)) * (u_h(i,k+1)-u_h(i,k)) + &
                     (v_h(i,k+1)-v_h(i,k)) * (v_h(i,k+1)-v_h(i,k))
    if (k == kb(i,j)) then
      duvert2(i,k) = (u_h(i,k)-u_h(i,kmb)) * (u_h(i,k)-u_h(i,kmb)) + &
                     (v_h(i,k)-v_h(i,kmb)) * (v_h(i,k)-v_h(i,kmb))
!  The buffer layer is usually thin or has about the same density
!  as the mixed layer, so gkbp1_gml is used instead of gkbp1_gkb.
!      sh2 = (duvert2(i,k+1) + gkbp1_gml(i,CS%nkml)*duvert2(i,k))
      sh2 = (duvert2(i,k+1) + gkbp1_gkb(i)*duvert2(i,k))
    else
      sh2 = (duvert2(i,k+1) + gkp1_gk(k)*duvert2(i,k))
    endif

    if (use_EOS) then
      h_Ri = G%m_to_H * (sh2 * I2gprime(i,k))
    else
      h_Ri = G%m_to_H * (0.5 * sh2 / G%g_prime(k+1))
    endif

    if (h(i,j,k) > RiNo_crit(i)*h_Ri) then
      FRi0(i,k) = 2.0*(RiNo_crit(i)*h_Ri - h(i,j,k))
    else
      if (k == kb(i,j)) then
        entmag = G%m_to_H * (shear_rate(i) * gratskbdt(i) * &
                             sqrt((gkb_gkbp1(i)+1.0)*sh2))
      else
        entmag = G%m_to_H * (shear_rate(i) * gratsdt(k) * &
                             sqrt((gk_gkp1(k)+1.0)*sh2))
      endif
      fm2 = 0.5*h(i,j,k) - 0.1*h_Ri - entmag
      fk2 = h_Ri*(0.2*h(i,j,k) + 2.0*RiNo_crit(i)*entmag)
      fr2 = sqrt(fm2*fm2 + fk2)

      if (fm2>=0.0) then
        FRi0(i,k) = 2.0*(-h(i,j,k) + (fm2 + fr2))
      else
        FRi0(i,k) = 2.0*(-h(i,j,k) + (fk2 / (-fm2+fr2)))
      endif
    endif
  endif ; enddo ; enddo

end subroutine Calculate_Rino_flux


subroutine Estimate_u_h(u, v, h, j, Ent_a, Ent_b, do_i, u_h, v_h, G)
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(in) :: u
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(in) :: v
  real, dimension(NXMEM_,NYMEM_,NZ_),  intent(in) :: h
  integer,                             intent(in)  :: j
  real, dimension(NXMEM_,NZ_),         intent(in)  :: Ent_a, Ent_b
  logical, dimension(NXMEM_),          intent(in)  :: do_i
  real, dimension(NXMEM_,NZ_),         intent(out) :: u_h, v_h
  type(ocean_grid_type), intent(in) :: G
!   This subroutine calculates u_h and v_h (velocities at thickness
! points), using the entrainments (in m) passed in as arguments.

! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      j - The meridional index upon which to work.
!  (in)      Ent_a - The amount of fluid entrained from the layer
!                    above within this time step, in units of m or kg m-2.
!  (in)      Ent_b - The amount of fluid entrained from the layer
!                    below within this time step, in units of m or kg m-2.
!  (out)     u_h - The zonal velocity at thickness points after
!                  entrainment, in m s-1.
!  (out)     v_h - The meridional velocity at thickness points after
!                  entrainment, in m s-1.
!  (in)      G - The ocean's grid structure.
  real :: b1(SZ1_(h))         ! b1 is used in the tridiagonal solver, in m-1.
  real :: c1(SZ1_(h),SZ3_(h)) ! c1 is used in the tridiagonal solver, nondim.
  real :: a_n(SZ1_(h)), a_s(SZ1_(h))  ! Fractional weights of the neighboring
  real :: a_e(SZ1_(h)), a_w(SZ1_(h))  ! velocity points, ~1/2 in the open
                                      ! ocean, nondimensional.
  real :: s                   ! The sum of neighboring velocity cell areas, m2.
  real :: Idenom              ! A scaling factor for the relative weights of
                              ! neighboring velocity points, in m-2.
  integer :: i, k, is, ie, nz
  is = G%isc ; ie = G%iec ; nz = G%ke

  do i=is,ie ; if (do_i(i)) then
    s = G%dxdy_u(I-1,j)+G%dxdy_u(I,j)
    if (s>0.0) then
      Idenom = sqrt(0.5*G%IDXDYh(i,j)/s)
      a_w(i) = G%dxdy_u(I-1,j)*Idenom
      a_e(i) = G%dxdy_u(I,j)*Idenom
    else
      a_w(i) = 0.0 ; a_e(i) = 0.0
    endif

    s = G%dxdy_v(i,J-1)+G%dxdy_v(i,J)
    if (s>0.0) then
      Idenom = sqrt(0.5*G%IDXDYh(i,j)/s)
      a_s(i) = G%dxdy_v(i,J-1)*Idenom
      a_n(i) = G%dxdy_v(i,J)*Idenom
    else
      a_s(i) = 0.0 ; a_n(i) = 0.0
    endif

    b1(i) = 1.0/(h(i,j,1) + Ent_b(i,1))
    u_h(i,1) = h(i,j,1)*b1(i)*(a_e(i)*u(I,j,1)+a_w(i)*u(I-1,j,1))
    v_h(i,1) = h(i,j,1)*b1(i)*(a_n(i)*v(i,J,1)+a_s(i)*v(i,J-1,1))
  endif ; enddo
  do k=2,nz ; do i=is,ie ; if (do_i(i)) then
    c1(i,k) = Ent_b(i,k-1) * b1(i)
    b1(i) = 1.0/((h(i,j,k) + Ent_b(i,k)) + (1.0-c1(i,k))*Ent_a(i,k))
    u_h(i,k)=(h(i,j,k)*(a_e(i)*u(I,j,k)+a_w(i)*u(I-1,j,k)) + &
               Ent_a(i,k)*u_h(i,k-1))*b1(i)
    v_h(i,k)=(h(i,j,k)*(a_n(i)*v(i,J,k)+a_s(i)*v(i,J-1,k)) + &
               Ent_a(i,k)*v_h(i,k-1))*b1(i)
  endif ; enddo ; enddo
  do k=nz-1,1,-1 ; do i=is,ie ; if (do_i(i)) then
    u_h(i,k) = u_h(i,k) + c1(i,k+1)*u_h(i,k+1)
    v_h(i,k) = v_h(i,k) + c1(i,k+1)*v_h(i,k+1)
  endif ; enddo ; enddo

end subroutine Estimate_u_h


subroutine entrain_H2000_init(Time, G, param_file, diag, CS)
  type(time_type),         intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ptrs), target, intent(inout) :: diag
  type(entrain_H2000_CS),  pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
  logical :: use_temperature
  character(len=128) :: version = '$Id: GOLD_entrain_H2000.F90,v 1.1.2.14 2010/09/07 13:53:17 rwh Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod  = "GOLD_entrain_H2000" ! This module's name.

  if (associated(CS)) then
    call GOLD_error(WARNING, "entrain_H2000_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag

  CS%bulkmixedlayer = .false. ; CS%RiNo_mix = .false.
  call read_param(param_file,"BULKMIXEDLAYER",CS%bulkmixedlayer)
  if (G%ke > 2) call read_param(param_file,"RINOMIX",CS%RiNo_mix)

  call read_param(param_file,"TEMPERATURE",use_temperature,.true.)
  CS%correct_density = use_temperature
  if (use_temperature) &
    call read_param(param_file,"CORRECT_DENSITY",CS%correct_density)

  call read_param(param_file,"KD",CS%Kd,.true.)
  call read_param(param_file,"KV",CS%Kv,.true.)
  CS%max_ent_it = 5
  call read_param(param_file,"MAX_ENT_IT",CS%max_ent_it)

  CS%Shearmix_rate = 0.1 ; CS%max_RiNo_it = 50 ; CS%RiNo_crit = 0.0
  if (CS%RiNo_mix) then
    call read_param(param_file,"SHEARMIX_RATE",CS%Shearmix_rate)
    call read_param(param_file,"MAX_RINO_IT",CS%max_RiNo_it)
    call read_param(param_file,"RINO_CRIT",CS%RiNo_crit,.true.)
  endif

  CS%Shearmix_rate_eq = CS%Shearmix_rate ; CS%Shearmix_lat_eq = 0.0
  CS%RiNo_crit_eq = CS%RiNo_crit
  if (CS%RiNo_mix) then
    call read_param(param_file,"SHEARMIX_RATE_EQ",CS%Shearmix_rate_eq)
    call read_param(param_file,"RINO_CRIT_EQ",CS%RiNo_crit_eq)
    call read_param(param_file,"SHEARMIX_LAT_EQ",CS%Shearmix_lat_eq)
  endif

  CS%id_Kd = register_diag_field('ocean_model', 'Kd', G%axeshl, Time, &
      'Diapycnal diffusivity', 'meter2 second-1')
  if (CS%id_Kd > 0) call safe_alloc_ptr(diag%Kd,G%isd,G%ied,G%jsd,G%jed,G%ke)

  ! Write all relevant parameters to the model log.
  call log_version(param_file, mod, version, tagname)
  call log_param(param_file, mod, "BULKMIXEDLAYER", CS%bulkmixedlayer, &
                 "If defined, use a refined Kraus-Turner-like bulk mixed \n"//&
                 "layer with transitional buffer layers.")
  call log_param(param_file, mod, "RINOMIX", CS%RiNo_mix, &
                 "Use Richardson number dependent mixing. The mixing rate \n"//&
                 "is proportional to the velocity shears when the shear \n"//&
                 "Richardson number drops below RINO_CRIT.  When \n"//&
                 "USE_H2000_SHEAR_MIXING is true, the scheme documented \n"//&
                 "in Hallberg (MWR, 2000) is used.", default=.false.)
  call log_param(param_file, mod, "TEMPERATURE", use_temperature, &
                 "If true, Temperature and salinity are used as state \n"//&
                 "variables.")
  if (use_temperature) &
    call log_param(param_file, mod, "CORRECT_DENSITY", CS%correct_density, &
                 "If true, the layer densities are restored toward their \n"//&
                 "target values by the diapycnal mixing, as described in \n"//&
                 "Hallberg (MWR, 2000).", default=.true.)
  call log_param(param_file, mod, "KD", CS%Kd, &
                 "The background diapycnal diffusivity of density. \n"//&
                 "Zero or the molecular value, ~1e-7 m2 s-1, may be used.", &
                 units="m2 s-1")
  call log_param(param_file, mod, "KV", CS%Kv, &
                 "The background kinematic viscosity. The molecular \n"//&
                 "value, ~1e-6 m2 s-1, may be used.", units="m2 s-1")
  call log_param(param_file, mod, "MAX_ENT_IT", CS%max_ent_it, &
                 "The maximum number of iterations that may be used to \n"//&
                 "calculate the interior diapycnal entrainment.", default=5)
  if (CS%RiNo_mix) then
    call log_param(param_file, mod, "SHEARMIX_RATE", CS%Shearmix_rate, &
                 "A nondimensional rate scale for shear-driven \n"//&
                 "entrainment, as described in Hallberg (MWR 2000).", &
                 units="nondim", default=0.1)
    call log_param(param_file, mod, "MAX_RINO_IT", CS%max_RiNo_it, &
                 "The maximum number of iterations that may be used to \n"//&
                 "estimate the Richardson number driven mixing.", &
                 units="nondim", default=50)
    call log_param(param_file, mod, "RINO_CRIT", CS%RiNo_crit, &
                 "The critical shear Richardson number for shear-driven \n"//&
                 "entrainment.  The theoretical value here is 1 using \n"//&
                 "the convention for bulk Richardson numbers, as in \n"//&
                 "Hallberg (MWR 2000).", units="nondim")
    call log_param(param_file, mod, "SHEARMIX_RATE_EQ", CS%Shearmix_rate_eq, &
                 "A potentially different value of SHEARMIX_RATE that is \n"//&
                 "used near the equator.  The default is to use the same \n"//&
                 "value as SHEARMIX_RATE.", units="nondim", &
                 default=CS%Shearmix_rate)
    call log_param(param_file, mod, "RINO_CRIT_EQ", CS%RiNo_crit_eq, &
                 "A potentially different value of RINO_CRIT that is \n"//&
                 "used near the equator.  The default is to use the same \n"//&
                 "value as RINO_CRIT.", units="nondim", default=CS%RiNo_crit)
    call log_param(param_file, mod, "SHEARMIX_LAT_EQ", CS%Shearmix_lat_eq, &
                 "The latitude for transitions from the equatorial values \n"//&
                 "of shear mixing to midlatitude values.", units="degrees", &
                 default=0.0)
  endif

end subroutine entrain_H2000_init

subroutine entrain_H2000_end(CS)
  type(entrain_H2000_CS), pointer :: CS

  if (associated(CS)) deallocate(CS)

end subroutine entrain_H2000_end

end module GOLD_entrain_H2000
