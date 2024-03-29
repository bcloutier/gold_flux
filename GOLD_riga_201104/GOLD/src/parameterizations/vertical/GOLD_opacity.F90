module GOLD_opacity
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
!*  Macros written all in capital letters are defined in GOLD_memory.h.*
!*                                                                     *
!*   This module will house the routines used to calcualte the opacity *
!* of the ocean.                                                       *
!*                                                                     *
!* CHL_from_file:                                                      *
!*   In this routine, the Morel (modified) and Manizza (modified)      *
!* schemes use the "blue" band in the paramterizations to determine    *
!* the e-folding depth of the incoming shortwave attenuation. The red  *
!* portion is lumped into the net heating at the surface.              *
!*                                                                     *
!* Morel, A., 1988: Optical modeling of the upper ocean in relation    *
!*   to itsbiogenous matter content (case-i waters)., J. Geo. Res.,    *
!*   93, 10,749-10,768.                                                *
!*                                                                     *
!* Manizza, M., C. LeQuere, A. J. Watson, and E. T. Buitenhuis, 2005:  *
!*  Bio-optical feedbacks amoung phytoplankton, upper ocean physics    *
!*  and sea-ice in a global model, Geophys. Res. Let., 32, L05603,     *
!*  doi:10.1029/2004GL020778.                                          *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, buoy, Rml, eaml, ebml, etc.           *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**
use GOLD_diag_mediator, only : time_type, diag_ptrs, safe_alloc_ptr, post_data
use GOLD_diag_mediator, only : query_averaging_enabled, register_diag_field
use GOLD_time_manager, only :  get_time
use GOLD_error_handler, only : GOLD_error, GOLD_mesg, FATAL, WARNING
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_file_parser, only : uppercase
use GOLD_grid, only : ocean_grid_type
use GOLD_io, only : slasher
use GOLD_tracer_flow_control, only : get_chl_from_model, tracer_flow_control_CS
use GOLD_variables, only : forcing, thermo_var_ptrs, optics_type
use time_interp_external_mod, only : init_external_field, time_interp_external
use time_interp_external_mod, only : time_interp_external_init
implicit none ; private

#include <GOLD_memory.h>

public set_opacity, opacity_init, opacity_end, opacity_manizza, opacity_morel

type, public :: opacity_CS ; private
  logical :: var_pen_sw      !   If true, use one of the CHL_A schemes
                             ! (specified below) to determine the e-folding
                             ! depth of incoming short wave radiation.
                             ! The default is false.
  integer :: opacity_scheme  !   An integer indicating which scheme should be
                             ! used to translate water properties into the
                             ! opacity (i.e., the e-folding depth) and (perhaps)
                             ! the number of bands of penetrating shortwave
                             ! radiation to use.
  real :: pen_sw_scale       !   The vertical absorption e-folding depth of the
                             ! penetrating shortwave radiation, in m.
  real :: pen_sw_frac        !   The fraction of shortwave radiation that is
                             ! penetrating with a constant e-folding approach.
  real :: blue_frac          !   The fraction of the penetrating shortwave
                             ! radiation that is in the blue band, ND.
  real :: opacity_land_value ! The value to use for opacity over land, in m-1.
                             ! The default is 10 m-1 - a value for muddy water.
  integer :: sbc_chl         ! An integer handle used in time interpolation of
                             ! chlorophyll read from a file.
  character(len=128) :: chl_file ! Data containing chl_a concentrations. Used
                             ! when var_pen_sw is defined and reading from file.
  logical ::  chl_from_file  !   If true, chl_a is read from a file.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
                             ! ocean diagnostic fields.
  type(tracer_flow_control_CS), pointer  :: tracer_flow_CSp => NULL() 
                    ! A pointer to the control structure of the tracer modules.

  integer :: id_sw_pen = -1, id_sw_vis_pen = -1, id_chl = -1
  integer, pointer :: id_opacity(:) => NULL()
end type opacity_CS

integer, parameter :: NO_SCHEME = 0, MANIZZA_05 = 1, MOREL_88 = 2

character*(10), parameter :: MANIZZA_05_STRING = "MANIZZA_05"
character*(10), parameter :: MOREL_88_STRING = "MOREL_88"

contains

subroutine set_opacity(optics, fluxes, G, CS)
  type(optics_type),                   intent(inout) :: optics
  type(forcing),                       intent(in)    :: fluxes  
  type(ocean_grid_type),               intent(in)    :: G
  type(opacity_CS),                    pointer       :: CS
! Arguments: (inout) opacity - The inverse of the vertical absorption decay
!                     scale for penetrating shortwave radiation, in m-1.
!            (inout) fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!            (in)    G - The ocean's grid structure.
!            (in)    CS - The control structure earlier set up by opacity_init.    

! local variables 
  integer :: i, j, k, n, is, ie, js, je, nz
  real :: inv_sw_pen_scale  ! The inverse of the e-folding scale, in m-1.
  real :: Inv_nbands        ! The inverse of the number of bands of penetrating
                            ! shortwave radiation.
  logical :: call_for_surface  ! if horizontal slice is the surface layer
  real :: tmp(SZI_(G),SZJ_(G),SZK_(G))  ! A 3-d temporary array.
  real :: chl(SZI_(G),SZJ_(G),SZK_(G))  ! The concentration of chlorophyll-A
                                        ! in mg m-3.
  real :: Pen_SW_tot(SZI_(G),SZJ_(G))   ! The penetrating shortwave radiation
                                        ! summed across all bands, in W m-2.
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not. associated(CS)) call GOLD_error(FATAL, "set_opacity: "// &
         "Module must be initialized via opacity_init before it is used.")

  if (CS%var_pen_sw) then
    if (CS%chl_from_file) then 
      call opacity_from_chl(optics, fluxes, G, CS)
    else 
      call get_chl_from_model(chl, G, CS%tracer_flow_CSp)
      call opacity_from_chl(optics, fluxes, G, CS, chl)
    endif
  else ! Use sw e-folding scale set by GOLD_input
    if (optics%nbands <= 1) then ; Inv_nbands = 1.0
    else ; Inv_nbands = 1.0 / real(optics%nbands) ; endif

    ! Make sure there is no division by 0.
    inv_sw_pen_scale = 1.0 / max(CS%pen_sw_scale, 0.1*G%Angstrom)

    do k=1,nz ; do j=js,je ; do i=is,ie  ; do n=1,optics%nbands
      optics%opacity_band(n,i,j,k) = inv_sw_pen_scale
    enddo ; enddo ; enddo ; enddo
    if (.not.associated(fluxes%sw) .or. (CS%pen_SW_scale <= 0.0)) then
      do j=js,je ; do i=is,ie ; do n=1,optics%nbands
        optics%sw_pen_band(n,i,j) = 0.0
      enddo ; enddo ; enddo
    else
      do j=js,je ; do i=is,ie ; do n=1,optics%nbands
        optics%sw_pen_band(n,i,j) = CS%pen_SW_frac * Inv_nbands * fluxes%sw(i,j)
      enddo ; enddo ; enddo
    endif
  endif
  
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_sw_pen > 0) then
      do j=js,je ; do i=is,ie
        Pen_SW_tot(i,j) = 0.0
        do n=1,optics%nbands
          Pen_SW_tot(i,j) = Pen_SW_tot(i,j) + optics%sw_pen_band(n,i,j)
        enddo
      enddo ; enddo
      call post_data(CS%id_sw_pen, Pen_SW_tot, CS%diag)
    endif
    if (CS%id_sw_vis_pen > 0) then
      if (CS%opacity_scheme == MANIZZA_05) then
        do j=js,je ; do i=is,ie
          Pen_SW_tot(i,j) = 0.0
          do n=1,min(optics%nbands,2)
            Pen_SW_tot(i,j) = Pen_SW_tot(i,j) + optics%sw_pen_band(n,i,j)
          enddo
        enddo ; enddo
      else
        do j=js,je ; do i=is,ie
          Pen_SW_tot(i,j) = 0.0
          do n=1,optics%nbands
            Pen_SW_tot(i,j) = Pen_SW_tot(i,j) + optics%sw_pen_band(n,i,j)
          enddo
        enddo ; enddo
      endif
      call post_data(CS%id_sw_vis_pen, Pen_SW_tot, CS%diag)
    endif
    do n=1,optics%nbands ; if (CS%id_opacity(n) > 0) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        tmp(i,j,k) = optics%opacity_band(n,i,j,k)
      enddo ; enddo ; enddo
      call post_data(CS%id_opacity(n), tmp, CS%diag)
    endif ; enddo
  endif

end subroutine set_opacity


subroutine opacity_from_chl(optics, fluxes, G, CS, chl_in)
  type(optics_type),              intent(inout)  :: optics
  type(forcing),                  intent(in)     :: fluxes  
  type(ocean_grid_type),          intent(in)     :: G
  type(opacity_CS),               pointer        :: CS
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in), optional :: chl_in
! Arguments: fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (out)     opacity - The inverse of the vertical absorption decay
!                           scale for penetrating shortwave radiation, in m-1.
!  (in)      G - The ocean's grid structure.
!  (in)      chl_in - A 3-d field of chlorophyll A, in mg m-3.

  real :: chl_data(SZI_(G),SZJ_(G)) ! The chlorophyll A concentrations in
                                    ! a layer, in mg/m^3.
  real :: Inv_nbands        ! The inverse of the number of bands of penetrating
                            ! shortwave radiation.
  real :: Inv_nbands_nir    ! The inverse of the number of bands of penetrating
                            ! near-infrafed radiation.
  real :: SW_pen_tot        ! The sum across the bands of the penetrating
                            ! shortwave radiation, in W m-2.
  real :: SW_vis_tot        ! The sum across the visible bands of shortwave
                            ! radiation, in W m-2.
  real :: SW_nir_tot        ! The sum across the near infrared bands of shortwave
                            ! radiation, in W m-2.
  type(time_type) :: day
  character(len=128) :: mesg
  integer :: days, seconds
  integer :: i, j, k, n, is, ie, js, je, nz, nbands
  logical :: multiband_vis_input, multiband_nir_input
  
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

!   In this model, the Morel (modified) and Manizza (modified) schemes
! use the "blue" band in the paramterizations to determine the e-folding
! depth of the incoming shortwave attenuation. The red portion is lumped
! into the net heating at the surface.
!
! Morel, A., Optical modeling of the upper ocean in relation to its biogenous
!   matter content (case-i waters).,J. Geo. Res., {93}, 10,749--10,768, 1988.
!

! Manizza, M., C.~L. Quere, A.~Watson, and E.~T. Buitenhuis, Bio-optical
!   feedbacks amoung phytoplankton, upper ocean physics and sea-ice in a
!   global model, Geophys. Res. Let., , L05,603, 2005.

  nbands = optics%nbands

  if (nbands <= 1) then ; Inv_nbands = 1.0
  else ; Inv_nbands = 1.0 / real(nbands) ; endif
  
  if (nbands <= 2) then ; Inv_nbands_nir = 0.0
  else ; Inv_nbands_nir = 1.0 / real(nbands - 2.0) ; endif

  multiband_vis_input = (associated(fluxes%sw_vis_dir) .and. &
                         associated(fluxes%sw_vis_dif))
  multiband_nir_input = (associated(fluxes%sw_nir_dir) .and. &
                         associated(fluxes%sw_nir_dif))

  do k=1,nz

    if (present(chl_in)) then
      do j=js,je ; do i=is,ie ; chl_data(i,j) = chl_in(i,j,k) ; enddo ; enddo
      do j=js,je ; do i=is,ie
        if ((G%hmask(i,j) > 0.5) .and. (chl_data(i,j) < 0.0)) then
          write(mesg,'(" Negative chl_in of ",(1pe12.4)," found at i,j,k = ", &
                    & 3(1x,i3), " lon/lat = ",(1pe12.4)," E ", (1pe12.4), " N.")') &
                     chl_data(i,j), i, j, k, G%geolonh(i,j), G%geolath(i,j)
          call GOLD_error(FATAL,"GOLD_opacity opacity_from_chl: "//trim(mesg))
        endif
      enddo ; enddo
    elseif (k==1) then
      ! Only the 2-d surface chlorophyll can be read in from a file.  The
      ! same value is assumed for all layers.
      call get_time(CS%Time,seconds,days)
      do j=js,je ; do i=is,ie ; chl_data(i,j) = 0.0 ; enddo ; enddo
      call time_interp_external(CS%sbc_chl, CS%Time, chl_data)
      do j=js,je ; do i=is,ie
        if ((G%hmask(i,j) > 0.5) .and. (chl_data(i,j) < 0.0)) then
          write(mesg,'(" Time_interp negative chl of ",(1pe12.4)," at i,j = ",&
                    & 2(i3), "lon/lat = ",(1pe12.4)," E ", (1pe12.4), " N.")') &
                     chl_data(i,j), i, j, G%geolonh(i,j), G%geolath(i,j)
          call GOLD_error(FATAL,"GOLD_opacity opacity_from_chl: "//trim(mesg))
        endif
      enddo ; enddo
    endif

    select case (CS%opacity_scheme)
      case (MANIZZA_05)
        if (k==1) then ; do j=js,je ; do i=is,ie
          SW_vis_tot = 0.0 ; SW_nir_tot = 0.0
          if (G%hmask(i,j) > 0.5) then
            if (multiband_vis_input) then
              SW_vis_tot = fluxes%sw_vis_dir(i,j) + fluxes%sw_vis_dif(i,j)
            else  ! Follow Manizza 05 in assuming that 42% of SW is visible.
              SW_vis_tot = 0.42 * fluxes%sw(i,j)
            endif
            if (multiband_nir_input) then
              SW_nir_tot = fluxes%sw_nir_dir(i,j) + fluxes%sw_nir_dif(i,j)
            else
              SW_nir_tot = fluxes%sw(i,j) - SW_vis_tot
            endif
          endif

          ! Band 1 is Manizza blue.
          optics%sw_pen_band(1,i,j) = CS%blue_frac*SW_vis_tot
          ! Band 2 (if used) is Manizza red.
          if (nbands > 1) &
            optics%sw_pen_band(2,i,j) = (1.0-CS%blue_frac)*SW_vis_tot
          ! All remaining bands are NIR, for lack of something better to do.
          do n=3,nbands
            optics%sw_pen_band(n,i,j) = Inv_nbands_nir * SW_nir_tot
          enddo
        enddo ; enddo ; endif

        do j=js,je ; do i=is,ie
          if (G%hmask(i,j) <= 0.5) then
            do n=1,optics%nbands
              optics%opacity_band(n,i,j,k) = CS%opacity_land_value
            enddo
          else
            ! Band 1 is Manizza blue.
            optics%opacity_band(1,i,j,k) = 0.0232 + 0.074*chl_data(i,j)**0.674
            if (nbands >= 2) &  !  Band 2 is Manizza red.
              optics%opacity_band(2,i,j,k) = 0.225 + 0.037*chl_data(i,j)**0.629
            ! All remaining bands are NIR, for lack of something better to do.
            do n=3,nbands ; optics%opacity_band(n,i,j,k) = 2.86 ; enddo
          endif
        enddo ; enddo

      case (MOREL_88)
        if (k==1) then ! Set up the surface fluxes.
          do j=js,je ; do i=is,ie
            SW_pen_tot = 0.0
            if (G%hmask(i,j) > 0.5) then ; if (multiband_vis_input) then
                SW_pen_tot = SW_pen_frac_morel(chl_data(i,j)) * &
                    (fluxes%sw_vis_dir(i,j) + fluxes%sw_vis_dif(i,j))
              else
                SW_pen_tot = SW_pen_frac_morel(chl_data(i,j)) * &
                    0.5*fluxes%sw(i,j)
            endif ; endif

            do n=1,optics%nbands
              optics%sw_pen_band(n,i,j) = Inv_nbands*SW_pen_tot
            enddo
          enddo ; enddo
        endif

        do j=js,je ; do i=is,ie
          optics%opacity_band(1,i,j,k) = CS%opacity_land_value
          if (G%hmask(i,j) > 0.5) &
            optics%opacity_band(1,i,j,k) = opacity_morel(chl_data(i,j))

          do n=2,optics%nbands
            optics%opacity_band(n,i,j,k) = optics%opacity_band(1,i,j,k)
          enddo
        enddo ; enddo

      case default
        call GOLD_error(FATAL, "opacity_from_chl: CS%opacity_scheme is not valid.")
    end select

    if ((k==1) .and. (CS%id_chl > 0)) then
      call post_data(CS%id_chl, chl_data, CS%diag)
    endif
  enddo
end subroutine opacity_from_chl

function opacity_morel(chl_data)
  real, intent(in)  :: chl_data
  real :: opacity_morel
! Argument : chl_data - The chlorophyll-A concentration in mg m-3.
!   The following are coefficients for the optical model taken from Morel and
! Antoine (1994). These coeficients represent a non uniform distribution of
! chlorophyll-a through the water column.  Other approaches may be more
! appropriate when using an interactive ecosystem model that predicts
! three-dimensional chl-a values.
  real, dimension(6), parameter :: &
       Z2_coef=(/7.925, -6.644, 3.662, -1.815, -0.218,  0.502/)
  real :: Chl, Chl2 ! The log10 of chl_data (in mg m-3), and Chl^2.

  Chl = log10(min(max(chl_data,0.02),60.0)) ; Chl2 = Chl*Chl
  opacity_morel = 1.0 / ( (Z2_coef(1) + Z2_coef(2)*Chl) + Chl2 * &
      ((Z2_coef(3) + Chl*Z2_coef(4)) + Chl2*(Z2_coef(5) + Chl*Z2_coef(6))) )
end function                       

function SW_pen_frac_morel(chl_data)
  real, intent(in)  :: chl_data
  real :: SW_pen_frac_morel
! Argument : chl_data - The chlorophyll-A concentration in mg m-3.
!   The following are coefficients for the optical model taken from Morel and
! Antoine (1994). These coeficients represent a non uniform distribution of
! chlorophyll-a through the water column.  Other approaches may be more
! appropriate when using an interactive ecosystem model that predicts
! three-dimensional chl-a values.
  real :: Chl, Chl2         ! The log10 of chl_data in mg m-3, and Chl^2.
  real, dimension(6), parameter :: & 
       V1_coef=(/0.321,  0.008, 0.132,  0.038, -0.017, -0.007/)

  Chl = log10(min(max(chl_data,0.02),60.0)) ; Chl2 = Chl*Chl
  SW_pen_frac_morel = 1.0 - ( (V1_coef(1) + V1_coef(2)*Chl) + Chl2 * &
       ((V1_coef(3) + Chl*V1_coef(4)) + Chl2*(V1_coef(5) + Chl*V1_coef(6))) )
end function SW_pen_frac_morel

function opacity_manizza(chl_data)
  real, intent(in)  :: chl_data
  real :: opacity_manizza
! Argument : chl_data - The chlorophyll-A concentration in mg m-3.
!   This sets the blue-wavelength opacity according to the scheme proposed by
! Manizza, M. et al, 2005.
  
  opacity_manizza = 0.0232 + 0.074*chl_data**0.674
end function  
   
subroutine opacity_init(Time, G, param_file, diag, tracer_flow, CS, optics)
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ptrs), target, intent(inout) :: diag
  type(tracer_flow_control_CS), target, intent(in) :: tracer_flow 
  type(opacity_CS),        pointer       :: CS  
  type(optics_type),       pointer       :: optics
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  
!  (in/out)  CS - A pointer that is set to point to the control structure
!                  for this module
  character(len=128) :: version = '$Id: GOLD_opacity.F90,v 1.1.2.2.2.20 2011/09/19 16:10:36 Robert.Hallberg Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=200) :: inputdir   ! The directory where NetCDF input files
  character(len=200) :: filename
  character(len=200) :: tmpstr
  character(len=40)  :: mod
  character(len=40)  :: bandnum, shortname
  character(len=200) :: longname
  character(len=40)  :: scheme_string
  logical :: use_scheme
  integer :: isd, ied, jsd, jed, nz, n
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke

  if (associated(CS)) then
    call GOLD_error(WARNING, "opacity_init called with an associated"// &
                             "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag
  CS%Time => Time
  CS%tracer_flow_CSp => tracer_flow

! parameters for CHL_A routines
  CS%var_pen_sw = .false.
  call read_param(param_file,"VAR_PEN_SW",CS%var_pen_sw)
  CS%opacity_scheme = NO_SCHEME ; scheme_string = ""
  if (CS%var_pen_sw) then
    use_scheme = .false. ; call read_param(param_file,"MOREL_PEN_SW",use_scheme)
    if (use_scheme) then
      call GOLD_mesg('The opacity scheme has been set by "#define MOREL_PEN_SW."'&
             //'  This is depricated.  Use "#define OPACITY_SCHEME MOREL_88" instead.')
      CS%opacity_scheme = MOREL_88 ; scheme_string = MOREL_88_STRING
    endif
    use_scheme = .false. ; call read_param(param_file,"MANIZZA_PEN_SW",use_scheme)
    if (use_scheme) then
      call GOLD_mesg('The opacity scheme has been set by "#define MANIZZA_PEN_SW."'&
             //'  This is depricated.  Use "#define OPACITY_SCHEME MANIZZA_05" instead.')
      CS%opacity_scheme = MANIZZA_05 ; scheme_string = MANIZZA_05_STRING
    endif

    tmpstr = "" ; call read_param(param_file,"OPACITY_SCHEME",tmpstr)
    if (len_trim(tmpstr) > 0) then
      tmpstr = uppercase(tmpstr)
      select case (tmpstr)
        case (MANIZZA_05_STRING)
          if ((CS%opacity_scheme /= MANIZZA_05) .and. &
              (CS%opacity_scheme /= NO_SCHEME)) call GOLD_error(FATAL, &
            "opacity_init: The opacity scheme has been set inconsistently "//&
               "in multiple places, last as MANIZZA_05 via OPACITY_SCHEME.")
          CS%opacity_scheme = MANIZZA_05 ; scheme_string = MANIZZA_05_STRING
        case (MOREL_88_STRING)
          if ((CS%opacity_scheme /= MOREL_88) .and. &
              (CS%opacity_scheme /= NO_SCHEME)) call GOLD_error(FATAL, &
            "opacity_init: The opacity scheme has been set inconsistently "//&
               "in multiple places, last as MOREL_88 via OPACITY_SCHEME.")
         CS%opacity_scheme = MOREL_88 ; scheme_string = MOREL_88_STRING
        case default
          call GOLD_error(FATAL, "opacity_init: #DEFINE OPACITY_SCHEME "//&
                                  trim(tmpstr) // "in input file is invalid.")
      end select
      call GOLD_mesg('opacity_init: opacity scheme set to "'//trim(tmpstr)//'".', 5)
    endif
    if (CS%opacity_scheme == NO_SCHEME) then
      call GOLD_error(WARNING, "opacity_init: No scheme has successfully "//&
               "been specified for the opacity.  Using the default MANIZZA_05.")
      CS%opacity_scheme = MANIZZA_05 ; scheme_string = MANIZZA_05_STRING
    endif

    CS%chl_from_file = .true.
    call read_param(param_file,"CHL_FROM_FILE",CS%chl_from_file)
    if (CS%chl_from_file) then 
      call time_interp_external_init()
      inputdir = "." ;
      call read_param(param_file,"INPUTDIR",inputdir)
      call read_param(param_file,"CHL_FILE",CS%chl_file,.true.)
      filename = trim(slasher(inputdir))//trim(CS%chl_file)
      CS%sbc_chl = init_external_field(filename,'CHL_A',domain=G%Domain%mpp_domain)
    endif

    CS%blue_frac = 0.5    
    call read_param(param_file, "BLUE_FRAC_SW", CS%blue_frac)
  else
    CS%pen_SW_scale = 0.0
    call read_param(param_file,"PEN_SW_SCALE",CS%pen_sw_scale)
    CS%pen_SW_frac = 0.0
    call read_param(param_file,"PEN_SW_FRAC",CS%pen_sw_frac)
  endif
  optics%nbands = 1 ; call read_param(param_file,"PEN_SW_NBANDS",optics%nbands)
  if (.not.ASSOCIATED(optics%min_wavelength_band)) &
    allocate(optics%min_wavelength_band(optics%nbands))
  if (.not.ASSOCIATED(optics%max_wavelength_band)) &
    allocate(optics%max_wavelength_band(optics%nbands))

  if (CS%opacity_scheme == MANIZZA_05) then
    optics%min_wavelength_band(1) =0
    optics%max_wavelength_band(1) =550
    if (optics%nbands >= 2) then
      optics%min_wavelength_band(2)=550
      optics%max_wavelength_band(2)=700
    endif
    if (optics%nbands > 2) then
      do n=3,optics%nbands 
        optics%min_wavelength_band(n) =700
        optics%max_wavelength_band(n) =2800
      enddo 
    endif
  endif       
 

  CS%opacity_land_value = 10.0
  call read_param(param_file,"OPACITY_LAND_VALUE",CS%opacity_land_value)

  if (.not.ASSOCIATED(optics%opacity_band)) &
    allocate(optics%opacity_band(optics%nbands,isd:ied,jsd:jed,nz))
  if (.not.ASSOCIATED(optics%sw_pen_band)) &
    allocate(optics%sw_pen_band(optics%nbands,isd:ied,jsd:jed))
  allocate(CS%id_opacity(optics%nbands)) ; CS%id_opacity(:) = -1

  CS%id_sw_pen = register_diag_field('ocean_model', 'SW_pen', G%axesh1, Time, &
      'Penetrating shortwave radiation flux into ocean', 'Watt meter-2')
  CS%id_sw_vis_pen = register_diag_field('ocean_model', 'SW_vis_pen', G%axesh1, Time, &
      'Visible penetrating shortwave radiation flux into ocean', 'Watt meter-2')
  do n=1,optics%nbands
    write(bandnum,'(i3)') n
    shortname = 'opac_'//trim(adjustl(bandnum))
    longname = 'Opacity for shortwave radiation in band '//trim(adjustl(bandnum))
    CS%id_opacity(n) = register_diag_field('ocean_model', shortname, G%axeshL, Time, &
      longname, 'meter-1')
  enddo
  if (CS%var_pen_sw) &
    CS%id_chl = register_diag_field('ocean_model', 'Chl_opac', G%axesh1, Time, &
        'Surface chlorophyll A concentration used to find opacity', 'mg meter-3')

  ! Write all relevant parameters to the model log.
  mod = "GOLD_opacity"
  call log_version(param_file, mod, version, tagname, "")
  call log_param(param_file, mod, "VAR_PEN_SW", CS%var_pen_sw, &
                 "If true, use one of the CHL_A schemes specified by \n"//&
                 "OPACITY_SCHEME to determine the e-folding depth of \n"//&
                 "incoming short wave radiation.", default=.false.)
  if (CS%var_pen_sw) then
    call log_param(param_file, mod, "CHL_FROM_FILE", CS%chl_from_file, &
                 "If true, chl_a is read from a file.", default=.true.)
    if (CS%chl_from_file) then
      call log_param(param_file, mod, "CHL_FILE", CS%chl_file, &
                 "CHL_FILE is the file containing chl_a concentrations in \n"//&
                 "the variable CHL_A. It is used when VAR_PEN_SW and \n"//&
                 "CHL_FROM_FILE are true.")
      call log_param(param_file, mod, "INPUTDIR/CHL_FILE", filename)
    endif
    call log_param(param_file, mod, "OPACITY_SCHEME", scheme_string, &
                 "This character string specifies how chlorophyll \n"//&
                 "concentrations are translated into opacities. Currently \n"//&
                 "valid options include:\n"//&
                 " \t\t  MANIZZA_05 - Use Manizza et al., GRL, 2005. \n"//&
                 " \t\t  MOREL_88 - Use Morel, JGR, 1988.", &
                 default=MANIZZA_05_STRING)
    if (CS%opacity_scheme == MANIZZA_05) then
      call log_param(param_file, mod, "BLUE_FRAC_SW", CS%blue_frac, &
                   "The fraction of the penetrating shortwave radiation \n"//&
                   "that is in the blue band.", default=0.5, units="nondim")
    endif
  else
    call log_param(param_file, mod, "PEN_SW_SCALE", CS%pen_sw_scale, &
                 "The vertical absorption e-folding depth of the \n"//&
                 "penetrating shortwave radiation.", units="m", default=0.0)
    call log_param(param_file, mod, "PEN_SW_FRAC", CS%pen_sw_frac, &
                 "The fraction of the shortwave radiation that penetrates \n"//&
                 "below the surface.", units="nondim", default=0.0)
  endif
  call log_param(param_file, mod, "PEN_SW_NBANDS", optics%nbands, &
                 "The number of bands of penetrating shortwave radiation.", &
                 default=1)
  call log_param(param_file, mod, "OPACITY_LAND_VALUE", CS%opacity_land_value, &
                 "The value to use for opacity over land. The default is \n"//&
                 "10 m-1 - a value for muddy water.", units="m-1", default=10.0)

end subroutine opacity_init


subroutine opacity_end(CS, optics)
  type(opacity_CS),  pointer           :: CS  
  type(optics_type), pointer, optional :: optics

  if (associated(CS%id_opacity)) deallocate(CS%id_opacity)
  if (associated(CS)) deallocate(CS)

  if (present(optics)) then ; if (associated(optics)) then
    if (ASSOCIATED(optics%opacity_band)) deallocate(optics%opacity_band)
    if (ASSOCIATED(optics%sw_pen_band)) deallocate(optics%sw_pen_band)
  endif ; endif

end subroutine opacity_end

end module GOLD_opacity
