
module ocean_rough_mod

!-----------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
use          mpp_mod, only: input_nml_file
#else
use          fms_mod, only: open_namelist_file
#endif

use       fms_mod, only: error_mesg, FATAL, file_exist,  &
                         check_nml_error, mpp_pe, mpp_root_pe, close_file, &
                         write_version_number, stdlog
use constants_mod, only: grav

implicit none
private

public :: compute_ocean_roughness, fixed_ocean_roughness

!-----------------------------------------------------------------------
character(len=256) :: version = '$Id: ocean_rough.F90,v 18.0.4.1 2010/09/14 19:29:20 pjp Exp $'
character(len=256) :: tagname = '$Name: riga_201104 $'
!-----------------------------------------------------------------------
!----- namelist -----

  real    :: roughness_init = 0.00044   ! not used in this version
  real    :: roughness_min  = 1.e-6
  real    :: charnock       = 0.032
  
  real    :: roughness_mom   = 5.8e-5
  real    :: roughness_heat  = 5.8e-5   ! was 4.00e-4
  real    :: roughness_moist = 5.8e-5
!  real, parameter :: zcoh1 = 0.0       ! Beljaars 1994 values
!  real, parameter :: zcoq1 = 0.0
! real, parameter :: zcoh1 = 1.4e-5
! real, parameter :: zcoq1 = 1.3e-4
  real            :: zcoh1 = 0.0 !miz
  real            :: zcoq1 = 0.0 !miz
  logical :: do_highwind     = .false.
  logical :: do_cap40        = .false.

  character(len=32) :: rough_scheme = 'fixed'   ! possible values:
                                                !   'fixed'
                                                !   'charnock'
                                                !   'beljaars'

namelist /ocean_rough_nml/ roughness_init, roughness_heat,  &
                           roughness_mom,  roughness_moist, &
                           roughness_min,                   &
                           charnock,                        &
                           rough_scheme, do_highwind,       &!miz
                           do_cap40, zcoh1, zcoq1            !sjl


!-----------------------------------------------------------------------

  logical :: do_init = .true.

!-----------------------------------------------------------------------
! ---- constants ----

! ..... high wind speed - rough sea
  real, parameter :: zcom1 = 1.8e-2    ! Charnock's constant
! ..... low wind speed - smooth sea
  real, parameter :: gnu   = 1.5e-5
  real, parameter :: zcom2 = 0.11
  real, parameter :: zcoh2 = 0.40
  real, parameter :: zcoq2 = 0.62


contains

!#######################################################################

 subroutine compute_ocean_roughness ( ocean, u_star,  &
                                      rough_mom, rough_heat, rough_moist )

 logical, intent(in)  :: ocean(:,:)
 real,    intent(in)  :: u_star(:,:)
 real,    intent(out) :: rough_mom(:,:), rough_heat(:,:), rough_moist(:,:)

!-----------------------------------------------------------------------
!  computes ocean roughness for momentum using wind stress
!  and sets roughness for heat/moisture using namelist value
!-----------------------------------------------------------------------

   real, dimension(size(ocean,1),size(ocean,2)) :: ustar2, xx1, xx2, w10 !miz
   real ::  a=0.001, b=0.028 !miz

   if (do_init) call ocean_rough_init


   if (trim(rough_scheme) == 'fixed') then

!  --- set roughness for momentum and heat/moisture ---

      call fixed_ocean_roughness ( ocean, rough_mom, rough_heat, &
                                          rough_moist )


!  --- compute roughness for momentum, heat, moisture ---

   else if (trim(rough_scheme) == 'beljaars' .or. &
            trim(rough_scheme) == 'charnock') then

      where (ocean)
          ustar2(:,:) = max(gnu*gnu,u_star(:,:)*u_star(:,:))          
          xx1(:,:) = gnu / sqrt(ustar2(:,:))
          xx2(:,:) = ustar2(:,:) / grav
      elsewhere
          rough_mom   = 0.0
          rough_heat  = 0.0
          rough_moist = 0.0
      endwhere

      if (trim(rough_scheme) == 'charnock') then
          where (ocean)
              rough_mom  (:,:) = charnock * xx2(:,:)
              rough_mom  (:,:) = max( rough_mom(:,:), roughness_min )
              rough_heat (:,:) = rough_mom  (:,:)
              rough_moist(:,:) = rough_mom  (:,:)
          endwhere
      else if (trim(rough_scheme) == 'beljaars') then
!     --- Beljaars scheme ---

! SJL*** High Wind correction following Moon et al 2007 ***
          if (do_highwind) then

              if ( do_cap40 ) then
     
              where (ocean)
                  w10(:,:) = 2.458 + u_star(:,:)*(20.255-0.56*u_star(:,:))  ! Eq(7) Moon et al.
                  where ( w10(:,:) > 12.5 )
! SJL mods: cap the growth of z0 with w10 up to 40 m/s
!                     rough_mom(:,:) = 0.001*(0.085*min(w10(:,:), 40.) - 0.58)    ! capped Eq(8b) Moon et al.
! z0 (w10=40) = 2.82E-3
                      rough_mom(:,:) = 0.001*(0.085*w10(:,:) - 0.58)    ! Eq(8b) Moon et al.
                      rough_mom(:,:) = min( rough_mom(:,:), 2.82E-3)
                  elsewhere     
                      rough_mom(:,:) = 0.0185/grav*u_star(:,:)**2  ! (8a) Moon et al.
                  endwhere
                  rough_heat (:,:) = zcoh1 * xx2(:,:) + zcoh2 * xx1(:,:)
                  rough_moist(:,:) = zcoq1 * xx2(:,:) + zcoq2 * xx1(:,:)
!             --- lower limit on roughness? ---
                  rough_mom  (:,:) = max( rough_mom  (:,:), roughness_min )
                  rough_heat (:,:) = max( rough_heat (:,:), roughness_min )
                  rough_moist(:,:) = max( rough_moist(:,:), roughness_min )
              endwhere

              else

              where (ocean)
                  w10(:,:) = 2.458 + u_star(:,:)*(20.255-0.56*u_star(:,:))  ! Eq(7) Moon et al.
                  where ( w10(:,:) > 12.5 )
                      rough_mom(:,:) = 0.001*(0.085*w10(:,:) - 0.58)    ! Eq(8b) Moon et al.
                  elsewhere     
                      rough_mom(:,:) = 0.0185/grav*u_star(:,:)**2  ! (8a) Moon et al.
                  endwhere
                  rough_heat (:,:) = zcoh1 * xx2(:,:) + zcoh2 * xx1(:,:)
                  rough_moist(:,:) = zcoq1 * xx2(:,:) + zcoq2 * xx1(:,:)
!             --- lower limit on roughness? ---
                  rough_mom  (:,:) = max( rough_mom  (:,:), roughness_min )
                  rough_heat (:,:) = max( rough_heat (:,:), roughness_min )
                  rough_moist(:,:) = max( rough_moist(:,:), roughness_min )
              endwhere

              endif
! SJL ----------------------------------------------------------------------------------------

          else
          where (ocean)
              rough_mom  (:,:) = zcom1 * xx2(:,:) + zcom2 * xx1(:,:)
              rough_heat (:,:) = zcoh1 * xx2(:,:) + zcoh2 * xx1(:,:)
              rough_moist(:,:) = zcoq1 * xx2(:,:) + zcoq2 * xx1(:,:)
!             --- lower limit on roughness? ---
              rough_mom  (:,:) = max( rough_mom  (:,:), roughness_min )
              rough_heat (:,:) = max( rough_heat (:,:), roughness_min )
              rough_moist(:,:) = max( rough_moist(:,:), roughness_min )
          endwhere
          endif
      endif
   endif

!-----------------------------------------------------------------------

 end subroutine compute_ocean_roughness

!#######################################################################

 subroutine fixed_ocean_roughness ( ocean, rough_mom, rough_heat, rough_moist )

 logical, intent(in)  :: ocean(:,:)
 real,    intent(out) :: rough_mom(:,:), rough_heat(:,:), rough_moist(:,:)

   if (do_init) call ocean_rough_init

    where (ocean)
       rough_mom   = roughness_mom
       rough_heat  = roughness_heat
       rough_moist = roughness_moist
    endwhere

 end subroutine fixed_ocean_roughness

!#######################################################################

 subroutine ocean_rough_init

   integer :: unit, ierr, io

!   ----- read and write namelist -----

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_rough_nml, iostat=io)
  ierr = check_nml_error(io, 'ocean_rough_nml')
#else
    if ( file_exist('input.nml')) then
          unit = open_namelist_file ('input.nml')
        ierr=1; do while (ierr /= 0)
           read  (unit, nml=ocean_rough_nml, iostat=io, end=10)
           ierr = check_nml_error(io,'ocean_rough_nml')
        enddo
 10     call close_file (unit)
    endif
#endif

!------- write version number and namelist ---------

    if ( mpp_pe() == mpp_root_pe() ) then
         call write_version_number(version, tagname)
         unit = stdlog()
         write (unit,nml=ocean_rough_nml)
         write (unit,11)
    endif

!------ constants -----

    roughness_moist = max (roughness_moist, roughness_min)
    roughness_heat  = max (roughness_heat , roughness_min)
    roughness_mom   = max (roughness_mom  , roughness_min)

    do_init = .false.

11 format (/,'namelist option USE_FIXED_ROUGH is no longer supported', &
           /,'use variable ROUGH_SCHEME instead')

 end subroutine ocean_rough_init

!#######################################################################

end module ocean_rough_mod

