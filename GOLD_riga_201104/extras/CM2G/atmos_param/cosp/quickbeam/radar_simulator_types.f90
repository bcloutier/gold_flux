 
!---------------------------------------------------------------------
!------------ FMS version number and tagname for this file -----------
        
! $Id: radar_simulator_types.f90,v 1.1.2.1.2.1.2.1 2009/10/06 17:54:20 rsh Exp $
! $Name: riga_201104 $

  module radar_simulator_types

   public radar_simulator_types_init

! Collection of common variables and types
! Part of QuickBeam v1.03 by John Haynes
! http://reef.atmos.colostate.edu/haynes/radarsim

  integer, parameter ::       &
  maxhclass = 20 	     ,& ! max number of hydrometeor classes
  nd = 85		     ,& ! number of discrete particles  
  nRe_types = 250		! number or Re size bins allowed in N and Z_scaled look up table

  real*8, parameter ::        &
  dmin = 0.1                 ,& ! min size of discrete particle
  dmax = 10000.                	! max size of discrete particle
   
  integer, parameter :: &
  mt_nfreq = 5              , &
  mt_ntt = 39               , &	! num temperatures in table
  mt_nf	= 14		    , &	! number of ice fractions in table  
  mt_nd = 85                   ! num discrete mode-p drop sizes in table


! ---- hydrometeor class type -----  
  
  type class_param
    real*8,  dimension(maxhclass) :: p1,p2,p3,dmin,dmax,apm,bpm,rho
    integer, dimension(maxhclass) :: dtype,col,cp,phase
    logical, dimension(maxhclass,nRe_types) :: scaled
    logical, dimension(maxhclass,mt_ntt,nRe_types) :: z_flag
    real*8,  dimension(maxhclass,mt_ntt,nRe_types) :: Ze_scaled,Zr_scaled,kr_scaled
    real*8,  dimension(maxhclass,nd,nRe_types) :: fc, rho_eff
    integer, dimension(maxhclass,nd,nRe_types) :: ifc
    integer, dimension(maxhclass) :: idd
  end type class_param

! ----- mie table structure -----
  
  type mie
    real*8 :: freq(mt_nfreq), tt(mt_ntt), f(mt_nf), D(mt_nd)
    real*8, dimension(mt_nd,mt_ntt,mt_nf,mt_nfreq) :: qext, qbsca
    integer :: phase(mt_ntt)
  end type mie

  real*8, dimension(:), allocatable :: &
    mt_ttl, &			! liquid temperatures (C)
    mt_tti, &			! ice temperatures (C)
    mt_qext, mt_qbsca		! extincion/backscatter efficiency

  integer*4 :: &
    cnt_liq, &			! liquid temperature count
      cnt_ice			! ice temperature count

  contains

subroutine radar_simulator_types_init

    
    integer :: i

! otherwise we still need to initialize temperature arrays for Ze 
! scaling (which is only done when not using mie table)
           
           cnt_ice=19
           cnt_liq=20
!      if (.not.(allocated(mt_ttl).and.allocated(mt_tti))) then
          allocate(mt_ttl(cnt_liq),mt_tti(cnt_ice))  ! note needed as th        is is global array ... 
                                                      ! which should be c        hanged in the future 
!     endif
                  
        do i=1,cnt_ice
           mt_tti(i)=(i-1)*5-90
        enddo 

        do i=1,cnt_liq
          mt_ttl(i)=(i-1)*5 - 60
       enddo

     
end subroutine radar_simulator_types_init


  end module radar_simulator_types
