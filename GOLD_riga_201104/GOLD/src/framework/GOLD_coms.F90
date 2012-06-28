module GOLD_coms

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

use GOLD_error_handler, only : GOLD_error, FATAL
use fms_mod, only : fms_end
use memutils_mod, only : print_memuse_stats
use mpp_mod, only : PE_here => mpp_pe, root_PE => mpp_root_pe, num_PEs => mpp_npes
use mpp_mod, only : broadcast => mpp_broadcast
use mpp_mod, only : sum_across_PEs => mpp_sum, min_across_PEs => mpp_min
use mpp_mod, only : max_across_PEs => mpp_max, GOLD_infra_init => mpp_init

implicit none ; private

public :: PE_here, root_PE, num_PEs, GOLD_infra_init, GOLD_infra_end
public :: broadcast, sum_across_PEs, min_across_PEs, max_across_PEs
public :: reproducing_sum

!   This module provides interfaces to the non-domain-oriented communication
! subroutines.

integer, parameter :: ni=6    ! The number of long integers to use to represent
                              ! a real number.
integer(kind=8), parameter :: prec=2_8**46 ! The precision of each integer.
real, parameter :: r_prec=2.0**46  ! A real version of prec.
real, parameter :: I_prec=1.0/(2.0**46) ! The inverse of prec.
integer, parameter :: max_count_prec=2**(63-46)-1 
                              ! The number of values that can be added together
                              ! with the current value of prec before there will
                              ! be roundoff problems.
real, parameter, dimension(ni) :: &
  pr = (/ r_prec**2, r_prec, 1.0, 1.0/r_prec, 1.0/r_prec**2, 1.0/r_prec**3 /)
real, parameter, dimension(ni) :: &
  I_pr = (/ 1.0/r_prec**2, 1.0/r_prec, 1.0, r_prec, r_prec**2, r_prec**3 /)

logical :: overflow_error = .false.

interface reproducing_sum
  module procedure reproducing_sum_2d, reproducing_sum_3d
end interface reproducing_sum

contains

function reproducing_sum_2d(array, isr, ier, jsr, jer) result(sum)
  real, dimension(:,:), intent(in) :: array
  integer,    optional, intent(in) :: isr, ier, jsr, jer
  real                             :: sum
  
  !   This subroutine uses a conversion to an integer representation
  ! of real numbers to give order-invariant sums that will reproduce
  ! across PE count.  This idea comes from R. Hallberg and A. Adcroft.

  integer(kind=8), dimension(ni)  :: ints_sum
  integer(kind=8) :: prec_error
  integer :: i, j, is, ie, js, je

  if (num_PEs() > max_count_prec) call GOLD_error(FATAL, &
    "reproducing_sum: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/num_PEs.")

  prec_error = (2_8**62 + (2_8**62 - 1)) / num_PEs()

  is = 1 ; ie = size(array,1) ; js = 1 ; je = size(array,2 )
  if (present(isr)) then
    if (isr < is) call GOLD_error(FATAL, &
      "Value of isr too small in reproducing_sum_2d.")
    is = isr
  endif
  if (present(ier)) then
    if (ier > ie) call GOLD_error(FATAL, &
      "Value of ier too large in reproducing_sum_2d.")
    ie = ier
  endif
  if (present(jsr)) then
    if (jsr < js) call GOLD_error(FATAL, &
      "Value of jsr too small in reproducing_sum_2d.")
    js = jsr
  endif
  if (present(jer)) then
    if (jer > je) call GOLD_error(FATAL, &
      "Value of jer too large in reproducing_sum_2d.")
    je = jer
  endif

  overflow_error = .false.
  ints_sum(:) = 0
  if ((je+1-js)*(ie+1-is) < max_count_prec) then
    do j=js,je ; do i=is,ie
      call increment_ints_fast(ints_sum, real_to_ints(array(i,j), prec_error));
    enddo ; enddo
    call carry_overflow(ints_sum, prec_error)
  elseif ((ie+1-is) < max_count_prec) then
    do j=js,je
      do i=is,ie
        call increment_ints_fast(ints_sum, real_to_ints(array(i,j), prec_error));
      enddo
      call carry_overflow(ints_sum, prec_error)
    enddo
  else
    do j=js,je ; do i=is,ie
      call increment_ints(ints_sum, real_to_ints(array(i,j), prec_error), &
                          prec_error);
    enddo ; enddo
  endif
  if (overflow_error) call GOLD_error(FATAL, "Overflow in reproducing_sum(_2d).")

  call sum_across_PEs(ints_sum, ni)
  
  sum = ints_to_real(ints_sum)   

end function reproducing_sum_2d

function reproducing_sum_3d(array, isr, ier, jsr, jer, sums) result(sum)
  real, dimension(:,:,:),        intent(in) :: array
  integer,    optional,          intent(in) :: isr, ier, jsr, jer
  real, dimension(:), optional, intent(out) :: sums
  real                             :: sum
  
  !   This subroutine uses a conversion to an integer representation
  ! of real numbers to give order-invariant sums that will reproduce
  ! across PE count.  This idea comes from R. Hallberg and A. Adcroft.

  integer(kind=8), dimension(ni)  :: ints_sum
  integer(kind=8), dimension(ni,size(array,3))  :: ints_sums
  integer(kind=8) :: prec_error
  integer :: i, j, k, is, ie, js, je, ke, isz, jsz
  
  if (num_PEs() > max_count_prec) call GOLD_error(FATAL, &
    "reproducing_sum: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/num_PEs.")

  prec_error = (2_8**62 + (2_8**62 - 1)) / num_PEs()

  is = 1 ; ie = size(array,1) ; js = 1 ; je = size(array,2) ; ke = size(array,3)
  if (present(isr)) then
    if (isr < is) call GOLD_error(FATAL, &
      "Value of isr too small in reproducing_sum(_3d).")
    is = isr
  endif
  if (present(ier)) then
    if (ier > ie) call GOLD_error(FATAL, &
      "Value of ier too large in reproducing_sum(_3d).")
    ie = ier
  endif
  if (present(jsr)) then
    if (jsr < js) call GOLD_error(FATAL, &
      "Value of jsr too small in reproducing_sum(_3d).")
    js = jsr
  endif
  if (present(jer)) then
    if (jer > je) call GOLD_error(FATAL, &
      "Value of jer too large in reproducing_sum(_3d).")
    je = jer
  endif
  jsz = je+1-js; isz = ie+1-is
  
  if (present(sums)) then
    if (size(sums) > ke) call GOLD_error(FATAL, "Sums is smaller than "//&
      "the vertical extent of array in reproducing_sum(_3d).")
    ints_sums(:,:) = 0
    overflow_error = .false.
    if (jsz*isz < max_count_prec) then
      do k=1,ke
        do j=js,je ; do i=is,ie
          call increment_ints_fast(ints_sums(:,k), &
                                   real_to_ints(array(i,j,k), prec_error));
        enddo ; enddo
        call carry_overflow(ints_sums(:,k), prec_error)
      enddo
    elseif (isz < max_count_prec) then
      do k=1,ke ; do j=js,je
        do i=is,ie
          call increment_ints_fast(ints_sums(:,k), &
                                   real_to_ints(array(i,j,k), prec_error));
        enddo
        call carry_overflow(ints_sums(:,k), prec_error)
      enddo ; enddo
    else
      do k=1,ke ; do j=js,je ; do i=is,ie
        call increment_ints(ints_sums(:,k), &
                            real_to_ints(array(i,j,k), prec_error), prec_error);
      enddo ; enddo ; enddo
    endif
    if (overflow_error) call GOLD_error(FATAL, "Overflow in reproducing_sum(_3d).")

    call sum_across_PEs(ints_sums(:,1:ke), ni*ke)

    sum = 0.0
    do k=1,ke
      sums(k) = ints_to_real(ints_sums(:,k))
      sum = sum + sums(k)
    enddo
  else
    ints_sum(:) = 0
    overflow_error = .false.
    if (jsz*isz < max_count_prec) then
      do k=1,ke
        do j=js,je ; do i=is,ie
          call increment_ints_fast(ints_sum, &
                                   real_to_ints(array(i,j,k), prec_error));
        enddo ; enddo
        call carry_overflow(ints_sum, prec_error)
      enddo
    elseif (isz < max_count_prec) then
      do k=1,ke ; do j=js,je
        do i=is,ie
          call increment_ints_fast(ints_sum, &
                                   real_to_ints(array(i,j,k), prec_error));
        enddo
        call carry_overflow(ints_sum, prec_error)
      enddo ; enddo
    else
      do k=1,ke ; do j=js,je ; do i=is,ie
        call increment_ints(ints_sum, real_to_ints(array(i,j,k), prec_error), &
                            prec_error);
      enddo ; enddo ; enddo
    endif
    if (overflow_error) call GOLD_error(FATAL, "Overflow in reproducing_sum(_3d).")

    call sum_across_PEs(ints_sum, ni)

    sum = ints_to_real(ints_sum)
  endif   

end function reproducing_sum_3d

function real_to_ints(r, prec_error) result(ints)
  real,                      intent(in) :: r
  integer(kind=8), optional, intent(in) :: prec_error
  integer(kind=8), dimension(ni)  :: ints
  !   This subroutine converts a real number to an equivalent representation
  ! using several long integers.

  real :: rs
  character(len=80) :: mesg
  integer(kind=8) :: ival, prec_err
  integer :: sgn, i

  prec_err = prec ; if (present(prec_error)) prec_err = prec_error
  
  sgn = 1 ; if (r<0.0) sgn = -1
  rs = abs(r)
  ints(:) = 0_8

  if (.not.(rs < prec_err*pr(1))) then
    write(mesg, '(ES13.5)') r
    call GOLD_error(FATAL,"Overflow in real_to_ints conversion of "//trim(mesg))
  endif
  do i=1,ni
    ival = int(rs*I_pr(i), 8)
    rs = rs - ival*pr(i)
    ints(i) = sgn*ival
  enddo

end function real_to_ints

function ints_to_real(ints) result(r)
  integer(kind=8), dimension(ni), intent(in) :: ints
  real :: r
  ! This subroutine reverses the conversion in real_to_ints.

  integer :: i
  
  r = 0.0
  do i=1,ni ; r = r + pr(i)*ints(i) ; enddo
end function ints_to_real

subroutine increment_ints(int_sum, int2, prec_error)
  integer(kind=8), dimension(ni), intent(inout) :: int_sum
  integer(kind=8), dimension(ni), intent(in)    :: int2
  integer(kind=8), optional,      intent(in)    :: prec_error
  
  ! This subroutine increments a number with another, both using the integer
  ! representation in real_to_ints.
  integer :: i

  do i=ni,2,-1
    int_sum(i) = int_sum(i) + int2(i)
    ! Carry the local overflow.
    if (int_sum(i) > prec) then
      int_sum(i) = int_sum(i) - prec
      int_sum(i-1) = int_sum(i-1) + 1
    elseif (int_sum(i) < -prec) then
      int_sum(i) = int_sum(i) + prec
      int_sum(i-1) = int_sum(i-1) - 1
    endif
  enddo
  int_sum(1) = int_sum(1) + int2(1)
  if (present(prec_error)) then
    if (abs(int_sum(1)) > prec_error) overflow_error = .true.
  else
    if (abs(int_sum(1)) > prec) overflow_error = .true.
  endif
   
end subroutine increment_ints

subroutine increment_ints_fast(int_sum, int2)
  integer(kind=8), dimension(ni), intent(inout) :: int_sum
  integer(kind=8), dimension(ni), intent(in)    :: int2
  
  ! This subroutine increments a number with another, both using the integer
  ! representation in real_to_ints, but without doing any carrying of overflow.
  integer :: i

  do i=1,ni ; int_sum(i) = int_sum(i) + int2(i) ; enddo
   
end subroutine increment_ints_fast

subroutine carry_overflow(int_sum, prec_error)
  integer(kind=8), dimension(ni), intent(inout) :: int_sum
  integer(kind=8),                intent(in)    :: prec_error

  ! This subroutine handles carrying of the overflow.
  integer :: i, num_carry

  do i=ni,2,-1 ; if (abs(int_sum(i)) > prec) then
    num_carry = int_sum(i) * I_prec
    int_sum(i) = int_sum(i) - num_carry*prec
    int_sum(i-1) = int_sum(i-1) + num_carry
  endif ; enddo
  if (abs(int_sum(1)) > prec_error) then
    overflow_error = .true.
  endif

end subroutine carry_overflow

subroutine GOLD_infra_end
  ! This subroutine should contain all of the calls that are required
  ! to close out the infrastructure cleanly.  This should only be called
  ! in ocean-only runs, as the coupler takes care of this in coupled runs.
  call print_memuse_stats( 'Memory HiWaterMark', always=.TRUE. )
  call fms_end
end subroutine GOLD_infra_end

end module GOLD_coms
