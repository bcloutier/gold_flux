module GOLD_fit_compress
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
!*  By Robert Hallberg, January 2001                                   *
!*                                                                     *
!*    This file contains the subroutine fit_compressibility, which     *
!*  determines optimal coefficients for a function that (1) compen-    *
!*  sates for typical compressibility when in situ density is divided  *
!*  by this function, (2) is analytically integrable, and (3) is       *
!*  linear in the fitted coefficients.                                 *
!*                                                                     *
!*    The form of this function is:                                    *
!*  B(p) = 1 + a[1]*p + a[2]*p^2 + sum_n=3toN(a[n]*p*exp(b[n]*p)       *
!*  The units of p are Pa.  The units of a[1], a[3], a[4], ... are     *
!*  therefore Pa-1, while the units of a[2] are Pa-2.                  *
!*                                                                     *
!*    Each input point is assumed to have the same degree of           *
!*  accuracy, so the vertical weighting that determines the profile    *
!*  can be specified by the the location of the input points.  There   *
!*  Is no need for the input points to coincide with the locations     *
!*  of grid points, density surfaces, or anything else, or even that   *
!*  the points are in order!  Warning messages are given if negative   *
!*  pressures or pressures exceeding 10^8 Pa (10,000 dbar) are input.  *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use GOLD_error_handler, only : GOLD_error, FATAL
use GOLD_EOS, only : calculate_compress, EOS_type

implicit none ; private

public fit_compressibility

integer, parameter :: NC = 10

contains

subroutine fit_compressibility(T, S, p, npts, nterm, a, b, eqn_of_state)
  real, dimension(:), intent(in)    :: T, S, p
  integer,            intent(in)    :: npts, nterm
  real,               intent(inout) :: a(:)
  real,               intent(in)    :: b(:)
  type(EOS_type),     pointer       :: eqn_of_state
!   This subroutine finds a function B(p) such that B(0) = 1 and
! 1/B dB/dp is a good fit to 1/rho drho/dp|T,s = 1 / rho*C_s^2.
!  The a[] are coefficients for a curve of the form:
! B(p) = 1 + a[2]*p + a[3]*p^2 + sum_n(a[n]*p*exp(b[n]*p))
! The units of p are Pa.  The units of a[2] and a[n>=4] are
!  therefore Pa-1, while the units of a[3] are Pa-2.

! Arguments: T - potential temperature relative to the surface in C.
!  (in)      S - salinity in PSU.
!  (in)      p - pressure in Pa.
!  (in)      npts - the number of points to fit.
!  (in)      nterm - the number of terms in the expression.
!  (in/out)  a - the coefficients of the curve that fits the
!                compressibility, may also be the first guess.
!  (in)      b - the exponential decay scales associated with the a
!                in Pa-1.
!  (in)      eqn_of_state - The type that selects the equation of state.

  integer :: i, it, maxit
  real :: chi2, lamda, chi2_prev
  real :: comp(npts), rho(npts), sigma(npts)
  character(len=256) :: mesg

  chi2 = 1.0 ; lamda = -1.0 ; maxit = 500

  if (nterm > NC) then
    write(mesg,'("Increase NC to ",I," in GOLD_fit_compress.F90 to enable the &
       &number of coefficients being requested via fit_compressibility.")') nterm
    call GOLD_error(FATAL, "GOLD_fit_compress: "//mesg)
  endif

  do i=1,npts
    if (p(i) < 0.0) &
      write(*,*) "Warning: A negative pressure, p[] = , has been included &
                 &in the input to fit_compress."
    if (p(i) > 1.0e8) &
      write (*,*) "Warning: A pressure, p[] = , exceeding 10^8 Pa (= 10,000 dbar)&
                  &has been included in the input to fit_compress."
!     This value of sigma (10^-12 - effectively 1 m s-1 RMS sound speed
!   error) is reasonable with a smooth reference profile of T & S, but
!   may be too aggressive with a more jagged curve.
!   The pressure gradient errors due to mismatches should go as p^-2
!   - one power to compensate the p^1 in the pressure gradient term,
!   the other to compensate for the fact that the dynamically active
!   thermal with shears decrease with depth.
    sigma(i) = 1.0e-12 / (1.0 + p(i)/(5.0e6))
  enddo

!   The a are coefficients for a curve of the form:
!   B(p) = 1 + a[1]*p + a[2]*p^2 + a[3]*p*exp(a[4]*p).
!   The following numbers are a good enough first guess.
  if (a(1) == 0.0 ) then
    a(1) = 1.0                !  a(1) is always 1.
    a(2) = 4.5e-10            !  a(2) is slope in Pa-1.
    a(3) = -5e-19             !  a(3) is the curvature in Pa-2.
    do i=4,nterm+1
    ! a[m] are slopes associated with exponential functiions in Pa-1.
      a(i) = -1e-12
    enddo
  endif
!  Calculate the density and sound speed squared = drho/dp in Pa^-1.
  call calculate_compress(T,S,p,rho,comp,1,npts,eqn_of_state)
  do i=1,npts
    comp(i) = comp(i) /rho(i)
  enddo

  maxit = 500
  it = 0
  do while (it <= maxit)
    if (it==1) chi2_prev = chi2
    if (it==maxit) lamda = 0.0

    call iterate_fit(p, comp, sigma, npts, nterm, a, b, chi2, lamda)

!   If chi2 is not changing quickly, do 5 more iterations to
! polish the parameters to the extent possible.
    if ((MOD(it,4)==3) ) then
      if (chi2 >= 0.99999*chi2_prev) then
        maxit = it + 3
      else
        chi2_prev = chi2
      endif
    endif
  it = it + 1
  enddo ! end of it loop

end subroutine fit_compressibility

subroutine eval_fn(p, a, b_arr, F_fit, dF_da, na)
  real,   intent(in)  :: p
  real,   intent(in)  :: a(:),b_arr(:)
  real,   intent(out) :: F_fit
  real,   intent(out) :: dF_da(:)
  integer,intent(in)  :: na
! This subroutine evaluates F = 1/B dB/dp and dF_da at the given
! pressure and values of a.
  real :: B, I_B, I_B2, ap, bp2
  real :: dp(NC+1), exp_gp(NC+1), gp(NC+1)
  integer :: m,n

  ap  = a(2)*p ; bp2 = a(3)*p*p

  B = 1.0 + ap + bp2
  F_fit = a(2) + 2.0*a(3)*p
  do m=4,na+1
    dp(m) = a(m)*p
    gp(m) = b_arr(m)*p
    exp_gp(m) = exp(gp(m))
    B = B + dp(m)*exp(gp(m))
    F_fit = F_fit + a(m)*(1.0+gp(m))*exp_gp(m)
  enddo
  I_B = 1.0 / B ; I_B2 = I_B*I_B
  F_fit = F_fit * I_B

  dF_da(2) = 1.0 - bp2
  dF_da(3) = 2.0 + ap
  do m=4,na+1
    dF_da(2) = dF_da(2) - dp(m)*gp(m)*exp_gp(m)
    dF_da(3) = dF_da(3) + dp(m)*(1.0-gp(m))*exp_gp(m)
    dF_da(m) = ((1.0 + gp(m)) + ap*gp(m)) + bp2*(gp(m)-1.0)
    do n=4,na+1
      dF_da(m) = dF_da(m) + dp(n)*exp_gp(n)*(gp(m) - gp(n))
    enddo
    dF_da(m) = dF_da(m) * (exp_gp(m) * I_B2)
  enddo
  dF_da(2) = dF_da(2) * I_B2
  dF_da(3) = dF_da(3) * (p * I_B2)

end subroutine eval_fn


subroutine SWAP(a,b)
  real, intent(inout) :: a,b
  real :: temp

  temp=a
  a=b
  b=temp
end subroutine SWAP


subroutine iterate_fit(p, F, sigma, npts, na, a, b, chi2, lamda)
  real, intent(in) :: p(:), F(:), sigma(:), b(:)
  integer, intent(in) :: npts, na
  real, intent(inout) :: a(:), chi2, lamda

! This function does one iteration of a linearized fit optimization.
  integer :: i, j, k
  real :: a_guess(NC+1), chi2_guess, da(NC+1), covar(NC+1,NC+1)
  real :: F_fit, wt, Isigma2, dF, dF_da(NC+1)
  integer :: piv_col, piv_row, used_as_piv(NC+1)
  real :: Ipivot_val, max_val, temp, t

  if (lamda < 0.0) then
    ! This is the first iteration with this subroutine.
    lamda=0.001
    call determine_chi2(p,F,sigma,npts,na,a,b,chi2)
  endif

!   The following block sets up the symmetric linearized curvature
! matrix and sets up the right hand side for the equation that is
! subsequently solved for the changes to a.

    do j=2,na+1
      do k=2,j+1
        covar(k,j) = 0.0
      enddo
      da(j) = 0.0
    enddo
    do i=1,npts
      call eval_fn(p(i),a,b,F_fit,dF_da,na)

      Isigma2 = 1.0 / (sigma(i)*sigma(i))
      dF = F(i) - F_fit
      do j=2,na+1
        wt = dF_da(j) * Isigma2
        do k=2,j!+1
          covar(k,j) = covar(k,j) + wt * dF_da(k)
        enddo
        da(j) = da(j) + dF * wt
      enddo
    enddo
    do j=3,na+1
      do k=2,j
        covar(j,k) = covar(k,j)
      enddo
    enddo
    do j=2,na+1
      covar(j,j) = covar(j,j) * (1.0+(lamda))
    enddo


! The following block does Gaussian elimination on the upper-left
! corner of covar, starting from 1,1, and determines da[] in place.


    do j=2,na+1
      used_as_piv(j)=0
    enddo
    do i=2,na+1
      ! Find the largest remaining value to use as the pivot.
      max_val=0.0
      do j=2,na+1
      if (used_as_piv(j) == 0) then
        do k=2,na+1
          if (used_as_piv(k) == 0) then
            if (ABS(covar(k,j)) >= max_val) then
              max_val = ABS(covar(k,j))
              piv_row = j ; piv_col = k
            endif
          endif
        enddo
      endif
      enddo
      used_as_piv(piv_col) = 1
      if (covar(piv_col,piv_row) == 0.0) then
        ! This should never happen for a reasonable function.
        call GOLD_error(FATAL, "GOLD_fit_compress: Fit Compress Singular Matrix.")
      endif

      ! Rearrange rows to put the pivot value on the diagonal.
      if (piv_row /= piv_col) then
        do k=2,na+1
          t = covar(k,piv_row) ; covar(k,piv_row) = covar(k,piv_col)
          covar(k,piv_col) = t
        enddo
        t = da(piv_row) ; da(piv_row) = da(piv_col) ; da(piv_col) = t
      endif

      ! Do one round of Gaussian elimination on the matrix.
      Ipivot_val=1.0/covar(piv_col,piv_col)
      covar(piv_col,piv_col)=1.0
      do k=2,na+1
        covar(k,piv_col) = covar(k,piv_col) * Ipivot_val
      enddo
      da(piv_col) = da(piv_col) * Ipivot_val
      do j=2,na+1
      if (j /= piv_col) then
        temp = covar(piv_col,j)
        covar(piv_col,j) = 0.0
        do k=2,na+1
          covar(k,j) = covar(k,j) - covar(k,piv_col)*temp
        enddo
        da(j) = da(j) - da(piv_col)*temp
      endif
      enddo
    enddo


  if (lamda == 0.0) return

  do j=2,na+1
    a_guess(j) = a(j) + da(j)
  enddo
  call determine_chi2(p,F,sigma,npts,na,a_guess,b,chi2_guess)
  if (chi2_guess < chi2) then
    ! This is a better estimate, so keep it.
    lamda = lamda * 0.1
    do j=2,na+1
      a(j) = a_guess(j)
    enddo
    chi2 = chi2_guess
  else
    ! Use a smaller step next time.
    lamda = lamda * 10.0
  endif

end subroutine iterate_fit

subroutine determine_chi2(p, F, sigma, npts, na, a, b, chi2)
real,    intent(in)  :: p(:), F(:), sigma(:), a(:), b(:)
integer, intent(in)  :: npts, na
real,    intent(out) :: chi2
! This subroutine evaluates chi squared.
  integer :: i
  real :: F_fit, dF_da(NC+1)

  chi2=0.0
  do i=1,npts
    call eval_fn(p(i),a,b,F_fit,dF_da,na)
    chi2 = chi2 +(F(i)-F_fit) * (F(i)-F_fit) / (sigma(i)*sigma(i))
  enddo

end subroutine determine_chi2

end module GOLD_fit_compress
