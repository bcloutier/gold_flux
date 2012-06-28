module GOLD_intrinsic_functions

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*   This module holds intrinsic functions which are not.              *
!*   supported by some compilers.                                      *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

  implicit none
  private

  public :: invcosh

  contains

  function invcosh(x)
    real, intent(in) :: x
    real :: invcosh
  
#ifdef __INTEL_COMPILER
    invcosh=acosh(x)
#else
    invcosh=log(x+sqrt(x*x-1))
#endif

  end function invcosh

end module GOLD_intrinsic_functions
