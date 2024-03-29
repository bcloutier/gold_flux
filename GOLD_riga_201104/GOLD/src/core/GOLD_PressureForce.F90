module GOLD_PressureForce
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
!*  By Robert Hallberg, April 1994 - June 2008                         *
!*                                                                     *
!*    This file contains the subroutine that directs the model to use  *
!*  the code that determines the horizontal accelerations due to       *
!*  pressure gradients.  The two options currently available are a     *
!*  traditional Montgomery potential form, and the analytic finite     *
!*  volume form described in Adcroft, Hallberg and Harrison, 2008,     *
!*  Ocean Modelling, 22, 106-113.                                      *
!*                                                                     *
!*    PressureForce takes 9 arguments, which are described below. If   *
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

use GOLD_CompressComp, only : Compress_CS
use GOLD_diag_mediator, only : diag_ptrs, time_type
use GOLD_error_handler, only : GOLD_error, GOLD_mesg, FATAL, WARNING, is_root_pe
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_grid, only : ocean_grid_type
use GOLD_PressureForce_AFV, only : PressureForce_AFV_Bouss, PressureForce_AFV_nonBouss
use GOLD_PressureForce_AFV, only : PressureForce_AFV_init, PressureForce_AFV_CS
use GOLD_PressureForce_Mont, only : PressureForce_Mont_Bouss, PressureForce_Mont_nonBouss
use GOLD_PressureForce_Mont, only : PressureForce_Mont_init, PressureForce_Mont_CS
use GOLD_tidal_forcing, only : calc_tidal_forcing, tidal_forcing_CS
use GOLD_variables, only : thermo_var_ptrs
implicit none ; private

#include <GOLD_memory.h>

public PressureForce, PressureForce_init, PressureForce_end

type, public :: PressureForce_CS ; private
  logical :: Analytic_FV_PGF ! If true, use the analytic finite volume form
                            ! (Adcroft et al., Ocean Mod. 2008) of the PGF.
  type(PressureForce_AFV_CS), pointer :: PressureForce_AFV_CSp => NULL()
  type(PressureForce_Mont_CS), pointer :: PressureForce_Mont_CSp => NULL()
end type PressureForce_CS

contains

subroutine PressureForce(h, tv, PFu, PFv, G, CS, p_atm, pbce, eta)
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in)   :: h
  type(thermo_var_ptrs), intent(inout)             :: tv
  real, dimension(NXMEMQ_,NYMEM_,NZ_), intent(out) :: PFu
  real, dimension(NXMEM_,NYMEMQ_,NZ_), intent(out) :: PFv
  type(ocean_grid_type),               intent(in)  :: G
  type(PressureForce_CS),              pointer     :: CS
  real, dimension(:,:),               optional, pointer     :: p_atm
  real, dimension(NXMEM_,NYMEM_,NZ_), optional, intent(out) :: pbce
  real, dimension(NXMEM_,NYMEM_),     optional, intent(out) :: eta

!    This subroutine works as a temporary interface between the model and the
! Boussinesq and non-Boussinesq pressure force routines.
! Descriptions of the variables are in each of the routines called in the
! following conditional block.

  if (CS%Analytic_FV_PGF) then
    if (G%Boussinesq) then
      call PressureForce_AFV_Bouss(h, tv, PFu, PFv, G, CS%PressureForce_AFV_CSp, p_atm, pbce, eta)
    else
      call PressureForce_AFV_nonBouss(h, tv, PFu, PFv, G, CS%PressureForce_AFV_CSp, p_atm, pbce, eta)
    endif
  else
    if (G%Boussinesq) then
      call PressureForce_Mont_Bouss(h, tv, PFu, PFv, G, CS%PressureForce_Mont_CSp, p_atm, pbce, eta)
    else
      call PressureForce_Mont_nonBouss(h, tv, PFu, PFv, G, CS%PressureForce_Mont_CSp, p_atm, pbce, eta)
    endif
  endif

end subroutine Pressureforce

subroutine PressureForce_init(Time, G, param_file, diag, CS, compress_CSp, &
                              tides_CSp)
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ptrs), target, intent(inout) :: diag
  type(PressureForce_CS),  pointer       :: CS
  type(Compress_CS), optional, pointer   :: compress_CSp
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
  character(len=128) :: mesg  ! A string for writing an output error message.
  character(len=128) :: version = '$Id: GOLD_PressureForce.F90,v 13.0.2.10.2.11 2010/08/17 15:06:19 rwh Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "GOLD_PressureForce" ! This module's name.

  if (associated(CS)) then
    call GOLD_error(WARNING, "PressureForce_init called with an associated "// &
                            "control structure.")
    return
  else ; allocate(CS) ; endif

  CS%Analytic_FV_PGF = .false.
  call read_param(param_file,"ANALYTIC_FV_PGF",CS%Analytic_FV_PGF)

  ! Write all relevant parameters to the model log.
  call log_version(param_file, mod, version, tagname)
  call log_param(param_file, mod, "ANALYTIC_FV_PGF", CS%Analytic_FV_PGF, &
                 "If true the pressure gradient forces are calculated \n"//&
                 "with a finite volume form that analytically integrates \n"//&
                 "the equations of state in pressure to avoid any \n"//&
                 "possibility of numerical thermobaric instability, as \n"//&
                 "described in Adcroft et al., O. Mod. (2008).  The \n"//&
                 "default should be changed to true.", default=.false.)

  if (CS%Analytic_FV_PGF) then
    call PressureForce_AFV_init(Time, G, param_file, diag, &
             CS%PressureForce_AFV_CSp, tides_CSp)
  else
    call PressureForce_Mont_init(Time, G, param_file, diag, &
             CS%PressureForce_Mont_CSp, compress_CSp, tides_CSp)
  endif

end subroutine PressureForce_init


subroutine PressureForce_end(CS)
  type(PressureForce_CS), pointer :: CS
  if (associated(CS)) deallocate(CS)
end subroutine PressureForce_end

end module GOLD_PressureForce
