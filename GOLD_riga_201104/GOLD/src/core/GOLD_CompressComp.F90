module GOLD_CompressComp
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
!*  By Robert Hallberg, July 2003                                      *
!*                                                                     *
!*    This file contains the subroutines that compensate for the       *
!*  compressibility of seawater.  read_compress or register_compress   *
!*  are used to specify the profiles whose compressibilties will be    *
!*  fit to an analytically tractible functional form and used as       *
!*  reference compressibilities.  uncompress_e_rho returns the         *
!*  compressibility compensated density and height, the non-Boussinesq *
!*  equivalent of which is documented in Hallberg (Ocean Mod., sub.).  *
!*  Grad_e_zstar calculates the gradients on geopotentials of the      *
!*  compensated heights - these are only nonzero if horizontally       *
!*  varying reference compressibilities are used.                      *
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

use GOLD_domains, only : pass_var, min_across_PEs
use GOLD_error_handler, only : GOLD_error, FATAL, WARNING, NOTE, is_root_pe
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_fit_compress, only : fit_compressibility
use GOLD_grid, only : ocean_grid_type
use GOLD_io, only : close_file, create_file, field_size, fieldtype
use GOLD_io, only : file_exists, read_data, SINGLE_FILE, vardesc, write_field
use GOLD_variables, only : thermo_var_ptrs
use GOLD_EOS, only : calculate_density, select_eqn_of_state, EOS_type

implicit none ; private

#include <GOLD_memory.h>

public read_compress, register_compress, compress_init
public uncompress_e_rho, Grad_z_estar, uncompress_p_alpha, Grad_p_pstar

integer, parameter :: NTERM = 6

#ifdef STATIC_MEMORY
# define NTERM_ NTERM
#else
# define NTERM_ :
#endif

type, public :: Compress_CS ; private
  real PTR_ :: cz(NXMEM_,NYMEM_,NTERM_) ! The coefficients of the depth-space
                             ! terms, used in Boussinesq models.
  real :: cr_cz(NTERM)       !   The ratio of the coefficients of each term
                             ! for the fit to B(z) to the fit for B'(z).
  real :: ez(NTERM)          ! The e-folding scales of the depth dependence,
                             ! in m-1.
  real PTR_ :: cp(NXMEM_,NYMEM_,NTERM_) ! The coefficients of the pressure-space
                             ! terms, used in non-Boussinesq models.
  real :: ca_cp(NTERM)       !   The ratio of the coefficients of each term
                             ! for the fit to B(p) to the fit for B'(p)
  real :: ep(NTERM)          ! The e-folding scales of the pressure dependence,
                             ! in Pa-1.
  logical :: Boussinesq_set = .false.
  logical :: nonBoussinesq_set = .false.
  type(EOS_type), pointer :: eqn_of_state => NULL() ! Type that indicates the
                                                    ! equation of state to use.
end type Compress_CS

contains

subroutine read_compress(profile_filename, T_name, S_name, Z_name, distance, D_0, &
    write_compress, compressibility_filename, G, CS)
  character(len=*), intent(in) :: profile_filename, T_name, S_name, Z_name
  real, intent(in) :: distance, D_0
  logical, intent(in) ::  write_compress
  character(len=*), intent(in) :: compressibility_filename
  type(ocean_grid_type), intent(inout) :: G
  type(Compress_CS), pointer   :: CS
!  This subroutine determines the spatially varying coefficients of
! the function that compensates for the compressibility of sea water.
! This function is used in uncompress_e_rho.  The coefficients
! themselves are found by fitting individual profiles of (T,S,z) from
! profile_filename in fit_compressibility and then smoothing them.  The data
! in profile_filename must have the same _horizontal_ resolution as the model
! run, but will be a depth-space dataset.

! Arguments: profile_filename - the path to the file in which the profiles
!                               that are to be fit are found.
!  (in)      T_name - variable name in profile_filename for potential
!                     temperature relative to the surface in C.
!  (in)      S_name - variable name in profile_filename for salinity in PSU.
!  (in)      Z_name - variable name in profile_filename for depths in m.
!  (in)      distance - the horizontal distance over which to smooth
!                       the profiles when the seafloor is at D_0.
!  (in)      D_0 - the depth at which the specified smoothing
!                  distance is realized.  In smooth_cz, the smoothing
!                  often has fluxes proportional to (D/D_0)^n.
!  (in)      write_compress - if .true., the coefficients of the fitted
!                  compressibilities are written to compressibility_filename.
!  (in)      compressibility_filename - the path to the file to which the
!                  fitted compressibilities are written and read.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure that has been set by a
!                 previous call to compress_init.
  integer :: ncid, k, status, npts
  real, allocatable :: T_c(:,:,:), S_c(:,:,:), Z(:)
  integer :: hard_fail
  real :: a(NTERM)    ! a[] are the coefficients for a curve:
  real :: esp(NTERM)  ! esp is the exponential decay scale in Pa-1.
  !  B(p) = 1 + a[1]*p + a[2]*p^2 + a[3]*p*exp(a[4]*p).
  !  The units of p are Pa.  The units of a[1], a[3], and a[4] are
  !  therefore Pa-1, while the units of a[2] are Pa-2.
  real, pointer :: p(:), T(:), S(:)
  real :: Rho0xG
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed, m

  integer :: fieldsize(4), reread
  logical :: used

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) call GOLD_error(FATAL,"uncompress_e_rho : "// &
       "Module GOLD_CompressComp must be initialized before it is used.")
  Rho0xG = G%Rho0 * G%g_Earth
  reread = reread_compress(compressibility_filename, G, CS)
  if ( reread == 0 ) return

  fieldsize = 0
  call field_size(profile_filename, T_name, fieldsize)
  npts = fieldsize(3)-1

  allocate(T_c(isd:ied,jsd:jed,npts+1) ) ; T_c = 0.0
  allocate(S_c(isd:ied,jsd:jed,npts+1) ) ; S_c = 0.0
  allocate(Z(npts+1))

  call read_data(profile_filename,T_name, T_c(:,:,:), domain=G%Domain%mpp_domain)
  call read_data(profile_filename,S_name, S_c(:,:,:), domain=G%Domain%mpp_domain)
  call read_data(profile_filename,Z_name, Z(:))

!    B(p) = 1 + a[1]*p + a[2]*p^2 + a[3]*p*exp(a[4]*p).
!    The units of p are Pa.  The units of a[1], a[3], and a[4] are
!    therefore Pa-1, while the units of a[2] are Pa-2.

  allocate(p(npts)) ; allocate(T(npts)) ; allocate(S(npts))

! Convert the input depths to pressures.
  if (Z(npts) > 0) then
    do i=1,npts
      p(i) = Rho0xG*Z(i)
    enddo
  else
    do i=1,npts
      p(i) = -1.0*Rho0xG*Z(i)
    enddo
  endif

  do m=1,NTERM
    if (m< 4) then
      esp(m) = 0.0
    else
      esp(m) = (-1.0/(5.0e6* real(m-3)))
    endif
  enddo

  if (G%Boussinesq) then
    CS%cr_cz(1) = 1.0 ; CS%cr_cz(2) = 2.0 ; CS%cr_cz(3) = 3.0
    do m=4,NTERM
      CS%ez(m) = -1.0*esp(m)*Rho0xG
      CS%cr_cz(m) = CS%ez(m)*CS%ez(m)
    enddo

    CS%cz(:,:,1) = 1.0
   else
    CS%ca_cp(1) = 1.0 ; CS%ca_cp(2) = 2.0 ; CS%ca_cp(3) = 3.0
    do m=4,NTERM
      CS%ep(m) = esp(m)
      CS%ca_cp(m) = CS%ep(m)*CS%ep(m)
    enddo

    CS%cp(:,:,1) = 1.0
  endif

  do j=js,je+1 ; do i=is,ie+1
!    do j=1,je-js+2 ; do i=1,ie-is+2
    do m=1,NTERM
      a(m) = 0.0
    enddo
    do k=1,npts
      T(k) = T_c(i,j,k) ; S(k) = S_c(i,j,k)
    enddo

    call fit_compressibility(T, S, p, npts, NTERM-1, a, esp, CS%eqn_of_state)


    if (G%Boussinesq) then
  ! Convert the values of a to the coefficients of curves
  ! B(z) =  z + cz[2]*z^2 + cz[3]*z^3
  !         + sum_m ( cz[m]*(1+(ez[m]*z - 1)*exp(ez[m]*z)) ).
  ! B'(z) =  1 + cr[2]*z + cr[3]*z^2 + sum_m cr[m]*z*exp(ez[m]*z).

      CS%cz(i,j,1) = 1.0
      CS%cz(i,j,2) = -0.5*a(2)*Rho0xG
      CS%cz(i,j,3) = a(3)*Rho0xG*Rho0xG/3.0
      do m=4,NTERM
        CS%cz(i,j,m) = -1.0*a(m)*Rho0xG / (CS%ez(m)*CS%ez(m))
      enddo
    else
  ! Convert the values of a to the coefficients of curves
  ! B(p) =  p + cp[2]*cp^2 + cp[3]*p^3
  !         + sum_m ( cp[m]*(1+(ep[m]*p-1)*exp(ep[m]*p)) ).
  ! B'(p) =  1 + cr[2]*p + cr[3]*p^2 + sum_m cr[m]*p*exp(ep[m]*p).
      CS%cp(i,j,1) = 1.0
      CS%cp(i,j,2) = -0.5*a(2)
      CS%cp(i,j,3) = -a(3)/3.0
      do m=4,NTERM
        CS%cp(i,j,m) = -a(m) / (CS%ep(m)*CS%ep(m))
      enddo
    endif  
  enddo ; enddo

  if (G%Boussinesq) then
    call smooth_cz(CS%cz, distance, D_0, G)
    do m=1,NTERM ; call pass_var(CS%cz(:,:,m), G%Domain) ; enddo
  else
    call smooth_cz(CS%cp, distance, D_0, G)
    do m=1,NTERM ; call pass_var(CS%cp(:,:,m), G%Domain) ; enddo
  endif

  deallocate(Z) ; deallocate(T_c) ; deallocate(S_c)
  deallocate(p) ; deallocate(T) ; deallocate(S)

  if (write_compress) call write_compress_file(compressibility_filename, G, CS)

end subroutine read_compress

subroutine write_compress_file(filepath, G, CS)
  character(len=*), intent(in) :: filepath
  type(ocean_grid_type), intent(inout) :: G
  type(Compress_CS), pointer   :: CS
! This subroutine writes out the compressibilities for this run.

! Arguments: filepath - the path to the file to which the fitted
!                      compressibilities are written.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure that has been set by a
!                 previous call to compress_init.

  type(vardesc) :: vars(2*NTERM-4)
  type(fieldtype) :: fields(2*NTERM-4)
  character(len=1) :: num
  integer :: unit, m

  if (G%Boussinesq) then
    vars(1) = vardesc("cz1","Coefficient 1",'h','1','1',"meter-1", 'd')
    vars(2) = vardesc("cz2","Coefficient 2",'h','1','1',"meter-2", 'd')
    do m=4,NTERM
      write (num,'(I1)') m-1
      vars(2*m-5) = vardesc("cz"//num,"Coefficient "//num,'h','1','1',"meter",'d')
      vars(2*m-4) = vardesc("ez"//num,"Exponential "//num,'1','1','1',"meter-1",'d')
    enddo
  else
    vars(1) = vardesc("cp1","Coefficient 1",'h','1','1',"Pa-1", 'd')
    vars(2) = vardesc("cp2","Coefficient 2",'h','1','1',"Pa-2", 'd')
    do m=4,NTERM
      write (num,'(I1)') m-1
      vars(2*m-5) = vardesc("cp"//num,"Coefficient "//num,'h','1','1',"Pa",'d')
      vars(2*m-4) = vardesc("ep"//num,"Exponential "//num,'1','1','1',"Pa-1",'d')
    enddo
  endif

  call create_file(unit, trim(filepath), vars, (2*NTERM-4), G, fields, SINGLE_FILE)

  if (G%Boussinesq) then
    call write_field(unit, fields(1), G%Domain%mpp_domain, CS%cz(:,:,2))
    call write_field(unit, fields(2), G%Domain%mpp_domain, CS%cz(:,:,3))
    do m=4,NTERM
      call write_field(unit, fields(2*m-5), G%Domain%mpp_domain, CS%cz(:,:,m))
      call write_field(unit, fields(2*m-4), CS%ez(m))
    enddo
  else
    call write_field(unit, fields(1), G%Domain%mpp_domain, CS%cp(:,:,2))
    call write_field(unit, fields(2), G%Domain%mpp_domain, CS%cp(:,:,3))
    do m=4,NTERM
      call write_field(unit, fields(2*m-5), G%Domain%mpp_domain, CS%cp(:,:,m))
      call write_field(unit, fields(2*m-4), CS%ep(m))
    enddo
  endif

  call close_file(unit)

end subroutine write_compress_file


function reread_compress(filepath, G, CS)
  character(len=*),      intent(in)    :: filepath
  type(ocean_grid_type), intent(inout) :: G
  type(Compress_CS),     pointer       :: CS
  integer  ::  reread_compress
!  This subroutine rereads the compressibilities from an existing file.
! Arguments: filepath - the path to the file frow which the fitted
!                       compressibilities are read, if filepath exists.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure that has been set by a
!                 previous call to compress_init.
  integer :: i, j, is, ie, js, je, m, err, cdfid
  character(len=48)  :: name          ! The variable name to be read.

  if (.not. file_exists(filepath) ) then
!    if (is_root_pe()) &
!      write(*,*) 'GOLD Failed to find compressiblity file ',filepath
    reread_compress = -1
    return
  endif
!  if (is_root_pe()) &
!    write(*,*) 'GOLD Reading from compressiblity file ',filepath

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (G%Boussinesq) then
    CS%cz(:,:,1) = 1.0 ; CS%cz(:,:,2:size(CS%cz,3)) = 0.0

    call read_data(filepath,"cz1", CS%cz(is:ie,js:je,2), domain=G%Domain%mpp_domain)
    call read_data(filepath,"cz2", CS%cz(is:ie,js:je,3), domain=G%Domain%mpp_domain)
    do m=4,NTERM
      write(name,'("cz",I1)') m-1
      call read_data(filepath,name, CS%cz(is:ie,js:je,m), domain=G%Domain%mpp_domain)
      write(name,'("ez",I1)') m-1
      call read_data(filepath,name, CS%ez(m))
    enddo

    do m=4,NTERM ; CS%cr_cz(m) = (CS%ez(m)*CS%ez(m)) ; enddo

    do m=1,NTERM ; call pass_var(CS%cz(:,:,m), G%Domain) ; enddo
  else
    CS%cp(:,:,1) = 1.0 ; CS%cp(:,:,2:size(CS%cp,3)) = 0.0

    call read_data(filepath,"cp1", CS%cp(is:ie,js:je,2), domain=G%Domain%mpp_domain)
    call read_data(filepath,"cp2", CS%cp(is:ie,js:je,3), domain=G%Domain%mpp_domain)
    do m=4,NTERM
      write(name,'("cp",I1)') m-1
      call read_data(filepath,name, CS%cp(is:ie,js:je,m), domain=G%Domain%mpp_domain)
      write(name,'("ep",I1)') m-1
      call read_data(filepath,name, CS%ep(m))
    enddo

    do m=4,NTERM ; CS%ca_cp(m) = (CS%ep(m)*CS%ep(m)) ; enddo

    do m=1,NTERM ; call pass_var(CS%cp(:,:,m), G%Domain) ; enddo
  endif

  reread_compress = 0
end function reread_compress


subroutine smooth_cz(cz, len, D_0, G)
  real,                  intent(inout) :: cz(NXMEM_,NYMEM_,NTERM_)
  real,                  intent(in)    :: len, D_0
  type(ocean_grid_type), intent(inout) :: G

  real, dimension(SZ1_(cz),SZ2_(cz))   :: ad_x, ad_y, I_D
  real, dimension(SZ1_(cz),SZ2_(cz),2) :: ct
  real :: max_dt
  real :: A
  integer :: i, j, is, ie, js, je, m, max_it, it, l1, l2
  real :: D_x, D_y
  real :: I_max_dt_here

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  l2 = 0
  A = len*len/(D_0*D_0)
  ct(:,:,:) = 0.0

  do j=js,je ; do i=is-1,ie
    D_x = 0.5*(G%D(i,j)+G%D(i+1,j))
    ad_x(i,j) = A*G%dy_u(i,j)*G%IDXu(i,j)*(G%D(i,j)*G%D(i+1,j))*D_x
  enddo ; enddo
  do j=js-1,je ; do i=is,ie
    D_y = 0.5*(G%D(i,j)+G%D(i,j+1))
    ad_y(i,j) = A*G%dx_v(i,j)*G%IDYv(i,j)*(G%D(i,j+1)*G%D(i,j))*D_y
  enddo ; enddo
  do j=js,je ; do i=is,ie
    if (G%D(i,j) > G%Angstrom) then
      I_D(i,j) =  (G%hmask(i,j)/(G%DXDYh(i,j)*G%D(i,j)))
    else
      I_D(i,j) = 0.0
    endif
  enddo ; enddo

! Find the maximum time step, take 95% of this.
  max_dt = 1.0e20
  do j=js,je ; do i=is,ie
!   Find the largest timestep that is nonoscilatory for diffusing the fields.
    I_max_dt_here = I_D(i,j) * (ad_x(i-1,j) + ad_x(i,j) + ad_y(i,j-1) + ad_y(i,j))
    if (0.95 < max_dt*I_max_dt_here) max_dt = 0.95 / I_max_dt_here
  enddo ; enddo
  call min_across_PEs(max_dt)
  max_it = int(ceiling(25.0/max_dt))
  max_dt = 25.0 / real( max_it)

  do j=js,je ; do i=is-1,ie
    ad_x(i,j) = ad_x(i,j) * max_dt
  enddo ; enddo
  do j=js-1,je ; do i=is,ie
    ad_y(i,j) = ad_y(i,j) * max_dt
  enddo ; enddo

  do m=2,NTERM
    do j=js,je ; do i=is,ie
      ct(i,j,1) = cz(i,j,m)
    enddo ; enddo

    do it=0,max_it-1
      l1 = MOD(it,2)+1 ; l2 = 3-l1
      call pass_var(ct(:,:,l1), G%Domain)
      do j=js,je ; do i=is,ie
        ct(i,j,l2) = ct(i,j, l1) + &
            I_D(i,j)*((ad_x(i,j)*(ct(i+1,j,l1)-ct(i,j,l1)) - &
                       ad_x(i-1,j)*(ct(i,j,l1)-ct(i-1,j,l1))) + &
                       (ad_y(i,j)*(ct(i,j+1,l1)-ct(i,j,l1)) - &
                       ad_y(i,j-1)*(ct(i,j,l1)-ct(i,j-1,l1)))) + &
                       max_dt*(cz(i,j,m) - ct(i,j,l1))
      enddo ; enddo
    enddo
    do j=js,je ; do i=is-1,ie
      cz(i,j,m) = ct(i,j,l2)
    enddo ; enddo
  enddo
end subroutine smooth_cz

subroutine register_compress(T, S, e, npts, G, CS)
  real, dimension(:) :: T, S, e
  integer :: npts
  type(ocean_grid_type), intent(inout) :: G
  type(Compress_CS),     pointer    :: CS
!  This subroutine determines the coefficients of a curve that
! compensates for the compressibility of sea water.  This curve is
! used in uncompress_e_rho.  The coefficients themselves are
! a fit to the profiles that are passed into this subroutine as
! arguments, and are determined in fit_compressibility.
!
! Arguments: T - potential temperature relative to the surface in C.
!  (in)      S - salinity in PSU.
!  (in)      e - interface height (negative downward) in m.
!  (in)      npts - the number of points to fit.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure that has been set by a
!                 previous call to compress_init.

  real :: a(NTERM)  !a[] are the coefficients for a curve:
  real :: esp(NTERM) ! esp is the exponential decay scale in Pa-1.
!  B(p) = 1 + a[1]*p + a[2]*p^2 + a[3]*p*exp(a[4]*p).
!  The units of p are Pa.  The units of a[1], a[3], and a[4] are
!  therefore Pa-1, while the units of a[2] are Pa-2.
  real :: p(npts)
  real :: cz_0(NTERM), cp_0(NTERM)
  real :: Rho0xG
  integer :: i, j, is, ie, js, je, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (.not.associated(CS)) call GOLD_error(FATAL,"register_compress: "// &
       "Module GOLD_CompressComp must be initialized before it is used.")
  Rho0xG = G%Rho0 * G%g_Earth

! Convert the input depths to pressures.
  do i=1,npts
    p(i) = -1.0*Rho0xG*e(i)
  enddo

  do m=1,NTERM
    a(m) = 0.0
    if (m<4) then
      esp(m) = 0.0
    else
      esp(m) = (-1.0/(5.0e6*real(m-3)))
    endif
  enddo

  call fit_compressibility(T, S, p, npts, NTERM-1, a, esp, CS%eqn_of_state)

! Convert the values of a to the coefficients of curves
! B(z) =  z + cz[1]*z^2 + cz[1]*z^3
!         + sum_m ( cz[m]*(1+(ez[m]*z-1)*exp(ez[m]*z)) ).
! B'(z) =  1 + cr[1]*z + cr[1]*z^2 + sum_m cr[m]*z*exp(ez[m]*z).

  if (G%Boussinesq) then
    CS%cr_cz(1) = 1.0 ; CS%cr_cz(2) = 2.0 ; CS%cr_cz(3) = 3.0
    do m=4,NTERM
      CS%ez(m) = -1.0*esp(m)*Rho0xG
      CS%cr_cz(m) = (CS%ez(m)*CS%ez(m))
    enddo

    cz_0(1) = 1.0 ; cz_0(2) = -0.5*a(2)*Rho0xG ; cz_0(3) = a(3)*Rho0xG**2/3.0
    do m=4,NTERM ; cz_0(m) = -1.0*a(m)*Rho0xG / (CS%ez(m)*CS%ez(m)) ; enddo

    CS%cz(:,:,1) = 1.0
    do m=1,NTERM ; do j=js,je+1 ; do i=is,ie+1
      CS%cz(i,j,m) = cz_0(m)
    enddo ; enddo ; enddo

    do m=1,NTERM ; call pass_var(CS%cz(:,:,m), G%Domain) ; enddo
  else
    CS%ca_cp(1) = 1.0 ; CS%ca_cp(2) = 2.0 ; CS%ca_cp(3) = 3.0
    do m=4,NTERM
      CS%ep(m) = esp(m)
      CS%ca_cp(m) = (CS%ep(m)*CS%ep(m))
    enddo

    ! The negative signs are here because a(z) are chosen to fit
    ! B'(z) ~ 1/rho drho/dp, and 1/rho drho/dp = -1/alpha dalpha/dp.
    cp_0(1) = 1.0 ; cp_0(2) = -0.5*a(2) ; cp_0(3) = -a(3)/3.0
    do m=4,NTERM ; cp_0(m) = -a(m) / (CS%ep(m)*CS%ep(m)) ; enddo

    CS%cp(:,:,1) = 1.0
    do m=1,NTERM ; do j=js,je+1 ; do i=is,ie+1
      CS%cp(i,j,m) = cp_0(m)
    enddo ; enddo ; enddo

    do m=1,NTERM ; call pass_var(CS%cp(:,:,m), G%Domain) ; enddo
  endif

end subroutine register_compress

subroutine uncompress_e_rho(e, tv, rho_star, G, CS, rho_star_bot)
  real, intent(inout) :: e(NXMEM_,NYMEM_,NZp1_)
  type(thermo_var_ptrs), intent(in) :: tv
  real, intent(out) :: rho_star(NXMEM_,NYMEM_,NZ_)
  type(ocean_grid_type), intent(in) :: G
  type(Compress_CS),     pointer    :: CS
  real, intent(out), optional :: rho_star_bot(NXMEM_,NYMEM_,NZ_)
!   This subroutine calculates compression compensated interface
! heights and in situ densities based upon the reference profile sent
! to register_compress.  To the extent that the reference T/S profile
! matches the local profile and the fit is good, the compensated
! density (rho_star) is approximately a neutral density.
!
! Arguments: e - interface height in m, usually negative, unadjusted
!                on input but compression adjusted on output.(in/out)
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (out)     rho_star - compression adjusted in situ density times
!                       (G_Earth/Rho0) evaluated at the pressure at
!                       the top of each layer, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure that has been set by a
!                 previous call to compress_init.
!  (out)     rho_star_bot - if present, this will be the compression adjusted
!                       in situ density times (G_Earth/Rho0) evaluated at the
!                       pressure at the bottom of each layer, in m s-2.
!  (out)     rho_in_situ - In-situ density at the top of a layer.
!  (in)      bouss - flag for Boussinesq or non-Boussinesq.

  real :: rho_in_situ(SZ1_(e)) !In-situ density at the top of a layer.
  real :: rho_is_bot(SZ1_(e))  !In-situ density at the layer's bottom.
  real :: press(SZ1_(e))       !Pressure at the top of a layer in Pa.
  real :: z               !Short hand for e, in m, negative downward.
  real :: G_rho0
  real :: Rho0xG
  real :: B1, IB1, exp_gam_z
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, m, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq

  if (.not.associated(CS)) call GOLD_error(FATAL, "uncompress_e_rho: "// &
       "Module GOLD_CompressComp must be initialized before it is used.")
  G_rho0 = G%g_Earth / G%Rho0
  Rho0xG = G%Rho0 * G%g_Earth
  do k=1,nz+1 ; do j=Jsq,Jeq+1
    do i=Isq,Ieq+1
      press(i) = -Rho0xG*e(i,j,k)
    enddo
    if (k<=nz) &
      call calculate_density(tv%T(:,j,k),tv%S(:,j,k),press, &
                             rho_in_situ,Isq,Ieq-Isq+2,tv%eqn_of_state)
    if (k>1 .and. (PRESENT(rho_star_bot)) ) &
      call calculate_density(tv%T(:,j,k-1),tv%S(:,j,k-1),press, &
                             rho_is_bot,Isq,Ieq-Isq+2,tv%eqn_of_state)
    
    do i=Isq,Ieq+1
      z = e(i,j,k)
      B1 = 2.0*CS%cz(i,j,2) + 3.0*CS%cz(i,j,3)*z
      e(i,j,k) = z + (CS%cz(i,j,2) + CS%cz(i,j,3)*z)*z*z
      do m=4,NTERM
        exp_gam_z = exp(CS%ez(m)*z)
        e(i,j,k) = e(i,j,k) + CS%cz(i,j,m)*(1.0 + exp_gam_z*(CS%ez(m)*z - 1.0))
        B1 = B1 + CS%cr_cz(m)*CS%cz(i,j,m)*exp_gam_z
      enddo
      IB1 = G_rho0 / (1.0 + z*B1)  
      if (k<=nz) rho_star(i,j,k) = IB1 * rho_in_situ(i)
      if ((k>1) .and. PRESENT(rho_star_bot)) &
        rho_star_bot(i,j,k-1) = IB1 * rho_is_bot(i)
    enddo
  enddo ; enddo

end subroutine uncompress_e_rho

subroutine uncompress_p_alpha(p, tv, alpha_star, G, CS, alpha_star_bot)
  real, dimension(NXMEM_,NYMEM_,NZp1_), intent(inout) :: p
  type(thermo_var_ptrs),                intent(in)  :: tv
  real, dimension(NXMEM_,NYMEM_,NZ_),   intent(out) :: alpha_star
  type(ocean_grid_type),                intent(in)  :: G
  type(Compress_CS),                    pointer     :: CS
  real, dimension(NXMEM_,NYMEM_,NZ_), optional, intent(out) :: alpha_star_bot
!   This subroutine calculates compression compensated interface
! heights and in situ densities based upon the reference profile sent
! to register_compress.  To the extent that the reference T/S profile
! matches the local profile and the fit is good, the compensated
! density (rho_star) is approximately a neutral density.
!
! Arguments: p - interface pressure in Pa, usually positive, unadjusted
!                on input but compression adjusted on output.(in/out)
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (out)     alpha_star - compression adjusted in situ specific volume
!                         evaluated at the pressure at the top
!                         of each layer, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure that has been set by a
!                 previous call to compress_init.
!  (out)     alpha_star_bot - if present, this will be the compression adjusted
!                       in situ specific volume evaluated at the pressure at
!                       the bottom of each layer, in m s-2.

  real :: rho_in_situ(SZ1_(p)) !In-situ density at the top of a layer.
  real :: rho_is_bot(SZ1_(p))  !In-situ density at the layer's bottom.
  real :: alpha_in_situ(SZ1_(p)) !In-situ density at the top of a layer.
  real :: alpha_is_bot(SZ1_(p))  !In-situ density at the layer's bottom.
  real :: press(SZ1_(p))       !Pressure at the top of a layer in Pa.
  real :: r                    ! Input pressure, in Pa.
  real :: I_Rho0xG
  real :: B1, IB1, exp_gam_z
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, m, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq

  if (.not.associated(CS)) call GOLD_error(FATAL, "uncompress_p_alpha: "// &
       "Module GOLD_CompressComp must be initialized before it is used.")
  I_Rho0xG = 1.0 / (G%Rho0*G%g_Earth)
  do k=1,nz+1 ; do j=Jsq,Jeq+1
    if (k<=nz) then
      call calculate_density(tv%T(:,j,k),tv%S(:,j,k),p(:,j,k), &
                             rho_in_situ,Isq,Ieq-Isq+2,tv%eqn_of_state)
      do i=Isq,Ieq+1 ; alpha_in_situ(i) = 1.0 / rho_in_situ(i) ; enddo
    endif
    if (k>1 .and. (PRESENT(alpha_star_bot)) ) then
      call calculate_density(tv%T(:,j,k-1),tv%S(:,j,k-1),p(:,j,k), &
                             rho_is_bot,Isq,Ieq-Isq+2,tv%eqn_of_state)
      do i=Isq,Ieq+1 ; alpha_is_bot(i) = 1.0 / rho_is_bot(i) ; enddo
    endif
    
    do i=Isq,Ieq+1
      r = p(i,j,k)
      B1 = 2.0*CS%cp(i,j,2) + 3.0*CS%cp(i,j,3)*r
      p(i,j,k) = r + (CS%cp(i,j,2) + CS%cp(i,j,3)*r)*r*r
      do m=4,NTERM
        exp_gam_z = exp(CS%ep(m)*r)
        p(i,j,k) = p(i,j,k) + CS%cp(i,j,m)*(1.0 + exp_gam_z*(CS%ep(m)*r - 1.0))
        B1 = B1 + CS%ca_cp(m)*CS%cp(i,j,m)*exp_gam_z
      enddo
      IB1 = 1.0 / (1.0 + r*B1)  
      if (k<=nz) alpha_star(i,j,k) = IB1 * alpha_in_situ(i)
      if ((k>1) .and. PRESENT(alpha_star_bot)) &
        alpha_star_bot(i,j,k-1) = IB1 * alpha_is_bot(i)
    enddo
  enddo ; enddo

end subroutine uncompress_p_alpha

subroutine Grad_z_estar(e, gx_e, gy_e, nz, G, CS)
  real, dimension(NXMEM_,NYMEM_,NZp1_) :: e
  real, dimension(NXMEMQ_,NYMEM_,NZp1_) :: gx_e
  real, dimension(NXMEM_,NYMEMQ_,NZp1_) :: gy_e
  integer :: nz
  type(ocean_grid_type), intent(in) :: G
  type(Compress_CS),     pointer :: CS
!   This subroutine determines the slopes of compensated heights (z*)
! relative to geopotentials.
! Arguments: e - interface height in m, usually negative, unadjusted
!                for compressibility, intent in.
!  (out)     gx_e - Zonal slopes (nondimensional).
!  (out)     gy_e - Meridional slopes (nondimensional).
!  (in)      nz - the number of layers to work on.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure that has been set by a
!                 previous call to compress_init.

  real :: zx, zy
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq

  if (.not.associated(CS)) call GOLD_error(FATAL,"Grad_z_estar: "// &
       "Module GOLD_CompressComp must be initialized before it is used.")

  do k=1,nz+1
    do j=js,je ; do I=Isq,Ieq
      zx = 0.5*(e(i,j,k) + e(i+1,j,k))
      gx_e(I,j,k) = ((CS%cz(i+1,j,2)-CS%cz(i,j,2)) + &
          (CS%cz(i+1,j,3)-CS%cz(i,j,3))*zx)*zx*zx
      do m=4,NTERM
        gx_e(I,j,k) = gx_e(I,j,k) + (CS%cz(i+1,j,m)-CS%cz(i,j,m)) * &
            (1.0 + exp(CS%ez(m)*zx)*(CS%ez(m)*zx - 1.0))
      enddo
      gx_e(I,j,k) = gx_e(I,j,k) * G%IDXu(I,j)
    enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie
      zy = 0.5*(e(i,j,k) + e(i,j+1,k))
      gy_e(i,J,k) = ((CS%cz(i,j+1,2)-CS%cz(i,j,2)) + &
          (CS%cz(i,j+1,3)-CS%cz(i,j,3))*zy)*zy*zy
      do m=4,NTERM
        gy_e(i,J,k) = gy_e(i,J,k) + (CS%cz(i,j+1,m)-CS%cz(i,j,m)) * &
            (1.0 + exp(CS%ez(m)*zy)*(CS%ez(m)*zy - 1.0))
      enddo
      gy_e(i,J,k) = gy_e(i,J,k) * G%IDYv(i,J)
    enddo ; enddo
  enddo
end subroutine Grad_z_estar

subroutine Grad_p_pstar(p, gx_p, gy_p, nz, G, CS)
  real, dimension(NXMEM_,NYMEM_,NZp1_),  intent(in)  :: p
  real, dimension(NXMEMQ_,NYMEM_,NZp1_), intent(out) :: gx_p
  real, dimension(NXMEM_,NYMEMQ_,NZp1_), intent(out) :: gy_p
  integer,                               intent(in)  :: nz
  type(ocean_grid_type),                 intent(in)  :: G
  type(Compress_CS),                     pointer     :: CS
!   This subroutine determines the gradients of compensated pressures (p*)
! along pressure surfaces.
! Arguments: p - pressure in Pa, unadjusted for compressibility, intent in.
!  (out)     gx_p - Zonal p* gradients along p surfaces (Pa m-1).
!  (out)     gy_p - Meridional p* gradients along p surfaces (Pa m-1).
!  (in)      nz - the number of layers to work on.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure that has been set by a
!                 previous call to compress_init.

  real :: p_u, p_v
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq

  if (.not.associated(CS)) call GOLD_error(FATAL,"Grad_z_estar: "// &
       "Module GOLD_CompressComp must be initialized before it is used.")

  do k=1,nz+1
    do j=js,je ; do I=Isq,Ieq
      p_u = 0.5*(p(i,j,k) + p(i+1,j,k))
      gx_p(I,j,k) = ((CS%cp(i+1,j,2)-CS%cp(i,j,2)) + &
          (CS%cp(i+1,j,3)-CS%cp(i,j,3))*p_u)*p_u*p_u
      do m=4,NTERM
        gx_p(I,j,k) = gx_p(I,j,k) + (CS%cp(i+1,j,m)-CS%cp(i,j,m)) * &
            (1.0 + exp(CS%ep(m)*p_u)*(CS%ep(m)*p_u - 1.0))
      enddo
      gx_p(I,j,k) = gx_p(I,j,k) * G%IDXu(I,j)
    enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie
      p_v = 0.5*(p(i,j,k) + p(i,j+1,k))
      gy_p(i,J,k) = ((CS%cp(i,j+1,2)-CS%cp(i,j,2)) + &
          (CS%cp(i,j+1,3)-CS%cp(i,j,3))*p_v)*p_v*p_v
      do m=4,NTERM
        gy_p(i,J,k) = gy_p(i,J,k) + (CS%cp(i,j+1,m)-CS%cp(i,j,m)) * &
            (1.0 + exp(CS%ep(m)*p_v)*(CS%ep(m)*p_v - 1.0))
      enddo
      gy_p(i,J,k) = gy_p(i,J,k) * G%IDYv(i,J)
    enddo ; enddo
  enddo
end subroutine Grad_p_pstar


subroutine compress_init(G, param_file, CS)
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  type(Compress_CS), pointer   :: CS
! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module.
  character(len=128) :: version = '$Id: GOLD_CompressComp.F90,v 13.0.2.9.2.9 2010/08/26 18:19:04 rwh Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'

  if (associated(CS)) then
    call GOLD_error(WARNING, "compress_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  call select_eqn_of_state(param_file, CS%eqn_of_state)
  if (G%Boussinesq) then
    ALLOC(CS%cz(G%isd:G%ied,G%jsd:G%jed,NTERM)) ; CS%cz(:,:,:) = 0.0
  else
    ALLOC(CS%cp(G%isd:G%ied,G%jsd:G%jed,NTERM)) ; CS%cp(:,:,:) = 0.0
  endif

  call log_version(param_file, "GOLD_CompressComp", version, tagname)

end subroutine compress_init

end module GOLD_CompressComp
