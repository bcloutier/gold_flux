module  GOLD_OCEAN_TOPAZ_MOD  !{

! ----------------------------------------------------------------
!                   GNU General Public License                        
! This file is a part of GOLD.                                             
!                                                                      
! GOLD is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! GOLD is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
! 
!<CONTACT EMAIL="John.Dunne@noaa.gov"> John P. Dunne
!</CONTACT>
!
!<REVIEWER EMAIL="Robert.Hallberg@noaa.gov"> Robert Hallberg
!</REVIEWER>
!
!<REVIEWER EMAIL="William.Cooke@noaa.gov"> William Cooke
!</REVIEWER>
!
!<OVERVIEW>
! GFDL Ocean Biogeochemistry Version 1 (TOPAZ1) Model
!</OVERVIEW>
!
!<DESCRIPTION>
!       Phytoplankton Biogeochemistry: Includes an explicit ecological model
!       including three phytoplankton groups (small, large/diatoms and
!       diazotrophs), growth limitation by light, temperature and a suite of
!       nutrients including nitrate, ammonia, phosphate, iron and silicate,
!       dissolved inorganic carbon, alkalinity, two kinds of dissolved organic
!       material, O2, nitrogen fixation and denitrification. CO2 gas exchange
!       is function of the biologically and physically forced solubility.
!       Additionally, changes in the vertical distribution of phytoplankton
!       affect heat absorption with climate feedbacks. Food web processing in
!       the euphotic zone and remineralization/dissolution through the ocean
!       interior are handled as in Dunne et al. (in prep).  CO2 and O2
!       equilibria and gas exchange follow OCMIP2 protocols.
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! http://www.ipsl.jussieu.fr/OCMIP/phase2/simulations/Biotic/HOWTO-Biotic.html
! </REFERENCE>
!
! <REFERENCE>
! Press, W. H., S. A. Teukosky, W. T. Vetterling, B. P. Flannery, 1992. 
! Numerical Recipes in FORTRAN, Second Edition, Cambridge University Press. 
! </REFERENCE>
!
! <REFERENCE>
! Enting, I.G., T. M. L. Wigley, M. Heimann, 1994. Future Emissions 
! and concentrations of carbon dioxide: key ocean / atmosphere / 
! land analyses, CSIRO Aust. Div. Atmos. Res. Tech. Pap. No. 31, 
! 118 pp.
! </REFERENCE>
! </INFO>
!
!------------------------------------------------------------------
!
!       Module ocean_topaz_mod
!
!       Implementation of routines to solve the GFDL Ocean Biogeochemistry
!       simulations.
!
!------------------------------------------------------------------
!

!
!------------------------------------------------------------------
!
!       Global definitions
!
!------------------------------------------------------------------
!

!
!----------------------------------------------------------------------
!
!       Modules
!
!----------------------------------------------------------------------
!

use field_manager_mod,  only: fm_field_name_len, fm_string_len
use fms_mod,            only: field_exist
use mpp_mod,            only: stdout, stdlog, mpp_error, mpp_sum, FATAL, NOTE
use mpp_domains_mod,    only: domain2d
use constants_mod,      only: WTMCO2, WTMO2

use coupler_util, only : extract_coupler_values, set_coupler_values
use coupler_types_mod,  only: ind_alpha, ind_csurf, ind_flux, coupler_2d_bc_type
use atmos_ocean_fluxes_mod, only: aof_set_coupler_flux
use GOLD_ocmip2_co2calc_mod, only : GOLD_ocmip2_co2calc, CO2_dope_vector

use GOLD_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use GOLD_diag_mediator, only : diag_ptrs, ocean_register_diag
use GOLD_diag_to_Z, only : diag_to_Z_CS, ocean_register_diag_with_z
use GOLD_checksums, only : hchksum
use GOLD_error_handler, only : GOLD_error, FATAL, WARNING, is_root_pe
use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
use GOLD_grid, only : ocean_grid_type
use GOLD_io, only : file_exists, read_data, slasher, vardesc
use GOLD_restart, only : register_restart_field, GOLD_restart_CS
use GOLD_sponge, only : set_up_sponge_field, sponge_CS
use GOLD_time_manager, only : time_type, get_time
use GOLD_tracer, only : register_tracer, advect_tracer_CS
use GOLD_tracer, only : add_tracer_diagnostics, add_tracer_OBC_values
use GOLD_tracer, only : tracer_vertdiff
use GOLD_variables, only : forcing, surface, ocean_OBC_type, thermo_var_ptrs
use GOLD_variables, only : optics_type

!
!----------------------------------------------------------------------
!
!       force all variables to be "typed"
!
!----------------------------------------------------------------------
!
implicit none
!
!----------------------------------------------------------------------
!
!       Make all routines and variables private by default
!
!----------------------------------------------------------------------
!
private

#include <GOLD_memory.h>
#ifdef _USE_TOPAZ
#include <fms_platform.h>

!
!----------------------------------------------------------------------
!
!       Public routines
!
!----------------------------------------------------------------------
!

public  :: register_TOPAZ
public  :: initialize_TOPAZ
public  :: TOPAZ_column_physics
public  :: TOPAZ_surface_state
public  :: TOPAZ_end
public  :: TOPAZ_CS
public  :: get_chl_from_TOPAZ
public  :: TOPAZ_coupler_flux_init
!
!----------------------------------------------------------------------
!
!       Private routines
!
!----------------------------------------------------------------------
!
private :: allocate_arrays
private :: register_diagnostics
!
!----------------------------------------------------------------------
!
!       Private parameters
!
!----------------------------------------------------------------------
!
character(len=32), parameter            :: package_name = 'ocean_topaz'
character(len=48), parameter            :: mod_name = 'ocean_topaz_mod'
character(len=fm_string_len), parameter :: &
                  default_restart_file       =  'ocean_topaz.res.nc',   &
                  default_ice_restart_file   =  'ice_topaz.res.nc',     &
                  default_ocean_restart_file =  'ocean_topaz.res.nc'
!
!       coefficients for O2 saturation
!
real, parameter :: a_0 = 2.00907
real, parameter :: a_1 = 3.22014
real, parameter :: a_2 = 4.05010
real, parameter :: a_3 = 4.94457
real, parameter :: a_4 = -2.56847e-01
real, parameter :: a_5 = 3.88767
real, parameter :: b_0 = -6.24523e-03
real, parameter :: b_1 = -7.37614e-03
real, parameter :: b_2 = -1.03410e-02
real, parameter :: b_3 = -8.17083e-03
real, parameter :: c_0 = -4.88682e-07
!
!       time conversions
!
real, parameter :: sperd = 24.0 * 3600.0
real, parameter :: spery = 365.25 * 24.0 * 3600.0
!
!----------------------------------------------------------------------
!
!       Private types
!
!----------------------------------------------------------------------
!
type phytoplankton
  real :: alpha,   &
    fdet0,         &
    fe_2_n_max,    &
    fe_2_n_static, &
    k_fe_2_n,      &
    k_fed,         &
    k_nh4,         &
    k_no3,         &
    k_p_2_n,       &
    k_po4,         &
    k_sio4,        &
    p_2_n_assem,   &
    p_2_n_static,  &
    P_C_max,       &
    plast_2_chl,   &
    r_nfix,        &
    r_other_min,   &
    r_uptake_max,  &
    si_2_n_max,    &
    si_2_n_static, &
    thetamax     
  real, _ALLOCATABLE, dimension(:,:,:)  :: &
    def_fe      _NULL , & ! Fe Deficiency
    def_p       _NULL , & ! P Deficiency
    felim       _NULL , & ! Fed Limitation
    irrlim      _NULL , & ! Light Limitation
    jgraz_fe    _NULL , & ! Fe grazing layer integral
    jgraz_n     _NULL , & ! Nitrogen grazing layer integral
    jgraz_sio2  _NULL , & ! Silicon grazing layer integral
    jprod_n2    _NULL , & ! Nitrogen fixation layer integral
    jprod_fe    _NULL , & ! Fe production layer integral
    jprod_nh4   _NULL , & ! NH4 production layer integral
    jprod_no3   _NULL , & ! NO3 production layer integral
    jprod_po4   _NULL , & ! PO4 production layer integral
    jprod_sio4  _NULL , & ! Silicon production layer integral
    liebig_lim  _NULL , & ! Overall nutrient limitation
    mu          _NULL , & ! Overall growth rate
    nh4lim      _NULL , & ! Ammonia limitation
    no3lim      _NULL , & ! Nitrate limitation
    po4lim      _NULL , & ! Phosphate limitation
    q_fe_2_n    _NULL , & ! Fe:N ratio
    q_p_2_n     _NULL , & ! P:N ratio
    q_p_2_n_opt _NULL , & ! Optimal P:N ratio
    silim       _NULL , & ! SiO4 limitation
    theta       _NULL     ! Chl:C ratio
  integer ::            &
    id_def_fe     = -1, & ! Diag id for Phyto. Fe Deficiency
    id_def_p      = -1, & ! Diag id for Phyto. P Deficiency
    id_felim      = -1, & ! Diag id for Phyto. Fed Limitation
    id_irrlim     = -1, & ! Diag id for Phyto. Light Limitation
    id_jgraz_fe   = -1, & ! Diag id for iron grazing layer integral
    id_jgraz_n    = -1, & ! Diag id for nitrogen grazing layer integral
    id_jgraz_sio2 = -1, & ! Diag id for silicon grazing layer integral
    id_jprod_n2   = -1, & ! Diag id for Nitrogen fixation layer integral
    id_jprod_fe   = -1, & ! Diag id for phyto. Fed production layer integral
    id_jprod_nh4  = -1, & ! Diag id for phyto. NH4 production layer integral
    id_jprod_no3  = -1, & ! Diag id for phyto. NO3 production layer integral
    id_jprod_po4  = -1, & ! Diag id for phyto. PO4 production layer integral
    id_jprod_sio4 = -1, & ! Diag id for phyto. SiO4 production layer integral
    id_liebig_lim = -1, & ! Diag id for Overall nutrient limitation
    id_mu         = -1, & ! Diag id for Overall growth rate
    id_nh4lim     = -1, & ! Diag id for Ammonia Limitation of Phyto
    id_no3lim     = -1, & ! Diag id for Nitrate Limitation of Phyto
    id_po4lim     = -1, & ! Diag id for Phosphate Limitation of Phyto
    id_q_fe_2_n   = -1, & ! Diag id for Fe:N ratio
    id_q_p_2_n    = -1, & ! Diag id for P:N ratio
    id_q_p_2_n_opt= -1, & ! Diag id for Optimal P:N ratio
    id_silim      = -1, & ! Diag id for SiO4 Limitation of Phyto
    id_theta      = -1
end type phytoplankton

integer, parameter :: NUM_PHYTO  = 3
!
! Array allocations and flux calculations assume that phyto(1) is the
! only phytoplankton group cabable of nitrogen uptake by N2 fixation while phyto(2:NUM_PHYTO) 
! are only cabable of nitrgen uptake by NH4 and NO3 uptake
!
integer, parameter :: DIAZO      = 1
integer, parameter :: LARGE      = 2
integer, parameter :: SMALL      = 3
type(phytoplankton), dimension(NUM_PHYTO) :: phyto

type tracer_type
! Tracer name
  character(len=64)                     :: name
! Tracer field
  real, _ALLOCATABLE, dimension(:,:,:)  :: field   _NULL
! Surface tracer flux
  real, _ALLOCATABLE, dimension(:,:)    :: stf     _NULL 
! Bottom tracer flux
  real, _ALLOCATABLE, dimension(:,:)    :: btf     _NULL 
! Tracer concentration in river runoff
  real, _ALLOCATABLE, dimension(:,:)    :: btm_reservoir     _NULL 
! Tracer concentration in river runoff
  real, _ALLOCATABLE, dimension(:,:)    :: triver  _NULL 
! Sinking rate
  real                                  :: sink_dist  
end type tracer_type

integer, parameter :: NUM_TRACERS = 38
integer, parameter :: NUM_PROG_TRACERS = 29
type(tracer_type), dimension(:), pointer :: Tr

type TOPAZ_CS  !{
  character(len=fm_string_len)          :: restart_file
  character(len=fm_field_name_len)      :: name

  logical  ::       &
  ca_2_n_static,    &                  ! If Ca:N is fixed in heterotrophs
  fe_2_n_static,    &                  ! If Fe:N is fixed in phytoplankton
  fe_ballast_assoc, &                  ! If iron scavenging is associated with ballast
  init,             &                  ! If tracers should be initializated
  p_2_n_static,     &                  ! If P:N is fixed in phytoplankton
  si_2_n_static,    &                  ! If Si:N is fixed in phytoplankton
  tracer_debug,     &
  tracer_debug_verbose

  real  ::          &
  atm_co2_flux,     &
  bio_tau,          &                  ! Miscellaneous
  c_2_n,            &                  ! Stoichiometry
  ca_2_n_het,       &                  ! CaCO3 formation in heterotrophis
  ca_2_n_het_static,&                  ! CaCO3 formation if fixed
  caco3_sat_max,    &                  ! Calcite maximum saturation state
  fe_2_n_sed,       &                  ! Iron
  fe_coast,         &                  ! Iron
  felig_2_don,      &                  ! Iron
  felig_bkg ,       &                  ! Iron
  gamma_cadet,      &                  ! Grazing
  gamma_cased_dis,  &                  ! Grazing
  gamma_irr_mem,    &                  ! Photosynthesis
  gamma_ldon,       &                  ! Dissolved Organic Material
  gamma_ndet,       &                  ! Grazing
  gamma_nhet,       &                  ! Grazing
  gamma_nitrif,     &                  ! Miscellaneous
  gamma_sidet,      &                  ! Grazing
  gamma_sdon,       &                  ! Dissolved Organic Material
  gamma_sdop,       &                  ! Dissolved Organic Material
  irr_inhibit,      &                  ! Miscellaneous  ! Photosynthesis?
  k_caco3_pres,     &                  ! Grazing
  k_diss_sio2,      &                  ! Grazing
  k_n_inhib_di,     &                  ! Photosynthesis
  k_o2,             &                  ! Stoichiometry
  kappa_eppley,     &                  ! Photosynthesis
  kappa_remin,      &                  ! Grazing
  kfe_2nd_order,    &                  ! Iron
  kfe_bal,          &
  kfe_des,          &                  ! Iron
  kfe_eq_lig,       &                  ! Iron
  kfe_org,          &
  k_lith,           &                  ! Miscellaneous
  ksp_caco3,        &                  ! Grazing
  lambda0,          &                  ! Grazing
  mass_2_n,         &                  ! Stoichiometry
  n_2_n_denit,      &                  ! Stoichiometry
  o2_min,           &                  ! Grazing
  o2_2_c,           &                  ! Stoichiometry
  o2_2_nfix,        &
  o2_2_nh4,         &                  ! Stoichiometry
  o2_2_no3,         &                  ! Stoichiometry
  o2_2_nitrif,      &                  ! Stoichiometry
  o2_inhib_di_pow,  &                  ! Photosynthesis
  o2_inhib_di_sat,  &                  ! Photosynthesis
  p_2_n_photo,      &                  ! P:N of photosynthesis machinery (Chloroplasts)
  p_2_n_RKR,        &                  ! P:N of standard Stoichiometry
  p_2_n_uptake,     &                  ! P:N of uptake machinery
  P_C_max_assem,    &                  ! Maximum assembly rate
  P_star,           &                  ! Grazing
  phi_nhet,         &                  ! Grazing
  phi_sdon,         &                  ! Dissolved Organic Material
  phi_sdop,         &                  ! Dissolved Organic Material
  phi_ldon,         &                  ! Dissolved Organic Material
  phi_lith,         &                  ! Miscellaneous
  phyto_min,        &                  ! Grazing
  plast_2_chl,      &                  ! Chloroplast:Chlorophyll conversion
  q_si_2_n_diss,    &                  ! Stoichiometric role in SiO2 dissolution
  r_bio_tau,        &                  ! Miscellaneous
  rpcaco3,          &                  ! Grazing
  rplith,           &                  ! Grazing
  rpsio2,           &                  ! Grazing
  thetamin,         &                  ! Photosynthesis
  wsink,            &                  ! Sinking
  zeta                                 ! Photosynthesis

  real, dimension(3)                    :: total_atm_co2
  real, _ALLOCATABLE, dimension(:,:)    :: alpha  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: csurf  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: co3_solubility  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: fcaco3_burial  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: fcaco3_redis  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: fcaco3_sed  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: fdenit_sed  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: ffe_sed  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: ffedet_btm  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: flithdet_btm  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jcadet  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jdenit_wc  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jdiss_sio2  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jfe_ads  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jfe_des  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jfe_graz  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jfe_coast  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jfedet  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jldon  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jndet  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jnh4  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jnh4_graz  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jnhet  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jnitrif  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jno3  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jo2  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jpdet  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jpo4  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jpo4_graz  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_cadet  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_lithdet  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_nhet  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_fedet  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_ndet  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jprod_pdet  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jsdon  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jsdop  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: jsidet  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: k1_co2  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: k2_co2  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: mask_coast  _NULL
  real, _ALLOCATABLE, dimension(:,:,:)  :: nLg_diatoms  _NULL
  real, _ALLOCATABLE, dimension(:,:)    :: pco2surf  _NULL
  integer                   :: &
  id_alpha         = -1,       & ! 
  id_btm_flux_alk  = -1,       & ! Bottom Flux - Alk
  id_btm_flux_dic  = -1,       & ! Bottom Flux - DIC
  id_btm_flux_ndet = -1,       & ! Bottom Flux - ndet
  id_btm_flux_no3  = -1,       & ! Bottom Flux - no3
  id_btm_flux_pdet = -1,       & ! Bottom Flux - pdet
  id_btm_flux_po4  = -1,       & ! Bottom Flux - po4
  id_co3_ion       = -1,       & ! Carbonate Ion
  id_co3_solubility = -1,      & ! Carbonate Ion Solubility
  id_csurf         = -1,       & ! CO2* water
  id_fcaco3        = -1,       & ! CaCO3 sinking flux
  id_fcaco3_burial = -1,       & ! CaCO3 sinking flux permanent burial
  id_fcaco3_redis  = -1,       & ! CaCO3 redissolution after sinking flux burial
  id_fcaco3_sed    = -1,       & ! CaCO3 sinking flux to sediment layer
  id_fdenit_sed    = -1,       & ! Sediment Denitrification flux
  id_ffe_sed       = -1,       & ! Sediment iron efflux
  id_ffedet_btm    = -1,       & ! Fe sinking flux burial
  id_flith         = -1,       & ! Lith sinking flux
  id_flithdet_btm  = -1,       & ! Lithogenic sinking flux burial
  id_fpofe         = -1,       & ! POFe sinking flux
  id_fpon          = -1,       & ! PON sinking flux
  id_fpop          = -1,       & ! POP sinking flux
  id_fsio2         = -1,       & ! Si sinking flux
  id_htotal        = -1,       & ! H+ ion concentration
  id_jcadet        = -1,       & ! CaCO3 change layer integral
  id_jdenit_wc     = -1,       & ! Water column Denitrification layer integral
  id_jdiss_sio2    = -1,       & ! SiO2 Dissolution during grazing layer integral
  id_jfe_ads       = -1,       & ! Iron adsorption layer integral
  id_jfe_des       = -1,       & ! Iron desorption layer integral
  id_jfe_graz      = -1,       & ! Dissolved iron grazing source layer integral
  id_jfe_coast     = -1,       & ! Coastal iron efflux layer integral
  id_jfedet        = -1,       & ! Loss of sinking iron layer integral
  id_jldon         = -1,       & ! Labile DON source layer integral
  id_jndet         = -1,       & ! Loss of sinking nitrogen layer integral
  id_jnh4          = -1,       & ! NH4 source layer integral
  id_jnh4_graz     = -1,       & ! NH4 grazing source layer integral
  id_jnhet         = -1,       & ! Heterotrophic N remineralization layer integral
  id_jnitrif       = -1,       & ! Nitrification layer integral
  id_jno3          = -1,       & ! NO3 source layer integral
  id_jo2           = -1,       & ! O2 source layer integral
  id_jpdet         = -1,       & ! Loss of sinking phosphorus layer integral
  id_jpo4          = -1,       & ! PO4 source layer integral
  id_jpo4_graz     = -1,       & ! PO4 source from grazing layer integral
  id_jprod_cadet   = -1,       & ! CaCO3 production layer integral  
  id_jprod_lithdet = -1,       & ! Lithogenic removal to sinking layer integral
  id_jprod_nhet    = -1,       & ! Heterotrophic N Production layer integral
  id_jprod_fedet   = -1,       & ! Detrital Fe production layer integral
  id_jprod_ndet    = -1,       & ! Detrital N production layer integral
  id_jprod_pdet    = -1,       & ! Detrital P production layer integral
  id_jsdon         = -1,       & ! Semilabile DON source layer integral
  id_jsdop         = -1,       & ! Semilabile DOP source layer integral
  id_jsidet         = -1,       & ! SiO4 source layer integral
  id_nLg_diatoms   = -1,       & ! Large diatom nitrogen
  id_pco2surf      = -1,       & ! Oceanic pCO2
  id_runoff_alk    = -1,       & ! Runoff - Alkalinity
  id_runoff_dic    = -1,       & ! Runoff - Dissolved Inorganic Carbon
  id_runoff_fed    = -1,       & ! Runoff - Dissolved Fe
  id_runoff_flux_alk = -1,     & ! Runoff flux - Alkalinity
  id_runoff_flux_dic = -1,     & ! Runoff flux - Dissolved Inorganic Carbon
  id_runoff_flux_fed = -1,     & ! Runoff flux - Dissolved Fe
  id_runoff_flux_ldon = -1,    & ! Runoff flux - LDON
  id_runoff_flux_lith = -1,    & ! Runoff flux - Lith
  id_runoff_flux_nh4 = -1,     & ! Runoff flux - NH4
  id_runoff_flux_no3 = -1,     & ! Runoff flux - NO3
  id_runoff_ldon   = -1,       & ! Runoff - LDON
  id_runoff_lith   = -1,       & ! Runoff - Lith
  id_runoff_nh4    = -1,       & ! Runoff - NH4
  id_runoff_no3    = -1,       & ! Runoff - NO3
  id_sfc_chl       = -1,       & ! Surface Chl
  id_sfc_no3       = -1,       & ! Surface NO3
  id_sfc_flux_co2  = -1,       & ! Surface Flux - CO2
  id_sfc_flux_fed  = -1,       & ! Surface Flux - Fed
  id_sfc_flux_lith = -1,       & ! Surface Flux - Lith
  id_sfc_flux_no3  = -1,       & ! Surface Flux - NO3
  id_sfc_flux_nh4  = -1,       & ! Surface Flux - NH4
  id_sfc_flux_o2   = -1,       & ! Surface Flux - O2
  id_tot_layer_int_c = -1,     & ! Total Carbon (DIC+OC+IC) boxwise layer integral
  id_tot_layer_int_fe = -1,    & ! Total Phosphorus (Fed+OFe) boxwise layer integral
  id_tot_layer_int_n = -1,     & ! Total Nitrogen (NO3+NH4+ON) boxwise layer integral
  id_tot_layer_int_p = -1,     & ! Total Phosphorus (PO4+OP) boxwise layer integral
  id_tot_layer_int_si = -1,    & ! Total Silicon (SiO4+SiO2) boxwise layer integral
  id_alk           = -1,       & ! Alkalinity Prognostic tracer
  id_cadet         = -1,       & ! Particulate Detrital CaCO3 Prognostic tracer
  id_dic           = -1,       & ! DIC Prognostic tracer
  id_fed           = -1,       & ! Dissolved Iron Prognostic tracer
  id_fedi          = -1,       & ! Diaz Iron Prognostic tracer
  id_felg          = -1,       & ! Large Iron Prognostic tracer
  id_fedet         = -1,       & ! Particulate Detrital Iron Prognostic tracer
  id_fesm          = -1,       & ! Small Iron Prognostic tracer
  id_ldon          = -1,       & ! Labile DON Prognostic tracer
  id_lith          = -1,       & ! Lithogenic Mineral Prognostic tracer
  id_lithdet       = -1,       & ! Lithogenic Mineral Prognostic tracer
  id_ndet          = -1,       & ! Particulate Detrital Nitrogen Prognostic tracer
  id_ndi           = -1,       & ! Diaz. Nitrogen Prognostic tracer
  id_nh4           = -1,       & ! Ammonium Prognostic tracer
  id_nhet          = -1,       & ! Heterotrophic N Prognostic tracer
  id_nlg           = -1,       & ! Large Nitrogen Prognostic tracer
  id_no3           = -1,       & ! Nitrate Prognostic tracer
  id_nsm           = -1,       & ! Small Nitrogen Prognostic tracer
  id_o2            = -1,       & ! Oxygen Prognostic tracer
  id_pdet          = -1,       & ! Particulate Detrital Phosphorus Prognostic tracer
  id_pdi           = -1,       & ! Diaz. Phosphorus Prognostic tracer
  id_plg           = -1,       & ! Large Phosphorus Prognostic tracer
  id_po4           = -1,       & ! Phosphate Prognostic tracer
  id_psm           = -1,       & ! Small Phosphorus Prognostic tracer
  id_sdon          = -1,       & ! Semilabile DON Prognostic tracer
  id_sdop          = -1,       & ! Semilabile DOP Prognostic tracer
  id_sidet         = -1,       & ! Particulate Detrital Silicon Prognostic tracer
  id_silg          = -1,       & ! Large Silicon Prognostic tracer
  id_sio4          = -1,       & ! Silicic Acid Prognostic tracer
  id_cased         = -1,       & ! Sediment CaCO3 Diagnostic tracer
  id_chl           = -1,       & ! Chlorophyll Diagnostic tracer
  id_fcadet_btm    = -1,       & ! CaCO3 sinking flux at bottom
  id_irr           = -1,       & ! Irradiance Diagnostic tracer
  id_irr_mem       = -1          ! Irradiance Memory Diagnostic tracer

  integer                   :: &
  ! Indices for various tracers follow.
  ind_alk          = 1,        & ! Alkalinity Prognostic tracer
  ind_cadet        = 2,        & ! DIC Prognostic tracer
  ind_dic          = 3,        & ! Particulate Detrital CaCO3 Prognostic tracer
  ind_fed          = 4,        & ! Dissolved Iron Prognostic tracer 
  ind_fedi         = 5,        & ! Diaz Iron Prognostic tracer   
  ind_felg         = 6,        & ! Large Iron Prognostic tracer   
  ind_fedet        = 7,        & ! Particulate Detrital Iron Prognostic tracer 
  ind_fesm         = 8,        & ! Small Iron Prognostic tracer 
  ind_ldon         = 9,        & ! Labile DON Prognostic tracer
  ind_lith         = 10,       & ! Lithogenic Mineral Prognostic tracer
  ind_lithdet      = 11,       & ! Particulate Detrital Lithogenic Prognostic tracer
  ind_ndet         = 12,       & ! Particulate Detrital Nitrogen Prognostic tracer
  ind_ndi          = 13,       & ! Diaz. Nitrogen Prognostic tracer
  ind_nh4          = 14,       & ! Ammonium Prognostic tracer 
  ind_nhet         = 15,       & ! Heterotrophic N Prognostic tracer
  ind_nlg          = 16,       & ! Large Nitrogen Prognostic tracer
  ind_no3          = 17,       & ! Nitrate Prognostic tracer
  ind_nsm          = 18,       & ! Small Nitrogen Prognostic tracer
  ind_o2           = 19,       & ! Oxygen Prognostic tracer
  ind_pdet         = 20,       & ! Particulate Detrital Phosphorus Prognostic tracer
  ind_pdi          = 21,       & ! Diaz. Phosphorus Prognostic tracer
  ind_plg          = 22,       & ! Large Phosphorus Prognostic tracer
  ind_po4          = 23,       & ! Phosphate Prognostic tracer
  ind_psm          = 24,       & ! Small Phosphorus Prognostic tracer
  ind_sdon         = 25,       & ! Semilabile DON Prognostic tracer
  ind_sdop         = 26,       & ! Semilabile DOP Prognostic tracer
  ind_sidet        = 27,       & ! Particulate Detrital Silicon Prognostic tracer
  ind_silg         = 28,       & ! Large Silicon Prognostic tracer
  ind_sio4         = 29,       & ! Silicic Acid Prognostic tracer
  ind_cased        = 30,       & ! Sediment CaCO3 Diagnostic tracer
  ind_chl          = 31,       & ! Chlorophyll Diagnostic tracer
  ind_fcadet_btm   = 32,       & ! CaCO3 sinking flux at bottom
  ind_irr          = 33,       & ! Irradiance Diagnostic tracer
  ind_irr_mem      = 34,       & ! Irradiance Memory Diagnostic tracer
  ind_htotal       = 35,       & ! H+ ion concentration Diagnostic tracer
  ind_co3_ion      = 36,       & ! Carbonate ion Diagnostic tracer
  ind_alpha        = 37,       & ! alpha Diagnostic tracer
  ind_csurf        = 38,       & ! CO2_csurf Diagnostic tracer
  ind_co2_flux     = -1,       & ! air_sea_gas_flux ocmip2
  ind_dry_dep_fed  = -1,       & ! air_sea_deposition dry
  ind_dry_dep_lith = -1,       & ! air_sea_deposition dry
  ind_dry_dep_nh4  = -1,       & ! air_sea_deposition dry
  ind_dry_dep_no3  = -1,       & ! air_sea_deposition dry
  ind_o2_flux      = -1,       & ! air_sea_gas_flux ocmip2
  ind_runoff_alk   = -1,       & ! land_sea_runoff river
  ind_runoff_dic   = -1,       & ! land_sea_runoff river
  ind_runoff_fed   = -1,       & ! land_sea_runoff river
  ind_runoff_ldon  = -1,       & ! land_sea_runoff river
  ind_runoff_lith  = -1,       & ! land_sea_runoff river
  ind_runoff_nh4   = -1,       & ! land_sea_runoff river
  ind_runoff_no3   = -1,       & ! land_sea_runoff river
  ind_wet_dep_fed  = -1,       & ! air_sea_deposition wet
  ind_wet_dep_lith = -1,       & ! air_sea_deposition wet
  ind_wet_dep_nh4  = -1,       & ! air_sea_deposition wet
  ind_wet_dep_no3  = -1          ! air_sea_deposition wet
  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
                             ! ocean diagnostic fields and control variables.

  integer :: nkml            ! The number of layers in the mixed layer.
  real  :: Rho_0
end type TOPAZ_CS  !}

type(CO2_dope_vector) :: CO2_dope_vec

integer :: &
  ind_co2_flux     = -1,       & ! air_sea_gas_flux ocmip2
  ind_dry_dep_fed  = -1,       & ! air_sea_deposition dry
  ind_dry_dep_lith = -1,       & ! air_sea_deposition dry
  ind_dry_dep_nh4  = -1,       & ! air_sea_deposition dry
  ind_dry_dep_no3  = -1,       & ! air_sea_deposition dry
  ind_o2_flux      = -1,       & ! air_sea_gas_flux ocmip2
  ind_runoff_alk   = -1,       & ! land_sea_runoff river
  ind_runoff_dic   = -1,       & ! land_sea_runoff river
  ind_runoff_fed   = -1,       & ! land_sea_runoff river
  ind_runoff_ldon  = -1,       & ! land_sea_runoff river
  ind_runoff_lith  = -1,       & ! land_sea_runoff river
  ind_runoff_nh4   = -1,       & ! land_sea_runoff river
  ind_runoff_no3   = -1,       & ! land_sea_runoff river
  ind_wet_dep_fed  = -1,       & ! air_sea_deposition wet
  ind_wet_dep_lith = -1,       & ! air_sea_deposition wet
  ind_wet_dep_nh4  = -1,       & ! air_sea_deposition wet
  ind_wet_dep_no3  = -1          ! air_sea_deposition wet
!
!----------------------------------------------------------------------
!
!       Public variables
!
!----------------------------------------------------------------------
!

!
!----------------------------------------------------------------------
!
!       Private variables
!
!----------------------------------------------------------------------
!

logical                                 :: module_initialized = .false.
logical                                 :: topaz_fluxes_initialized = .false.

character(len=128) :: version = '$Id: GOLD_OCEAN_TOPAZ.F90,v 1.1.4.17 2010/07/28 14:22:43 aja Exp $'
character(len=128) :: tagname = '$Name: GOLD_ogrp $'

!
!----------------------------------------------------------------------
!
!       Input parameters:
!
!  htotal_in            = default value for htotal for an initial run
!  htotal_scale_lo      = scaling parameter to chose htotallo
!  htotal_scale_hi      = scaling parameter to chose htotalhi
!
!----------------------------------------------------------------------
!

real                                    :: htotal_in
real                                    :: htotal_scale_hi
real                                    :: htotal_scale_lo

!
!----------------------------------------------------------------------
!
!       Calculated parameters (with possible initial input values):
!
!  global_wrk_duration  = total time during calculation of global
!                         variables
!
!----------------------------------------------------------------------
!

integer                                         :: id_o2_sat = -1
integer                                         :: id_sc_co2 = -1
integer                                         :: id_sc_o2 = -1
real, allocatable, dimension(:,:)               :: sc_no_term
real, allocatable, dimension(:,:)               :: sc_co2
real, allocatable, dimension(:,:)               :: sc_o2
real, allocatable, dimension(:,:,:)             :: htotalhi
real, allocatable, dimension(:,:,:)             :: htotallo
integer                                         :: instances
real, allocatable, dimension(:,:)               :: o2_saturation
real, allocatable, dimension(:)                 :: tk
real, allocatable, dimension(:)                 :: ts
real, allocatable, dimension(:)                 :: ts2
real, allocatable, dimension(:)                 :: ts3
real, allocatable, dimension(:)                 :: ts4
real, allocatable, dimension(:)                 :: ts5
real, allocatable, dimension(:)                 :: tt

!
!-----------------------------------------------------------------------
!
!       Subroutine and function definitions
!
!-----------------------------------------------------------------------
!

contains

! <DESCRIPTION>
!   Routine to simply pass back the values of the Chlorophyll tracer field.
! </DESCRIPTION>
subroutine get_chl_from_TOPAZ(Chl_array, topaz)
real, dimension(:,:,:), intent(out) :: Chl_array
type(TOPAZ_CS),         intent(in)  :: topaz

if ( _ALLOCATED(Tr(topaz%ind_chl)%field)) then 
  Chl_array = Tr(topaz%ind_chl)%field
else
  Chl_array = 0.0
endif

end subroutine get_chl_from_TOPAZ

! <DESCRIPTION>
!       Set up any extra fields needed by the ocean-atmosphere gas fluxes
! Initialize variables, read in namelists, calculate constants for a given run
! and allocate diagnostic arrays
! </DESCRIPTION>

function register_TOPAZ(G, param_file, topaz, diag, &
                                     tr_adv_CSp, restart_CS)
  type(ocean_grid_type), intent(in)   :: G
  type(param_file_type), intent(in)   :: param_file
  type(TOPAZ_CS),   pointer      :: topaz
  type(diag_ptrs), target, intent(in) :: diag
  type(advect_tracer_CS), pointer     :: tr_adv_CSp
  type(GOLD_restart_CS),   pointer     :: restart_CS
! This subroutine is used to register tracer fields and subroutines
! to be used with GOLD.
! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  tr_adv_CSp - A pointer that is set to point to the control structure
!                  for the tracer advection and diffusion module.
!  (in)      restart_CS - A pointer to the restart control structure.

logical :: register_TOPAZ
integer :: i, j, k
integer :: isc, iec, jsc, jec, nz
type(vardesc)   :: temp_desc
character(len=40)  :: mod = "GOLD_OCEAN_TOPAZ" ! This module's name.
character(len=48) :: name

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nz = G%ke
  CO2_dope_vec%isc = G%isc ; CO2_dope_vec%iec = G%iec 
  CO2_dope_vec%jsc = G%jsc ; CO2_dope_vec%jec = G%jec
  CO2_dope_vec%isd = G%isd ; CO2_dope_vec%ied = G%ied
  CO2_dope_vec%jsd = G%jsd ; CO2_dope_vec%jed = G%jed

!
!       allocate storage for topaz array
!

allocate ( topaz )
topaz%diag => diag
topaz%name = " "

!
!-----------------------------------------------------------------------
!     dynamically allocate the global topaz arrays
!-----------------------------------------------------------------------
!

call allocate_arrays(G, topaz, phyto)

!
!       Set up the ocean-atmosphere gas flux fields
!

!
!       Coupler fluxes
!

  call TOPAZ_coupler_flux_init
  topaz%ind_co2_flux     = ind_co2_flux    
  topaz%ind_o2_flux      = ind_o2_flux     
  topaz%ind_runoff_alk   = ind_runoff_alk  
  topaz%ind_runoff_dic   = ind_runoff_dic  
  topaz%ind_runoff_fed   = ind_runoff_fed  
  topaz%ind_runoff_ldon  = ind_runoff_ldon 
  topaz%ind_runoff_lith  = ind_runoff_lith 
  topaz%ind_runoff_nh4   = ind_runoff_nh4  
  topaz%ind_runoff_no3   = ind_runoff_no3  
  topaz%ind_dry_dep_fed  = ind_dry_dep_fed 
  topaz%ind_wet_dep_fed  = ind_wet_dep_fed 
  topaz%ind_dry_dep_lith = ind_dry_dep_lith
  topaz%ind_wet_dep_lith = ind_wet_dep_lith
  topaz%ind_dry_dep_nh4  = ind_dry_dep_nh4 
  topaz%ind_wet_dep_nh4  = ind_wet_dep_nh4 
  topaz%ind_dry_dep_no3  = ind_dry_dep_no3 
  topaz%ind_wet_dep_no3  = ind_wet_dep_no3 

!
!       Set up the field input
!


  allocate(Tr(NUM_TRACERS))
  name = ""
  if (topaz%name .ne. " " ) name = '('//trim(topaz%name)//')'

!   The following vardesc types contain a package of metadata about each tracer,
! including, in order, the following elements: name; longname; horizontal
! staggering ('h') for collocation with thickness points ; vertical staggering
! ('L') for a layer variable ; temporal staggering ('s' for snapshot) ; units ;
! and precision in non-restart output files ('f' for 32-bit float or 'd' for
! 64-bit doubles). For most tracers, only the name, longname and units should
! be changed.  See GOLD_variables for the full type description.
!
!       ALK (Total carbonate alkalinity)
!
  temp_desc = vardesc('alk','Alkalinity ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_alk, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       Cadet (Sinking detrital/particulate CaCO3)
!
  temp_desc = vardesc('cadet','Detrital CaCO3 ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_cadet, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       DIC (Dissolved inorganic carbon)
!
  temp_desc = vardesc('dic','Dissolved Inorganic Carbon', &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_dic, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       Fed (assumed to be all available to phytoplankton)
!
  temp_desc = vardesc('fed','Dissolved Iron ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_fed, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       Fedet (Sinking detrital/particulate iron)
!
  temp_desc = vardesc('fedet','Detrital Iron ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_fedet, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       FeDi (Iron in N2-fixing phytoplankton for variable Fe:N ratios)
!
  temp_desc = vardesc('fedi','Diazotroph Iron ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_fedi, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       FeLg (Iron in large phytoplankton to allow for variable Fe:N ratios)
!
  temp_desc = vardesc('felg','Large Iron ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_felg, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       FeSm(Iron in small phytoplankton to allow for variable Fe:N ratios)
!
  temp_desc = vardesc('fesm','Small Iron ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_fesm, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       LDON (Labile organic nitrogen; assumed to have a fixed, Redfield, 
!       Ketchum and Richards (1963) C:N:P ratio)
!
  temp_desc = vardesc('ldon','labile DON ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_ldon, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       LITH (Lithogenic aluminosilicate particles)
!
  temp_desc = vardesc('lith','Lithogenic Aluminosilicate ' // trim(name), &
                      'h','L','s',"g kg-1",'f')
  call set_prog_tracer(topaz%ind_lith, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       LITHdet (Detrital Lithogenic aluminosilicate particles)
!
  temp_desc = vardesc('lithdet','Detrital Lithogenic Aluminosilicate ' // trim(name), &
                      'h','L','s',"g kg-1",'f')
  call set_prog_tracer(topaz%ind_lithdet, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       Ndet (Sinking detrital/particulate Nitrogen)
!
  temp_desc = vardesc('ndet','Detrital Nitrogen ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_ndet, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       NDi (assumed to be facultative N2-fixers, with a variable N:P ratio
!
  temp_desc = vardesc('ndi','Diazotroph Nitrogen ' // trim(name) , &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_ndi, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       NH4
!
  temp_desc = vardesc('nh4','Ammonia ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_nh4, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       Heterotrophic N (assumed to be the sum of heterotrphic bacteria,
!            microzooplankton and mesozooplankton - effectively a biomass
!            storage reservoir for nutrients and carbon; assumed to have 
!            a fixed C:N:P ratio)
!
  temp_desc = vardesc('nhet','Heterotrophic Nitrogen ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_nhet, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       NLg (assumed to be a dynamic combination of diatoms and other 
!            eukaryotes all effectively greater than 5 um in diameter,
!            and having a fixed C:N ratio)
!
  temp_desc = vardesc('nlg','Large Nitrogen ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_nlg, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       NO3
!
  temp_desc = vardesc("no3","Nitrate",'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_no3, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)

!
!       NSm (Nitrogen in picoplankton and nanoplankton - effectively less
!            than 5 um ; assumed to have a fixed C:N:P ratio))
!
  temp_desc = vardesc('nsm','Small Nitrogen ' // trim(name) , &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_nsm, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       O2
!
  temp_desc = vardesc('o2','Oxygen ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_o2, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       Pdet (Sinking detrital/particulate Phosphorus)
!
  temp_desc = vardesc('pdet','Detrital Phosphorus ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_pdet, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       PDi (Phosphorus in diaz. phytoplankton for variable N:P ratios)
!
  temp_desc = vardesc('pdi','Diaz Phosphorus ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_pdi, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       PLg (Phosphorus in large phytoplankton for variable N:P ratios)
!
  temp_desc = vardesc('plg','Large Phosphorus ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_plg, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       PO4
!
  temp_desc = vardesc('po4','Phosphate ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_po4, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       PSm (Phosphorus in small phytoplankton for variable N:P ratios)
!
  temp_desc = vardesc('psm','Small Phosphorus ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_psm, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       Sidet (Sinking detrital/particulate Silicon)
!
  temp_desc = vardesc('sidet','Detrital Silicon ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_sidet, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       SiLg (Silicon in large phytoplankton for variable Si:N ratios
!
  temp_desc = vardesc('silg','Large Silicon ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_silg, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       SDON (Semilabile dissolved organic nitrogen for variable N:P ratios)
!
  temp_desc = vardesc('sdon','Semilabile DON ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_sdon, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       SDOP (Semilabile dissolved organic phosphorus for variable N:P ratios)
!
  temp_desc = vardesc('sdop','Semilabile SDOP ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_sdop, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       SiO4
!
  temp_desc = vardesc('sio4','Silicate ' // trim(name), &
                      'h','L','s',"mol kg-1",'f')
  call set_prog_tracer(topaz%ind_sio4, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       Cased (CaCO3 in sediments)
!
  temp_desc = vardesc('cased','Sediment CaCO3 ' // trim(name), &
                      'h','1','s',"mmol m-2",'f')
  call set_diag_tracer(topaz%ind_cased, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       Chl (Chlorophyll)
!
  temp_desc = vardesc('chl','Chlorophyll ' // trim(name), &
                      'h','L','s',"ug kg-1",'f')
  call set_diag_tracer(topaz%ind_chl, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       fcadet_btm (CaCO3 flux to sediments)
!
  temp_desc = vardesc('fcadet_btm','CaCO3 flux to Sediments ' // trim(name), &
                      'h','1','s',"mmol m-2 s-1",'f')
  call set_diag_tracer(topaz%ind_fcadet_btm, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       Irr (Irradiance)
!
  temp_desc = vardesc('irr','Irradiance','h','L','s',"Watts m-2",'f')
  call set_diag_tracer(topaz%ind_irr, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       Irr_mem (Irradiance Memory)
!
  temp_desc = vardesc('irr_mem','Irradiance memory','h','L','s',"Watts m-2",'f')
  call set_diag_tracer(topaz%ind_irr_mem, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       htotal (H+ ion concentration)
!
  temp_desc = vardesc('htotal','H+ ion concentration','h','L','s'," ",'f')
  call set_diag_tracer(topaz%ind_htotal, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       CO3_ion (Carbonate ion)
!
  temp_desc = vardesc('co3_ion','Carbonate ion','h','L','s',"mol kg-1 ",'f')
  call set_diag_tracer(topaz%ind_co3_ion, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       Alpha (alpha for CO2)
!
  temp_desc = vardesc('alpha','alpha for CO2','h','1','s'," ",'f')
  call set_diag_tracer(topaz%ind_alpha, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)
!
!       CO2_csurf (csurf for CO2)
!
  temp_desc = vardesc('co2_csurf','csurf for CO2','h','1','s',"mol kg-1 ",'f')
  call set_diag_tracer(topaz%ind_csurf, Tr, temp_desc, param_file, G, &
                       restart_CS, tr_adv_CSp)

  topaz%restart_file  = default_restart_file
  topaz%init     = .false.
  call read_param(param_file, "restart_file", topaz%restart_file )
  call read_param(param_file, "init",    topaz%init    )

  htotal_scale_lo = 0.01
  htotal_scale_hi = 100.0 
  htotal_in          = 1.0e-08
  call read_param(param_file, "htotal_scale_lo", htotal_scale_lo)
  call read_param(param_file, "htotal_scale_hi", htotal_scale_hi)
  call read_param(param_file, "htotal_in",htotal_in)

  if (topaz%init) then
  write(*,*) "Initializing htotal"
  
!    topaz%htotal(:,:,:) = htotal_in
    Tr(topaz%ind_htotal)%field(:,:,:) = htotal_in
  endif  !}
!
!-----------------------------------------------------------------------
! Initialize TOPAZ parameters
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! Stoichiometry
!-----------------------------------------------------------------------
!
! Values taken from OCMIP-II Biotic protocols after Anderson
! and Sarmiento (1994)
!
  topaz%c_2_n                = 106.0 / 16                     ! mol C mol N-1
  topaz%mass_2_n             = 106.0 / 16.0 * 12.0 * 1.87     ! g mol N-1
  topaz%n_2_n_denit          = 6.5                            ! dimensionless
  topaz%o2_2_c               = 150.0 / 106                    ! mol O2 mol C-1
  topaz%o2_2_nfix            = (118.0 +3/5*(150-118)) / 16    ! mol O2 mol N-1
  topaz%o2_2_no3             = 150.0 / 16                     ! mol O2 mol N-1
  topaz%o2_2_nh4             = 118.0 / 16                     ! mol O2 mol N-1
  topaz%o2_2_nitrif          = 2.0                            ! mol O2 mol N-1
  call read_param(param_file, "c_2_n",                   topaz%c_2_n)
  call read_param(param_file, "mass_2_n",                topaz%mass_2_n)
  call read_param(param_file, "n_2_n_denit",             topaz%n_2_n_denit)
  call read_param(param_file, "o2_2_c",                  topaz%o2_2_c)
  call read_param(param_file, "o2_2_no3",                topaz%o2_2_no3)
  call read_param(param_file, "o2_2_nfix",               topaz%o2_2_nfix)
  call read_param(param_file, "o2_2_nh4",                topaz%o2_2_nh4)
  call read_param(param_file, "o2_2_nitrif",             topaz%o2_2_nitrif)
!
! CaCO3 to nitrogen uptake ratio set in order to obtain a global CaCO3 export flux out of
! the euphotic zone of approximately 0.6 Pg C in CaCO3 after Sarmiento et al (2002)
!
  topaz%ca_2_n_het    = 0.015 * 106.0 / 16                    ! mol Ca mol N-1
  call read_param(param_file, "ca_2_n_het",              topaz%ca_2_n_het)
!
! Maximum limit on calcite CaCO3 saturation state for uptake
!
  topaz%caco3_sat_max    = 10.0                               ! dimensionless
  call read_param(param_file, "caco3_sat_max",           topaz%caco3_sat_max)
!
! Static small phytoplankton/coccholith/foram-grazed/pteropod-grazed CaCO3 to
! nitrogen uptake ratio for a single, globally constant value instead of a
! dynamic value.
!
  topaz%ca_2_n_static        = .false.
  call read_param(param_file, "ca_2_n_static",           topaz%ca_2_n_static)
!
! Static small phytoplankton/coccholith/foram-grazed/pteropod-grazed CaCO3 to
! nitrogen uptake ratio for a single, globally constant value instead of a
! dynamic value.
!
  topaz%ca_2_n_het_static = 0.012 * 106.0 / 16                ! mol Ca mol N-1
  call read_param(param_file, "ca_2_n_static_het",       topaz%ca_2_n_static) 
!
! P:N limitation of phytoplankton growth, is taken from the fraction of the
! cell dedicated to assembly after Klausmeier et al (2004, Optimal nitrogen-to-phosphorus
! stoichiometry of phytoplankton, Nature, 429, 171-174), but also consistent with other
! work on the limitation of growth and optimal stoichiometry by Geider et al (1998, A 
! dynamic regulatory model of phytoplanktonic acclimation to light, nutrients,
! and temperature, Limnol. Oceanogr., 43, 679-694), Christian (2005, 
! Biogeochemical cycling in the oligotrophic ocean: Redfield and non-Redfield  
! models, Limnol. Oceanogr., 50, 646-657) and Hall et al. (2005,
! Constratints on primary producer N:P stoichiometry along N:P supply ratio gradients
! Ecology, 86, 1894-1904).  p_2_n_assem is based on the N:P of ribosomes of 7.2 for
! eukaryotes (large) and 6.0 for prokaryotes (small, diaz) from Elser et al. (1996, Organism
! size, life history, and N:P stoichiometry, BioScience, 46, 674-684).
!
  phyto(DIAZO)%p_2_n_assem   = 1.0 / 6.0                      ! mol P mol N-1
  phyto(LARGE)%p_2_n_assem   = 1.0 / 7.2                      ! mol P mol N-1
  phyto(SMALL)%p_2_n_assem   = 1.0 / 6.0                      ! mol P mol N-1
  topaz%p_2_n_photo          = 0.0128                         ! mol P mol N-1
  call read_param(param_file, "p_2_n_assem_Di",          phyto(DIAZO)%p_2_n_assem)
  call read_param(param_file, "p_2_n_assem_Lg",          phyto(LARGE)%p_2_n_assem)
  call read_param(param_file, "p_2_n_assem_Sm",          phyto(SMALL)%p_2_n_assem)
  call read_param(param_file, "p_2_n_photo",             topaz%p_2_n_photo)
!
! Standard P:N value of 16 from Redfield, Ketchum and Richards (1963; The influence 
! of organisms on the composition of sea-water; The Sea, 2, 26-77).  Alternatively,
! some studies suggest a value of 15 from Goldman (1980) as 
! reprinted in Broecker and Peng (1982; Tracers in the Sea), but this was found to
! lead to excessive P limitation in preliminary tests.
!
  topaz%p_2_n_RKR            = 1.0 / 16.0                     ! mol P mol N-1
  call read_param(param_file, "p_2_n_RKR",               topaz%p_2_n_RKR)
!
! In contrast to Klausmeier et al, a finite value of
! phosphorus in uptake machinery (p_2_n_uptake) is assumed in order to simulate
! the need for ATP in active transport.  This is set to p_2_n_photo.
!
  topaz%p_2_n_uptake         = 0.0128                         ! mol P mol N-1
  call read_param(param_file, "p_2_n_uptake",            topaz%p_2_n_uptake)
!
! Calculation of the optimal P:N ratio is taken from Klausmeier et al (2004).  Allocation
! of phytoplankton cell structures to "other" is set to a minimum of 0.2 for large phytoplankton
! and 0.3 for small phytoplankton to be consistent with their different sizes and
! corresponding thetamax values.  Allocation to photosynthesis is taken from the realized
! Chlorophyll to Carbon ratio (theta) using the plastids to chlorophyll ratio.
! Note that thetamax should never be set such that r_other_min+thetamax*plastid_2_chl
! is greater than 1.
! For diazotrophs, an additional allocation is specified for nitrogenase and other
! structures necessary for nitrogen fixation, and r_other_min is correspondingly lowered
! in order to match the observed range in N:P ratios of 14-182 in White, A. E., Y. H 
! Spitz, D. M. Karl and R. M. Letelier (2006, Flexible elemental stoichiometry in 
! Trichodesmium spp. and its ecological implications. Limnol Oceanogr, 51, 1777-1790.
!
! The relationship between chlorophyll and carbon is taken from J. T. O. Kirk and R. A. E.
! Tilney-Bassett (1978, The Plastids: Their chemistry, structure, growth and inheritance.
! Revised second edition, Elsevier/North-Holland Biomedical Press, New York.)
! for eukaryotes (large) that have chloroplasts and 6.0 for prokaryotes (small, diaz) which
! only have the thylakoid.

!  topaz%plast_2_chl          = 1.0 / 1.87 / 0.07              ! g C g Chl-1
  phyto(DIAZO)%plast_2_chl   = 1.0 / 1.87 /0.08               ! g C g Chl-1
  phyto(LARGE)%plast_2_chl   = 1.0 / 1.87 /0.043              ! g C g Chl-1
  phyto(SMALL)%plast_2_chl   = 1.0 / 1.87 /0.08               ! g C g Chl-1
  phyto(DIAZO)%r_nfix        = 0.5                            ! dimensionless
  phyto(LARGE)%r_nfix        = 0.0                            ! dimensionless
  phyto(SMALL)%r_nfix        = 0.0                            ! dimensionless
  phyto(DIAZO)%r_other_min   = 0.1                            ! dimensionless
  phyto(LARGE)%r_other_min   = 0.2                            ! dimensionless
  phyto(SMALL)%r_other_min   = 0.2                            ! dimensionless
  phyto(DIAZO)%r_uptake_max  = 0.2                            ! dimensionless
  phyto(LARGE)%r_uptake_max  = 0.2                            ! dimensionless
  phyto(SMALL)%r_uptake_max  = 0.4                            ! dimensionless
!  call read_param(param_file, "plast_2_chl",             topaz%plast_2_chl)
  call read_param(param_file, "plast_2_chl_Di",   phyto(DIAZO)%plast_2_chl)
  call read_param(param_file, "plast_2_chl_Lg",   phyto(LARGE)%plast_2_chl)
  call read_param(param_file, "plast_2_chl_Sm",   phyto(SMALL)%plast_2_chl)
  call read_param(param_file, "r_nfix_Di",        phyto(DIAZO)%r_nfix)
  call read_param(param_file, "r_nfix_Lg",        phyto(LARGE)%r_nfix)
  call read_param(param_file, "r_nfix_Sm",        phyto(SMALL)%r_nfix)
  call read_param(param_file, "r_other_min_Di",   phyto(DIAZO)%r_other_min)
  call read_param(param_file, "r_other_min_Lg",   phyto(LARGE)%r_other_min)
  call read_param(param_file, "r_other_min_Sm",   phyto(SMALL)%r_other_min)
  call read_param(param_file, "r_uptake_max_Di",  phyto(DIAZO)%r_uptake_max)
  call read_param(param_file, "r_uptake_max_Lg",  phyto(LARGE)%r_uptake_max)
  call read_param(param_file, "r_uptake_max_Sm",  phyto(SMALL)%r_uptake_max)
!
! Static phosphorus to nitrogen uptake ratio for a single, globally constant value for each
! phytoplankton group instead of dynamic values for each.
!
  topaz%p_2_n_static         = .false.
  call read_param(param_file, "p_2_n_static",            topaz%p_2_n_static)
!
! If static, Diazotrophs are assumed to have low P:N after Letelier, R. M. & Karl, 
! D. M. Role of Trichodesmium spp. in the productivity of the subtropical North
! Pacific Ocean. Mar. Ecol. Prog. Ser. 133, 263-273 (1996).
!
  phyto(DIAZO)%p_2_n_static  = 1.0 / 40.0                     ! mol P mol N-1
  call read_param(param_file, "p_2_n_static_Di",  phyto(DIAZO)%p_2_n_static)
!
! If static, Large and Small are assumed to have RKR P:N.
!
  phyto(LARGE)%p_2_n_static  = 1.0 / 16.0                     ! mol P mol N-1
  phyto(SMALL)%p_2_n_static  = 1.0 / 16.0                     ! mol P mol N-1
  call read_param(param_file, "p_2_n_static_Lg",  phyto(LARGE)%p_2_n_static)
  call read_param(param_file, "p_2_n_static_Sm",  phyto(SMALL)%p_2_n_static)
!
! Maximum diatom silicon to nitrogen uptake ratio after Mongin, M., D. M. Nelson,
! P. Pondaven, M. A. Brzezinski and P. Treguer (2003): Simulation of upper-ocean
! biogeochemistry with a flexible-composition phytoplankton model: C, N and Si
! cycling in the western Sargasso Sea. Deep-Sea Res. I, 50, 1445-1480.
!
  phyto(LARGE)%si_2_n_max    = 5.0                            ! mol Si mol N-1
  call read_param(param_file, "si_2_n_max_Lg",    phyto(LARGE)%si_2_n_max)
!
! Static large phytoplankton/diatom silicon to nitrogen uptake ratio for a single,
! globally constant value instead of a dynamic value.
!
  topaz%si_2_n_static        = .false.
  call read_param(param_file, "si_2_n_static",           topaz%si_2_n_static)
!
! Value of the static large phytoplankton/diatom silicon to nitrogen uptake ratio for
! a single, globally constant value instead of a dynamic value.  Rather than the
! canonical value of 1.0 after culture work of Bzezinski (1985), this 
! value is set at the iron-stressed value of 2.0 after Hutchins, D. A. and K. W.
! Bruland (1998): Iron-limited diatom growth and Si:N uptake ratios in a coastal
! upwelling regime. Nature, 393, 561-564 and Takeda, S. (1998): Influence of iron
! variability on nutrient consumption ratio of diatoms in oceanic waters. Nature, 393,
! 774-777.
!
  phyto(LARGE)%si_2_n_static = 2.0                            ! mol Si mol N-1
  call read_param(param_file, "si_2_n_static_Lg", phyto(LARGE)%si_2_n_static)
!
!-----------------------------------------------------------------------
! Monod half saturation coefficients.  k_PO4 for small phytoplankton was taken
! from the minimum of 5 nM observed in Cotner JB, Ammerman JA, Peele ER, Bentzen 
! E (1997)  Phosphorus-limited bacterioplankton growth in the Sargasso Sea. 
! Aquat Microb Ecol 13:141.), and was assumed to follow the effective size
! relationship of Lg/Sm = 3.  k_nh4 was assumed to be equal to k_po4.  k_no3 was
! assumed to be a factor of ten larger than k_nh4 to allow nh4 inhibition of
! no3 uptake after the formulation of Sharada and Yajnik (2005;
! Evaluation of six relations of the kinetics of uptake by phytoplankton in
! multi-nutrient environment using JGOFS experimental results. DSR II, 52,
! 1892-1909.) and Frost and 
! Franzen (1992).  k_SiO4 taken from Dugdale, R.C. and Wilkerson, F.P. (1998, 
! Silicate regulation of new production in the equatorial Pacific upwelling.
! Nature, 391, 270-273).  
!-----------------------------------------------------------------------
!
  phyto(DIAZO)%k_fed         = 15.0e-9                        ! mol Fed/kg
  phyto(LARGE)%k_fed         = 15.0e-9                        ! mol Fed/kg
  phyto(SMALL)%k_fed         = 5.0e-9                         ! mol Fed/kg
  phyto(LARGE)%k_nh4         = 6.0e-7                         ! mol NH4/kg
  phyto(SMALL)%k_nh4         = 2.0e-7                         ! mol NH4/kg
  phyto(LARGE)%k_no3         = 6.0e-6                         ! mol NO3/kg
  phyto(SMALL)%k_no3         = 2.0e-6                         ! mol NO3/kg
  phyto(DIAZO)%k_po4         = 6.0e-7                         ! mol PO4/kg
  phyto(LARGE)%k_po4         = 6.0e-7                         ! mol PO4/kg
  phyto(SMALL)%k_po4         = 2.0e-7                         ! mol PO4/kg
!  phyto(LARGE)%k_sio4        = 3.0e-6                         ! mol SiO4/kg
!## jgj NOTE: k_sio4_Sm = 1.0e-6 in mom4p1
  phyto(LARGE)%k_sio4        = 1.0e-6                         ! mol SiO4/kg
  call read_param(param_file, "k_fed_Di",         phyto(DIAZO)%k_fed)
  call read_param(param_file, "k_fed_Lg",         phyto(LARGE)%k_fed)
  call read_param(param_file, "k_fed_Sm",         phyto(SMALL)%k_fed)
  call read_param(param_file, "k_nh4_Lg",         phyto(LARGE)%k_nh4)
  call read_param(param_file, "k_nh4_Sm",         phyto(SMALL)%k_nh4)
  call read_param(param_file, "k_no3_Lg",         phyto(LARGE)%k_no3)
  call read_param(param_file, "k_no3_Sm",         phyto(SMALL)%k_no3)
  call read_param(param_file, "k_po4_Di",         phyto(DIAZO)%k_po4)
  call read_param(param_file, "k_po4_Lg",         phyto(LARGE)%k_po4)
  call read_param(param_file, "k_po4_Sm",         phyto(SMALL)%k_po4)
  call read_param(param_file, "k_sio4_Lg",        phyto(LARGE)%k_sio4)
!
!-----------------------------------------------------------------------
! Iron
!-----------------------------------------------------------------------
!
! Whether or not to allow mineral ballast dissolution to
! return iron to the dissolved phase - a "false" value
! assumes that all iron is associated with organic material.
! A true value assumes that iron is distributed between
! mineral and organic matter by mass leading to a deeper
! regeneration length scale.
!
  topaz%fe_ballast_assoc     = .true. 
  call read_param(param_file, "fe_ballast_assoc",        topaz%fe_ballast_assoc)
!
! Background concentration of iron ligand of 2.0e-9 taken from Rue, E. L. and 
! K. W. Bruland (1995) Mar. Chem., 50, 117-138.  Alternatively, a value of 
! 1.0e-9 can be taken from Parekh, P., M. J. Follows and E. A. Boyle (2005) 
! Decoupling of iron and phosphate in the global ocean. Glob. Biogeochem. 
! Cycles, 19,  doi: 10.1029/2004GB002280.
!
  topaz%felig_bkg            = 2.0e-9                         ! mol Fe kg-1
  call read_param(param_file, "felig_bkg", topaz%felig_bkg)
!
! Ratio of iron ligand to semilabile and labile don taken from the ratio of
! ligand to dissolved organic carbon in deep water.
!
  topaz%felig_2_don          = 2.0e-3 / 40.0 * 106.0 / 16.0   ! mol Fe mol N-1
  call read_param(param_file, "felig_2_don",             topaz%felig_2_don)
!
! Iron limitation of the Chl:C, through the def_fe factor, to allow
! iron to modulate diazotrophic phytoplankton light utilization
! efficiency. This value is set to a relatively high value after Raven et al.
!
! Iron limitation of the Chl:C, through the def_fe factor, to allow
! iron to modulate small phytoplankton light utilization efficiency.
! This value is set to a very low value after Sunda and Huntsman (1995, 
! Mar. Chem, 50, 189-206) for Small phytoplankton and is assumed to
! have an implicit allometric relationship of eff_size_Lg_2_Sm.
!
! values raised to stimulate expression of iron limitation.

  phyto(DIAZO)%k_fe_2_n      = 60.0e-6 * 106.0 / 16.0         ! mol Fe mol N-1
  phyto(LARGE)%k_fe_2_n      = 30.0e-6 * 106.0 / 16.0         ! mol Fe mol N-1
  phyto(SMALL)%k_fe_2_n      = 10.0e-6 * 106.0 / 16.0         ! mol Fe mol N-1
  call read_param(param_file, "k_fe_2_n_Di",      phyto(DIAZO)%k_fe_2_n)
  call read_param(param_file, "k_fe_2_n_Lg",      phyto(LARGE)%k_fe_2_n)
  call read_param(param_file, "k_fe_2_n_Sm",      phyto(SMALL)%k_fe_2_n)
!
! Maximum Fe:N level where uptake ceases for Small Phytoplankton...
! that is, where the phytoplankton get "full" of iron.  This maximum
! value was set using the "fe replete" value of 23 umol Fe / mol C for
! Synechococcus observed by Porta et al (2003; J. hycol, 39, 64-73)
! Since the Droop function diminishes to zero at half the maximum, the
! Droop maximum is set to twice the observed maximum.
!
  phyto(SMALL)%fe_2_n_max    = 46.e-6 * 106.0 / 16.0          ! mol Fe mol N-1
  call read_param(param_file, "fe_2_n_max_Sm",    phyto(SMALL)%fe_2_n_max)
!
! Maximum Fe:N level where uptake ceases for Large Phytoplankton...
! that is, where the phytoplankton get "full" of iron.  This maximum
! value was set from a "fe replete" value of 333 umol Fe / mol C for the
! coast diatom T weissflogii observed by Sunda and Huntsman (1995, 
! Mar. Chem, 50, 189-206).
!
  phyto(LARGE)%fe_2_n_max    = 666.0e-6 * 106.0 / 16.0        ! mol Fe mol N-1
  phyto(DIAZO)%fe_2_n_max    = 666.0e-6 * 106.0 / 16.0        ! mol Fe mol N-1
  call read_param(param_file, "fe_2_n_max_Lg",    phyto(LARGE)%fe_2_n_max)
  call read_param(param_file, "fe_2_n_max_Di",    phyto(DIAZO)%fe_2_n_max)
!
! Ratio of iron influx from bottom sediment boundaries based on nitrogen flux
!
!  topaz%fe_2_n_sed           = 5.0e-6 * 117/16               ! mol Fe mol N-1
!## jgj NOTE: to match mom4p1
  topaz%fe_2_n_sed           = 100.0e-6 * 106.0 / 16.0         ! mol Fe mol N-1
  call read_param(param_file, "fe_2_n_sed",              topaz%fe_2_n_sed)
!
! Static iron to nitrogen uptake ratio for a single, globally constant value for each
! phytoplankton group instead of dynamic values for each.
!
  topaz%fe_2_n_static = .false.
  call read_param(param_file, "fe_2_n_static",           topaz%fe_2_n_static)
!
! If static, Diazotrophs are assumed to have a Fe:N after the value of the Fe:C ratio
! necessary to achieve half the maximum growth rate of Trichodesmium in
! Berman-Frank, I., J. T. Cullen, Y. Shaked, R. M. Sherrell, and P. G. Falkowski
! (2001) Iron availability, cellular iron quotas, and nitrogen fixation in
! Trichodesmium. Limnol. Oceanogr., 46, 1249-1260.
!
!  call fm_util_set_value('fe_2_n_static_Di', 30.0e-6 * 106.0 / 16.0) ! mol Fe mol N-1
!
! If static, large and small phytoplankton are assumed to have an Fe:N after the value of
! the Fe:C ratio necessary to achieve half the maximum growth rate of T. oceanica in 
! Sunda and Huntsman (1995) Iron uptake and growth limitation in oceanic and coastal 
! phytoplankton. Mar. Chem., 50, 189-196.
!
!  call fm_util_set_value('fe_2_n_static_Lg', 3.0e-6 * 106.0 / 16.0 ) ! mol Fe mol N-1
!  call fm_util_set_value('fe_2_n_static_Sm', 3.0e-6 * 106.0 / 16.0 ) ! mol Fe mol N-1
!
! values raised to stimulate expression of iron limitation.
!
  phyto(DIAZO)%fe_2_n_static  = 60.0e-6 * 106.0 / 16.0         ! mol Fe mol N-1
  phyto(LARGE)%fe_2_n_static  = 45.0e-6 * 106.0 / 16.0         ! mol Fe mol N-1
  phyto(SMALL)%fe_2_n_static  = 15.0e-6 * 106.0 / 16.0         ! mol Fe mol N-1
  call read_param(param_file, "fe_2_n_static_Di",  phyto(DIAZO)%fe_2_n_static)
  call read_param(param_file, "fe_2_n_static_Lg",  phyto(LARGE)%fe_2_n_static)
  call read_param(param_file, "fe_2_n_static_Sm",  phyto(SMALL)%fe_2_n_static)
!
! Rate kinetics of iron influx from side boundaries
!
!  topaz%fe_coast             = 1.5e-12                         ! mol Fe m kg-1 s-1
!## jgj NOTE: to match mom4p1
  topaz%fe_coast             = 1.0e-11                         ! mol Fe m kg-1 s-1
  call read_param(param_file, "fe_coast",                topaz%fe_coast)
!
! Second-order iron scavenging in order to prevent high iron
! accumulations in high deposition regions (like the tropical
! Atlantic).
!
  topaz%kfe_2nd_order        = 1.0e10/sperd                    ! mol Fe-1 kg s-1
  call read_param(param_file, "kfe_2nd_order",           topaz%kfe_2nd_order)
!
! Adsorption rate coefficient for ballast.  This can be set to a low value
! to prevent iron from accumulating in the deep ocean and keep a no3-like
! profile instead of a sio4-like profile.
!
  topaz%kfe_bal             = 0.0/sperd                        ! g ballast-1 m3 s-1
  call read_param(param_file, "kfe_bal",                topaz%kfe_bal)
!
! Desorption rate coefficient.  Set to 0.0068 d-1 after Bacon and Anderson 
! (1982) J. Geophys. Res., 87, No. C3, 2045-2056. and Clegg and Whittfield 
! (1993) Deep-Sea Res., 40, 1529-1545 from Thorium.
!
  topaz%kfe_des              = 0.0068/sperd                    ! s-1
  call read_param(param_file, "kfe_des",                 topaz%kfe_des)
!
! Equilibrium constant for (free and inorganically bound) iron binding with 
! organic ligands taken from Parekh, P., M. J. Follows and E. A. Boyle (2005) 
! Decoupling of iron and phosphate in the global ocean. Glob. Biogeochem. 
! Cycles, 19, doi: 10.1029/2004GB002280.
!
  topaz%kfe_eq_lig           = 1.0e11                          ! mol lig-1 kg
  call read_param(param_file, "kfe_eq_lig",              topaz%kfe_eq_lig)
!
! Adsorption rate coefficient for detrital organic material.  This was set
! to a low value to prevent iron from accumulating in the deep ocean and
! keep a no3-like profile instead of a sio4-like profile.
  topaz%kfe_org             = 0.0/sperd                       ! g org-1 m3 s-1
  call read_param(param_file, "kfe_org",                topaz%kfe_org)
!
!-----------------------------------------------------------------------
! Photosynthesis
!-----------------------------------------------------------------------
!
! Phytoplankton growth altered from Geider et al (1997)
! and Moore et al (2002).  Thetamax values
! are at the high end in order to account for the additional
! iron limitation term.  The factor of 6.022e17 is to convert
! from umol to quanta and 2.77e18 to convert from quanta/sec
! to Watts given the average energy spectrum for underwater
! PAR from the Seabird sensor.  Values of P_C_max are decreased
! relative to Geider et al. (1997) by a factor of 4 to account
! the difference in reference temperatures and, for Small and Large
! increased by a factor of 5 to account for the addition of po4 limitation.
! Values of thetamax are increased relative to Geider et al. (1997)
! by a factor of 2 to account for the def_fe factor.
!
! alpha is assumed to not have a size relationship
!
! theta_max is assumed to have an implicit allometric relationship of 
! eff_size**(2/3) from the surface area to volume relationship.
!
  phyto(DIAZO)%alpha         = 1.0e-5 * 2.77e18 / 6.022e17    ! g C g Chl-1 m2 W-1 s-1
  phyto(LARGE)%alpha         = 2.0e-5 * 2.77e18 / 6.022e17    ! g C g Chl-1 m2 W-1 s-1
  phyto(SMALL)%alpha         = 2.0e-5 * 2.77e18 / 6.022e17    ! g C g Chl-1 m2 W-1 s-1
  topaz%kappa_eppley         = 0.063                          ! deg C-1
  phyto(DIAZO)%P_C_max       = 0.6e-5                         ! s-1
  phyto(LARGE)%P_C_max       = 1.5e-5                         ! s-1
  phyto(SMALL)%P_C_max       = 2.0e-5                         ! s-1
  phyto(DIAZO)%thetamax      = 0.040                          ! g Chl g C-1
  phyto(LARGE)%thetamax      = 0.060                          ! g Chl g C-1
  phyto(SMALL)%thetamax      = 0.040                          ! g Chl g C-1
  topaz%thetamin             = 0.005                          ! g Chl g C-1
  topaz%zeta                 = 0.1                            ! dimensionless
  call read_param(param_file, "alpha_Di",         phyto(DIAZO)%alpha)
  call read_param(param_file, "alpha_Lg",         phyto(LARGE)%alpha)
  call read_param(param_file, "alpha_Sm",         phyto(SMALL)%alpha)
  call read_param(param_file, "kappa_eppley",            topaz%kappa_eppley)
  call read_param(param_file, "P_C_max_Di",       phyto(DIAZO)%P_C_max)
  call read_param(param_file, "P_C_max_Lg",       phyto(LARGE)%P_C_max)
  call read_param(param_file, "P_C_max_Sm",       phyto(SMALL)%P_C_max)
  call read_param(param_file, "thetamax_Di",      phyto(DIAZO)%thetamax)
  call read_param(param_file, "thetamax_Lg",      phyto(LARGE)%thetamax)
  call read_param(param_file, "thetamax_Sm",      phyto(SMALL)%thetamax)
  call read_param(param_file, "thetamin",                topaz%thetamin)
  call read_param(param_file, "zeta",                    topaz%zeta)
!
! Diazotrophs are assumed to be inhibited by nitrate after Holl and Montoya 
! (2005;  Interactions between nitrate uptake and nitrogen fixation in 
! continuous cultures of the marine diazotroph trichodesmium cyanobacteria); 
! J. Phycol, 41, 1178-1183).
!
  topaz%k_n_inhib_Di     = 7.0e-6                            ! mol NO3 kg-1
  call read_param(param_file, "k_n_inhib_Di",        topaz%k_n_inhib_Di)
!

! Diazotrophs are also assumed to be inhibited by oxygen after Stewart and
! Pearson (1970; Effects of aerobic and anaerobic conditions on growth and
! metabolism of blue-green algae. Proc. Soc. Lond. B., 175, 293-311) and
! Berman-Frank et al (2005; Inhibition of nitrogenase by oxygen in marine
! cyanobacteria controls the global nitrogen and  oxygen cycles. Biogeosciences
! Discussions, 2, 261-273).
!
  topaz%o2_inhib_Di_pow      = 4.0                            ! dimensionless
  topaz%o2_inhib_Di_sat      = 3.0e-4                         ! mol O2 kg-1
  call read_param(param_file, "o2_inhib_Di_pow",         topaz%o2_inhib_Di_pow)
  call read_param(param_file, "o2_inhib_Di_sat",         topaz%o2_inhib_Di_sat)
!
! Chl:C response rate constant for phytoplankton calibrated to 1 d-1
! after Owens et al (1980, Diel Periodicity in cellular Chlorophyll
! content of marine diatoms, Mar. Biol, 59, 71-77).
!
  topaz%gamma_irr_mem        = 1.0 / sperd                    ! s-1
  call read_param(param_file, "gamma_irr_mem",           topaz%gamma_irr_mem)
!
!-----------------------------------------------------------------------
! Grazing
!-----------------------------------------------------------------------
!
! Values of fractional detritus production from the global
!  synthesis of Dunne et al. (submitted)
!
  phyto(LARGE)%fdet0         = 0.93                           ! dimensionless
  phyto(SMALL)%fdet0         = 0.18                           ! dimensionless
  call read_param(param_file, "fdet0_Lg",         phyto(LARGE)%fdet0)
  call read_param(param_file, "fdet0_Sm",         phyto(SMALL)%fdet0)
!
! Rate constant for remineralization of heterotrophic biomass
! 
  topaz%gamma_nhet           = 1.0 / (30.0 * sperd)           ! s-1
  call read_param(param_file, "gamma_nhet",              topaz%gamma_nhet)
!
! Dissolution of sio2 was set as a temperature-dependent
! fraction of grazed material to be roughly in line with
! the work of Kamatani (1982)
!
  topaz%k_diss_sio2          = 3.0                            ! s-1
  call read_param(param_file, "k_diss_sio2",             topaz%k_diss_sio2)
!
! Half saturation of oxic remineralization rate.
!
  topaz%k_o2                 = 20.0e-6                        ! mol O2/kg
  call read_param(param_file, "k_o2",                    topaz%k_o2)
!
! Temperature-dependence of fractional detritus production
! from the global synthesis of Dunne et al. (submitted)
!
  topaz%kappa_remin          = -0.032                         ! deg C-1
  call read_param(param_file, "kappa_remin",             topaz%kappa_remin)
!
! T=0 phytoplankton specific grazing rate from the global
! synthesis of Dunne et al. (2005)
!
  topaz%lambda0              = 0.19 / sperd                   ! s-1
  call read_param(param_file, "lambda0",                 topaz%lambda0)
!
! Minimum oxygen concentration for oxic remineralization.
! this is necessary for both numerical stability and to
! qeue the switch to denitrification
!
!  topaz%o2_min               = 2.0 * 1.0e-06                  ! mol O2 kg-1
!## jgj NOTE: to match mom4p1
  topaz%o2_min               = 1.0 * 1.0e-06                  ! mol O2 kg-1
  call read_param(param_file, "o2_min",                  topaz%o2_min)
!
! Pivot phytoplankton concentration for grazing-based
! variation in ecosystem structure from the global
! synthesis of Dunne et al. (submitted)
!
  topaz%P_star               = 1.9e-6 * 16.0 / 106.0          ! mol N kg-1
  call read_param(param_file, "P_star",                  topaz%P_star)
!
! Minimum phytoplankton concentration or grazing.  This is
! necessary for numerical stability.
! 
  topaz%phyto_min            = 1.0e-10                        ! mol N kg-1
  call read_param(param_file, "phyto_min",               topaz%phyto_min)
!
! Fraction of phytoplankton grazing and dom consumption eventually going to NH4
! that is temporarily stored in heterotrophic biomass
! 
  topaz%phi_nhet             = 0.75                           ! dimensionless
  call read_param(param_file, "phi_nhet",                topaz%phi_nhet)
!
! SiO2 dissolution is set to globally dissolve 50% after Nelson et al. (1995) 
! through grazing.  The temperature functionality is set to a combination of
! stoichiometry and Eppley temperature formulation to give roughly the range of 
! observations in Kamatani (1982) with respect to frustrule thickness and
! temperature by utilizing the inverse of the Eppley temperature
! functionality and a normalization to stoichiometry ( q_si_2_n_diss).
! The value of q_si_2_n_diss was set so as to simultaneously reproduce the low
! silicon export efficiencies (~0.1) observed in the equatorial Pacific by
! Blain et al. (1997, DSR I; Dunne et al., 1999, GBC) and high export efficiencies
! of ~0.64 observed in the Southern Ocean by Bzrezinski et al., 2001, retaining
! a ~0.5 global average after Nelson et al. (1995).
!
  topaz%q_si_2_n_diss        = 3.0                            ! mol mol Si mol N
  call read_param(param_file, "q_si_2_n_diss",                 topaz%rpcaco3)
!
! Organic matter protection by mineral - after Klaas and
! Archer (2002)
!
  topaz%rpcaco3              = 0.070/12*16/106.0*100          ! mol N mol Ca-1
  topaz%rplith               = 0.065/12*16/106.0              ! mol N g lith-1
  topaz%rpsio2               = 0.026/12*16/106.0*60           ! mol N mol Si-1
  call read_param(param_file, "rpcaco3",                 topaz%rpcaco3)
  call read_param(param_file, "rplith",                  topaz%rplith )
  call read_param(param_file, "rpsio2",                  topaz%rpsio2 )
!
!-----------------------------------------------------------------------
! Remineralization length scales
!-----------------------------------------------------------------------
!
! Sinking velocity of detritus - 20 m d-1 consistent with a characteristic sinking
! velocity of 100 m d-1 of marine aggregates and a disaggregation rate constant
! of 5 d-1 (Clegg and Whitfield, 1992; Dunne, 1999)
!
!## jgj NOTE: gamma_ndet, gamma_cadet, gamma_sidet are tagged to wsink (currently 20/sperd)
!## jgj NOTE: mom4p1 default is currently 10.0 / sperd
  topaz%wsink                = 20.0 / sperd                  ! m s-1
  call read_param(param_file, "wsink",                   topaz%wsink)

  Tr(:)%sink_dist = 0.0;
  Tr(topaz%ind_cadet)%sink_dist = topaz%wsink                !20.0/sperd
  Tr(topaz%ind_fedet)%sink_dist = topaz%wsink                !20.0/sperd
  Tr(topaz%ind_lithdet)%sink_dist = topaz%wsink              !20.0/sperd
  Tr(topaz%ind_ndet)%sink_dist = topaz%wsink                 !20.0/sperd
  Tr(topaz%ind_pdet)%sink_dist = topaz%wsink                 !20.0/sperd
  Tr(topaz%ind_sidet)%sink_dist = topaz%wsink                !20.0/sperd
!
! Value of gamma_ndet to approximate upper e-folding of the globally-tuned
! "Martin curve" used in the OCMIP-II Biotic configuration of (z/75)^-0.9
! that gives a value of exp(-1) at 228 m from 75 m for an e-folding scale
! of 188 m.
!
!!  topaz%gamma_ndet      = 20.0 / sperd / 188.0               ! s-1
  topaz%gamma_ndet      = topaz%wsink / 188.0                ! s-1
  call read_param(param_file, "gamma_ndet",         topaz%gamma_ndet)

! Cadet dissolution rate constant: Most deep traps show little correlation
! with depth, suggesting little water column dissolution in much of the
! ocean.  Values were calibrated to get the observed 67% transfer efficiency at
! Ocean Weather Station PAPA in the North Pacific between 1000-3800 m
! equivalent to a (1-omega)-modulated dissolution length scale of 1343 m
! based on the Honjo data on the JGOFS database assuming an omega of 0.81 from
! GLODAP and converted to a dissolution rate constant assuming a sinking
! velocity of 20 m d-1.
!
!!  topaz%gamma_cadet       = 20.0 / sperd / 1343.0            ! s-1
  topaz%gamma_cadet       = topaz%wsink / 1343.0             ! s-1
!
! Sidet dissolution rate constant assuming a dissolution length scale of 2000 m
! consistent with Gnanadesikan (2000) and assuming a sinking velocity of 20 m d-1.
!
!!  topaz%gamma_sidet       = 20.0 / sperd / 2000.0            ! s-1
  topaz%gamma_sidet       = topaz%wsink / 2000.0             ! s-1
  call read_param(param_file, "gamma_cadet",          topaz%gamma_cadet)
  call read_param(param_file, "gamma_sidet",          topaz%gamma_sidet)
!
!-----------------------------------------------------------------------
! CaCO3 dissolution
!-----------------------------------------------------------------------
!
! Coefficients for calcite solubility taken from Sayles, F. L. (1985, CaCO3
! solubility in marine sediments: evidence for equilibrium and non-equilibrium
! behavior, Geochim. Cosmochim. Acta, 49, 877-888)
!
  topaz%ksp_caco3            = 4.95e-7                       ! mol2 kg-2
  topaz%k_caco3_pres         = 1.0e-17                       ! mol2 m-4 s-2
!
! Redissolution of previously-deposited CaCO3 sediments
!
  topaz%gamma_cased_dis       = 1.0e-3 / spery               ! s-1
  call read_param(param_file, "ksp_caco3",               topaz%ksp_caco3     )
  call read_param(param_file, "k_caco3_pres",            topaz%k_caco3_pres  )
  call read_param(param_file, "gamma_cased_dis",         topaz%gamma_cased_dis)
!
!-----------------------------------------------------------------------
! Dissolved Organic Material
!-----------------------------------------------------------------------
!
!
! Dissolved Organic Material remineralization rate constants
! and fractional production ratios, all to be consistent
! with the work of Abell et al. (2000, Distributions of TOP, TON and TOC
! in the North pacific subtropical gyre: Implications for nutrient supply
! in the surface ocean and remineralization in the upper thermocline,
! J. Mar. Res., 58, 203-222) and DOC and DON data provided by
! Dennis Hansell (personal communication).  Assuming refractory/deep
! concentrations of DOC=42 uM, DON=1.8 uM, and DOP=0.0 uM, we allow the
! semilabile and labile pools to have RKR C:N and reproduce observed
! surface expression in north Pacific subtropical gyre (DOC=72 uM, 
! DON=6 uM, and DOP=0.2 uM) and depth penetration. 
!
! Warning: phi_sdon + phi_ldon should be < 1.0.
!
!
!  topaz%gamma_sdon           = 1.0 / (18.0 *365.0 * sperd)    ! s-1
!  topaz%gamma_sdop           = 1.0 / (4.0 *365.0 * sperd)     ! s-1
!## jgj NOTE: to match mom4p1
  topaz%gamma_sdon           = 1.0 / (18.0 * spery)            ! s-1
  topaz%gamma_sdop           = 1.0 / (4.0 * spery)             ! s-1
  topaz%phi_sdon             = 0.02                           ! dimensionless
  topaz%phi_sdop             = 0.04                           ! dimensionless
  call read_param(param_file, "gamma_sdon",              topaz%gamma_sdon)
  call read_param(param_file, "gamma_sdop",              topaz%gamma_sdop)
  call read_param(param_file, "phi_sdon",                topaz%phi_sdon)
  call read_param(param_file, "phi_sdop",                topaz%phi_sdop)
!
! The remineralization rate constant for labile DOP (bio_tau_ldon)
! was set to 3 months after Archer et al. (1997, GBC, 11, 435-452).
! The fraction going to labile DOC was inspired by data-model
! comparisons to Libby and Wheeler (1997, Deep-Sea Res. I, 44, 345-361)
!
  topaz%gamma_ldon           = 1.0 / (90.0 * sperd)           ! s-1
  topaz%phi_ldon             = 0.06                           ! dimensionless
  call read_param(param_file, "gamma_ldon",              topaz%gamma_ldon)
  call read_param(param_file, "phi_ldon",                topaz%phi_ldon)
!
!-----------------------------------------------------------------------
! Miscellaneous
!-----------------------------------------------------------------------
!
! Debug flag to calculate global integrals for tracers
!
  topaz%tracer_debug         = .false.
  call read_param(param_file, "tracer_debug",            topaz%tracer_debug)
  topaz%tracer_debug_verbose = .false.
  call read_param(param_file, "tracer_debug_verbose",    topaz%tracer_debug_verbose)
!
! Nitrification rate constant assumed to be light-limited with an inhibition
! factor.  gamma_nitrif was tuned to reproduce the scaling observed in Ward et
! al. (1982; Microbial nitrification rates in the  primary nitrite maximum off
! southern California, Deep-Sea Res., 29, 247-255), and irr_inhibit was tuned to
! reproduce Olson (1981; Differential photoinhibition of marine nitrifying
! bacteria: a possible mechanism for the formulation of the primary nitrite
! maximum, J. Mar. Res., 39, 227-238).
!
  topaz%gamma_nitrif         = 1.0 / (30.0 * sperd)           ! s-1
  topaz%irr_inhibit          = 2.0                            ! m2 W-1
  call read_param(param_file, "gamma_nitrif",            topaz%gamma_nitrif)
  call read_param(param_file, "irr_inhibit",             topaz%irr_inhibit )
!
! Scavenging rate coefficient for lithogenic material relative to large
! phytoplankton concentration via large phytoplankton grazing.
!
  topaz%phi_lith             = 0.002                          ! dimensionless
  call read_param(param_file, "phi_lith",                topaz%phi_lith)
!
! Adsorption rate coefficient for lithogenic material onto sinking material.
! This was set to a small but non-zero to prevent lithogenic
! material from accumulating in the deep ocean.
!
!  topaz%k_lith               = 10.0/365.0/sperd               ! s-1
!## jgj NOTE: to match mom4p1
  topaz%k_lith               = 1e-6/sperd                      ! s-1
  call read_param(param_file, "k_lith",                  topaz%k_lith)

!
! Constant for productivity mask
!
  topaz%bio_tau              = 30.0 * sperd                   ! s
  call read_param(param_file, "bio_tau", topaz%bio_tau)

  topaz%r_bio_tau = 1.0 / topaz%bio_tau

  topaz%nkml = 1 ; call read_param(param_file,"NKML",topaz%nkml)
  call read_param(param_file,"RHO_0",topaz%Rho_0,.true.)

  mod = "GOLD_OCEAN_TOPAZ"
  call log_version(param_file, mod, version, tagname)
  call log_param(param_file, mod, "restart_file", topaz%restart_file )
  call log_param(param_file, mod, "init", topaz%init )
  call log_param(param_file, mod, "htotal_scale_lo", htotal_scale_lo)
  call log_param(param_file, mod, "htotal_scale_hi", htotal_scale_hi)
  call log_param(param_file, mod, "htotal_in",htotal_in)
  call log_param(param_file, mod, "c_2_n", topaz%c_2_n)
  call log_param(param_file, mod, "mass_2_n", topaz%mass_2_n)
  call log_param(param_file, mod, "n_2_n_denit", topaz%n_2_n_denit)
  call log_param(param_file, mod, "o2_2_c", topaz%o2_2_c)
  call log_param(param_file, mod, "o2_2_no3", topaz%o2_2_no3)
  call log_param(param_file, mod, "o2_2_nfix", topaz%o2_2_nfix)
  call log_param(param_file, mod, "o2_2_nh4", topaz%o2_2_nh4)
  call log_param(param_file, mod, "o2_2_nitrif", topaz%o2_2_nitrif)
  call log_param(param_file, mod, "ca_2_n_het", topaz%ca_2_n_het)
  call log_param(param_file, mod, "caco3_sat_max", topaz%caco3_sat_max)
  call log_param(param_file, mod, "ca_2_n_static", topaz%ca_2_n_static)
  call log_param(param_file, mod, "ca_2_n_static_het", topaz%ca_2_n_static) 
  call log_param(param_file, mod, "p_2_n_assem_Di", phyto(DIAZO)%p_2_n_assem)
  call log_param(param_file, mod, "p_2_n_assem_Lg", phyto(LARGE)%p_2_n_assem)
  call log_param(param_file, mod, "p_2_n_assem_Sm", phyto(SMALL)%p_2_n_assem)
  call log_param(param_file, mod, "p_2_n_photo", topaz%p_2_n_photo)
  call log_param(param_file, mod, "p_2_n_RKR", topaz%p_2_n_RKR)
  call log_param(param_file, mod, "p_2_n_uptake", topaz%p_2_n_uptake)
!  call log_param(param_file, mod, "plast_2_chl", topaz%plast_2_chl)
  call log_param(param_file, mod, "plast_2_chl_Di", phyto(DIAZO)%plast_2_chl)
  call log_param(param_file, mod, "plast_2_chl_Lg", phyto(LARGE)%plast_2_chl)
  call log_param(param_file, mod, "plast_2_chl_Sm", phyto(SMALL)%plast_2_chl)
  call log_param(param_file, mod, "r_nfix_Di", phyto(DIAZO)%r_nfix)
  call log_param(param_file, mod, "r_nfix_Lg", phyto(LARGE)%r_nfix)
  call log_param(param_file, mod, "r_nfix_Sm", phyto(SMALL)%r_nfix)
  call log_param(param_file, mod, "r_other_min_Di", phyto(DIAZO)%r_other_min)
  call log_param(param_file, mod, "r_other_min_Lg", phyto(LARGE)%r_other_min)
  call log_param(param_file, mod, "r_other_min_Sm", phyto(SMALL)%r_other_min)
  call log_param(param_file, mod, "r_uptake_max_Di", phyto(DIAZO)%r_uptake_max)
  call log_param(param_file, mod, "r_uptake_max_Lg", phyto(LARGE)%r_uptake_max)
  call log_param(param_file, mod, "r_uptake_max_Sm", phyto(SMALL)%r_uptake_max)
  call log_param(param_file, mod, "p_2_n_static", topaz%p_2_n_static)
  call log_param(param_file, mod, "p_2_n_static_Di", phyto(DIAZO)%p_2_n_static)
  call log_param(param_file, mod, "p_2_n_static_Lg", phyto(LARGE)%p_2_n_static)
  call log_param(param_file, mod, "p_2_n_static_Sm", phyto(SMALL)%p_2_n_static)
  call log_param(param_file, mod, "si_2_n_max_Lg", phyto(LARGE)%si_2_n_max)
  call log_param(param_file, mod, "si_2_n_static", topaz%si_2_n_static)
  call log_param(param_file, mod, "si_2_n_static_Lg", phyto(LARGE)%si_2_n_static)
  call log_param(param_file, mod, "k_fed_Di", phyto(DIAZO)%k_fed)
  call log_param(param_file, mod, "k_fed_Lg", phyto(LARGE)%k_fed)
  call log_param(param_file, mod, "k_fed_Sm", phyto(SMALL)%k_fed)
  call log_param(param_file, mod, "k_nh4_Lg", phyto(LARGE)%k_nh4)
  call log_param(param_file, mod, "k_nh4_Sm", phyto(SMALL)%k_nh4)
  call log_param(param_file, mod, "k_no3_Lg", phyto(LARGE)%k_no3)
  call log_param(param_file, mod, "k_no3_Sm", phyto(SMALL)%k_no3)
  call log_param(param_file, mod, "k_po4_Di", phyto(DIAZO)%k_po4)
  call log_param(param_file, mod, "k_po4_Lg", phyto(LARGE)%k_po4)
  call log_param(param_file, mod, "k_po4_Sm", phyto(SMALL)%k_po4)
  call log_param(param_file, mod, "k_sio4_Lg", phyto(LARGE)%k_sio4)
  call log_param(param_file, mod, "fe_ballast_assoc", topaz%fe_ballast_assoc)
  call log_param(param_file, mod, "felig_bkg", topaz%felig_bkg)
  call log_param(param_file, mod, "felig_2_don", topaz%felig_2_don)
  call log_param(param_file, mod, "k_fe_2_n_Di", phyto(DIAZO)%k_fe_2_n)
  call log_param(param_file, mod, "k_fe_2_n_Lg", phyto(LARGE)%k_fe_2_n)
  call log_param(param_file, mod, "k_fe_2_n_Sm", phyto(SMALL)%k_fe_2_n)
  call log_param(param_file, mod, "fe_2_n_max_Sm", phyto(SMALL)%fe_2_n_max)
  call log_param(param_file, mod, "fe_2_n_max_Lg", phyto(LARGE)%fe_2_n_max)
  call log_param(param_file, mod, "fe_2_n_max_Di", phyto(DIAZO)%fe_2_n_max)
  call log_param(param_file, mod, "fe_2_n_sed", topaz%fe_2_n_sed)
  call log_param(param_file, mod, "fe_2_n_static", topaz%fe_2_n_static)
  call log_param(param_file, mod, "fe_2_n_static_Di", phyto(DIAZO)%fe_2_n_static)
  call log_param(param_file, mod, "fe_2_n_static_Lg", phyto(LARGE)%fe_2_n_static)
  call log_param(param_file, mod, "fe_2_n_static_Sm", phyto(SMALL)%fe_2_n_static)
  call log_param(param_file, mod, "fe_coast", topaz%fe_coast)
  call log_param(param_file, mod, "kfe_2nd_order", topaz%kfe_2nd_order)
  call log_param(param_file, mod, "kfe_bal", topaz%kfe_bal)
  call log_param(param_file, mod, "kfe_des", topaz%kfe_des)
  call log_param(param_file, mod, "kfe_eq_lig", topaz%kfe_eq_lig)
  call log_param(param_file, mod, "kfe_org", topaz%kfe_org)
  call log_param(param_file, mod, "alpha_Di", phyto(DIAZO)%alpha)
  call log_param(param_file, mod, "alpha_Lg", phyto(LARGE)%alpha)
  call log_param(param_file, mod, "alpha_Sm", phyto(SMALL)%alpha)
  call log_param(param_file, mod, "kappa_eppley", topaz%kappa_eppley)
  call log_param(param_file, mod, "P_C_max_Di", phyto(DIAZO)%P_C_max)
  call log_param(param_file, mod, "P_C_max_Lg", phyto(LARGE)%P_C_max)
  call log_param(param_file, mod, "P_C_max_Sm", phyto(SMALL)%P_C_max)
  call log_param(param_file, mod, "thetamax_Di", phyto(DIAZO)%thetamax)
  call log_param(param_file, mod, "thetamax_Lg", phyto(LARGE)%thetamax)
  call log_param(param_file, mod, "thetamax_Sm", phyto(SMALL)%thetamax)
  call log_param(param_file, mod, "thetamin", topaz%thetamin)
  call log_param(param_file, mod, "zeta", topaz%zeta)
  call log_param(param_file, mod, "k_n_inhib_Di", topaz%k_n_inhib_Di)
  call log_param(param_file, mod, "o2_inhib_Di_pow", topaz%o2_inhib_Di_pow)
  call log_param(param_file, mod, "o2_inhib_Di_sat", topaz%o2_inhib_Di_sat)
  call log_param(param_file, mod, "gamma_irr_mem", topaz%gamma_irr_mem)
  call log_param(param_file, mod, "fdet0_Lg", phyto(LARGE)%fdet0)
  call log_param(param_file, mod, "fdet0_Sm", phyto(SMALL)%fdet0)
  call log_param(param_file, mod, "gamma_nhet", topaz%gamma_nhet)
  call log_param(param_file, mod, "k_diss_sio2", topaz%k_diss_sio2)
  call log_param(param_file, mod, "k_o2", topaz%k_o2)
  call log_param(param_file, mod, "kappa_remin", topaz%kappa_remin)
  call log_param(param_file, mod, "lambda0", topaz%lambda0)
  call log_param(param_file, mod, "o2_min", topaz%o2_min)
  call log_param(param_file, mod, "P_star", topaz%P_star)
  call log_param(param_file, mod, "phyto_min", topaz%phyto_min)
  call log_param(param_file, mod, "phi_nhet", topaz%phi_nhet)
  call log_param(param_file, mod, "q_si_2_n_diss", topaz%rpcaco3)
  call log_param(param_file, mod, "rpcaco3", topaz%rpcaco3)
  call log_param(param_file, mod, "rplith", topaz%rplith )
  call log_param(param_file, mod, "rpsio2", topaz%rpsio2 )
  call log_param(param_file, mod, "wsink", topaz%wsink)
  call log_param(param_file, mod, "gamma_ndet", topaz%gamma_ndet)
  call log_param(param_file, mod, "gamma_cadet", topaz%gamma_cadet)
  call log_param(param_file, mod, "gamma_sidet", topaz%gamma_sidet)
  call log_param(param_file, mod, "ksp_caco3", topaz%ksp_caco3 )
  call log_param(param_file, mod, "k_caco3_pres", topaz%k_caco3_pres )
  call log_param(param_file, mod, "gamma_cased_dis", topaz%gamma_cased_dis)
  call log_param(param_file, mod, "gamma_sdon", topaz%gamma_sdon)
  call log_param(param_file, mod, "gamma_sdop", topaz%gamma_sdop)
  call log_param(param_file, mod, "phi_sdon", topaz%phi_sdon)
  call log_param(param_file, mod, "phi_sdop", topaz%phi_sdop)
  call log_param(param_file, mod, "gamma_ldon", topaz%gamma_ldon)
  call log_param(param_file, mod, "phi_ldon", topaz%phi_ldon)
  call log_param(param_file, mod, "tracer_debug", topaz%tracer_debug)
  call log_param(param_file, mod, "tracer_debug_verbose", topaz%tracer_debug_verbose)
  call log_param(param_file, mod, "gamma_nitrif", topaz%gamma_nitrif)
  call log_param(param_file, mod, "irr_inhibit", topaz%irr_inhibit )
  call log_param(param_file, mod, "k_lith", topaz%k_lith)
  call log_param(param_file, mod, "bio_tau", topaz%bio_tau)
  call log_param(param_file, mod, "NKML", topaz%nkml)
  call log_param(param_file, mod, "RHO_0", topaz%Rho_0)

  module_initialized = .true.
  register_TOPAZ = .true.
end function register_TOPAZ

!#######################################################################
! <SUBROUTINE NAME="initialize_TOPAZ">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
!
!       Note: this subroutine should be merged into ocean_topaz_start
! </DESCRIPTION>
subroutine initialize_TOPAZ(restart, day, G, h, OBC, topaz, &
                                      sponge_CSp, diag_to_Z_CSp)
  logical,                            intent(in) :: restart
  type(time_type), target,            intent(in) :: day
  type(ocean_grid_type),              intent(in) :: G
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in) :: h
  type(ocean_OBC_type),               pointer    :: OBC
  type(TOPAZ_CS),                     pointer    :: topaz
  type(sponge_CS),                    pointer    :: sponge_CSp
  type(diag_to_Z_CS),                 pointer    :: diag_to_Z_CSp


! Arguments: restart - .true. if the fields have already been read from
!                     a restart file.
!  (in)      day - Time of the start of the run.
!  (in)      G - The ocean's grid structure.
!  (in)      h - Layer thickness, in m.
!  (in)      OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in/out)  CS - The control structure returned by a previous call to
!                 GOLD_register_TOPAZ.
!  (in/out)  sponge_CSp - A pointer to the control structure for the sponges, if
!                         they are in use.  Otherwise this may be unassociated.
!  (in/out)  diag_to_Z_Csp - A pointer to the control structure for diagnostics
!                            in depth space.


!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

integer         :: i_bnd_off
integer         :: j_bnd_off
integer         :: i
integer         :: j
integer         :: k
integer         :: n
integer         :: isc, iec, jsc, jec, nk
!
!---------------------------------------------------------------------
! Use shortened naming convention
!---------------------------------------------------------------------
!
  call register_diagnostics(G, day, topaz, diag_to_Z_CSp)


end subroutine initialize_TOPAZ

!#######################################################################
! <SUBROUTINE NAME="TOPAZ_column_physics">
!
! <DESCRIPTION>
!     Calculate the surface boundary conditions
! </DESCRIPTION>

subroutine TOPAZ_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                          G, topaz, tv, optics)
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in) :: h_old, h_new, ea, eb
  type(forcing),                      intent(in) :: fluxes
  real,                               intent(in) :: dt
  type(ocean_grid_type),              intent(in) :: G
  type(TOPAZ_CS),                     pointer    :: topaz
  type(thermo_var_ptrs),              intent(in) :: tv
  type(optics_type),                  intent(in) :: optics
!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! CFCs are relatively simple, as they are passive tracers. with only a surface
! flux as a source.

! Arguments: h_old -  Layer thickness before entrainment, in m or kg m-2.
!  (in)      h_new -  Layer thickness after entrainment, in m or kg m-2.
!  (in)      ea - an array to which the amount of fluid entrained
!                 from the layer above during this call will be
!                 added, in m or kg m-2.
!  (in)      eb - an array to which the amount of fluid entrained
!                 from the layer below during this call will be
!                 added, in m or kg m-2.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.
!  (in)      topaz - The control structure returned by a previous call to
!                 GOLD_register_TOPAZ.
!  (in)      tv - The structure containing thermodynamic variables.
!  (in)      optics - The structure containing optical properties.
!
! The arguments to this subroutine are redundant in that
!     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]

real :: hold0(SZI_(G))        ! The original topmost layer thickness,
                              ! with surface mass fluxes added back, m.
real :: b1(SZI_(G))           ! b1 and c1 are variables used by the
real :: c1(SZI_(G),SZK_(G))   ! tridiagonal solver.

real                                              :: temp_array(SZI_(G),SZJ_(G))
real, pointer, dimension(:,:,:)                   :: tracer, tracer_th
real, pointer, dimension(:,:)                     :: tracer_sfc_flux, &
                                                     tracer_btm_flux, &
                                                     btm_reservoir
real, dimension(G%isd:G%ied, G%jsd:G%jed,G%ke)    :: grid_tmask, &
                                                     rho_dzt, &
                                                     dzt
real, dimension(G%isd:G%ied, G%jsd:G%jed)         :: hblt_depth
integer, dimension(G%isd:G%ied,G%jsd:G%jed)       :: grid_kmt
!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer :: i, j, k, n
integer :: ind_no3,       ind_nh4,        ind_fed,        ind_ldon, &
           ind_dic,       ind_alk,        ind_o2,         ind_lith, &
           ind_irr
integer :: isc, iec, jsc, jec, nk

! =====================================================================
!     begin executable code
! =====================================================================
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nk = G%ke

  grid_tmask(:,:,:) = 0.0
  do j = jsc, jec ; do i = isc, iec  !{
    if (G%hmask(i,j) .gt. 0) &
    grid_tmask(i,j,:) = 1.0
  enddo ; enddo !} i,j

  rho_dzt(:,:,:) = G%H_to_kg_m2 * G%Angstrom
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
    rho_dzt(i,j,k) = G%H_to_kg_m2 * h_old(i,j,k)
  enddo ; enddo ; enddo !}

  dzt(:,:,:) = 1.0
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
    dzt(i,j,k) = G%H_to_m * h_old(i,j,k)
  enddo ; enddo ; enddo !}

! Boussinesq model 
  hblt_depth(:,:) = G%H_to_m * G%Angstrom
  do j=jsc,jec ; do i=isc,iec ; 
    hblt_depth(i,j) = G%H_to_m * h_old(i,j,1)
  enddo ; enddo
  do k=2,topaz%nkml ; do j=jsc,jec ; do i=isc,iec
    hblt_depth(i,j) = hblt_depth(i,j) + G%H_to_m * h_old(i,j,k)
  enddo ; enddo ; enddo

  grid_kmt(:,:) = 0
  do j = jsc, jec ; do i = isc, iec   !{
! Tell the code that a layer thicker than 1m is the bottom layer.
    if (G%hmask(i,j) .gt. 0) then
      grid_kmt(i,j) = nk
    endif
  enddo ; enddo !} i,j

  do j =jsc, jec ; do i = isc, iec  !{
    topaz%mask_coast(i,j) = 0.0
    if (G%hmask(i,j) .gt. 0) then !{
      if (G%hmask(i-1,j) .eq. 0 .or. G%hmask(i,j-1) .eq. 0 .or. &
          G%hmask(i+1,j) .eq. 0 .or. G%hmask(i,j+1) .eq. 0) then !{
        topaz%mask_coast(i,j) = 1.0
      endif !}
    endif !}
  enddo ; enddo !} i,j
!
!---------------------------------------------------------------------
! Use shortened naming convention
!---------------------------------------------------------------------
!

  ind_dic     = topaz%ind_dic      ; ind_o2      = topaz%ind_o2       
  ind_alk     = topaz%ind_alk      ; ind_fed     = topaz%ind_fed
  ind_ldon    = topaz%ind_ldon     ; ind_lith    = topaz%ind_lith
  ind_nh4     = topaz%ind_nh4      ; ind_no3     = topaz%ind_no3
  ind_irr     = topaz%ind_irr

!
!---------------------------------------------------------------------
!     use the surface fluxes from the coupler
!       stf is in mol/m^2/s, flux from coupler is positive upwards
!       deposition is in mol/m^2/s, flux from coupler is positive upwards
!       river flux is in m s-1
!       river runoff tracer concentration (triver) is in mol/kg
!---------------------------------------------------------------------
!
! Alk has river flux and negative contribution from NO3 deposition flux
!
  Tr(ind_alk)%stf(:,:) = 0.0
  if (topaz%ind_runoff_alk .gt. 0 ) then 
  call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_runoff_alk, ind_flux, &
                              Tr(ind_alk)%triver, isc, iec, jsc, jec, +1.0)
    Tr(ind_alk)%stf = Tr(ind_alk)%stf + Tr(ind_alk)%triver * (fluxes%liq_runoff + fluxes%froz_runoff)
  endif
!
! DIC has CO2 has gas exchange and river flux
!
  call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_co2_flux, ind_flux, &
                              Tr(ind_dic)%stf, isc, iec, jsc, jec, -1.0)
  if (topaz%ind_runoff_dic .gt. 0 ) then 
  call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_runoff_dic, ind_flux, &
                              Tr(ind_dic)%triver, isc, iec, jsc, jec, +1.0)
    Tr(ind_dic)%stf = Tr(ind_dic)%stf + Tr(ind_dic)%triver * (fluxes%liq_runoff + fluxes%froz_runoff)
  endif
!
! Fed has deposition and river flux
!
  if (topaz%ind_dry_dep_fed .gt. 0 ) then 
    call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_dry_dep_fed, ind_flux,&
                                Tr(ind_fed)%stf, isc, iec, jsc, jec, -1.0)
  else
    Tr(ind_fed)%stf(:,:) = 0.0
  endif
  if (topaz%ind_wet_dep_fed .gt. 0 ) then 
    call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_wet_dep_fed, ind_flux,&
                                temp_array, isc, iec, jsc, jec, -1.0)
    Tr(ind_fed)%stf = Tr(ind_fed)%stf + temp_array
  endif
  if (topaz%ind_runoff_fed .gt. 0 ) then 
  call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_runoff_fed, ind_flux, &
                            Tr(ind_fed)%triver, isc, iec, jsc, jec, +1.0)
    Tr(ind_fed)%stf = Tr(ind_fed)%stf + Tr(ind_fed)%triver * (fluxes%liq_runoff + fluxes%froz_runoff)
  endif
!
! Lith has deposition and river flux
!
  if (topaz%ind_dry_dep_lith .gt. 0 ) then 
    call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_dry_dep_lith,ind_flux,&
                                Tr(ind_lith)%stf, isc, iec, jsc, jec, -1.0)
  else
    Tr(ind_lith)%stf(:,:) = 0.0
  endif
  if (topaz%ind_wet_dep_lith .gt. 0 ) then 
    call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_wet_dep_lith,ind_flux,&
                                temp_array, isc, iec, jsc, jec, -1.0)
    Tr(ind_lith)%stf = Tr(ind_lith)%stf + temp_array
  endif
  if (topaz%ind_runoff_lith .gt. 0 ) then 
    call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_runoff_lith, ind_flux, &
                                Tr(ind_lith)%triver, isc, iec, jsc, jec, +1.0)
    Tr(ind_lith)%stf = Tr(ind_lith)%stf + Tr(ind_lith)%triver * (fluxes%liq_runoff + fluxes%froz_runoff)
  endif
!
! LDON has river flux only
!
    if (topaz%ind_runoff_ldon .gt. 0 ) then 
    call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_runoff_ldon, ind_flux, &
                            Tr(ind_ldon)%triver, isc, iec, jsc, jec, +1.0)
    Tr(ind_ldon)%stf = Tr(ind_ldon)%stf + Tr(ind_ldon)%triver * (fluxes%liq_runoff + fluxes%froz_runoff)
  endif
!
! NH4 has deposition and river flux
!
  if (topaz%ind_dry_dep_nh4 .gt. 0 ) then 
    call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_dry_dep_nh4, ind_flux,&
                                Tr(ind_nh4)%stf, isc, iec, jsc, jec, -1.0)
  else
    Tr(ind_nh4)%stf(:,:) = 0.0
  endif
  if (topaz%ind_wet_dep_nh4 .gt. 0 ) then 
    call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_wet_dep_nh4, ind_flux,&
                                temp_array, isc, iec, jsc, jec, -1.0)
    Tr(ind_nh4)%stf = Tr(ind_nh4)%stf + temp_array
  endif
  if (topaz%ind_runoff_nh4 .gt. 0 ) then 
    call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_runoff_nh4, ind_flux, &
                                Tr(ind_nh4)%triver, isc, iec, jsc, jec, +1.0)
    Tr(ind_nh4)%stf = Tr(ind_nh4)%stf + Tr(ind_nh4)%triver * (fluxes%liq_runoff + fluxes%froz_runoff)
  endif
!
! NO3 has deposition, river flux, and negative deposition contribution to alkalinity
!
  if (topaz%ind_dry_dep_no3 .gt. 0 ) then 
    call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_dry_dep_no3, ind_flux,&
                                Tr(ind_no3)%stf, isc, iec, jsc, jec, -1.0)
    Tr(ind_alk)%stf = Tr(ind_alk)%stf - Tr(ind_no3)%stf
  else
    Tr(ind_no3)%stf(:,:) = 0.0
  endif
  if (topaz%ind_wet_dep_no3 .gt. 0 ) then 
    call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_wet_dep_no3, ind_flux,&
                                temp_array, isc, iec, jsc, jec, -1.0)
    Tr(ind_no3)%stf = Tr(ind_no3)%stf + temp_array
    Tr(ind_alk)%stf = Tr(ind_alk)%stf - temp_array
  endif
  if (topaz%ind_runoff_no3 .gt. 0 ) then 
    call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_runoff_no3, ind_flux, &
                                Tr(ind_no3)%triver, isc, iec, jsc, jec, +1.0)
    Tr(ind_no3)%stf = Tr(ind_no3)%stf + Tr(ind_no3)%triver * (fluxes%liq_runoff + fluxes%froz_runoff)
  endif
!
! O2 has gas exchange only
!
  call extract_coupler_values(fluxes%tr_fluxes, topaz%ind_o2_flux, ind_flux, &
                              Tr(ind_o2)%stf, isc, iec, jsc, jec, -1.0)
!
! Global integral calc
!
  if(topaz%tracer_debug) &
    call global_tracer_integrals(G, topaz, Tr, h_old, dt, 'Before Source')


  call ocean_topaz_source(G, topaz, G%dxdyh, grid_tmask, grid_kmt, &
                          rho_dzt, dzt, hblt_depth, fluxes, dt, tv, optics)
!
! Global integral calc
!
  if(topaz%tracer_debug) &
    call global_tracer_integrals(G, topaz, Tr, h_old, dt, 'After Source')


  ! Use a tridiagonal solver to determine the concentrations after the
  ! surface source is applied and diapycnal advection and diffusion occurs.
  do n = 1, NUM_PROG_TRACERS
    tracer      => Tr(n)%field
    tracer_sfc_flux => Tr(n)%stf
    tracer_btm_flux => Tr(n)%btf
!
!---------------------------------------------------------------------
! Cadet, Fedet and Lithdet all have sinking fluxes through the bottom
!---------------------------------------------------------------------
!
    if (n==topaz%ind_fedet .or. n==topaz%ind_cadet .or. n==topaz%ind_lithdet) then
      btm_reservoir => Tr(n)%btm_reservoir
      call tracer_vertdiff(h_old, ea, eb, dt, tracer, G, btm_reservoir=btm_reservoir,     &
                           sfc_flux=tracer_sfc_flux,                          &
                           btm_flux=tracer_btm_flux, sink_rate=Tr(n)%sink_dist)
!
!---------------------------------------------------------------------
! Everything else has a no sinking flux through the bottom
!---------------------------------------------------------------------
!
    else
      call tracer_vertdiff(h_old, ea, eb, dt, tracer, G,                                  &
                           sfc_flux=tracer_sfc_flux,&
                           btm_flux=tracer_btm_flux, sink_rate=Tr(n)%sink_dist)
    endif
  enddo

!
! Global integral calc
!
  if(topaz%tracer_debug) &
    call global_tracer_integrals(G, topaz, Tr, h_new, dt, 'After vertdiff')
!
!---------------------------------------------------------------------
! Get bottom flux of Cadet, Fedet and Lithdet and reset bottom flux boundary condition
!---------------------------------------------------------------------
!
  do j =jsc, jec ; do i = isc, iec  !{
    Tr(topaz%ind_fcadet_btm)%field(i,j,1) = Tr(topaz%ind_cadet)%btm_reservoir(i,j) / dt
    topaz%ffedet_btm(i,j) = Tr(topaz%ind_fedet)%btm_reservoir(i,j) / dt
    topaz%flithdet_btm(i,j) = Tr(topaz%ind_lithdet)%btm_reservoir(i,j) / dt
    Tr(topaz%ind_cadet)%btm_reservoir(i,j) = 0.0
    Tr(topaz%ind_fedet)%btm_reservoir(i,j) = 0.0
    Tr(topaz%ind_lithdet)%btm_reservoir(i,j) = 0.0
  enddo ; enddo !} i,j
!
! Global integral calc
!
  if(topaz%tracer_debug) &
    call global_tracer_integrals(G, topaz, Tr, h_new, dt, 'After btm discard')
!
!-----------------------------------------------------------------------
!       Save variables for diagnostics
!-----------------------------------------------------------------------
!
  if (topaz%id_alk .gt. 0) &
    call post_data(topaz%id_alk, Tr(topaz%ind_alk)%field, topaz%diag)
  if (topaz%id_cadet .gt. 0) &
    call post_data(topaz%id_cadet, Tr(topaz%ind_cadet)%field, topaz%diag)
  if (topaz%id_dic .gt. 0) &
    call post_data(topaz%id_dic, Tr(topaz%ind_dic)%field, topaz%diag)
  if (topaz%id_fed .gt. 0) &
    call post_data(topaz%id_fed, Tr(topaz%ind_fed)%field, topaz%diag)
  if (topaz%id_fedet .gt. 0) &
    call post_data(topaz%id_fedet, Tr(topaz%ind_fedet)%field, topaz%diag)
  if (topaz%id_fedi .gt. 0) &
    call post_data(topaz%id_fedi, Tr(topaz%ind_fedi)%field, topaz%diag)
  if (topaz%id_felg .gt. 0) &
    call post_data(topaz%id_felg, Tr(topaz%ind_felg)%field, topaz%diag)
  if (topaz%id_fesm .gt. 0) &
    call post_data(topaz%id_fesm, Tr(topaz%ind_fesm)%field, topaz%diag)
  if (topaz%id_ldon .gt. 0) &
    call post_data(topaz%id_ldon, Tr(topaz%ind_ldon)%field, topaz%diag)
  if (topaz%id_lith .gt. 0) &
    call post_data(topaz%id_lith, Tr(topaz%ind_lith)%field, topaz%diag)
  if (topaz%id_lithdet .gt. 0) &
    call post_data(topaz%id_lithdet, Tr(topaz%ind_lithdet)%field, topaz%diag)
  if (topaz%id_ndet .gt. 0) &
    call post_data(topaz%id_ndet, Tr(topaz%ind_ndet)%field, topaz%diag)
  if (topaz%id_ndi .gt. 0) &
    call post_data(topaz%id_ndi, Tr(topaz%ind_ndi)%field, topaz%diag)
  if (topaz%id_nh4 .gt. 0) &
    call post_data(topaz%id_nh4, Tr(topaz%ind_nh4)%field, topaz%diag)
  if (topaz%id_nhet .gt. 0) &
    call post_data(topaz%id_nhet, Tr(topaz%ind_nhet)%field, topaz%diag)
  if (topaz%id_nlg .gt. 0) &
    call post_data(topaz%id_nlg, Tr(topaz%ind_nlg)%field, topaz%diag)
  if (topaz%id_no3 .gt. 0) &
    call post_data(topaz%id_no3, Tr(topaz%ind_no3)%field, topaz%diag)
  if (topaz%id_nsm .gt. 0) &
    call post_data(topaz%id_nsm, Tr(topaz%ind_nsm)%field, topaz%diag)
  if (topaz%id_o2 .gt. 0) &
    call post_data(topaz%id_o2, Tr(topaz%ind_o2)%field, topaz%diag)
  if (topaz%id_pdet .gt. 0) &
    call post_data(topaz%id_pdet, Tr(topaz%ind_pdet)%field, topaz%diag)
  if (topaz%id_pdi .gt. 0) &
    call post_data(topaz%id_pdi, Tr(topaz%ind_pdi)%field, topaz%diag)
  if (topaz%id_plg .gt. 0) &
    call post_data(topaz%id_plg, Tr(topaz%ind_plg)%field, topaz%diag)
  if (topaz%id_psm .gt. 0) &
    call post_data(topaz%id_psm, Tr(topaz%ind_psm)%field, topaz%diag)
  if (topaz%id_po4 .gt. 0) &
    call post_data(topaz%id_po4, Tr(topaz%ind_po4)%field, topaz%diag)
  if (topaz%id_sdon .gt. 0) &
    call post_data(topaz%id_sdon, Tr(topaz%ind_sdon)%field, topaz%diag)
  if (topaz%id_sdop .gt. 0) &
    call post_data(topaz%id_sdop, Tr(topaz%ind_sdop)%field, topaz%diag)
  if (topaz%id_sidet .gt. 0) &
    call post_data(topaz%id_sidet, Tr(topaz%ind_sidet)%field, topaz%diag)
  if (topaz%id_silg .gt. 0) &
    call post_data(topaz%id_silg, Tr(topaz%ind_silg)%field, topaz%diag)
  if (topaz%id_sio4 .gt. 0) &
    call post_data(topaz%id_sio4, Tr(topaz%ind_sio4)%field, topaz%diag)
  if (topaz%id_cased .gt. 0) &
    call post_data(topaz%id_cased, Tr(topaz%ind_cased)%field, topaz%diag)
  if (topaz%id_chl .gt. 0) &
    call post_data(topaz%id_chl, Tr(topaz%ind_chl)%field, topaz%diag)
  if (topaz%id_fcadet_btm .gt. 0) &
    call post_data(topaz%id_fcadet_btm, Tr(topaz%ind_fcadet_btm)%field, topaz%diag)
  if (topaz%id_irr .gt. 0) &
    call post_data(topaz%id_irr, Tr(topaz%ind_irr)%field, topaz%diag)
  if (topaz%id_irr_mem .gt. 0) &
    call post_data(topaz%id_irr_mem, Tr(topaz%ind_irr_mem)%field, topaz%diag)
  if (topaz%id_htotal .gt. 0)               &
    call post_data(topaz%id_htotal, Tr(topaz%ind_htotal)%field,  topaz%diag)
  if (topaz%id_co3_ion .gt. 0)              &
    call post_data(topaz%id_co3_ion, Tr(topaz%ind_co3_ion)%field, topaz%diag)
  if (topaz%id_alpha .gt. 0)              &
    call post_data(topaz%id_alpha, Tr(topaz%ind_alpha)%field, topaz%diag)
  if (topaz%id_csurf .gt. 0)              &
    call post_data(topaz%id_csurf, Tr(topaz%ind_csurf)%field, topaz%diag)


  if (topaz%id_runoff_alk .gt. 0) &
    call post_data(topaz%id_runoff_alk, Tr(ind_alk)%triver, topaz%diag)
  if (topaz%id_runoff_flux_alk .gt. 0) &
    call post_data(topaz%id_runoff_flux_alk, Tr(ind_alk)%triver *            &
              (fluxes%liq_runoff(:,:) + fluxes%froz_runoff(:,:)), topaz%diag)

  if (topaz%id_runoff_dic .gt. 0) &
    call post_data(topaz%id_runoff_dic, Tr(ind_dic)%triver, topaz%diag)
  if (topaz%id_runoff_flux_dic .gt. 0) &
    call post_data(topaz%id_runoff_flux_dic, Tr(ind_dic)%triver *            &
              (fluxes%liq_runoff(:,:) + fluxes%froz_runoff(:,:)), topaz%diag)

  if (topaz%id_runoff_fed .gt. 0) &
    call post_data(topaz%id_runoff_fed, Tr(ind_fed)%triver, topaz%diag)
  if (topaz%id_runoff_flux_fed .gt. 0) &
    call post_data(topaz%id_runoff_flux_fed, Tr(ind_fed)%triver *            &
              (fluxes%liq_runoff(:,:) + fluxes%froz_runoff(:,:)), topaz%diag)

  if (topaz%id_runoff_ldon .gt. 0) &
    call post_data(topaz%id_runoff_ldon, Tr(ind_ldon)%triver, topaz%diag)
  if (topaz%id_runoff_flux_ldon .gt. 0) &
    call post_data(topaz%id_runoff_flux_ldon, Tr(ind_ldon)%triver *          &
              (fluxes%liq_runoff(:,:) + fluxes%froz_runoff(:,:)), topaz%diag)

  if (topaz%id_runoff_lith .gt. 0) &
    call post_data(topaz%id_runoff_lith, Tr(ind_lith)%triver, topaz%diag)
  if (topaz%id_runoff_flux_lith .gt. 0) &
    call post_data(topaz%id_runoff_flux_lith, Tr(ind_lith)%triver *          &
              (fluxes%liq_runoff(:,:) + fluxes%froz_runoff(:,:)), topaz%diag)

  if (topaz%id_runoff_nh4 .gt. 0) &
    call post_data(topaz%id_runoff_nh4, Tr(ind_nh4)%triver, topaz%diag)
  if (topaz%id_runoff_flux_nh4 .gt. 0) &
    call post_data(topaz%id_runoff_flux_nh4, Tr(ind_nh4)%triver *            &
              (fluxes%liq_runoff(:,:) + fluxes%froz_runoff(:,:)), topaz%diag)

  if (topaz%id_runoff_no3 .gt. 0) &
    call post_data(topaz%id_runoff_no3, Tr(ind_no3)%triver, topaz%diag)
  if (topaz%id_runoff_flux_no3 .gt. 0) &
    call post_data(topaz%id_runoff_flux_no3, Tr(ind_no3)%triver *            &
              (fluxes%liq_runoff(:,:) + fluxes%froz_runoff(:,:)), topaz%diag)

  if (topaz%id_sfc_flux_co2 .gt. 0) &
    call post_data(topaz%id_sfc_flux_co2, Tr(ind_dic)%stf, topaz%diag)

  if (topaz%id_sfc_flux_o2 .gt. 0) &
    call post_data(topaz%id_sfc_flux_o2, Tr(ind_O2)%stf, topaz%diag)

  if (topaz%id_sfc_flux_fed .gt. 0) &
    call post_data(topaz%id_sfc_flux_fed, Tr(ind_fed)%stf, topaz%diag)

  if (topaz%id_sfc_flux_lith .gt. 0) &
    call post_data(topaz%id_sfc_flux_lith, Tr(ind_lith)%stf, topaz%diag)

  if (topaz%id_sfc_flux_nh4 .gt. 0) &
    call post_data(topaz%id_sfc_flux_nh4, Tr(ind_nh4)%stf, topaz%diag)

  if (topaz%id_sfc_flux_no3 .gt. 0) &
    call post_data(topaz%id_sfc_flux_no3, Tr(ind_no3)%stf, topaz%diag)

  if (topaz%id_btm_flux_alk .gt. 0) &
    call post_data(topaz%id_btm_flux_alk, Tr(ind_alk)%btf, topaz%diag)

  if (topaz%id_btm_flux_dic .gt. 0) &
    call post_data(topaz%id_btm_flux_dic, Tr(ind_dic)%btf, topaz%diag)

  if (topaz%id_btm_flux_ndet .gt. 0) &
    call post_data(topaz%id_btm_flux_ndet, Tr(topaz%ind_ndet)%btf, topaz%diag)

  if (topaz%id_btm_flux_no3 .gt. 0) &
    call post_data(topaz%id_btm_flux_no3, Tr(ind_no3)%btf, topaz%diag)

  if (topaz%id_btm_flux_pdet .gt. 0) &
    call post_data(topaz%id_btm_flux_pdet, Tr(topaz%ind_pdet)%btf, topaz%diag)

  if (topaz%id_btm_flux_po4 .gt. 0) &
    call post_data(topaz%id_btm_flux_po4, Tr(topaz%ind_po4)%btf, topaz%diag)

  if (topaz%id_ffedet_btm .gt. 0)           &
      call post_data(topaz%id_ffedet_btm, topaz%ffedet_btm, topaz%diag)

  if (topaz%id_flithdet_btm .gt. 0)         &
      call post_data(topaz%id_flithdet_btm, topaz%flithdet_btm, topaz%diag)

end subroutine TOPAZ_column_physics


!#######################################################################
! <SUBROUTINE NAME="TOPAZ_surface_state">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>
subroutine TOPAZ_surface_state(state, h, G, topaz)
  type(surface),                      intent(inout) :: state
  real, dimension(NXMEM_,NYMEM_,NZ_), intent(in)    :: h
  type(ocean_grid_type),              intent(in)    :: G
  type(TOPAZ_CS),                     pointer       :: topaz
!   This subroutine sets up the fields that the coupler needs to calculate the
! CFC fluxes between the ocean and atmosphere.
! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      h - Layer thickness, in m.
!  (in)      G - The ocean's grid structure.
!  (in)      topaz - The control structure returned by a previous call to
!                 GOLD_register_TOPAZ.
!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

integer         :: i
integer         :: j
integer         :: k
integer         :: n
integer         :: ind, ind_dic, ind_alk, ind_po4, ind_o2, ind_sio4, ind_htotal
integer         :: isc, iec, jsc, jec, nk
real, dimension(G%isd:G%ied, G%jsd:G%jed)         :: alpha, Csurf

real            :: epsln=1.0e-30

!
!---------------------------------------------------------------------
! Use shortened naming convention
!---------------------------------------------------------------------
!
  ind_dic     = topaz%ind_dic   ; ind_alk     = topaz%ind_alk
  ind_po4     = topaz%ind_po4   ; ind_o2      = topaz%ind_o2
  ind_sio4    = topaz%ind_sio4  ; ind_htotal  = topaz%ind_htotal
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nk = G%ke

if (topaz%init ) then
  do j = jsc, jec ; do i = isc, iec  !{
    htotallo(i,j,1) = htotal_scale_lo * Tr(ind_htotal)%field(i,j,1)
    htotalhi(i,j,1) = htotal_scale_hi * Tr(ind_htotal)%field(i,j,1)
  enddo ; enddo ; !} i, j

  call GOLD_ocmip2_co2calc(CO2_dope_vec,  &
       G%hmask(:,:),                &
       state%SST(:,:),              &
       state%SSS(:,:),              &
       Tr(ind_dic)%field(:,:,1),    &
       Tr(ind_po4)%field(:,:,1),    &
       Tr(ind_sio4)%field(:,:,1),   &
       Tr(ind_alk)%field(:,:,1),    &
       htotallo(:,:,1),             &
       htotalhi(:,:,1),             &
       Tr(ind_htotal)%field(:,:,1), &
       co2star=Tr(topaz%ind_csurf)%field(:,:,1),  &
       alpha=Tr(topaz%ind_alpha)%field(:,:,1),      &
       pco2surf=topaz%pco2surf(:,:),&
       co3_ion=Tr(topaz%ind_co3_ion)%field(:,:,1))
       
       topaz%init = .false.
endif
!
!---------------------------------------------------------------------
!  Compute the Schmidt number of CO2 in seawater using the 
!  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
!  7373-7382).
!---------------------------------------------------------------------
!
    do j = jsc, jec ; do i = isc, iec  !{
      sc_co2(i,j) =                                                &
!Old Wanninkhof numbers
!       2073.1    + tv%T(i,j,1) *              &! or state%SST(i,j)
!      (-125.62   + tv%T(i,j,1) *            &! or state%SST(i,j)
!         (3.6276 + tv%T(i,j,1) *            &! or state%SST(i,j)
!        (-0.043219))) * G%hmask(i,j)
!New Wanninkhof numbers
       2068.9    + state%SST(i,j) *            &
      (-118.63   + state%SST(i,j) *            &
         (2.9311 + state%SST(i,j) *            &
        (-0.027))) * G%hmask(i,j)
      sc_no_term(i,j) = sqrt(660.0 / (sc_co2(i,j) + epsln))
!
! sc_no_term is applied to both variables here rather than in the coupler 
! in order to limit the number of variables sent to the coupler.
!
! Units of alpha and csurf are converted from mol/kg in GOLD_ocmip2_co2calc to
! coupler units of mol m-3
!
      alpha(i,j) = Tr(topaz%ind_alpha)%field(i,j,1) * topaz%Rho_0 * sc_no_term(i,j)
      Csurf(i,j) = Tr(topaz%ind_csurf)%field(i,j,1) * topaz%Rho_0 * sc_no_term(i,j)
    enddo ; enddo  !}
  call set_coupler_values(alpha, state%tr_fields, topaz%ind_co2_flux, &
                          ind_alpha, isc, iec, jsc, jec)
  call set_coupler_values(Csurf, state%tr_fields, topaz%ind_co2_flux, &
                          ind_csurf, isc, iec, jsc, jec)

!
!---------------------------------------------------------------------
!  Compute the oxygen saturation concentration at 1 atm total
!  pressure in mol/kg given the temperature (t, in deg C) and
!  the salinity (s, in permil)
!
!  From Garcia and Gordon (1992), Limnology and Oceonography.
!  The formula used is from page 1310, eq (8).
!
!  *** Note: the "a3*ts^2" term (in the paper) is incorrect. ***
!  *** It shouldn't be there.                                ***
!
!  o2_saturation is defined between T(freezing) <= T <= 40 deg C and
!                                   0 permil <= S <= 42 permil
!
! check value: T = 10 deg C, S = 35 permil,
!              o2_saturation = 0.282015 mol m-3
!---------------------------------------------------------------------
!
    do j = jsc, jec ; do i = isc, iec  !{
      tt(i) = 298.15 - state%SST(i,j)
      tk(i) = 273.15 + state%SST(i,j)
      ts(i) = log(tt(i) / tk(i))
      ts2(i) = ts(i) * ts(i)
      ts3(i) = ts2(i) * ts(i)
      ts4(i) = ts3(i) * ts(i)
      ts5(i) = ts4(i) * ts(i)
      o2_saturation(i,j) =                                     &
        exp(a_0 + a_1*ts(i) + a_2*ts2(i) +                     &
            a_3*ts3(i) + a_4*ts4(i) + a_5*ts5(i) +             &
            state%SSS(i,j) *  &
            (b_0 + b_1*ts(i) + b_2*ts2(i) + b_3*ts3(i) +       &
             c_0*state%SSS(i,j)))
!
!       convert from ml/l to mol m-3
!
        o2_saturation(i,j) = o2_saturation(i,j) * (1000.0/22391.6)
!
!---------------------------------------------------------------------
!  Compute the Schmidt number of O2 in seawater using the 
!  formulation proposed by Keeling et al. (1998, Global Biogeochem.
!  Cycles, 12, 141-163).
!---------------------------------------------------------------------
!
      sc_o2(i,j) =                              &
!Old Keeling et al numbers
!        1638.0   + state%SST(i,j) *             &
!        (-81.83  + state%SST(i,j) *             &
!          (1.483 + state%SST(i,j) *             &
!         (-0.008004))) * G%hmask(i,j)
!New Wanninkhof numbers
        1929.7   + state%SST(i,j) *             &
       (-117.46  + state%SST(i,j) *             &
          (3.116 + state%SST(i,j) *             &
         (-0.0306))) * G%hmask(i,j)
      sc_no_term(i,j) = sqrt(660.0 / (sc_o2(i,j) + epsln))
        alpha(i,j) = o2_saturation(i,j) * sc_no_term(i,j)
        Csurf(i,j) = Tr(ind_o2)%field(i,j,1) * topaz%Rho_0 * &
        sc_no_term(i,j)
    enddo ; enddo  !} j
  call set_coupler_values(alpha, state%tr_fields, topaz%ind_o2_flux, &
                          ind_alpha, isc, iec, jsc, jec)
  call set_coupler_values(Csurf, state%tr_fields, topaz%ind_o2_flux, &
                          ind_csurf, isc, iec, jsc, jec)


end subroutine TOPAZ_surface_state

!#######################################################################
! <SUBROUTINE NAME="TOPAZ_end">
!
! <DESCRIPTION>
!     Write out to restart file for various topaz quantities for this run.
! </DESCRIPTION>
subroutine TOPAZ_end(topaz)
type(TOPAZ_CS),                      intent(in)   :: topaz

end subroutine TOPAZ_end



subroutine register_diagnostics(G, model_time, topaz, CS_diag_Z)
type(ocean_grid_type), intent(in)    :: G
type(time_type),       intent(in)    :: model_time
type(TOPAZ_CS),        intent(inout) :: topaz
type(diag_to_Z_CS),    pointer       :: CS_diag_Z

real, dimension(:,:,:), pointer :: tr_ptr
type(vardesc)   :: vardesc_temp
!       register the fields
!
!   The following vardesc types contain a package of metadata about each tracer,
! including, in order, the following elements: name; longname; horizontal
! staggering ('h') for collocation with thickness points ; vertical staggering
! ('L') for a layer variable ; temporal staggering ('s' for snapshot) ; units ;
! and precision in non-restart output files ('f' for 32-bit float or 'd' for
! 64-bit doubles). For most tracers, only the name, longname and units should
! be changed.  See GOLD_variables for the full type description.
! ocean_register_diag is a wrapper for register_diag_field.

  vardesc_temp = vardesc("alk","Alkalinity",'h','L','s',"mol-eq kg-1",'f')
  tr_ptr => Tr(topaz%ind_alk)%field
  topaz%id_alk = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("cadet","Particulate Detrital CaCO3", 'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_cadet)%field
  topaz%id_cadet = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("dic","Dissolved Inorganic Carbon",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_dic)%field
  topaz%id_dic = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("fed","Dissolved Iron",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_fed)%field
  topaz%id_fed = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("fedet","Particulate Detrital Iron",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_nhet)%field
  topaz%id_fedet = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("fedi","Diazotroph Iron",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_fedi)%field
  topaz%id_fedi = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("felg","Large Iron",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_felg)%field
  topaz%id_felg = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("fesm","Small Iron",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_fesm)%field
  topaz%id_fesm = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("ldon","Labile Dissolved Organic Nitrogen",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_ldon)%field
  topaz%id_ldon = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("lith","Lithogenic Mineral",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_lith)%field
  topaz%id_lith = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("lithdet","Detrital Lithogenic Mineral",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_lithdet)%field
  topaz%id_lithdet = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  tr_ptr => Tr(topaz%ind_nhet)%field
  vardesc_temp = vardesc("ndet","Particulate Detrital Nitrogen",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_ndet)%field
  topaz%id_ndet = ocean_register_diag_with_z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("ndi","Diazotroph Nitrogen",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_ndi)%field
  topaz%id_ndi = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("nh4","Ammonia",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_nh4)%field
  topaz%id_nh4 = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("nhet","heterotrophic Nitrogen",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_nhet)%field
  topaz%id_nhet = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("nlg","Large Nitrogen",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_nlg)%field
  topaz%id_nlg = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("no3","Nitrate",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_no3)%field
  topaz%id_no3 = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("nsm","Small Nitrogen",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_nsm)%field
  topaz%id_nsm = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("o2","Oxygen",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_o2)%field
  topaz%id_o2 = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("pdet","Particulate Detrital Phosphorus",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_pdet)%field
  topaz%id_pdet = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("pdi","Diaz Phosphorus",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_pdi)%field
  topaz%id_pdi = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("plg","Large Phosphorus",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_plg)%field
  topaz%id_plg = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("psm","Small Phosphorus",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_psm)%field
  topaz%id_psm = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("po4","Phosphate",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_po4)%field
  topaz%id_po4 = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("sdon","Semilabile Dissolved Organic Nitrogen",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_sdon)%field
  topaz%id_sdon = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("sdop","Semilabile Dissolved Organic Phosphorus",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_sdop)%field
  topaz%id_sdop = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("sidet","Particulate Detrital Silicon",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_sidet)%field
  topaz%id_sidet = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("silg","Large Silicon",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_silg)%field
  topaz%id_silg = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("sio4","Silicic Acid",'h','L','s',"mol kg-1",'f')
  tr_ptr => Tr(topaz%ind_sio4)%field
  topaz%id_sio4 = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("cased","Sediment CaCO3",'h','1','s',"mol m-2",'f')
  topaz%id_cased = ocean_register_diag(vardesc_temp, G, model_time)
  vardesc_temp = vardesc("chl","Chlorophyll",'h','L','s',"ug kg-1",'f')
  tr_ptr => Tr(topaz%ind_chl)%field
  topaz%id_chl = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("fcadet_btm","CaCO3 flux at bottom",'h','1','s','mol m-2 s-1','f')
  topaz%id_fcadet_btm = ocean_register_diag(vardesc_temp, G, model_time)
  vardesc_temp = vardesc("irr","Irradiance",'h','L','s',"W m-2",'f')
  tr_ptr => Tr(topaz%ind_irr)%field
  topaz%id_irr = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)
  vardesc_temp = vardesc("irr_mem","Irradiance Memory",'h','L','s',"W m-2",'f')
  tr_ptr => Tr(topaz%ind_irr_mem)%field
  topaz%id_irr_mem = ocean_register_diag_with_Z(tr_ptr, vardesc_temp, G, model_time, CS_diag_Z)

  vardesc_temp = vardesc("sc_co2","Schmidt number - CO2",'h','1','s',' ','f')
  id_sc_co2 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("sc_o2","Schmidt number - O2",'h','1','s',' ','f')
  id_sc_o2 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("o2_saturation","O2 saturation",'h','1','s','mol kg-1','f')
  id_o2_sat = ocean_register_diag(vardesc_temp, G, model_time)

! Register the diagnostics for the various phytoplankton 
!
! Register Limitation Diagnostics
!
  vardesc_temp = vardesc("def_fe_Di","Diaz. Phyto. Fe Deficiency",'h','L','s',' ','f')
  phyto(DIAZO)%id_def_fe = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("def_fe_Lg","Large Phyto. Fe Deficiency",'h','L','s',' ','f')
  phyto(LARGE)%id_def_fe = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("def_fe_Sm","Small Phyto. Fe Deficiency",'h','L','s',' ','f')
  phyto(SMALL)%id_def_fe = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("felim_Di","Diaz. Phyto. Fed Limitation",'h','L','s',' ','f')
  phyto(DIAZO)%id_felim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("felim_Lg","Large Phyto. Fed Limitation",'h','L','s',' ','f')
  phyto(LARGE)%id_felim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("felim_Sm","Small Phyto. Fed Limitation",'h','L','s',' ','f')
  phyto(SMALL)%id_felim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("irrlim_Di","Diaz. Phyto. Light Limitation",'h','L','s',' ','f')
  phyto(DIAZO)%id_irrlim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("irrlim_Lg","Large Phyto. Light Limitation",'h','L','s',' ','f')
  phyto(LARGE)%id_irrlim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("irrlim_Sm","Small Phyto. Light Limitation",'h','L','s',' ','f')
  phyto(SMALL)%id_irrlim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("liebig_lim_Di","Diaz. Phyto. Overall Nutrient Limitation",'h','L','s',' ','f')
  phyto(DIAZO)%id_liebig_lim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("liebig_lim_Lg","Large Phyto. Overall Nutrient Limitation",'h','L','s',' ','f')
  phyto(LARGE)%id_liebig_lim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("liebig_lim_Sm","Small Phyto. Overall Nutrient Limitation",'h','L','s',' ','f')
  phyto(SMALL)%id_liebig_lim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("mu_Di","Diaz. Phyto. Overall Growth Rate",'h','L','s','s-1','f')
  phyto(DIAZO)%id_mu = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("mu_Lg","Large Phyto. Overall Growth Rate",'h','L','s','s-1','f')
  phyto(LARGE)%id_mu = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("mu_Sm","Small Phyto. Overall Growth Rate",'h','L','s','s-1','f')
  phyto(SMALL)%id_mu = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("nh4lim_Lg","Ammonia Limitation of Large Phyto",'h','L','s',' ','f')
  phyto(LARGE)%id_nh4lim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("nh4lim_Sm","Ammonia Limitation of Small Phyto",'h','L','s',' ','f')
  phyto(SMALL)%id_nh4lim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("no3lim_Lg","Nitrate Limitation of Large Phyto",'h','L','s',' ','f')
  phyto(LARGE)%id_no3lim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("no3lim_Sm","Nitrate Limitation of Small Phyto",'h','L','s',' ','f')
  phyto(SMALL)%id_no3lim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("po4lim_Di","Phosphate Limitation of Diaz. Phyto",'h','L','s',' ','f')
  phyto(DIAZO)%id_po4lim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("po4lim_Lg","Phosphate Limitation of Large Phyto",'h','L','s',' ','f')
  phyto(LARGE)%id_po4lim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("po4lim_Sm","Phosphate Limitation of Small Phyto",'h','L',' ',' ','f')
  phyto(SMALL)%id_po4lim = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("q_fe_2_n_Di","Fe:N ratio of Diaz. Phyto",'h','L','s','mol/mol','f')
  phyto(DIAZO)%id_q_fe_2_n = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("q_fe_2_n_Lg","Fe:N ratio of Large Phyto",'h','L','s','mol/mol','f')
  phyto(LARGE)%id_q_fe_2_n = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("q_fe_2_n_Sm","Fe:N ratio of Small Phyto",'h','L','s','mol/mol','f')
  phyto(SMALL)%id_q_fe_2_n = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("q_p_2_n_Di","P:N ratio of Diaz. Phyto",'h','L','s','mol/mol','f')
  phyto(DIAZO)%id_q_p_2_n = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("q_p_2_n_Lg","P:N ratio of Large Phyto",'h','L','s','mol/mol','f')
  phyto(LARGE)%id_q_p_2_n = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("q_p_2_n_Sm","P:N ratio of Small Phyto",'h','L','s','mol/mol','f')
  phyto(SMALL)%id_q_p_2_n = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("q_p_2_n_opt_Di","Optimal P:N ratio of Diaz. Phyto",'h','L','s','mol/mol','f')
  phyto(DIAZO)%id_q_p_2_n_opt = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("q_p_2_n_opt_Lg","Optimal P:N ratio of Large Phyto",'h','L','s','mol/mol','f')
  phyto(LARGE)%id_q_p_2_n_opt = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("q_p_2_n_opt_Sm","Optimal P:N ratio of Small Phyto",'h','L','s','mol/mol','f')
  phyto(SMALL)%id_q_p_2_n_opt = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("silim_Lg","SiO4 Limitation of Large Phyto",'h','L',' ',' ','f')
  phyto(LARGE)%id_silim = ocean_register_diag(vardesc_temp, G, model_time)
!
! Register Grazing Diagnostics
!
  vardesc_temp = vardesc("jgraz_fe_Di","Diazotroph iron grazing layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(DIAZO)%id_jgraz_fe = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jgraz_fe_Lg","Large iron grazing layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(LARGE)%id_jgraz_fe = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jgraz_fe_Sm","Small iron grazing layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(SMALL)%id_jgraz_fe = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jgraz_n_Di","Diazotroph nitrogen grazing layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(DIAZO)%id_jgraz_n = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jgraz_n_Lg","Large nitrogen grazing layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(LARGE)%id_jgraz_n = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jgraz_n_Sm","Small nitrogen grazing layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(SMALL)%id_jgraz_n = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jgraz_sio2_Lg","Large silicon grazing layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(LARGE)%id_jgraz_sio2 = ocean_register_diag(vardesc_temp, G, model_time)
!
! Register Production Diagnostics
!
  vardesc_temp = vardesc("jprod_n2_Di","Nitrogen fixation layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(DIAZO)%id_jprod_n2 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_fe_Di","Diaz. phyto. Fed production layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(DIAZO)%id_jprod_fe = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_fe_Lg","Large phyto. Fed production layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(LARGE)%id_jprod_fe = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_fe_Sm","Small phyto. Fed production layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(SMALL)%id_jprod_fe = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_nh4_Di","Diaz. phyto. NH4 production layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(DIAZO)%id_jprod_nh4 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_nh4_Lg","Large phyto. NH4 production layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(LARGE)%id_jprod_nh4 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_nh4_Sm","Small phyto. NH4 production layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(SMALL)%id_jprod_nh4 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_no3_Di","Diaz. phyto. NO3 production layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(DIAZO)%id_jprod_no3 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_no3_Lg","Large phyto. NO3 production layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(LARGE)%id_jprod_no3 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_no3_Sm","Small phyto. NO3 Production layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(SMALL)%id_jprod_no3 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_po4_Di","Diaz. phyto. PO4 production layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(DIAZO)%id_jprod_po4 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_po4_Lg","Large phyto. PO4 production layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(LARGE)%id_jprod_po4 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_po4_Sm","Small phyto. PO4 production layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(SMALL)%id_jprod_po4 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_sio4_Lg","Large phyto. SiO4 production layer integral",'h','L','s','mol m-2 s-1','f')
  phyto(LARGE)%id_jprod_sio4 = ocean_register_diag(vardesc_temp, G, model_time)
!
! Register non Phytoplankton Diagnostics
!
  vardesc_temp = vardesc("runoff_alk","Runoff - Alkalinity",'h','1','s','mol-eq kg-1','f')
  topaz%id_runoff_alk = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("runoff_flux_alk","Runoff flux - Alkalinity",'h','1','s','mol-eq m-2 s-1','f')
  topaz%id_runoff_flux_alk = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("runoff_dic","Runoff - Dissolved Inorganic Carbon",'h','1','s','mol kg-1','f')
  topaz%id_runoff_dic = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("runoff_flux_dic","Runoff flux - Dissolved Inorganic Carbon",'h','1','s','mol m-2 s-1','f')
  topaz%id_runoff_flux_dic = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("runoff_fed","Runoff - Dissolved Fe",'h','1','s','mol kg-1','f')
  topaz%id_runoff_fed = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("runoff_flux_fed","Runoff flux - Dissolved Fe",'h','1','s','mol m-2 s-1','f')
  topaz%id_runoff_flux_fed = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("runoff_ldon","Runoff - LDON",'h','1','s','mol kg-1','f')
  topaz%id_runoff_ldon = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("runoff_flux_ldon","Runoff flux - LDON",'h','1','s','mol m-2 s-1','f')
  topaz%id_runoff_flux_ldon = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("runoff_lith","Runoff - Lith",'h','1','s','g kg-1','f')
  topaz%id_runoff_lith = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("runoff_flux_lith","Runoff flux - Lith",'h','1','s','g m-2 s-1','f')
  topaz%id_runoff_flux_lith = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("runoff_nh4","Runoff - NH4",'h','1','s','mol kg-1','f')
  topaz%id_runoff_nh4 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("runoff_flux_nh4","Runoff flux - NH4",'h','1','s','mol m-2 s-1','f')
  topaz%id_runoff_flux_nh4 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("runoff_no3","Runoff - NO3",'h','1','s','mol kg-1','f')
  topaz%id_runoff_no3 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("runoff_flux_no3","Runoff flux - NO3",'h','1','s','mol m-2 s-1','f')
  topaz%id_runoff_flux_no3 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("sfc_flux_co2","Surface Flux - CO2",'h','1','s','mol m-2 s-1','f')
  topaz%id_sfc_flux_co2 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("sfc_flux_o2","Surface Flux - O2",'h','1','s','mol m-2 s-1','f')
  topaz%id_sfc_flux_o2 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("sfc_flux_fed","Surface Flux - Fed",'h','1','s','mol m-2 s-1','f')
  topaz%id_sfc_flux_fed = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("sfc_flux_lith","Surface Flux - Lith",'h','1','s','g m-2 s-1','f')
  topaz%id_sfc_flux_lith = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("sfc_flux_no3","Surface Flux - NO3",'h','1','s','mol m-2 s-1','f')
  topaz%id_sfc_flux_no3 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("sfc_flux_nh4","Surface Flux - NH4",'h','1','s','mol m-2 s-1','f')
  topaz%id_sfc_flux_nh4 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("btm_flux_alk","Bottom Flux - Alk",'h','1','s','mol m-2 s-1','f')
  topaz%id_btm_flux_alk = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("btm_flux_ndet","Bottom Flux - Ndet",'h','1','s','mol m-2 s-1','f')
  topaz%id_btm_flux_ndet = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("btm_flux_no3","Bottom Flux - NO3",'h','1','s','mol m-2 s-1','f')
  topaz%id_btm_flux_no3 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("btm_flux_pdet","Bottom Flux - Pdet",'h','1','s','mol m-2 s-1','f')
  topaz%id_btm_flux_pdet = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("btm_flux_po4","Bottom Flux - PO4",'h','1','s','mol m-2 s-1','f')
  topaz%id_btm_flux_po4 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("csurf","CO2* water",'h','L','s','mol kg-1','f')
  topaz%id_csurf = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("co3_ion","Carbonate Ion",'h','L','s','mol kg-1','f')
  topaz%id_co3_ion = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("co3_solubility","Carbonate Ion Solubility",'h','L','s','mol kg-1','f')
  topaz%id_co3_solubility = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("pco2surf","Oceanic pCO2",'h','1','s','uatm','f')
  topaz%id_pco2surf = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("sfc_chl","Surface Chl",'h','1','s','ug kg-1','f')
  topaz%id_sfc_chl = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("sfc_no3","Surface NO3",'h','1','s','mol kg-1','f')
  topaz%id_sfc_no3 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("htotal","h+ ion concentration",'h','L','s','mol kg-1','f')
  topaz%id_htotal = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("tot_layer_int_c","Total Carbon (DIC+OC+IC) boxwise layer integral",'h','L','s','mol m-2','f')
  topaz%id_tot_layer_int_c = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("tot_layer_int_fe","Total Iron (Fed_OFe) boxwise layer integral",'h','L','s','mol m-2','f')
  topaz%id_tot_layer_int_fe = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("tot_layer_int_n","Total Nitrogen (NO3+NH4+ON) boxwise layer integral",'h','L','s','mol m-2','f')
  topaz%id_tot_layer_int_n = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("tot_layer_int_p","Total Phosphorus (PO4+OP) boxwise layer integral",'h','L','s','mol m-2','f')
  topaz%id_tot_layer_int_p = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("tot_layer_int_si","Total Silicon (SiO4+SiO2) boxwise layer integral",'h','L','s','mol m-2','f')
  topaz%id_tot_layer_int_si = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("nLg_diatoms","Large diatom nitrogen",'h','L','s',' ','f')
  topaz%id_nLg_diatoms = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_cadet","CaCO3 production layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jprod_cadet = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_lithdet","Lithogenic removal to sinking layer integral",'h','L','s','g m-2 s-1','f')
  topaz%id_jprod_lithdet = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_nhet","heterotrophic N Production layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jprod_nhet = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_fedet","Detrital Fedet production layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jprod_fedet = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_ndet","Detrital PON production layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jprod_ndet = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jprod_pdet","Detrital phosphorus production layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jprod_pdet = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jcadet","CaCO3 change layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jcadet = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jfe_ads","Iron adsorption layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jfe_ads = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jfe_des","Iron desorption layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jfe_des = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jfe_graz","Dissolved iron grazing source layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jfe_graz = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jfe_coast","Coastal iron efflux layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jfe_coast = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jfedet","Loss of sinking iron layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jfedet = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jldon","Labile DON source layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jldon = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jndet","Loss of sinking nitrogen layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jndet = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jnh4","NH4 source layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jnh4 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jnh4_graz","NH4 grazing source layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jnh4_graz = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jnhet","heterotrophic N remineralization layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jnhet = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jnitrif","Nitrification layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jnitrif = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jno3","NO3 source layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jno3 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jo2","O2 source layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jo2 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jpdet","Loss of sinking phosphorus layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jpdet = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jpo4","PO4 source layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jpo4 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jpo4_graz","PO4 source from grazing layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jpo4_graz = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jdenit_wc","Water column Denitrification layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jdenit_wc = ocean_register_diag(vardesc_temp, G, model_time)
       
  vardesc_temp = vardesc("jdiss_sio2","SiO2 Dissolution during grazing layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jdiss_sio2 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jsdon","Semilabile DON source layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jsdon = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jsdop","Semilabile DOP source layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jsdop = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("jsidet","SiO4 source layer integral",'h','L','s','mol m-2 s-1','f')
  topaz%id_jsidet = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("fcaco3","CaCO3 detrital sinking flux",'h','L','s','mol m-2 s-1','f')
  topaz%id_fcaco3 = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("fcaco3_burial","CaCO3 permanent burial flux",'h','1','s','mol m-2 s-1','f')
  topaz%id_fcaco3_burial = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("fcaco3_redis","CaCO3 redissolution from sediments",'h','1','s','mol m-2 s-1','f')
  topaz%id_fcaco3_redis = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("fcaco3_sed","CaCO3 flux to sediment layer",'h','1','s','mol m-2 s-1','f')
  topaz%id_fcaco3_sed = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("fdenit_sed","Sediment Denitrification flux",'h','1','s','mol m-2 s-1','f')
  topaz%id_fdenit_sed = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("ffe_sed","Sediment iron efflux",'h','1','s','mol m-2 s-1','f')
  topaz%id_ffe_sed = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("ffedet_btm","Fe detrital sinking flux burial",'h','1','s','mol m-2 s-1','f')
  topaz%id_ffedet_btm = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("flith","Lith sinking flux",'h','L','s','g m-2 s-1','f')
  topaz%id_flith = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("flithdet_btm","Lithogenic detrital sinking flux burial",'h','1','s','g m-2 s-1','f')
  topaz%id_flithdet_btm = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("fpofe","POFe sinking flux",'h','L','s','mol m-2 s-1','f')
  topaz%id_fpofe = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("fpon","PON sinking flux",'h','L','s','mol m-2 s-1','f')
  topaz%id_fpon = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("fpop","POP sinking flux",'h','L','s','mol m-2 s-1','f')
  topaz%id_fpop = ocean_register_diag(vardesc_temp, G, model_time)

  vardesc_temp = vardesc("fsio2","Si sinking flux",'h','L','s','mol m-2 s-1','f')
  topaz%id_fsio2 = ocean_register_diag(vardesc_temp, G, model_time)

end subroutine register_diagnostics


subroutine set_prog_tracer(index, tracer_array, tracer_desc, param_file, G, &
                           restart_CS, tr_adv_CSp)
! Replaces call to otpm_set_prog_tracer
  integer,                         intent(in)  :: index
  type(tracer_type), dimension(:), pointer     :: tracer_array
  type(vardesc)        ,           intent(in)  :: tracer_desc
  type(param_file_type),           intent(in)  :: param_file
  type(ocean_grid_type),           intent(in)  :: G
  type(advect_tracer_CS),          pointer     :: tr_adv_CSp
  type(GOLD_restart_CS),            pointer     :: restart_CS
! Arguments: 
! (in)      G          - The ocean's grid structure.
! (in)      param_file - A structure indicating the open file to parse for
!                        model parameter values.
! (in/out)  tr_adv_CSp - A pointer that is set to point to the control structure
!                        for the tracer advection and diffusion module.
! (in)      restart_CS - A pointer to the restart control structure.

  integer :: isd, ied, jsd, jed, nz
  real, dimension(:,:,:), pointer :: tr_ptr

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke
  if ( index < 0 .or. index > size(tracer_array) ) &
     call mpp_error(FATAL, "set_prog_tracer : Tracer index is out of range")

 call mpp_error(NOTE, "Registering tracer "//trim(tracer_desc%name))

  tracer_array(index)%name = tracer_desc%name
  allocate(tracer_array(index)%field(isd:ied,jsd:jed,nz))
  tracer_array(index)%field(:,:,:) = 0.0
  allocate(tracer_array(index)%stf(isd:ied,jsd:jed))
  tracer_array(index)%stf  (:,:)   = 0.0
  allocate(tracer_array(index)%btf(isd:ied,jsd:jed))
  tracer_array(index)%btf  (:,:)   = 0.0
  allocate(tracer_array(index)%btm_reservoir(isd:ied,jsd:jed))
  tracer_array(index)%btm_reservoir(:,:) = 0.0
  allocate(tracer_array(index)%triver(isd:ied,jsd:jed))
  tracer_array(index)%triver(:,:)  = 0.0

  ! This pointer assignment is needed to force the compiler not to do a copy in
  ! the registration calls.
  tr_ptr => tracer_array(index)%field
  ! Register tracer for the restart file.
  call register_restart_field(tr_ptr, tr_ptr, tracer_desc, .false., restart_CS)
  ! Register tracer for horizontal advection & diffusion.
  call register_tracer(tr_ptr, tracer_desc%name, param_file, tr_adv_CSp)

end subroutine set_prog_tracer

subroutine set_diag_tracer(index, tracer_array, tracer_desc, param_file, G, &
                           restart_CS, tr_adv_CSp)
! Replaces call to otpm_set_diag_tracer
! Same as set_diag_tracer except for call to register_tracer
  integer,                         intent(in) :: index
  type(tracer_type), dimension(:), pointer     :: tracer_array
  type(vardesc)        ,           intent(in)  :: tracer_desc
  type(param_file_type),           intent(in)  :: param_file
  type(ocean_grid_type),           intent(in)  :: G
  type(advect_tracer_CS),          pointer     :: tr_adv_CSp
  type(GOLD_restart_CS),           pointer     :: restart_CS
! Arguments: 
! (in)      G          - The ocean's grid structure.
! (in)      param_file - A structure indicating the open file to parse for
!                        model parameter values.
! (in/out)  tr_adv_CSp - A pointer that is set to point to the control structure
!                        for the tracer advection and diffusion module.
! (in)      restart_CS - A pointer to the restart control structure.

  integer :: isd, ied, jsd, jed, nz
  real, dimension(:,:,:), pointer :: tr_ptr

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke
  if ( index < 0 .or. index > size(tracer_array) ) &
     call mpp_error(FATAL, "set_diag_tracer : Tracer index is out of range")
  tracer_array(index)%name = tracer_desc%name
  if (tracer_desc%z_grid .eq. '1') then
    allocate(tracer_array(index)%field(isd:ied,jsd:jed,1))
  elseif (tracer_desc%z_grid .eq. 'L') then
    allocate(tracer_array(index)%field(isd:ied,jsd:jed,nz))
  elseif (tracer_desc%z_grid .eq. 'i') then
    allocate(tracer_array(index)%field(isd:ied,jsd:jed,nz+1))
  else
    call mpp_error(FATAL, "set_diag_tracer : unknown vertical grid "//trim(tracer_desc%z_grid))
  endif
  tracer_array(index)%field(:,:,:) = 0.0

  ! This pointer assignment is needed to force the compiler not to do a copy in
  ! the registration calls.
  tr_ptr => tracer_array(index)%field
  ! Register tracer for the restart file.
  call register_restart_field(tr_ptr, tr_ptr, tracer_desc, .false., restart_CS)

end subroutine set_diag_tracer


!#######################################################################
! <SUBROUTINE NAME="allocate_arrays">
!
! <DESCRIPTION>
!     Dynamically allocate arrays
! </DESCRIPTION>
!

subroutine allocate_arrays(G, topaz, phyto)  !{
  type(ocean_grid_type),             intent(in)    :: G
  type(TOPAZ_CS),               intent(inout) :: topaz
  type(phytoplankton), dimension(:), intent(inout) :: phyto(:)

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!
integer :: isd, ied, jsd, jed, nk ,n

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nk = G%ke 

!-----------------------------------------------------------------------
!     
!       global variables
!
!
! allocate and initialize array elements for gas exchange
!
  allocate( sc_no_term(isd:ied,jsd:jed) )               ; sc_no_term(:,:)         = 0.0
  allocate( sc_co2(isd:ied,jsd:jed) )                   ; sc_co2(:,:)             = 0.0
  allocate( sc_o2(isd:ied,jsd:jed) )                    ; sc_o2(:,:)              = 0.0
  allocate( htotallo(isd:ied,jsd:jed,nk) )              ; htotallo(:,:,:)         = 0.0
  allocate( htotalhi(isd:ied,jsd:jed,nk) )              ; htotalhi(:,:,:)         = 0.0
  allocate( o2_saturation(isd:ied,jsd:jed) )            ; o2_saturation(:,:)      = 0.0
  allocate( tt(isd:ied) )                               ; tt(:)                   = 0.0
  allocate( tk(isd:ied) )                               ; tk(:)                   = 0.0
  allocate( ts(isd:ied) )                               ; ts(:)                   = 0.0
  allocate( ts2(isd:ied) )                              ; ts2(:)                  = 0.0
  allocate( ts3(isd:ied) )                              ; ts3(:)                  = 0.0
  allocate( ts4(isd:ied) )                              ; ts4(:)                  = 0.0
  allocate( ts5(isd:ied) )                              ; ts5(:)                  = 0.0
!
! allocate and initialize array elements of all phytoplankton groups
!
  do n = 1, NUM_PHYTO
    allocate(phyto(n)%def_fe(isd:ied,jsd:jed,nk))       ; phyto(n)%def_fe         = 0.0
    allocate(phyto(n)%def_p(isd:ied,jsd:jed,nk))        ; phyto(n)%def_p          = 0.0
    allocate(phyto(n)%felim(isd:ied,jsd:jed,nk))        ; phyto(n)%felim          = 0.0
    allocate(phyto(n)%irrlim(isd:ied,jsd:jed,nk))       ; phyto(n)%irrlim         = 0.0
    allocate(phyto(n)%jgraz_fe(isd:ied,jsd:jed,nk))     ; phyto(n)%jgraz_fe       = 0.0
    allocate(phyto(n)%jgraz_n(isd:ied,jsd:jed,nk))      ; phyto(n)%jgraz_n        = 0.0
    allocate(phyto(n)%jprod_fe(isd:ied,jsd:jed,nk))     ; phyto(n)%jprod_fe       = 0.0
    allocate(phyto(n)%jprod_nh4(isd:ied,jsd:jed,nk))    ; phyto(n)%jprod_nh4      = 0.0
    allocate(phyto(n)%jprod_no3(isd:ied,jsd:jed,nk))    ; phyto(n)%jprod_no3      = 0.0
    allocate(phyto(n)%jprod_po4(isd:ied,jsd:jed,nk))    ; phyto(n)%jprod_po4      = 0.0
    allocate(phyto(n)%liebig_lim(isd:ied,jsd:jed,nk))   ; phyto(n)%liebig_lim     = 0.0
    allocate(phyto(n)%mu(isd:ied,jsd:jed,nk))           ; phyto(n)%mu             = 0.0
    allocate(phyto(n)%po4lim(isd:ied,jsd:jed,nk))       ; phyto(n)%po4lim         = 0.0
    allocate(phyto(n)%q_fe_2_n(isd:ied,jsd:jed,nk))     ; phyto(n)%q_fe_2_n       = 0.0
    allocate(phyto(n)%q_p_2_n(isd:ied,jsd:jed,nk))      ; phyto(n)%q_p_2_n        = 0.0
    allocate(phyto(n)%q_p_2_n_opt(isd:ied,jsd:jed,nk))  ; phyto(n)%q_p_2_n_opt    = 0.0
    allocate(phyto(n)%theta(isd:ied,jsd:jed,nk))        ; phyto(n)%theta          = 0.0
  enddo
!
! allocate and initialize array elements of only some phytoplankton groups
!
  do n = 2, NUM_PHYTO
    allocate(phyto(n)%nh4lim(isd:ied,jsd:jed,nk))      ; phyto(n)%nh4lim          = 0.0
    allocate(phyto(n)%no3lim(isd:ied,jsd:jed,nk))      ; phyto(n)%no3lim          = 0.0
  enddo
!
! allocate and initialize array elements of only one phytoplankton group
!
  allocate(phyto(DIAZO)%jprod_n2(isd:ied,jsd:jed,nk))   ; phyto(DIAZO)%jprod_n2   = 0.0
  allocate(phyto(LARGE)%jgraz_sio2(isd:ied,jsd:jed,nk)) ; phyto(LARGE)%jgraz_sio2 = 0.0
  allocate(phyto(LARGE)%jprod_sio4(isd:ied,jsd:jed,nk)) ; phyto(LARGE)%jprod_sio4 = 0.0
  allocate(phyto(LARGE)%silim(isd:ied,jsd:jed,nk))      ; phyto(LARGE)%silim      = 0.0
!
! allocate and initialize array elements of topaz
!
  allocate(topaz%co3_solubility(isd:ied,jsd:jed,nk))    ; topaz%co3_solubility    = 0.0
  allocate(topaz%fcaco3_burial(isd:ied,jsd:jed))        ; topaz%fcaco3_burial     = 0.0
  allocate(topaz%fcaco3_redis(isd:ied,jsd:jed))         ; topaz%fcaco3_redis      = 0.0
  allocate(topaz%fcaco3_sed(isd:ied,jsd:jed))           ; topaz%fcaco3_sed        = 0.0
  allocate(topaz%fdenit_sed(isd:ied,jsd:jed))           ; topaz%fdenit_sed        = 0.0
  allocate(topaz%ffe_sed(isd:ied,jsd:jed))              ; topaz%ffe_sed           = 0.0
  allocate(topaz%ffedet_btm(isd:ied,jsd:jed))           ; topaz%ffedet_btm        = 0.0
  allocate(topaz%flithdet_btm(isd:ied,jsd:jed))         ; topaz%flithdet_btm      = 0.0
  allocate(topaz%jcadet(isd:ied,jsd:jed,nk))            ; topaz%jcadet            = 0.0
  allocate(topaz%jdenit_wc(isd:ied,jsd:jed,nk))         ; topaz%jdenit_wc         = 0.0
  allocate(topaz%jdiss_sio2(isd:ied,jsd:jed,nk))        ; topaz%jdiss_sio2        = 0.0
  allocate(topaz%jfe_ads(isd:ied,jsd:jed,nk))           ; topaz%jfe_ads           = 0.0
  allocate(topaz%jfe_des(isd:ied,jsd:jed,nk))           ; topaz%jfe_des           = 0.0
  allocate(topaz%jfe_graz(isd:ied,jsd:jed,nk))          ; topaz%jfe_graz          = 0.0
  allocate(topaz%jfe_coast(isd:ied,jsd:jed,nk))         ; topaz%jfe_coast         = 0.0
  allocate(topaz%jfedet(isd:ied,jsd:jed,nk))            ; topaz%jfedet            = 0.0
  allocate(topaz%jldon(isd:ied,jsd:jed,nk))             ; topaz%jldon             = 0.0
  allocate(topaz%jndet(isd:ied,jsd:jed,nk))             ; topaz%jndet             = 0.0
  allocate(topaz%jnh4(isd:ied,jsd:jed,nk))              ; topaz%jnh4              = 0.0
  allocate(topaz%jnh4_graz(isd:ied,jsd:jed,nk))         ; topaz%jnh4_graz         = 0.0
  allocate(topaz%jnhet(isd:ied,jsd:jed,nk))             ; topaz%jnhet             = 0.0
  allocate(topaz%jnitrif(isd:ied,jsd:jed,nk))           ; topaz%jnitrif           = 0.0
  allocate(topaz%jno3(isd:ied,jsd:jed,nk))              ; topaz%jno3              = 0.0
  allocate(topaz%jo2(isd:ied,jsd:jed,nk))               ; topaz%jo2               = 0.0
  allocate(topaz%jpdet(isd:ied,jsd:jed,nk))             ; topaz%jpdet             = 0.0
  allocate(topaz%jpo4(isd:ied,jsd:jed,nk))              ; topaz%jpo4              = 0.0
  allocate(topaz%jpo4_graz(isd:ied,jsd:jed,nk))         ; topaz%jpo4_graz         = 0.0
  allocate(topaz%jprod_cadet(isd:ied,jsd:jed,nk))       ; topaz%jprod_cadet       = 0.0
  allocate(topaz%jprod_lithdet(isd:ied,jsd:jed,nk))     ; topaz%jprod_lithdet     = 0.0
  allocate(topaz%jprod_nhet(isd:ied,jsd:jed,nk))        ; topaz%jprod_nhet        = 0.0
  allocate(topaz%jprod_fedet(isd:ied,jsd:jed,nk))       ; topaz%jprod_fedet       = 0.0
  allocate(topaz%jprod_ndet(isd:ied,jsd:jed,nk))        ; topaz%jprod_ndet        = 0.0
  allocate(topaz%jprod_pdet(isd:ied,jsd:jed,nk))        ; topaz%jprod_pdet        = 0.0
  allocate(topaz%jsdon(isd:ied,jsd:jed,nk))             ; topaz%jsdon             = 0.0
  allocate(topaz%jsdop(isd:ied,jsd:jed,nk))             ; topaz%jsdop             = 0.0
  allocate(topaz%jsidet(isd:ied,jsd:jed,nk))            ; topaz%jsidet            = 0.0
  allocate(topaz%k1_co2(isd:ied,jsd:jed,nk))            ; topaz%k1_co2            = 0.0
  allocate(topaz%k2_co2(isd:ied,jsd:jed,nk))            ; topaz%k2_co2            = 0.0
  allocate(topaz%mask_coast(isd:ied,jsd:jed))           ; topaz%mask_coast        = 0.0
  allocate(topaz%nLg_diatoms(isd:ied,jsd:jed,nk))       ; topaz%nLg_diatoms       = 0.0
  allocate(topaz%pco2surf(isd:ied,jsd:jed))             ; topaz%pco2surf          = 0.0
end subroutine  allocate_arrays  !}
! </SUBROUTINE> NAME="allocate_arrays"



!#######################################################################
! <SUBROUTINE NAME="ocean_topaz_source">
!
! <DESCRIPTION>
!     compute the source terms for the topaz tracers, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! </DESCRIPTION>

subroutine ocean_topaz_source(G, topaz,  &
      grid_dat, grid_tmask, grid_kmt,&
      rho_dzt, dzt, hblt_depth, fluxes, dt, tv, optics)  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!     Arguments
!-----------------------------------------------------------------------
!
type(ocean_grid_type),               intent(in)   :: G
type(TOPAZ_CS),                      pointer      :: topaz
real,    dimension(G%isd:,G%jsd:),   intent(in)   :: grid_dat
real,    dimension(G%isd:,G%jsd:,:), intent(in)   :: grid_tmask
integer, dimension(G%isd:,G%jsd:),   intent(in)   :: grid_kmt
real,    dimension(G%isd:,G%jsd:,:), intent(in)   :: rho_dzt
real,    dimension(G%isd:,G%jsd:,:), intent(in)   :: dzt
real,    dimension(G%isd:,G%jsd:),   intent(in)   :: hblt_depth
type(forcing),                       intent(in)   :: fluxes
real,                                intent(in)   :: dt
type(thermo_var_ptrs),               intent(in)   :: tv
type(optics_type),                   intent(in)   :: optics
!
! This subroutine calculates the biogeochemical source terms for the TOPAZ
! variables.
!
!  (in)      G - The ocean's grid structure.
!  (in)      topaz - The control structure returned by a previous call to
!                 GOLD_register_TOPAZ.
!  (in)      grid_dat - XY area of tracer cell (m2)
!  (in)      grid_tmask - tracer mask - true if ocean tracer cell, false if land
!  (in)      grid_kmt - bottom layer - variable in MOM, but either 0 or nk GOLD
!  (in)      rho_dzt - density multiplied by layer thickness
!  (in)      dzt - layer thicknesses (m)
!  (in)      hblt_depth - Depth of active mixing in mixed layer
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      tv - The structure containing thermodynamic variables.
!  (in)      optics - The structure containing optical properties.
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_topaz_source'
character(len=256), parameter   :: error_header =                            &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
character(len=256), parameter   :: warn_header =                             &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
character(len=256), parameter   :: note_header =                             &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '): '

integer    :: isc, iec, jsc, jec, nk
!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer :: i, j, k, n, nb
integer :: &
 ind_alk,     ind_cadet,    ind_dic,      ind_fed,      ind_fedet,    ind_fedi,  &
 ind_felg,    ind_fesm,     ind_ldon,     ind_lith,     ind_lithdet,  ind_ndet,  &
 ind_ndi,     ind_nh4,      ind_nhet,     ind_nlg,      ind_no3,      ind_nsm,   &
 ind_o2,      ind_pdet,     ind_pdi,      ind_plg,      ind_po4,      ind_psm,   &
 ind_sdon,    ind_sdop,     ind_sidet,    ind_silg,     ind_sio4,                &
 ind_cased,   ind_chl,      ind_co3_ion,  ind_fcadet_btm, ind_htotal,   ind_irr, &
 ind_irr_mem,  kblt
logical :: used

real :: epsln=1.0e-30
real :: feprime
real :: graz_lg_terms
real :: jgraz_n
real :: jgraz_p
real :: jgraz_fe
real :: jprod_di_tot2nterm
real :: log_btm_flx
real :: P_C_m
real :: p_lim_nhet
real :: r_assem,  r_photo,   r_uptake
real :: tmp_irrad_ML
real :: tmp_hblt
real, dimension(max(optics%nbands,1)) :: tmp_irr_band
real :: tmp_Irrad, tmp_opacity

real,dimension(NUM_PHYTO) :: p_2_n_max


real,dimension(G%isc:G%iec,G%jsc:G%jec,1:G%ke) :: &
   ccadet,   cfed,       cfedet,      cldon,      clith,     clithdet,   cndet,       &
   cnh4,     cnhet,      cno3,        c_o2,       cpdet,     cpo4,       csdon,       &
   csidet,   csiLg,      csio4,       expkT,      irr_mix,   frac_det_prod,           &
   q_si_2_n_Lg_Diatoms,  cco3_ion,    zt

real,dimension(G%isc:G%iec,G%jsc:G%jec,1:G%ke,NUM_PHYTO) :: cn

!
!---------------------------------------------------------------------
! Use shortened naming convention
!---------------------------------------------------------------------
!
  ind_alk     = topaz%ind_alk      ;  ind_cadet   = topaz%ind_cadet
  ind_dic     = topaz%ind_dic      ;  ind_fed     = topaz%ind_fed     
  ind_fedi    = topaz%ind_fedi     ;  ind_felg    = topaz%ind_felg    
  ind_fedet   = topaz%ind_fedet    ;  ind_fesm    = topaz%ind_fesm     
  ind_ldon    = topaz%ind_ldon     ;  ind_lith    = topaz%ind_lith    
  ind_lithdet = topaz%ind_lithdet  ;  ind_ndet    = topaz%ind_ndet    
  ind_ndi     = topaz%ind_ndi      ;  ind_nh4     = topaz%ind_nh4     
  ind_nhet    = topaz%ind_nhet     ;  ind_nlg     = topaz%ind_nlg     
  ind_no3     = topaz%ind_no3      ;  ind_nsm     = topaz%ind_nsm     
  ind_o2      = topaz%ind_o2       ;  ind_pdet    = topaz%ind_pdet    
  ind_pdi     = topaz%ind_pdi      ;  ind_plg     = topaz%ind_plg
  ind_po4     = topaz%ind_po4      ;  ind_psm     = topaz%ind_psm
  ind_sdon    = topaz%ind_sdon     ;  ind_sdop    = topaz%ind_sdop
  ind_sidet    = topaz%ind_sidet   ;  ind_silg    = topaz%ind_silg    
  ind_sio4    = topaz%ind_sio4     ;  ind_cased   = topaz%ind_cased
  ind_chl     = topaz%ind_chl      ;  ind_co3_ion = topaz%ind_co3_ion
  ind_htotal  = topaz%ind_htotal   ;  ind_irr     = topaz%ind_irr
  ind_irr_mem = topaz%ind_irr_mem  ;  ind_fcadet_btm  = topaz%ind_fcadet_btm

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; 
  nk = G%ke ; 
  
k=1
    do j = jsc, jec ; do i = isc, iec  !{
        htotallo(i,j,k) = htotal_scale_lo * Tr(ind_htotal)%field(i,j,k)
        htotalhi(i,j,k) = htotal_scale_hi * Tr(ind_htotal)%field(i,j,k)
    enddo ; enddo ; !} i, j
    call GOLD_ocmip2_co2calc(CO2_dope_vec,  &
       grid_tmask(:,:,1),                &
       tv%T(:,:,1),              &
       tv%S(:,:,1),              &
       Tr(ind_dic)%field(:,:,1),    &
       Tr(ind_po4)%field(:,:,1),    &
       Tr(ind_sio4)%field(:,:,1),   &
       Tr(ind_alk)%field(:,:,1),    &
       htotallo(:,:,1),             &
       htotalhi(:,:,1),             &
       Tr(ind_htotal)%field(:,:,1), &
       co2star=Tr(topaz%ind_csurf)%field(:,:,1),  &
       alpha=Tr(topaz%ind_alpha)%field(:,:,1),      &
       pco2surf=topaz%pco2surf(:,:),&
       co3_ion=Tr(ind_co3_ion)%field(:,:,1))


! The k=1 loop is done in surface_state
  do k = 2, nk ;
    do j = jsc, jec ; do i = isc, iec  !{
        htotallo(i,j,k) = htotal_scale_lo * Tr(ind_htotal)%field(i,j,k)
        htotalhi(i,j,k) = htotal_scale_hi * Tr(ind_htotal)%field(i,j,k)
    enddo ; enddo ; !} i, j

    call GOLD_ocmip2_co2calc(CO2_dope_vec,&
         grid_tmask(:,:,k),         &
         tv%T(:,:,k),               &
         tv%S(:,:,k),               &
         Tr(ind_dic)%field(:,:,k),  &
         Tr(ind_po4)%field(:,:,k),  &
         Tr(ind_sio4)%field(:,:,k), &
         Tr(ind_alk)%field(:,:,k),  &
         htotallo(:,:,k),           &
         htotalhi(:,:,k),           &
         Tr(ind_htotal)%field(:,:,k),&
         co3_ion=Tr(ind_co3_ion)%field(:,:,k))
  enddo  !} k

  n=DIAZO
  p_2_n_max(n) = phyto(DIAZO)%p_2_n_assem * (1.0 - phyto(n)%r_other_min -                 &
    phyto(n)%r_nfix) + topaz%p_2_n_RKR * phyto(n)%r_other_min
  do n = 2, NUM_PHYTO   !{
    p_2_n_max(n) = phyto(n)%p_2_n_assem * (1.0 - phyto(n)%r_other_min) +                  &
      topaz%p_2_n_RKR * phyto(n)%r_other_min
  enddo !} n
!
!---------------------------------------------------------------------
! Calculate positive tracer concentrations
!---------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
        c_o2(i,j,k)    = max(0.0,Tr(ind_o2)%field(i,j,k))
        ccadet(i,j,k)  = max(0.0,Tr(ind_cadet)%field(i,j,k))
        cco3_ion(i,j,k)= max(0.0,Tr(ind_co3_ion)%field(i,j,k))
        cfed(i,j,k)    = max(0.0,Tr(ind_fed)%field(i,j,k))
        cfedet(i,j,k)  = max(0.0,Tr(ind_fedet)%field(i,j,k))
        cldon(i,j,k)   = max(0.0,Tr(ind_ldon)%field(i,j,k))
        clith(i,j,k)   = max(0.0,Tr(ind_lith)%field(i,j,k))
        clithdet(i,j,k)= max(0.0,Tr(ind_lithdet)%field(i,j,k))
        cndet(i,j,k)   = max(0.0,Tr(ind_ndet)%field(i,j,k))
        cn(i,j,k,DIAZO)= max(0.0,Tr(ind_ndi)%field(i,j,k))
        cn(i,j,k,LARGE)= max(0.0,Tr(ind_nlg)%field(i,j,k))
        cn(i,j,k,SMALL)= max(0.0,Tr(ind_nsm)%field(i,j,k))
        cnh4(i,j,k)    = max(0.0,Tr(ind_nh4)%field(i,j,k))
        cnhet(i,j,k)   = max(0.0,Tr(ind_nhet)%field(i,j,k))
        cno3(i,j,k)    = max(0.0,Tr(ind_no3)%field(i,j,k))
        cpdet(i,j,k)   = max(0.0,Tr(ind_pdet)%field(i,j,k))
        cpo4(i,j,k)    = max(0.0,Tr(ind_po4)%field(i,j,k))
        csdon(i,j,k)   = max(0.0,Tr(ind_sdon)%field(i,j,k))
        csidet(i,j,k)  = max(0.0,Tr(ind_sidet)%field(i,j,k))
        csio4(i,j,k)   = max(0.0,Tr(ind_sio4)%field(i,j,k))
        csiLg(i,j,k)   = max(0.0,Tr(ind_silg)%field(i,j,k))

        zt(i,j,k) = 0.0

        phyto(DIAZO)%q_fe_2_n(i,j,k) = min(phyto(DIAZO)%fe_2_n_max, max(0.0,              &
          Tr(ind_fedi)%field(i,j,k) / max(epsln,Tr(ind_ndi)%field(i,j,k))))
        phyto(LARGE)%q_fe_2_n(i,j,k) = min(phyto(LARGE)%fe_2_n_max,max(0.0,               &
          Tr(ind_felg)%field(i,j,k) / max(epsln,Tr(ind_nlg)%field(i,j,k))))
        phyto(SMALL)%q_fe_2_n(i,j,k) = min(phyto(SMALL)%fe_2_n_max,max(0.0,               &
          Tr(ind_fesm)%field(i,j,k) / max(epsln,Tr(ind_nsm)%field(i,j,k))))
        phyto(DIAZO)%q_p_2_n(i,j,k) = min(p_2_n_max(DIAZO), max(0.0,                      &
          Tr(ind_pdi)%field(i,j,k) / max(epsln,Tr(ind_ndi)%field(i,j,k))))
        phyto(LARGE)%q_p_2_n(i,j,k) = min(p_2_n_max(LARGE), max(0.0,                      &
          Tr(ind_plg)%field(i,j,k) / max(epsln,Tr(ind_nlg)%field(i,j,k))))
        phyto(SMALL)%q_p_2_n(i,j,k) = min(p_2_n_max(SMALL), max(0.0,                      &
          Tr(ind_psm)%field(i,j,k) / max(epsln,Tr(ind_nsm)%field(i,j,k))))
  enddo ; enddo ; enddo !} i,j,k
!
!---------------------------------------------------------------------
!     Phytoplankton growth and grazing through the water column
!---------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! Calculate nutrient limitation terms
! 1d-30 added to avoid divide by zero where necessary
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
!
!-----------------------------------------------------------------------
! N limitation with NH4 inhibition after Frost and Franzen (1992)
!-----------------------------------------------------------------------
!
    do n= 2, NUM_PHYTO   !{
      phyto(n)%no3lim(i,j,k) = cno3(i,j,k) / ((phyto(SMALL)%k_no3 + cno3(i,j,k)) * (1.0 +     &
	  cnh4(i,j,k) / phyto(SMALL)%k_nh4))
      phyto(n)%nh4lim(i,j,k) = cnh4(i,j,k) / (phyto(SMALL)%k_nh4 + cnh4(i,j,k))
    enddo !} n
!
!-----------------------------------------------------------------------
! The rest are straight Michaelis Menten
!-----------------------------------------------------------------------
!
! phyto(LARGE)%k_sio4 was topaz(n)%k_sio4_Sm in MOM_TOPAZ code.
    phyto(LARGE)%silim(i,j,k) = csio4(i,j,k) / (phyto(LARGE)%k_sio4 + csio4(i,j,k))
    do n= 1, NUM_PHYTO   !{
      phyto(n)%po4lim(i,j,k) = cpo4(i,j,k) / (phyto(n)%k_po4 + cpo4(i,j,k))
      phyto(n)%felim(i,j,k) = cfed(i,j,k) / (phyto(n)%k_fed + cfed(i,j,k))
!
!-----------------------------------------------------------------------
! Calculate phosphorus and iron deficiency.  Phosphorus deficiency is approximated
! from the deviation from the maximum based on the amount of phosphorus needed for
! assembly.  Iron deficiency is approximated based on a sigmoidal half-saturation
! based on the observed relationship of Sunda and Huntsman (1997).
!-----------------------------------------------------------------------
!
      phyto(n)%def_p(i,j,k) = phyto(n)%q_p_2_n(i,j,k) / p_2_n_max(n)
      phyto(n)%def_fe(i,j,k) = phyto(n)%q_fe_2_n(i,j,k)**2.0 / (phyto(n)%k_fe_2_n**2.0 +  &
        phyto(n)%q_fe_2_n(i,j,k)**2.0)
    enddo !} n
  enddo ; enddo ;enddo !} i,j,k
!
! Calculate Leibig nutrient terms to for static and dynamic stoichiometry
!
  if (topaz%p_2_n_static) then  !{ p:n is static
    if (topaz%fe_2_n_static) then  !{ p:n and fe:n are static
      do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
        n=DIAZO
        phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%po4lim(i,j,k),phyto(n)%felim(i,j,k))
        do n= 2, NUM_PHYTO   !{
          phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%no3lim(i,j,k) +                       &
            phyto(n)%nh4lim(i,j,k), phyto(n)%po4lim(i,j,k), phyto(n)%felim(i,j,k))
        enddo !} n
      enddo ;  enddo ;  enddo !} i,j,k
    else  !{ p:n is static but fe:n is dynamic
      do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
        n=DIAZO
        phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%po4lim(i,j,k), phyto(n)%def_fe(i,j,k))
        do n= 2, NUM_PHYTO   !{
          phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k),&
            phyto(n)%po4lim(i,j,k), phyto(n)%def_fe(i,j,k))
        enddo !} n
      enddo ;  enddo ;  enddo !} i,j,k
    endif
  else
    if (topaz%fe_2_n_static) then  !{ fe:n is static but p:n is dynamic
      do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
        n=DIAZO
        phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%felim(i,j,k), phyto(n)%def_p(i,j,k))
        do n= 2, NUM_PHYTO   !{
          phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k),&
            phyto(n)%felim(i,j,k),phyto(n)%def_p(i,j,k))
        enddo !} n
      enddo ;  enddo ;  enddo !} i,j,k
    else  !{ pe:n and fe:n are dynamic
      do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
        n=DIAZO
        phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%def_fe(i,j,k), phyto(n)%def_p(i,j,k))
        do n= 2, NUM_PHYTO   !{
          phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k),&
            phyto(n)%def_fe(i,j,k),phyto(n)%def_p(i,j,k))
        enddo !} n
      enddo ;  enddo ;  enddo !} i,j,k
    endif
  endif
!
!-----------------------------------------------------------------------
! Calculate general ancillary terms
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! Phytoplankton in the actively mixed layer are assumed to be photoadapted to
! the mean light level in that layer as defined in the KPP routine plus one
! more vertical box to account for mixing directly below the boundary layer
!-----------------------------------------------------------------------
!
  do j = jsc, jec ; do i = isc, iec   !{
    do nb=1,optics%nbands !{
      if (optics%max_wavelength_band(nb) .lt. 710) then !{
        tmp_irr_band(nb) = optics%sw_pen_band(nb,i,j)
      else
        tmp_irr_band(nb) = 0.0
      endif !}
    enddo !}

    kblt = 0 ; tmp_irrad_ML = 0.0 ; tmp_hblt = 0.0
    do k = 1, nk !{
      tmp_Irrad = 0.0
      do nb=1,optics%nbands !{
        tmp_opacity = optics%opacity_band(nb,i,j,k)
        tmp_Irrad = tmp_Irrad + tmp_irr_band(nb) * exp(-tmp_opacity * dzt(i,j,k) * 0.5)
      !   Change tmp_irr_band from being the value atop layer k to the value
      ! at the bottom of layer k.
        tmp_irr_band(nb) = tmp_irr_band(nb) * exp(-tmp_opacity * dzt(i,j,k))
      enddo !}
      Tr(ind_irr)%field(i,j,k) = tmp_Irrad
      irr_mix(i,j,k) = tmp_Irrad
      if ((k == 1) .or. (tmp_hblt .lt. hblt_depth(i,j))) then !{
        kblt = kblt+1
        tmp_irrad_ML = tmp_irrad_ML + irr_mix(i,j,k) * dzt(i,j,k)
        tmp_hblt = tmp_hblt + dzt(i,j,k)
      endif !}
    enddo !} k-loop
    irr_mix(i,j,1:kblt) = tmp_irrad_ML / max(1.0e-6,tmp_hblt)
  enddo ; enddo !} i,j

  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
!
!-----------------------------------------------------------------------
! Temperature functionality of growth and grazing
!-----------------------------------------------------------------------
!
    expkT(i,j,k) = exp(topaz%kappa_eppley * tv%T(i,j,k))
!
!-----------------------------------------------------------------------
! Phytoplankton photoadaptation timescale
!-----------------------------------------------------------------------
!
    Tr(ind_irr_mem)%field(i,j,k) = Tr(ind_irr_mem)%field(i,j,k) + (irr_mix(i,j,k) - &
      Tr(ind_irr_mem)%field(i,j,k)) * min(1.0,topaz%gamma_irr_mem * dt)
  enddo ; enddo ; enddo !} i,j,k
!
!-----------------------------------------------------------------------
! This model predicts the Chl:N ratio at each time-step as an equilibrated
! phytoplankton response to the combined pressures of light, major nutrient
! and iron limitation.  Phytoplankton uptake is generally modelled after
! Geider et al. (1997) of steady state N and CO2 uptake, but also includes
! the following important modifications:
!
! 1) The temperature effect of Eppley (1972) is used instead
!    of that in Geider et al (1997) for both simplicity and
!    to incorporate combined effects on uptake, incorporation
!    into organic matter and photorespiration.  Values of PCmax
!    are normalized to 0C rather than 20C in Geider et al. (1997)
! 2) The Fe:N ratio is allowed to modulate the Chl:N ratio to be
!    consistent with Sunda and Huntsman (1997) through the "def_fe"
!    factor - the phytoplankton Fe:N ratio normalized to a saturated
!    value (k_fe_2_n) necessary to synthesize chlorophyll.
!    The def_fe calculation:
!
!    def_fe  = (Fe:N)**2 / [(Fe:N)irr**2 + (Fe:N)**2]
!
! 3) Values of the maximum Chl:C ratio (thetamax) are increased and
!    values of alpha decreased to account for the additional iron
!    term in the theta equation.
! 4) A thetamin value is also incorporated to set a minimum level of
!    chlorophyll per carbon
!
! While major nutrient limitation is handled through Michaelis Menten
! limitation of the phytoplankton specific growth prefactor (P_C_m), iron
! limitation is handled indirectly through modulatation of the Chl:N ratio. 
! This allows a compensatory relationship between irradiance and iron
! availability on phytoplankton specific growth, i.e. if they have a lot of
! light, they don't need a lot of iron and vice versa. Def_fe is assumed to be
! a quadratic function of the Fe:N ratio nomalized to vary between 0 and 1. 
! This relationship is a simple/crude representation of the complex
! physiological requirements and functionality of iron by separating
! phytoplankton iron into three components:
!     1) a "basal" requirement of iron for phytoplankton respiration and protein
!            synthesis (e.g. the electron transport chain)
!     2) Chlorophyll synthesis for photosynthesis
!     3) Luxury uptake
! While somewhat mathematically ad-hoc, this representation is grounded in the
! observed relationship between Chl:C, Fe:C, dissolved Fe and phytoplankton
! specific growth rates of Sunda and Huntsman (1997) as well our general
! understanding of the role of iron in phytoplankton physiology (e.g Geider and
! La Rocha; 1994; Photosynthesis Res., 39, 275-301).  P_C_m is calculated as a Liebig 
! nutrient- or stoichiometry- limited rate modulated by the Eppley-type Temperature
! function.
!
!   P_C_m = P_C_max * e**(k*T) * Liebig_lim
!
! Chl:C calculation after Geider et al (1997) but with addition of a minimum Chl term
! that is also modulated by nutrient limitation after optimal allocation work:
!
!   Chl:C = (Chl:C_max-Chl:C_min) / (1+ (Chl:C_max-Chl:C_min) * a * I / 2 / P_C_m)
!      + Chl:C_min * Liebig_lim
!
! Growth rate calculation after Geider et al (1997):
!
! mu = P_C_m / (1 + zeta) * (1 - e**( - alpha * I * Chl:C / P_C_m))
!
! Total Chlorophyll is calculated for use in the short-wave absorption module,
! ocean_shortwave_pen.F90.  Note: For this chlorophyll concentration to have
! an active effect on irradiance and shortwave absorption, the variable
! read_chl=.false. must be set in ocean_shortwave_pen_nml.  Otherwise,
! that module will use the data override values.
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
    Tr(ind_chl)%field(i,j,k) = 0.0
    do n = 1, NUM_PHYTO   !{
      P_C_m = phyto(n)%liebig_lim(i,j,k) * phyto(n)%P_C_max * expkT(i,j,k) + epsln
      phyto(n)%theta(i,j,k) = (phyto(n)%thetamax - topaz%thetamin) / (1.0 +               &
        (phyto(n)%thetamax - topaz%thetamin) * phyto(n)%alpha *                           &
        Tr(ind_irr_mem)%field(i,j,k) * 0.5 / P_C_m) + topaz%thetamin *                    &
        phyto(n)%liebig_lim(i,j,k)
      Tr(ind_chl)%field(i,j,k) = Tr(ind_chl)%field(i,j,k) + topaz%c_2_n * 12.0e6 *        &
        phyto(n)%theta(i,j,k) * cn(i,j,k,n)
      phyto(n)%irrlim(i,j,k) = (1.0 - exp(-phyto(n)%alpha * Tr(ind_irr)%field(i,j,k) *    &
        phyto(n)%theta(i,j,k) / P_C_m))
      phyto(n)%mu(i,j,k) = P_C_m / (1.0 + topaz%zeta) * phyto(n)%irrlim(i,j,k)
    enddo !} n
  enddo ;  enddo ; enddo !} i,j,k

!
!-----------------------------------------------------------------------
! Apply production terms:
!
! NO3 and NH4 uptake are calculated as fractions of total N
! uptake
!
! Diazotrophs produce organic N from N2
!
! PO4 production is assumed to be stoichiometric to N with
! same stoichiometry for Sm and Lg but higher for Di
!
! Large and Small phytoplankton are allowed to always take up as much Iron
! as they can after Sunda and Huntman (1997).  Small phytoplankton are forced
! to diminish their uptake at saturated levels of the Fe:C ratio in small
! phytoplankton (to mimic their general lack of luxury storage capacity).
!
! Si uptake is make to be consistent with the Si:N ratio
! synthesis of Martin-Jezequel et al (2000) and the Droop
! quota argument of Mongin et al. (submitted)
!
! CaCO3 formation is set to go directly to detritus as a constant
! fraction of Sm production after Moore et al (2002)
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
!
! Nitrogen fixation is assumed to be inhibited by high oxygen after Stewart and
! Pearson (1970; Effects of aerobic and anaerobic conditions on growth and
! metabolism of blue-green algae. Proc. Soc. Lond. B., 175, 293-311).  Nitrogen
! fixation is also assumed to be limited by nitrate after Holl and Montoya
! (2005;  Interactions between fixed nitrogen uptake and nitrogen fixation in
! continuous cultures of the marine diazotroph trichodesmium cyanobacteria); J.
! Phycol, 41, 1178-1183).
!
    n=DIAZO
    jprod_di_tot2nterm=phyto(n)%mu(i,j,k) * cn(i,j,k,n) * (1.0 -                          &
      c_o2(i,j,k)**topaz%o2_inhib_Di_pow / (c_o2(i,j,k)**topaz%o2_inhib_Di_pow +          &
      topaz%o2_inhib_Di_sat**topaz%o2_inhib_Di_pow)) / (cno3(i,j,k) + cnh4(i,j,k) +       &
      topaz%k_n_inhib_Di)
    phyto(n)%jprod_n2(i,j,k) = jprod_di_tot2nterm * topaz%k_n_inhib_Di      
    phyto(n)%jprod_nh4(i,j,k) = jprod_di_tot2nterm * cnh4(i,j,k)
    phyto(n)%jprod_no3(i,j,k) = jprod_di_tot2nterm * cno3(i,j,k)
    do n = 2, NUM_PHYTO
      phyto(n)%jprod_no3(i,j,k) = phyto(n)%mu(i,j,k) * cn(i,j,k,n) *                      &
        phyto(n)%no3lim(i,j,k) / (phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k) + epsln)
      phyto(n)%jprod_nh4(i,j,k) = phyto(n)%mu(i,j,k) * cn(i,j,k,n) *                      &
        phyto(n)%nh4lim(i,j,k) / (phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k) + epsln)
    enddo !} n
  enddo ; enddo ; enddo !} i,j,k

  if (topaz%p_2_n_static) then  !{
    do k = 1, nk  ;    do j = jsc, jec ;      do i = isc, iec   !{
      n=DIAZO
      phyto(n)%jprod_po4(i,j,k) = (phyto(n)%jprod_n2(i,j,k) + phyto(n)%jprod_nh4(i,j,k) + &
        phyto(n)%jprod_no3(i,j,k)) * phyto(n)%p_2_n_static
      do n = 1, NUM_PHYTO
        phyto(n)%jprod_po4(i,j,k) = (phyto(n)%jprod_no3(i,j,k) +                          &
            phyto(n)%jprod_nh4(i,j,k)) * phyto(n)%p_2_n_static
      enddo !} n
    enddo ; enddo ; enddo !} i,j,k
  else
    do k = 1, nk  ;    do j = jsc, jec ;      do i = isc, iec   !{
!
! Phosphate uptake patterned after the optimal allocation theory of Klausmeier et al.
! 2004 in which N:P is determined by the allocation towards photosynthetic machinery,
! nutrient uptake, assemby, and other.  Allocation for "other" is held fixed at 0.2 
! after Klausmeier et al. (2004).  Allocation towards photosynthetic machinery
! is taken from the instantaneous C:Chl ratio which is renormalized to the ratio of
! chlorophyll to chloroplasts in the cell, allocation towards uptake is taken from
! the degree of nutrient limitation normalized to the remaining space, and assembly is
! then taken from the difference, i.e.:
!    R_other=0.2
!    R_photo=1/1.87/0.07*Chl2C
!    R_uptake=(1-min(plim,nlim))*(1-R_other-R_photo)
!    R_assembly=1-R_other-R_photo-R_uptake
!    p2n_opt=p2n_assembly*R_assembly+p2n_photo*R_photo+p2n_other*R_other
!
! For efficiency, this is all condensed into a single calculation.  Phytoplankton are
! then allowed to accelerate their approach towards this optimum by taking up more or
! less phosphate than the optimum N:P linearly with respect to their relative N:P offset.
!
! In order to simulate the role of storage/vacuoles, Eukaryotes (Lg) are allowed to maximize
! r_assem while Prokaryotes (Sm, Di) are forced to maximize r_other
!
! Here, the maximum assembly rate is assumed to be equal to the maximum photosynthetic
! rate multiplied by 1/(1 - r_other_min).  However, there is some indication that the
! maximum assembly rate should be double the maximum photosynthetic rate
! based on comparisons of the envelope for maximum growth rates of bacteria of
! ~0.72e-5 * exp(0.063*T) [s-1] for whole organism phytoplankton after Eppley 
! (1972, Fish. Bull., 70, 1063-1084) and ~1.6e-5 * exp(0.063*T) [s-1] for whole
! organism heterotrophic bacteria after Ducklow and Hill (1985, Limnol. 
! Oceanogr., 30, 239-259).  Assuming that the r_assem of heterotrophic bacteria is 0.8
! at the maximum growth rate (e.g. r_other=0.2 after Klausmeier et al., 2004) and 
! r_assem of phytoplankton is 0.64 at the maximum growth rate (Klausmeier et al., 2004)
! gives P_C_max_assem = 1.6e-5 / 0.72e-5 * 0.64 / 0.8 * P_C_max.  The conversion is highly
! uncertain however, and leave room for future tuning.

      n=DIAZO
      r_photo = phyto(n)%plast_2_chl * phyto(n)%theta(i,j,k)
      r_uptake = (1.0 - min(phyto(n)%po4lim(i,j,k), phyto(n)%felim(i,j,k))) *             &
            phyto(n)%r_uptake_max
!
! diaz phytoplankton are assumed to take up only as much PO4 for assembly as they need
! based on their overall nutrient limited growth rate, thereby assuming that they cannot
! store extra PO4.
!
      r_assem = min(1.0 - phyto(n)%r_nfix - phyto(n)%r_other_min - r_photo - r_uptake,    &
        (1.0 - phyto(n)%r_other_min) * phyto(n)%liebig_lim(i,j,k))
      phyto(n)%q_p_2_n_opt(i,j,k) = phyto(n)%p_2_n_assem * r_assem + topaz%p_2_n_uptake * &
        r_uptake + topaz%p_2_n_photo * r_photo + topaz%p_2_n_RKR * (1.0 -                 &
        phyto(n)%r_nfix - r_photo - r_uptake - r_assem)
      phyto(n)%jprod_po4(i,j,k) = min(phyto(n)%po4lim(i,j,k) * phyto(n)%P_C_max *         &
        expkT(i,j,k) * cn(i,j,k,n), phyto(n)%jprod_n2(i,j,k) + phyto(n)%jprod_no3(i,j,k) +&
        phyto(n)%jprod_nh4(i,j,k)) * phyto(n)%q_p_2_n_opt(i,j,k)
      n=LARGE
      r_photo = phyto(n)%plast_2_chl * phyto(n)%theta(i,j,k)
      r_uptake = (1.0 - min(phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k),              &
        phyto(n)%po4lim(i,j,k), phyto(n)%felim(i,j,k))) * phyto(n)%r_uptake_max
!
! large phytoplankton are assumed to take up only as much PO4 for assembly as they need
! based on their overall PO4 nutrient limited growth rate, thereby assuming that they can
! store extra PO4.
!
      r_assem = min(1.0 - phyto(n)%r_other_min - r_photo - r_uptake,                      &
        (1.0 - phyto(n)%r_other_min) * phyto(n)%po4lim(i,j,k))
      phyto(n)%q_p_2_n_opt(i,j,k) = phyto(n)%p_2_n_assem * r_assem +                      &
        topaz%p_2_n_uptake * r_uptake + topaz%p_2_n_photo * r_photo + topaz%p_2_n_RKR *   &
        (1.0 -r_photo - r_uptake - r_assem)
      phyto(n)%jprod_po4(i,j,k) = min(phyto(n)%po4lim(i,j,k) * phyto(n)%P_C_max *         &
        expkT(i,j,k) * cn(i,j,k,n),phyto(n)%jprod_no3(i,j,k) +                            &
        phyto(n)%jprod_nh4(i,j,k)) * phyto(n)%q_p_2_n_opt(i,j,k)
      n=SMALL
      r_photo = phyto(n)%plast_2_chl * phyto(n)%theta(i,j,k)
      r_uptake = (1.0 - min(phyto(n)%no3lim(i,j,k) + phyto(n)%nh4lim(i,j,k),              &
        phyto(n)%po4lim(i,j,k), phyto(n)%felim(i,j,k))) * phyto(n)%r_uptake_max
!
! small phytoplankton are assumed to take up only as much PO4 for assembly as they need
! based on their overall nutrient limited growth rate, thereby assuming that they cannot
! store extra PO4.
!
      r_assem = min(1.0 - phyto(n)%r_other_min - r_photo - r_uptake,                      &
        (1.0 - phyto(n)%r_other_min) * phyto(n)%liebig_lim(i,j,k))
      phyto(n)%q_p_2_n_opt(i,j,k) = phyto(n)%p_2_n_assem * r_assem +                      &
        topaz%p_2_n_uptake * r_uptake + topaz%p_2_n_photo * r_photo + topaz%p_2_n_RKR *   &
        (1.0 - r_photo - r_uptake - r_assem)
      phyto(n)%jprod_po4(i,j,k) = min(phyto(n)%po4lim(i,j,k) * phyto(n)%P_C_max *         &
        expkT(i,j,k) * cn(i,j,k,n), phyto(n)%jprod_no3(i,j,k) +                           &
        phyto(n)%jprod_nh4(i,j,k)) * phyto(n)%q_p_2_n_opt(i,j,k)
    enddo ; enddo ; enddo !} i,j,k
  endif

  if (topaz%fe_2_n_static) then  !{
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
      phyto(DIAZO)%jprod_fe(i,j,k) = (phyto(DIAZO)%jprod_n2(i,j,k) +                      &
        phyto(DIAZO)%jprod_nh4(i,j,k) + phyto(DIAZO)%jprod_no3(i,j,k)) *                  &
        phyto(DIAZO)%fe_2_n_static
      do n = 2, NUM_PHYTO  !{
        phyto(n)%jprod_Fe(i,j,k) = (phyto(n)%jprod_no3(i,j,k) +                           &
          phyto(n)%jprod_nh4(i,j,k)) * phyto(n)%fe_2_n_static
      enddo
    enddo ; enddo ; enddo !} i,j,k
  else
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
      do n = 1, NUM_PHYTO  !{
        phyto(n)%jprod_Fe(i,j,k) = phyto(n)%P_C_max * expkT(i,j,k) * cn(i,j,k,n) *        &
          phyto(n)%felim(i,j,k) * (phyto(n)%fe_2_n_max - phyto(n)%q_fe_2_n(i,j,k))
      enddo   !} n
    enddo ; enddo ; enddo !} i,j,k
  endif

  if (topaz%si_2_n_static) then  !{
    do k = 1, nk  ; do j = jsc, jec ; do i = isc, iec   !{
      topaz%nLg_diatoms(i,j,k) = cn(i,j,k,LARGE) * phyto(LARGE)%silim(i,j,k)
      q_si_2_n_Lg_diatoms(i,j,k) = phyto(LARGE)%si_2_n_static
      phyto(LARGE)%jprod_sio4(i,j,k) = (phyto(LARGE)%jprod_no3(i,j,k) +                   &
        phyto(LARGE)%jprod_nh4(i,j,k)) * phyto(LARGE)%silim(i,j,k) *                      &
        phyto(LARGE)%si_2_n_static
    enddo ; enddo ; enddo !} i,j,k
  else
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
!
! The fraction of large phytoplankton that is diatoms is set equal to silim.  All
! other large phytoplankton are assumed to be non-diatom, e.g. green algae,
! dinoflagellates, phaeocystis, etc.  SiO4 uptake is limited by the maximum quota
! as for PO4 and Fed.
!
      topaz%nLg_diatoms(i,j,k) = cn(i,j,k,LARGE) * phyto(LARGE)%silim(i,j,k)
      q_si_2_n_Lg_diatoms(i,j,k) = min(phyto(LARGE)%si_2_n_max, csilg(i,j,k) /            &
        (topaz%nLg_diatoms(i,j,k) + epsln))
      phyto(LARGE)%jprod_sio4(i,j,k) = phyto(LARGE)%P_C_max * expkT(i,j,k) *              &
        topaz%nLg_diatoms(i,j,k) * phyto(LARGE)%silim(i,j,k) * max(0.0,                   &
        phyto(LARGE)%si_2_n_max - q_si_2_n_Lg_diatoms(i,j,k))
    enddo ; enddo ; enddo !} i,j,k
  endif

!
!-----------------------------------------------------------------------
!
!     Food Web Processing
!
! Phytoplankton loss:
! Sm loss is proportional to Sm**2,
! Lg loss is proportional to Lg**(4/3)
! after Dunne et al, 2005.
! for the purposes of grazing, Di and Lg are both considered part of the
! "Large" or traditional ecosystem pool. Prey switching/selectivity is assumed
! to factor into determining the overall grazing pressure on these two groups
! after Fasham et al. (1990; JMR, ) except that the effective size is used
! as a proxy for particle number for encounter rates with zooplankton predators.
! For this case, the sum of these species is used for the overall grazing
! pressure (nlg_tot=ndi+nlg) and the individual grazing pressure on each 
! component is determined by:
!
!  Pi / eff_size_i / sum((Pi / eff_size_i)^2)^(1/2)
!
! Additionally two criteria for numerical stability are added:
!    1) the absolute first order rate constant is
! never allowed to be greater than 1/dt for numerical stability.
!    2) a Michaelis-Menton type of threshold using a half
! saturation value of phyto_min is set to prevent phytoplankton
! from going extinct at low concentrations.
!
! Nitrification is set to be inhibited by light after Ward et al. (1982)
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
    graz_Lg_terms = topaz%lambda0 * expkT(i,j,k) * ((cn(i,j,k,DIAZO) + cn(i,j,k,LARGE)) / &
      topaz%P_star)**(1.0/3.0) * (cn(i,j,k,DIAZO) + cn(i,j,k,LARGE)) /                    &
      (cn(i,j,k,DIAZO) + cn(i,j,k,LARGE) + topaz%phyto_min) /                             &
      (cn(i,j,k,DIAZO)**2 + cn(i,j,k,LARGE)**2 + epsln)**(0.5)
    phyto(DIAZO)%jgraz_n(i,j,k) = min( 1.0 / dt , graz_Lg_terms * cn(i,j,k,DIAZO)) *    &
      cn(i,j,k,DIAZO)
    phyto(LARGE)%jgraz_n(i,j,k) = min( 1.0 / dt , graz_Lg_terms * cn(i,j,k,LARGE)) *    &
      cn(i,j,k,LARGE)
    phyto(SMALL)%jgraz_n(i,j,k) = min( 1.0 / dt , topaz%lambda0 * expkT(i,j,k) *        &
      cn(i,j,k,SMALL) ** 2.0 / (topaz%P_star * (cn(i,j,k,SMALL) + topaz%phyto_min))) *    &
      cn(i,j,k,SMALL)
  enddo ; enddo ; enddo !} i,j,k

  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
    jgraz_n = 0.0
    jgraz_p = 0.0
    do n = 1, NUM_PHYTO  !{
      jgraz_n = jgraz_n + phyto(n)%jgraz_n(i,j,k)
      jgraz_p = jgraz_p + phyto(n)%q_p_2_n(i,j,k) * phyto(n)%jgraz_n(i,j,k)
    enddo   !} n
    topaz%jprod_ndet(i,j,k) = (phyto(SMALL)%fdet0 * (phyto(SMALL)%jgraz_n(i,j,k) +         &
      phyto(DIAZO)%jgraz_n(i,j,k)) + phyto(LARGE)%fdet0 * phyto(LARGE)%jgraz_n(i,j,k)) *  &
      (1.0 - topaz%phi_sdon - topaz%phi_ldon) * exp(topaz%kappa_remin * tv%T(i,j,k))
    frac_det_prod(i,j,k) = topaz%jprod_ndet(i,j,k) / (jgraz_n + epsln)
    p_lim_nhet = min(1.0, jgraz_p / (jgraz_n + epsln) / topaz%p_2_n_RKR)
    topaz%jprod_pdet(i,j,k) = frac_det_prod(i,j,k) * jgraz_p * p_lim_nhet
    topaz%jsdon(i,j,k) = topaz%phi_sdon * jgraz_n
    topaz%jsdop(i,j,k) = topaz%phi_sdop * jgraz_p
    topaz%jldon(i,j,k) = topaz%phi_ldon * jgraz_n * p_lim_nhet
    topaz%jprod_nhet(i,j,k) = (jgraz_n - topaz%jprod_ndet(i,j,k) - topaz%jsdon(i,j,k) -    &
      topaz%jldon(i,j,k)) * topaz%phi_nhet * p_lim_nhet
    topaz%jnh4_graz(i,j,k) = jgraz_n - topaz%jprod_ndet(i,j,k) - topaz%jsdon(i,j,k) -      &
      topaz%jldon(i,j,k) - topaz%jprod_nhet(i,j,k)
    topaz%jpo4_graz(i,j,k) = jgraz_p - topaz%jprod_pdet(i,j,k) - topaz%jsdop(i,j,k) -      &
      (topaz%jldon(i,j,k) + topaz%jprod_nhet(i,j,k)) * topaz%p_2_n_RKR
    topaz%jnhet(i,j,k) = topaz%gamma_nhet * expkT(i,j,k) * cnhet(i,j,k) 
    topaz%jnitrif(i,j,k) = topaz%gamma_nitrif * expkT(i,j,k) * cnh4(i,j,k) *              &
      phyto(SMALL)%nh4lim(i,j,k) * (1.0 - Tr(ind_irr_mem)%field(i,j,k) /                  &
      (topaz%irr_inhibit + Tr(ind_irr_mem)%field(i,j,k)))
!
! Lithogenic material is assumed to get converted into the sinking particulate phase
! through incorporation into mesozooplankton fecal pellets with an efficiency of phi_lith.
!
    topaz%jprod_lithdet(i,j,k) = (phyto(LARGE)%jgraz_n(i,j,k) / (cn(i,j,k,LARGE) +        &
      epsln) * topaz%phi_lith + topaz%k_lith) * clith(i,j,k)
  enddo ; enddo ; enddo !} i,j,k

  do j = jsc, jec ;      do i = isc, iec   !{
    zt(i,j,1) = dzt(i,j,1)
  enddo ; enddo !} i,j

  do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
    zt(i,j,k) = zt(i,j,k-1) + dzt(i,j,k)
  enddo ; enddo ; enddo !} i,j,k

  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
!
!---------------------------------------------------------------------
! CaCO3 solubility taken from Sayles, F. L. (1985, CaCO3 solubility in marine 
! sediments: evidence for equilibrium and non-equilibrium behavior, Geochim. 
! Cosmochim. Acta, 49, 877-888) in which:
!
!   co3_solubility = ksp_caco3*exp(-deltaV / R / T * press / 10) / 
!                                              (Ca_ave * S / S_ave / rho / rho)
!
! where deltaV = -41.2 cm^3/mol, R = 82.057 cm^3 mol-1 degree K-1, Ca_ave = 
! 0.010233, S_ave = 35 and rho = 1.025 to give the accumulated constants: 
! -deltaV/(R*10) = 0.05021, S_ave/Ca_ave = 3420.31 and unit conversion for rho^2 of 1e-6.

    topaz%co3_solubility(i,j,k) = max(topaz%ksp_caco3 * exp ( 0.05021 / (tv%T(i,j,k) +    &
      273.15) * zt(i,j,k)) * 3.42031e-3 * topaz%Rho_0 * topaz%Rho_0 / max(epsln,          &
      tv%S(i,j,k)),epsln)
  enddo ; enddo ; enddo !} i,j,k

  if (topaz%ca_2_n_static) then  !{
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
      topaz%jprod_cadet(i,j,k) = (phyto(DIAZO)%jgraz_n(i,j,k) +                           &
        phyto(LARGE)%jgraz_n(i,j,k) + phyto(SMALL)%jgraz_n(i,j,k)) *                      &
        topaz%ca_2_n_het_static
    enddo ; enddo ; enddo !} i,j,k
  else
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
! CaCO3 production is assumed to be proportional to both calcite supersaturation and
! mesozooplankton grazing where a relative transfer efficiency of 0.1 is assumed for the 
! Diazotroph and Small to microzooplankton to mesozooplankton.  This was done to agree with
! observations of a tropical maximum in the Ca:Corg ratio (Sarmiento et al., 2002; 
! Jin et al., 2006)
      topaz%jprod_cadet(i,j,k) = (phyto(LARGE)%jgraz_n(i,j,k) + 0.01 *                    &
        (phyto(DIAZO)%jgraz_n(i,j,k) + phyto(SMALL)%jgraz_n(i,j,k))) * topaz%ca_2_n_het * &
        min(topaz%caco3_sat_max, max(0.0,cco3_ion(i,j,k) / topaz%co3_solubility(i,j,k) -  &
        1.0))
    enddo ; enddo ; enddo !} i,j,k
  endif
!
!-----------------------------------------------------------------------
! Iron and Silicon Processing:
!
! Iron is recycled with the same efficiency as nitrogen.
!
! SiO2 dissolution is set to globally dissolve 50% after Nelson et al. (1995) 
! through grazing.  The temperature functionality is set to a combination
! Michaelis  Menton and Eppley temperature formulation to give roughly the range
! of  observations in Nelson et al. (1995), Blain et al. (1999) and Bzrezenski's
! work... 
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
    jgraz_fe = 0.0
    do n = 1, NUM_PHYTO  !{
     phyto(n)%jgraz_fe(i,j,k) = phyto(n)%jgraz_n(i,j,k) * phyto(n)%q_fe_2_n(i,j,k)
     jgraz_fe = jgraz_fe + phyto(n)%jgraz_fe(i,j,k)
    enddo   !} n
    topaz%jprod_fedet(i,j,k) = frac_det_prod(i,j,k) * jgraz_fe
    topaz%jfe_graz(i,j,k) = jgraz_fe - topaz%jprod_fedet(i,j,k)
    phyto(LARGE)%jgraz_sio2(i,j,k) = phyto(LARGE)%jgraz_n(i,j,k) *                        &
      csiLg(i,j,k) / (cn(i,j,k,LARGE) + epsln)
    topaz%jdiss_sio2(i,j,k) = phyto(LARGE)%jgraz_sio2(i,j,k) *                            &
      exp(-q_si_2_n_Lg_diatoms(i,j,k) / (topaz%q_si_2_n_diss * expkT(i,j,k)))
  enddo ; enddo ; enddo !} i,j,k
!
!-----------------------------------------------------------------------
!   Ballast Protection Interior Remineralization Scheme and Iron scavenging
!-----------------------------------------------------------------------
!
!
  do k=1,nk ; do j=jsc,jec ; do i=isc,iec  !{
    topaz%jcadet(i,j,k) = topaz%gamma_cadet * (1.0 - min(1.0, cco3_ion(i,j,k) /           &
      topaz%co3_solubility(i,j,k))) * ccadet(i,j,k)
    topaz%jsidet(i,j,k) = topaz%gamma_sidet * csidet(i,j,k)
    topaz%jdenit_wc(i,j,k) = 0.0
!
!---------------------------------------------------------------------
! Remineralization of unprotected organic material and
! previously protected particulate organic material
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!   Under oxic conditions
!---------------------------------------------------------------------
!
    if (c_o2(i,j,k) .gt. topaz%o2_min) then  !{
      topaz%jndet(i,j,k) = topaz%gamma_ndet * c_o2(i,j,k) / (topaz%k_o2 + c_o2(i,j,k)) *  &
        max(0.0, cndet(i,j,k) - (topaz%rpcaco3 * ccadet(i,j,k) + topaz%rplith *           &
        clithdet(i,j,k) + topaz%rpsio2 * csidet(i,j,k)))
    else !}{
!
!---------------------------------------------------------------------
!   Under suboxic conditions
!---------------------------------------------------------------------
!
      topaz%jndet(i,j,k) = topaz%gamma_ndet * topaz%o2_min / (topaz%k_o2 + topaz%o2_min) *&
        max(0.0, cndet(i,j,k) - (topaz%rpcaco3 * ccadet(i,j,k) + topaz%rplith *           &
        clithdet(i,j,k) + topaz%rpsio2 * csidet(i,j,k)))
      topaz%jdenit_wc(i,j,k) = topaz%jndet(i,j,k) * topaz%n_2_n_denit
    endif !}
!
!---------------------------------------------------------------------
! Apply N change to P assuming equal partitioning between protected,
! previously protected and unprotected particulate organic material
!---------------------------------------------------------------------
!
    topaz%jpdet(i,j,k) = topaz%jndet(i,j,k) / (cndet(i,j,k) + epsln) * cpdet(i,j,k)
!
!---------------------------------------------------------------------
! Apply N change to Fe incorporating adsorption and desorption
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
! Calculate free and inorganically associated iron concentration for scavenging
!---------------------------------------------------------------------
!
    feprime = 1.0 + topaz%kfe_eq_lig * (topaz%felig_bkg + topaz%felig_2_don *             &
      (cldon(i,j,k) + csdon(i,j,k)) - cfed(i,j,k))
    feprime = (-feprime + (feprime * feprime + 4.0 * topaz%kfe_eq_lig *                   &
      cfed(i,j,k))**(0.5)) / (2.0 * topaz%kfe_eq_lig)
!
!---------------------------------------------------------------------
! The absolute first order rate constant is never allowed to be greater than
! 1/dt for numerical stability.
!---------------------------------------------------------------------
!
    topaz%jfe_ads(i,j,k) = min(1.0/dt, topaz%kfe_org * (cndet(i,j,k) * topaz%mass_2_n + &
      topaz%kfe_bal * (csidet(i,j,k) * 60.0 + ccadet(i,j,k) * 100.0 + clithdet(i,j,k))) * &
      topaz%Rho_0 * topaz%wsink + topaz%kfe_2nd_order * feprime) * feprime
    topaz%jfe_des(i,j,k)=topaz%kfe_des * cfedet(i,j,k)
!
!---------------------------------------------------------------------
! Choose between associating particulate Fe with ballast and organic matter, 
! or just with organic matter
!---------------------------------------------------------------------
!
    if (topaz%fe_ballast_assoc) then  !{
      topaz%jfedet(i,j,k) = (topaz%jndet(i,j,k) * topaz%mass_2_n + topaz%jsidet(i,j,k) *  &
        60.0 + topaz%jcadet(i,j,k) * 100.0) / (cndet(i,j,k) * topaz%mass_2_n +            &
        csidet(i,j,k) * 60.0 + ccadet(i,j,k) * 100.0 + clith(i,j,k) + epsln) *            &
        cfedet(i,j,k)
    else !}{
      topaz%jfedet(i,j,k) = topaz%jndet(i,j,k) / (cndet(i,j,k) + epsln) * cfedet(i,j,k)
    endif !}
!
!---------------------------------------------------------------------
! Calculate iron and lithogenic loss to sediments
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
! Apply sediment flux to all ocean cells adjacent or corner to land
!---------------------------------------------------------------------
!
    topaz%jfe_coast(i,j,k) = topaz%fe_coast * topaz%mask_coast(i,j) *                     &
      grid_tmask(i,j,k) / sqrt(grid_dat(i,j))
  enddo ; enddo ; enddo  !} i,j,k

!
!---------------------------------------------------------------------
! Account for remineralization/dissolution of sinking flux, and
! sediment processed in bottom box
!---------------------------------------------------------------------
!
  do j = jsc, jec ; do i = isc, iec  !{
    k = grid_kmt(i,j)
    if (k .gt. 0) then !{
!
!---------------------------------------------------------------------
! Subtract sedimentary denitrification after Middelburg et al. 1996,
! GBC, 10, 661-673
!---------------------------------------------------------------------
!
      if (cndet(i,j,k) .gt. 0.0 .and. topaz%wsink .gt. 0.0) then !{
! convert flux to umol C cm-2 d-1
        log_btm_flx = log10(cndet(i,j,k) * topaz%Rho_0 * topaz%wsink * topaz%c_2_n *      &
          sperd * 100.0)
        topaz%fdenit_sed(i,j) = min(cndet(i,j,k) * topaz%Rho_0 * topaz%wsink,             &
          10**(-0.9543 + 0.7662 * log_btm_flx - 0.235 * log_btm_flx**2.0) /               &
          (topaz%c_2_n * sperd * 100.0)) * topaz%n_2_n_denit * cno3(i,j,k) /              &
          (phyto(SMALL)%k_no3 + cno3(i,j,k))
      else
        topaz%fdenit_sed(i,j) = 0.0
      endif !}
!
!---------------------------------------------------------------------
! Calculate iron addition from sediments as a function of organic matter
! supply
!---------------------------------------------------------------------
!
      topaz%ffe_sed(i,j) = topaz%fe_2_n_sed * cndet(i,j,k) * topaz%Rho_0 * topaz%wsink
!
!---------------------------------------------------------------------
! Determine the flux of CaCO3 retained in sediment using the metamodel calibrated to the
! Hales (2003) steady state model of CaCO3 burial where spery*100 converts from
! mol m-2 s-1 to g m-2 y-1
!---------------------------------------------------------------------
!
      topaz%fcaco3_sed(i,j) = Tr(ind_fcadet_btm)%field(i,j,1) * min(1.0, max(0.01,        &
        (Tr(ind_fcadet_btm)%field(i,j,1) * spery * 100.0) / (16.0 +                       &
        (Tr(ind_fcadet_btm)%field(i,j,1) * spery * 100.0)) *                              &
        (cco3_ion(i,j,k) / topaz%co3_solubility(i,j,k) * 1.418)**3.96 + 0.027 *           &
        log(max(0.1, clithdet(i,j,k) * topaz%wsink * spery * 100.0)) - 0.072))
!
!---------------------------------------------------------------------
! Allow slow dissolution and ultimate burial of sediment CaCO3 assuming a 10 cm mixed
! layer advecting downward at lithogenic and CaCO3-based sediment accumulation rate
! assuming a density of 2.7 g cm-3 and a porosity of 0.7 - 2.7e6*(1-0.7)=8.1e5.
!---------------------------------------------------------------------
!
      topaz%fcaco3_redis(i,j) = Tr(ind_cased)%field(i,j,1) * max(0.0, 1.0 -               &
        cco3_ion(i,j,k) / topaz%co3_solubility(i,j,k)) * topaz%gamma_cased_dis * 0.1
      topaz%fcaco3_burial(i,j) = Tr(topaz%ind_cased)%field(i,j,1) *                       &
        (clithdet(i,j,k) * topaz%wsink + topaz%fcaco3_sed(i,j) * 100.0) / 8.1e5
      Tr(topaz%ind_cased)%field(i,j,1) = max(0.0,Tr(topaz%ind_cased)%field(i,j,1) +       &
        (topaz%fcaco3_redis(i,j) - topaz%fcaco3_burial(i,j)) * 10.0 * dt)
!
!-----------------------------------------------------------------------
!     Calculate external bottom fluxes for tracer_vertdiff
!-----------------------------------------------------------------------
!
      Tr(ind_alk)%btf(i,j) = - (Tr(ind_fcadet_btm)%field(i,j,1) - topaz%fcaco3_sed(i,j)) -  &
       topaz%fcaco3_redis(i,j) + topaz%fdenit_sed(i,j)
      Tr(ind_dic)%btf(i,j) = - (Tr(ind_fcadet_btm)%field(i,j,1) - topaz%fcaco3_sed(i,j)) -  &
       topaz%fcaco3_redis(i,j) - topaz%fdenit_sed(i,j) / topaz%n_2_n_denit * topaz%c_2_n
      Tr(ind_fed)%btf(i,j) = - topaz%ffe_sed(i,j)
      Tr(ind_ndet)%btf(i,j) = topaz%fdenit_sed(i,j) / topaz%n_2_n_denit
      Tr(ind_no3)%btf(i,j) = topaz%fdenit_sed(i,j)
      Tr(ind_pdet)%btf(i,j) = Tr(ind_ndet)%btf(i,j) * cpdet(i,j,k) / (cndet(i,j,k) + epsln)
      Tr(ind_po4)%btf(i,j) = - Tr(ind_pdet)%btf(i,j)
    endif !}
  enddo ; enddo  !} i, j

!
!-----------------------------------------------------------------------
!
!     CALCULATE SOURCE/SINK TERMS FOR EACH TRACER
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!     Phytoplankton Nitrogen and Phosphorus
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
!
! Diazotrophic Phytoplankton Nitrogen
!
    Tr(ind_ndi)%field(i,j,k) = Tr(ind_ndi)%field(i,j,k) + (phyto(DIAZO)%jprod_n2(i,j,k) + &
      phyto(DIAZO)%jprod_no3(i,j,k) + phyto(DIAZO)%jprod_nh4(i,j,k) -                     &
      phyto(DIAZO)%jgraz_n(i,j,k)) * dt
!
! Large Phytoplankton Nitrogen
!
    Tr(ind_nlg)%field(i,j,k) = Tr(ind_nlg)%field(i,j,k) + (phyto(LARGE)%jprod_no3(i,j,k) +&
      phyto(LARGE)%jprod_nh4(i,j,k) - phyto(LARGE)%jgraz_n(i,j,k)) * dt
!
! Small Phytoplankton Nitrogen
!
    Tr(ind_nsm)%field(i,j,k) = Tr(ind_nsm)%field(i,j,k) + (phyto(SMALL)%jprod_no3(i,j,k) +&
      phyto(SMALL)%jprod_nh4(i,j,k) - phyto(SMALL)%jgraz_n(i,j,k)) * dt
!
! Diazotrophic Phytoplankton Phosphorus
!
    Tr(ind_pdi)%field(i,j,k) = Tr(ind_pdi)%field(i,j,k) +                                 &
      (phyto(DIAZO)%jprod_po4(i,j,k) - phyto(DIAZO)%jgraz_n(i,j,k) *                      &
      phyto(DIAZO)%q_p_2_n(i,j,k)) * dt
!
! Large Phytoplankton Phosphorus
!
    Tr(ind_plg)%field(i,j,k) = Tr(ind_plg)%field(i,j,k) +                                 &
      (phyto(LARGE)%jprod_po4(i,j,k) - phyto(LARGE)%jgraz_n(i,j,k) *                      &
      phyto(LARGE)%q_p_2_n(i,j,k)) * dt
!
! Small Phytoplankton Phosphorus
!
    Tr(ind_psm)%field(i,j,k) = Tr(ind_psm)%field(i,j,k) +                                 &
      (phyto(SMALL)%jprod_po4(i,j,k) - phyto(SMALL)%jgraz_n(i,j,k) *                      &
      phyto(SMALL)%q_p_2_n(i,j,k)) * dt
  enddo ; enddo ; enddo  !} i,j,k
!
!-----------------------------------------------------------------------
!     Phytoplankton Silicon and Iron
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
!
! Large Phytoplankton Silicon
!
    Tr(ind_silg)%field(i,j,k) = Tr(ind_silg)%field(i,j,k) +                               &
     (phyto(LARGE)%jprod_sio4(i,j,k) - phyto(LARGE)%jgraz_sio2(i,j,k)) * dt
! Diazotrophic Phytoplankton Iron
!
    Tr(ind_fedi)%field(i,j,k) = Tr(ind_fedi)%field(i,j,k) +                               &
     (phyto(DIAZO)%jprod_fe(i,j,k) - phyto(DIAZO)%jgraz_fe(i,j,k)) * dt
!
! Large Phytoplankton Iron
!
    Tr(ind_felg)%field(i,j,k) = Tr(ind_felg)%field(i,j,k) +                               &
     (phyto(LARGE)%jprod_fe(i,j,k) - phyto(LARGE)%jgraz_fe(i,j,k)) * dt
!
! Small Phytoplankton Iron
!
    Tr(ind_fesm)%field(i,j,k) = Tr(ind_fesm)%field(i,j,k) +                               &
      (phyto(SMALL)%jprod_fe(i,j,k) - phyto(SMALL)%jgraz_fe(i,j,k)) * dt
  enddo ; enddo ; enddo  !} i,j,k
!
!-----------------------------------------------------------------------
!     Heterotrophic N
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
    Tr(ind_nhet)%field(i,j,k)=Tr(ind_nhet)%field(i,j,k) + (topaz%jprod_nhet(i,j,k) -     &
      topaz%jnhet(i,j,k)) * dt
  enddo ; enddo ; enddo  !} i,j,k
!
!-----------------------------------------------------------------------
!     NO3
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
    topaz%jno3(i,j,k) =  topaz%jnitrif(i,j,k) - phyto(DIAZO)%jprod_no3(i,j,k) -          &
      phyto(LARGE)%jprod_no3(i,j,k) - phyto(SMALL)%jprod_no3(i,j,k) -                    &
      topaz%jdenit_wc(i,j,k)
    Tr(ind_no3)%field(i,j,k) = Tr(ind_no3)%field(i,j,k) + topaz%jno3(i,j,k) * dt
  enddo ; enddo ; enddo  !} i,j,k
!
!-----------------------------------------------------------------------
!     Other nutrients
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
!
! NH4
!
    topaz%jnh4(i,j,k) = topaz%jnh4_graz(i,j,k) + topaz%jnhet(i,j,k) + topaz%gamma_ldon *  &
      cldon(i,j,k) + topaz%gamma_sdon * csdon(i,j,k) + topaz%jndet(i,j,k)-                &
      phyto(DIAZO)%jprod_nh4(i,j,k) - phyto(LARGE)%jprod_nh4(i,j,k) -                     &
      phyto(SMALL)%jprod_nh4(i,j,k) - topaz%jnitrif(i,j,k)
   Tr(ind_nh4)%field(i,j,k) = Tr(ind_nh4)%field(i,j,k) + topaz%jnh4(i,j,k) * dt
!
! PO4
!
    topaz%jpo4(i,j,k) =  topaz%jpo4_graz(i,j,k) + (topaz%jnhet(i,j,k) + topaz%gamma_ldon *&
      cldon(i,j,k)) * topaz%p_2_n_RKR + topaz%gamma_sdop * Tr(ind_sdop)%field(i,j,k) +    &
      topaz%jpdet(i,j,k) - phyto(DIAZO)%jprod_po4(i,j,k) - phyto(LARGE)%jprod_po4(i,j,k) -&
      phyto(SMALL)%jprod_po4(i,j,k)
    Tr(ind_po4)%field(i,j,k) = Tr(ind_po4)%field(i,j,k) + topaz%jpo4(i,j,k) * dt
!
! SiO4
!
    Tr(ind_sio4)%field(i,j,k) = Tr(ind_sio4)%field(i,j,k) + (topaz%jsidet(i,j,k) -        &
      phyto(LARGE)%jprod_sio4(i,j,k) + topaz%jdiss_sio2(i,j,k)) * dt
!
! Fed
!
    Tr(ind_fed)%field(i,j,k) = Tr(ind_fed)%field(i,j,k) + (topaz%jfe_graz(i,j,k) -        &
      phyto(DIAZO)%jprod_Fe(i,j,k) - phyto(LARGE)%jprod_Fe(i,j,k) -         &
      phyto(SMALL)%jprod_Fe(i,j,k) - topaz%jfe_ads(i,j,k) + topaz%jfe_des(i,j,k) +        &
      topaz%jfedet(i,j,k) + topaz%jfe_coast(i,j,k)) * dt
  enddo ; enddo ; enddo  !} i,j,k
!
!-----------------------------------------------------------------------
!     Detrital Components
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
!
! Cadet
!
    Tr(ind_cadet)%field(i,j,k) = Tr(ind_cadet)%field(i,j,k) + (topaz%jprod_cadet(i,j,k) - &
      topaz%jcadet(i,j,k)) * dt
!
! Fedet
!
    Tr(ind_fedet)%field(i,j,k) = Tr(ind_fedet)%field(i,j,k) + (topaz%jprod_fedet(i,j,k) + &
      topaz%jfe_ads(i,j,k) - topaz%jfe_des(i,j,k) - topaz%jfedet(i,j,k)) * dt
!
! Lithdet
!
    Tr(ind_lithdet)%field(i,j,k) = Tr(ind_lithdet)%field(i,j,k) +                         &
      topaz%jprod_lithdet(i,j,k) * dt
!
! Ndet
!
    Tr(ind_ndet)%field(i,j,k) = Tr(ind_ndet)%field(i,j,k) + (topaz%jprod_ndet(i,j,k) -    &
      topaz%jndet(i,j,k)) * dt
!
! Pdet
!
    Tr(ind_pdet)%field(i,j,k) = Tr(ind_pdet)%field(i,j,k) + (topaz%jprod_pdet(i,j,k) -    &
      topaz%jpdet(i,j,k)) * dt
!
! Sidet
!
    Tr(ind_sidet)%field(i,j,k) = Tr(ind_sidet)%field(i,j,k) +                             &
      (phyto(LARGE)%jgraz_sio2(i,j,k) - topaz%jdiss_sio2(i,j,k) - topaz%jsidet(i,j,k)) *  &
      dt
  enddo ; enddo ; enddo  !} i,j,k
!
!-----------------------------------------------------------------------
!     Dissolved Organic Matter
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
!
! Semilabile Dissolved Organic Nitrogen
!
    topaz%jsdon(i,j,k) = topaz%jsdon(i,j,k) -  topaz%gamma_sdon * csdon(i,j,k)
    Tr(ind_sdon)%field(i,j,k) = Tr(ind_sdon)%field(i,j,k) + topaz%jsdon(i,j,k) * dt
!
! Semilabile Dissolved Organic Phosphorus
!
    topaz%jsdop(i,j,k) = topaz%jsdop(i,j,k) - topaz%gamma_sdop * Tr(ind_sdop)%field(i,j,k)
    Tr(ind_sdop)%field(i,j,k) = Tr(ind_sdop)%field(i,j,k) + topaz%jsdop(i,j,k) * dt
!
! Labile Dissolved Organic Nitrogen
!
    topaz%jldon(i,j,k) = topaz%jldon(i,j,k) - topaz%gamma_ldon * cldon(i,j,k)
    Tr(ind_ldon)%field(i,j,k) = Tr(ind_ldon)%field(i,j,k) + topaz%jldon(i,j,k) * dt
  enddo ; enddo ; enddo  !} i,j,k
!
!-----------------------------------------------------------------------
!     O2
!
! O2 production from nitrate, ammonia and nitrogen fixation and
! O2 consumption from production of NH4 from non-sinking particles,
! sinking particles and DOM and O2 consumption from nitrification
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j =jsc, jec ; do i = isc, iec  !{
    topaz%jo2(i,j,k) = (topaz%o2_2_no3 * (phyto(DIAZO)%jprod_no3(i,j,k) +                &
      phyto(LARGE)%jprod_no3(i,j,k) + phyto(SMALL)%jprod_no3(i,j,k)) + topaz%o2_2_nh4 *  &
      (phyto(DIAZO)%jprod_nh4(i,j,k) + phyto(LARGE)%jprod_nh4(i,j,k) +                   &
      phyto(SMALL)%jprod_nh4(i,j,k) + phyto(DIAZO)%jprod_n2(i,j,k))) * grid_tmask(i,j,k)
!
!-----------------------------------------------------------------------
! If O2 is present
!-----------------------------------------------------------------------
!
    if (c_o2(i,j,k) .gt. topaz%o2_min) then  !{
      topaz%jo2(i,j,k) = topaz%jo2(i,j,k) - topaz%o2_2_nh4 * (topaz%jnh4_graz(i,j,k)      &
        + topaz%jnhet(i,j,k) + topaz%jndet(i,j,k) + topaz%gamma_sdon * csdon(i,j,k)       &
        + topaz%gamma_ldon * cldon(i,j,k)) - topaz%o2_2_nitrif * topaz%jnitrif(i,j,k)
    endif  !}
    Tr(ind_o2)%field(i,j,k) = Tr(ind_o2)%field(i,j,k) + topaz%jo2(i,j,k) * dt
  enddo ; enddo ; enddo  !} i,j,k
!
!-----------------------------------------------------------------------
!     The Carbon system
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
!
! Alkalinity
!
! Alkalinity + organic matter + 2 X NH4 is assumed to be conserved except for
! the effects Alk in river runoff and CaCO3 sedimentation and sediment erosion.
! 
    Tr(ind_alk)%field(i,j,k) = Tr(ind_alk)%field(i,j,k) + (2.0 * (topaz%jcadet(i,j,k) -   &
      topaz%jprod_cadet(i,j,k)) - topaz%jno3(i,j,k) + topaz%jnh4(i,j,k) -                 &
      topaz%jdenit_wc(i,j,k) + phyto(DIAZO)%jprod_n2(i,j,k)) * dt
!
! Dissolved Inorganic Carbon
!
    Tr(ind_dic)%field(i,j,k) = Tr(ind_dic)%field(i,j,k) + (topaz%c_2_n *                  &
      (topaz%jno3(i,j,k) + topaz%jnh4(i,j,k) + topaz%jdenit_wc(i,j,k) -                   &
      phyto(DIAZO)%jprod_n2(i,j,k)) + topaz%jcadet(i,j,k) - topaz%jprod_cadet(i,j,k)) * dt
  enddo ; enddo ; enddo !} i,j,k
!
!-----------------------------------------------------------------------
!     Lithogenic aluminosilicate particulates
!-----------------------------------------------------------------------
!
  do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
    Tr(ind_lith)%field(i,j,k) = Tr(ind_lith)%field(i,j,k) - topaz%jprod_lithdet(i,j,k) * dt
  enddo ; enddo ; enddo  !} i,j,k
!
!-----------------------------------------------------------------------
!       Save variables for diagnostics
!-----------------------------------------------------------------------
!

  do n= 1, NUM_PHYTO
  !Beware that some of these might not be legitimate output variables
  ! e.g. The original code did not have phyto(DIAZO)%jprod_nh4
  ! However the calls to ocean_register_diag should have been done properly.
    if (phyto(n)%id_def_fe .gt. 0)          &
      call post_data(phyto(n)%id_def_fe,     phyto(n)%def_fe,             topaz%diag)
    if (phyto(n)%id_felim .gt. 0)           &
      call post_data(phyto(n)%id_felim,      phyto(n)%felim,              topaz%diag)
    if (phyto(n)%id_irrlim .gt. 0)          &
      call post_data(phyto(n)%id_irrlim,     phyto(n)%irrlim,             topaz%diag)
    if (phyto(n)%id_jgraz_n .gt. 0)         &
      call post_data(phyto(n)%id_jgraz_n,    phyto(n)%jgraz_n*rho_dzt,    topaz%diag)
    if (phyto(n)%id_jgraz_fe .gt. 0)        &
      call post_data(phyto(n)%id_jgraz_fe,   phyto(n)%jgraz_fe*rho_dzt,   topaz%diag)
    if (phyto(n)%id_jprod_fe .gt. 0)        &
      call post_data(phyto(n)%id_jprod_fe,   phyto(n)%jprod_fe*rho_dzt,   topaz%diag)
    if (phyto(n)%id_jprod_nh4 .gt. 0)       &
      call post_data(phyto(n)%id_jprod_nh4,  phyto(n)%jprod_nh4*rho_dzt,  topaz%diag)
    if (phyto(n)%id_jprod_no3 .gt. 0)       &
      call post_data(phyto(n)%id_jprod_no3,  phyto(n)%jprod_no3*rho_dzt,  topaz%diag)
    if (phyto(n)%id_jprod_po4 .gt. 0)       &
      call post_data(phyto(n)%id_jprod_po4,  phyto(n)%jprod_po4*rho_dzt,  topaz%diag)
    if (phyto(n)%id_liebig_lim .gt. 0)      &
      call post_data(phyto(n)%id_liebig_lim,phyto(n)%liebig_lim,          topaz%diag)
    if (phyto(n)%id_mu .gt. 0)              &
      call post_data(phyto(n)%id_mu,        phyto(n)%mu,                  topaz%diag)
    if (phyto(n)%id_nh4lim .gt. 0)          &
      call post_data(phyto(n)%id_nh4lim,     phyto(n)%nh4lim,             topaz%diag)
    if (phyto(n)%id_no3lim .gt. 0)          &
      call post_data(phyto(n)%id_no3lim,     phyto(n)%no3lim,             topaz%diag)
    if (phyto(n)%id_po4lim .gt. 0)          &
      call post_data(phyto(n)%id_po4lim,     phyto(n)%po4lim,             topaz%diag)
    if (phyto(n)%id_q_fe_2_n .gt. 0)        &
      call post_data(phyto(n)%id_q_fe_2_n,   phyto(n)%q_fe_2_n,           topaz%diag)
    if (phyto(n)%id_q_p_2_n .gt. 0)         &
      call post_data(phyto(n)%id_q_p_2_n,    phyto(n)%q_p_2_n,            topaz%diag)
    if (phyto(n)%id_q_p_2_n_opt .gt. 0)     &
      call post_data(phyto(n)%id_q_p_2_n_opt, phyto(n)%q_p_2_n_opt,       topaz%diag)
    if (phyto(n)%id_theta .gt. 0)           &
      call post_data(phyto(n)%id_theta,      phyto(n)%theta,              topaz%diag)
  enddo

  if (phyto(DIAZO)%id_jprod_n2 .gt. 0)      &
      call post_data(phyto(DIAZO)%id_jprod_n2,phyto(DIAZO)%jprod_n2*rho_dzt,topaz%diag)
  if (phyto(LARGE)%id_jgraz_sio2 .gt. 0)    &
      call post_data(phyto(LARGE)%id_jgraz_sio2,phyto(LARGE)%jgraz_sio2*rho_dzt,topaz%diag)
  if (phyto(LARGE)%id_jprod_sio4 .gt. 0)    &
      call post_data(phyto(LARGE)%id_jprod_sio4,phyto(LARGE)%jprod_sio4*rho_dzt,topaz%diag)
  if (phyto(LARGE)%id_silim .gt. 0)         &
      call post_data(phyto(LARGE)%id_silim,  phyto(LARGE)%silim,          topaz%diag)

  if (id_sc_co2 .gt. 0)                     &
      call post_data(id_sc_co2,              sc_co2,                      topaz%diag)
  if (id_sc_o2 .gt. 0)                      &
      call post_data(id_sc_o2,               sc_o2,                       topaz%diag)
  if (id_o2_sat .gt. 0)                     &
      call post_data(id_o2_sat,              o2_saturation,               topaz%diag)
  if (topaz%id_alpha .gt. 0)                &
      call post_data(topaz%id_alpha,         Tr(topaz%ind_alpha)%field(:,:,1),topaz%diag)
  if (topaz%id_csurf .gt. 0)                &
      call post_data(topaz%id_csurf,         Tr(topaz%ind_csurf)%field(:,:,1),topaz%diag)
  if (topaz%id_pco2surf .gt. 0)             &
      call post_data(topaz%id_pco2surf,      topaz%pco2surf,              topaz%diag)
  if (topaz%id_sfc_chl .gt. 0)              &
      call post_data(topaz%id_sfc_chl,       Tr(ind_chl)%field(:,:,1),    topaz%diag)
  if (topaz%id_sfc_no3 .gt. 0)              &
      call post_data(topaz%id_sfc_no3,       Tr(ind_no3)%field(:,:,1),    topaz%diag)
  if (topaz%id_htotal .gt. 0)               &
      call post_data(topaz%id_htotal,        Tr(ind_htotal)%field,  topaz%diag)
  if (topaz%id_co3_ion .gt. 0)              &
      call post_data(topaz%id_co3_ion,       Tr(topaz%ind_co3_ion)%field, topaz%diag)
  if (topaz%id_co3_solubility .gt. 0)       &
      call post_data(topaz%id_co3_solubility, topaz%co3_solubility,       topaz%diag)
!
!---------------------------------------------------------------------
! Calculate total carbon  = Dissolved Inorganic Carbon + Phytoplankton Carbon
!   + Dissolved Organic Carbon (including refractory) + Heterotrophic Biomass
!   + Detrital Orgainc and Inorganic Carbon
! For the oceanic carbon budget, a constant 42 uM of dissolved organic
! carbon is added to represent the refractory component.
! For the oceanic nitrogen budget, a constant 2 uM of dissolved organic
! nitrogen is added to represent the refractory component.
!---------------------------------------------------------------------
!
  if (topaz%id_tot_layer_int_c .gt. 0)  &
    call post_data(topaz%id_tot_layer_int_c,(Tr(ind_dic)%field + 4.2e-5 +                 &
        Tr(ind_cadet)%field + topaz%c_2_n * (Tr(ind_ndi)%field + Tr(ind_nlg)%field +      &
        Tr(ind_nsm)%field + Tr(ind_ldon)%field + Tr(ind_sdon)%field + Tr(ind_nhet)%field +&
        Tr(ind_ndet)%field)) * rho_dzt,topaz%diag)
  if (topaz%id_tot_layer_int_fe .gt. 0)  &
    call post_data(topaz%id_tot_layer_int_fe,(Tr(ind_fed)%field + Tr(ind_fedi)%field +    &
      Tr(ind_felg)%field + Tr(ind_fesm)%field + Tr(ind_fedet)%field) * rho_dzt,topaz%diag)
  if (topaz%id_tot_layer_int_n .gt. 0)  &
    call post_data(topaz%id_tot_layer_int_n,(Tr(ind_no3)%field + 2.0e-6 +                 &
      Tr(ind_nh4)%field + Tr(ind_ndi)%field + Tr(ind_nlg)%field + Tr(ind_nsm)%field +     &
      Tr(ind_nsm)%field + Tr(ind_ldon)%field + Tr(ind_sdon)%field + Tr(ind_nhet)%field +  &
      Tr(ind_ndet)%field) * rho_dzt,topaz%diag)
  if (topaz%id_tot_layer_int_p .gt. 0)  &
    call post_data(topaz%id_tot_layer_int_p,(Tr(ind_po4)%field + Tr(ind_pdi)%field +      &
      Tr(ind_plg)%field + Tr(ind_psm)%field + Tr(ind_sdop)%field + Tr(ind_pdet)%field +   &
      topaz%p_2_n_RKR * (Tr(ind_nhet)%field + Tr(ind_ldon)%field)) * rho_dzt,topaz%diag)
  if (topaz%id_tot_layer_int_si .gt. 0)  &
    call post_data(topaz%id_tot_layer_int_si,(Tr(ind_sio4)%field + Tr(ind_silg)%field +   &
      Tr(ind_sidet)%field) * rho_dzt,topaz%diag)
  if (topaz%id_nLg_diatoms .gt. 0)          &
      call post_data(topaz%id_nLg_diatoms,   topaz%nLg_diatoms,           topaz%diag)
  if (topaz%id_jprod_cadet .gt. 0)          &
      call post_data(topaz%id_jprod_cadet,   topaz%jprod_cadet*rho_dzt,   topaz%diag)
  if (topaz%id_jprod_fedet .gt. 0)          &
      call post_data(topaz%id_jprod_fedet,   topaz%jprod_fedet*rho_dzt,   topaz%diag)
  if (topaz%id_jprod_lithdet .gt. 0)        &
      call post_data(topaz%id_jprod_lithdet, topaz%jprod_lithdet*rho_dzt, topaz%diag)
  if (topaz%id_jprod_ndet .gt. 0)           &
      call post_data(topaz%id_jprod_ndet,    topaz%jprod_ndet*rho_dzt,    topaz%diag)
  if (topaz%id_jprod_nhet .gt. 0)           &
      call post_data(topaz%id_jprod_nhet,    topaz%jprod_nhet*rho_dzt,    topaz%diag)
  if (topaz%id_jprod_pdet .gt. 0)           &
      call post_data(topaz%id_jprod_pdet,    topaz%jprod_pdet*rho_dzt,    topaz%diag)
  if (topaz%id_jcadet .gt. 0)               &
      call post_data(topaz%id_jcadet,        topaz%jcadet*rho_dzt,        topaz%diag)
  if (topaz%id_jfe_ads .gt. 0)              &
      call post_data(topaz%id_jfe_ads,       topaz%jfe_ads*rho_dzt,       topaz%diag)
  if (topaz%id_jfe_des .gt. 0)              &
      call post_data(topaz%id_jfe_des,       topaz%jfe_des*rho_dzt,       topaz%diag)
  if (topaz%id_jfe_graz .gt. 0)             &
      call post_data(topaz%id_jfe_graz,      topaz%jfe_graz*rho_dzt,      topaz%diag)
  if (topaz%id_jfe_coast .gt. 0)            &
      call post_data(topaz%id_jfe_coast, topaz%jfe_coast*rho_dzt,         topaz%diag)
  if (topaz%id_jfedet .gt. 0)               &
      call post_data(topaz%id_jfedet,        topaz%jfedet*rho_dzt,        topaz%diag)
  if (topaz%id_jldon .gt. 0)                &
      call post_data(topaz%id_jldon,         topaz%jldon*rho_dzt,         topaz%diag)
  if (topaz%id_jndet .gt. 0)                &
      call post_data(topaz%id_jndet,         topaz%jndet*rho_dzt,         topaz%diag)
  if (topaz%id_jnh4 .gt. 0)                 &
      call post_data(topaz%id_jnh4,          topaz%jnh4*rho_dzt,          topaz%diag)
  if (topaz%id_jnh4_graz .gt. 0)            &
      call post_data(topaz%id_jnh4_graz,     topaz%jnh4_graz*rho_dzt,     topaz%diag)
  if (topaz%id_jnhet .gt. 0)                &
      call post_data(topaz%id_jnhet,         topaz%jnhet*rho_dzt,         topaz%diag)
  if (topaz%id_jnitrif .gt. 0)              &
      call post_data(topaz%id_jnitrif,       topaz%jnitrif*rho_dzt,       topaz%diag)

  if (topaz%id_jno3 .gt. 0)                 &
      call post_data(topaz%id_jno3,          topaz%jno3*rho_dzt,          topaz%diag)
  if (topaz%id_jo2 .gt. 0)                  &
      call post_data(topaz%id_jo2,           topaz%jo2*rho_dzt,           topaz%diag)
  if (topaz%id_jpdet .gt. 0)                &
      call post_data(topaz%id_jpdet,         topaz%jpdet*rho_dzt,          topaz%diag)
  if (topaz%id_jpo4 .gt. 0)                 &
      call post_data(topaz%id_jpo4,          topaz%jpo4*rho_dzt,          topaz%diag)
  if (topaz%id_jpo4_graz .gt. 0)            &
      call post_data(topaz%id_jpo4_graz,     topaz%jpo4_graz*rho_dzt,     topaz%diag)
  if (topaz%id_jdenit_wc .gt. 0)            &
      call post_data(topaz%id_jdenit_wc,     topaz%jdenit_wc*rho_dzt,     topaz%diag)
  if (topaz%id_jdiss_sio2 .gt. 0)           &
      call post_data(topaz%id_jdiss_sio2,    topaz%jdiss_sio2*rho_dzt,    topaz%diag)
  if (topaz%id_jsdon .gt. 0)                &
      call post_data(topaz%id_jsdon,         topaz%jsdon*rho_dzt,         topaz%diag)
  if (topaz%id_jsdop .gt. 0)                &
      call post_data(topaz%id_jsdop,         topaz%jsdop*rho_dzt,         topaz%diag)
  if (topaz%id_jsidet .gt. 0)               &
      call post_data(topaz%id_jsidet,        topaz%jsidet*rho_dzt,        topaz%diag)
  if (topaz%id_fcaco3 .gt. 0)               &
      call post_data(topaz%id_fcaco3,        ccadet*topaz%wsink,          topaz%diag)
  if (topaz%id_fcaco3_burial .gt. 0)        &
      call post_data(topaz%id_fcaco3_burial, topaz%fcaco3_burial,         topaz%diag)
  if (topaz%id_fcaco3_redis .gt. 0)         &
      call post_data(topaz%id_fcaco3_redis,  topaz%fcaco3_redis,          topaz%diag)
  if (topaz%id_fcaco3_sed .gt. 0)           &
      call post_data(topaz%id_fcaco3_sed,    topaz%fcaco3_sed,            topaz%diag)
  if (topaz%id_fdenit_sed .gt. 0)           &
      call post_data(topaz%id_fdenit_sed,    topaz%fdenit_sed,            topaz%diag)
  if (topaz%id_ffe_sed .gt. 0)              &
      call post_data(topaz%id_ffe_sed,       topaz%ffe_sed,               topaz%diag)
  if (topaz%id_flith .gt. 0)                &
    call post_data(topaz%id_flith,         Tr(ind_lithdet)%field * topaz%Rho_0 *        &
      topaz%wsink * grid_tmask,topaz%diag)
  if (topaz%id_fpofe .gt. 0)                &
    call post_data(topaz%id_fpofe,         Tr(ind_fedet)%field * topaz%Rho_0 *          &
      topaz%wsink * grid_tmask,topaz%diag)
  if (topaz%id_fpon .gt. 0)                 &
    call post_data(topaz%id_fpon,          Tr(ind_ndet)%field * topaz%Rho_0 *           &
      topaz%wsink * grid_tmask,topaz%diag)
  if (topaz%id_fpop .gt. 0)                 &
    call post_data(topaz%id_fpop,          Tr(ind_pdet)%field * topaz%Rho_0 *           &
      topaz%wsink * grid_tmask,topaz%diag)
  if (topaz%id_fsio2 .gt. 0)                &
    call post_data(topaz%id_fsio2,         Tr(ind_sidet)%field  *topaz%Rho_0 *          &
      topaz%wsink  *grid_tmask,topaz%diag)

end subroutine  ocean_topaz_source  !}
! </SUBROUTINE> NAME="ocean_topaz_source"

subroutine global_tracer_integrals(G, topaz, Tr, h, dt, context)
type(ocean_grid_type),               intent(in)   :: G
type(TOPAZ_CS),                      intent(in)   :: topaz
type(tracer_type), dimension(:),     intent(in)   :: Tr
real,    dimension(G%isd:,G%jsd:,:), intent(in)   :: h
real,                                intent(in)   :: dt
character(len=*),                    intent(in)   :: context

  integer :: i, j, k, n, is, ie, js, je, nz, ntr
  real,    dimension(G%isd:G%ied,G%jsd:G%jed)   :: temp

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  ntr = NUM_PROG_TRACERS

  do n=1,ntr
    temp(is:ie,js:je) = 0.0
    do k = 1,nz ; do j=js,je ; do i=is,ie
      if (G%hmask(i,j) > 0.5) &
        temp(i,j) = temp(i,j) + Tr(n)%field(i,j,k)*h(i,j,k) * topaz%rho_0*G%DXDYh(i,j)
    enddo ; enddo ; enddo
    do j=js,je ; do i=is,ie
      if (G%hmask(i,j) > 0.5) &
        temp(i,j) = temp(i,j) + Tr(n)%btm_reservoir(i,j) * topaz%rho_0*G%DXDYh(i,j)
    enddo ; enddo
    call hchksum(temp, trim(context)//" Tr%field "//trim(Tr(n)%name), G, 0)
    do j=js,je ; do i=is,ie
      if (G%hmask(i,j) > 0.5) &
        temp(i,j) = Tr(n)%stf(i,j)*G%DXDYh(i,j)*dt
    enddo ; enddo
    call hchksum(temp, trim(context)//" Tr%stf "//trim(Tr(n)%name), G, 0)
    do j=js,je ; do i=is,ie
      if (G%hmask(i,j) > 0.5) &
        temp(i,j) = Tr(n)%btf(i,j)*G%DXDYh(i,j)*dt
    enddo ; enddo
    call hchksum(temp, trim(context)//" Tr%btf "//trim(Tr(n)%name), G, 0)
  enddo

  if(topaz%tracer_debug_verbose) then
! If verbose debugging is activated then this prints out the checksums of the actual tracer field
! as opposed to the total mass/flux of the tracer.
    do n=1,ntr
      call hchksum(Tr(n)%field, trim(context)//" Tr%field raw"//trim(Tr(n)%name), G, 0)
      call hchksum(Tr(n)%stf, trim(context)//" Tr%stf raw"//trim(Tr(n)%name), G, 0)
      call hchksum(Tr(n)%btf, trim(context)//" Tr%btf raw"//trim(Tr(n)%name), G, 0)
    enddo
  endif


end subroutine  global_tracer_integrals   !}

subroutine TOPAZ_coupler_flux_init

  if (topaz_fluxes_initialized ) return

  ind_co2_flux = aof_set_coupler_flux('co2_flux',        &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',      &
       mol_wt = WTMCO2, param = (/ 9.36e-07, 9.7561e-06 /),            &
       ice_restart_file = default_ice_restart_file,                         &
       ocean_restart_file = default_ocean_restart_file,                     &
       caller = "register_TOPAZ")

  ind_o2_flux = aof_set_coupler_flux('o2_flux',          &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',      &
       mol_wt = WTMO2, param = (/ 9.36e-07, 9.7561e-06 /),             &
       ice_restart_file = default_ice_restart_file,                         &
       ocean_restart_file = default_ocean_restart_file,                     &
       caller = "register_TOPAZ")

  ind_runoff_alk = aof_set_coupler_flux('runoff_alk',    &
       flux_type = 'land_sea_runoff', implementation = 'river',        &
       param = (/ 1.0e-03 /),                                          &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")

  ind_runoff_dic = aof_set_coupler_flux('runoff_dic',    &
       flux_type = 'land_sea_runoff', implementation = 'river',        &
       param = (/ 12.011e-03 /),                                       &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")

  ind_runoff_fed = aof_set_coupler_flux('runoff_fed',    &
       flux_type = 'land_sea_runoff', implementation = 'river',        &
       param = (/ 55.847e-03 /),                                       &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")

  ind_runoff_ldon = aof_set_coupler_flux('runoff_ldon',  &
       flux_type = 'land_sea_runoff', implementation = 'river',        &
       param = (/ 14.0067e-03 /),                                      &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")

  ind_runoff_lith = aof_set_coupler_flux('runoff_lith',  &
       flux_type = 'land_sea_runoff', implementation = 'river',        &
       param = (/ 1.0e-03 /),                                          &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")

  ind_runoff_nh4 = aof_set_coupler_flux('runoff_nh4',    &
       flux_type = 'land_sea_runoff', implementation = 'river',        &
       param = (/ 14.0067e-03 /),                                      &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")

  ind_runoff_no3 = aof_set_coupler_flux('runoff_no3',    &
       flux_type = 'land_sea_runoff', implementation = 'river',        &
       param = (/ 14.0067e-03 /),                                      &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")

  ind_dry_dep_fed = aof_set_coupler_flux('dry_dep_fed',  &
       flux_type = 'air_sea_deposition', implementation = 'dry',       &
       param = (/ 55.847e-03 /),                                       &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")

  ind_wet_dep_fed = aof_set_coupler_flux('wet_dep_fed',  &
       flux_type = 'air_sea_deposition', implementation = 'wet',       &
       param = (/ 55.847e-03 /),                                       &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")

  ind_dry_dep_lith = aof_set_coupler_flux('dry_dep_lith',&
       flux_type = 'air_sea_deposition', implementation = 'dry',       &
       param = (/ 1.0e-03 /),                                          &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")

  ind_wet_dep_lith = aof_set_coupler_flux('wet_dep_lith',&
       flux_type = 'air_sea_deposition', implementation = 'wet',       &
       param = (/ 1.0e-03 /),                                          &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")

  ind_dry_dep_nh4 = aof_set_coupler_flux('dry_dep_nh4',  &
       flux_type = 'air_sea_deposition', implementation = 'dry',       &
       param = (/ 14.0067e-03 /),                                      &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")

  ind_wet_dep_nh4 = aof_set_coupler_flux('wet_dep_nh4',  &
       flux_type = 'air_sea_deposition', implementation = 'wet',       &
       param = (/ 14.0067e-03 /),                                      &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")

  ind_dry_dep_no3 = aof_set_coupler_flux('dry_dep_no3',  &
       flux_type = 'air_sea_deposition', implementation = 'dry',       &
       param = (/ 14.0067e-03 /),                                      &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")

  ind_wet_dep_no3 = aof_set_coupler_flux('wet_dep_no3',  &
       flux_type = 'air_sea_deposition', implementation = 'wet',       &
       param = (/ 14.0067e-03 /),                                      &
       ice_restart_file = default_ice_restart_file,                         &
       caller = "register_TOPAZ")
  
  topaz_fluxes_initialized = .true.     

end subroutine TOPAZ_coupler_flux_init
#endif
end module  GOLD_OCEAN_TOPAZ_MOD  !}
