!============= DUMMY VERSION OF atmos_ocean_fluxes_mod =========================

module atmos_ocean_fluxes_mod

implicit none ; private

public :: aof_set_coupler_flux

contains

function aof_set_coupler_flux(name, flux_type, implementation, atm_tr_index, &
     param, flag, ice_restart_file, ocean_restart_file, &
     units, caller)  result (coupler_index)
  character(len=*), intent(in)  :: name, flux_type, implementation
  integer, intent(in), optional :: atm_tr_index
  real, intent(in), dimension(:), optional    :: param
  logical, intent(in), dimension(:), optional :: flag
  character(len=*), intent(in), optional :: &
    ice_restart_file, ocean_restart_file, units, caller
  integer :: coupler_index
  ! None of these arguments are used for anything.

  coupler_index = -1

end function aof_set_coupler_flux

end module atmos_ocean_fluxes_mod
