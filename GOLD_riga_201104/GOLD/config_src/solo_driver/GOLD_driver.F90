program GOLD_main
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
!*            The Generalized Ocean Layer Dynamics model               *
!*                               GOLD                                  *
!*                                                                     *
!*  By Robert Hallberg                                                 *
!*                                                                     *
!*    This file is the ocean-only driver for the Generalized Ocean     *
!*  Layered Dynamics (GOLD) ocean model.  A separate ocean interface   *
!*  for use with coupled models is provided in ocean_model_GOLD.F90.   *
!*  These two drivers are kept in separate directories for convenience *
!*  of code selection during compiling.  This file orchestrates the    *
!*  calls to the GOLD initialization routines, to the subroutine that  *
!*  steps the model, and coordinates the output and saving restarts.   *
!*  A description of all of the files that constitute GOLD is found in *
!*  the comments at the beginning of GOLD.F90.  The arguments of each  *
!*  subroutine are described where the subroutine is defined.          *
!*                                                                     *
!*  Macros written all in capital letters are defined in GOLD_memory.h *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

  use GOLD_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
  use GOLD_cpu_clock, only : CLOCK_COMPONENT
  use GOLD_diag_mediator, only : enable_averaging, disable_averaging, diag_mediator_end
  use GOLD, only : initialize_GOLD, step_GOLD, GOLD_control_struct
  use GOLD, only : calculate_surface_state, GOLD_end
  use GOLD_domains, only : GOLD_infra_init, GOLD_infra_end
  use GOLD_error_handler, only : GOLD_error, GOLD_mesg, WARNING, FATAL, is_root_pe
  use GOLD_file_parser, only : read_param, log_param, log_version, param_file_type
  use GOLD_file_parser, only : close_param_file
  use GOLD_grid, only : ocean_grid_type
  use GOLD_io, only : file_exists, open_file, close_file
  use GOLD_io, only : check_nml_error, io_infra_init, io_infra_end
  use GOLD_io, only : APPEND_FILE, ASCII_FILE, READONLY_FILE, SINGLE_FILE
  use GOLD_restart, only : save_restart
  use GOLD_sum_output, only : write_energy, accumulate_net_input
  use GOLD_sum_output, only : GOLD_sum_output_init, sum_output_CS
  use GOLD_surface_forcing, only : set_forcing, average_forcing, forcing_save_restart
  use GOLD_surface_forcing, only : surface_forcing_init, surface_forcing_CS
  use GOLD_time_manager, only : time_type, set_date, set_time, get_date, time_type_to_real
  use GOLD_time_manager, only : operator(+), operator(-), operator(*), operator(/)
  use GOLD_time_manager, only : operator(>), operator(<), operator(>=)
  use GOLD_time_manager, only : increment_date, set_calendar_type, month_name
  use GOLD_time_manager, only : JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use GOLD_variables, only : forcing, surface, directories
  use GOLD_write_cputime, only : write_cputime, GOLD_write_cputime_init
  use GOLD_write_cputime, only : write_cputime_start_clock, write_cputime_CS

  use GOLD_ice_shelf, only : initialize_ice_shelf, ice_shelf_end, ice_shelf_CS
  use GOLD_ice_shelf, only : shelf_calc_flux, ice_shelf_save_restart
! , add_shelf_flux_forcing, add_shelf_flux_IOB
  implicit none

#include <GOLD_memory.h>

  type(forcing)         :: fluxes ! A structure containing pointers to
                                  ! the ocean forcing fields.
  type(surface)         :: state  ! A structure containing pointers to
                                  ! the ocean surface state fields.
  type(ocean_grid_type), pointer :: grid ! A pointer to a structure containing
                                  ! metrics and related information.
  logical :: use_ice_shelf = .false. ! If .true., use the ice shelf model for
                                  ! part of the domain.
  logical :: permit_restart = .true. ! .true. if restart files are to be saved.
  integer :: m, n

  integer :: nmax=2000000000;   ! nmax is the number of iterations
                                ! after which to stop so that the
                                ! simulation does not exceed its CPU
                                ! time limit.  nmax is determined by
                                ! evaluating the CPU time used between
                                ! successive calls to write_energy.
                                ! Initially it is set to be very large.
  type(directories) :: dirs     ! A structure containing several relevant directory paths.

  type(time_type), target :: Time ! A copy of the ocean model's time.
                                ! Other modules can set pointers to this and
                                ! change it to manage diagnostics.
  type(time_type) :: Master_Time  ! The ocean model's master clock. No other
                                ! modules are ever given access to this.
  type(time_type) :: Time1      ! The value of the ocean model's time at the
                                ! start of a call to step_GOLD.
  type(time_type) :: Start_time ! The start time of the simulation.
  type(time_type) :: segment_start_time ! The start time of this run segment.
  type(time_type) :: Time_end     ! End time for the segment or experiment.
  type(time_type) :: write_energy_time ! The next time to write to the energy file.
  type(time_type) :: restart_time ! The next time to write restart files.
  type(time_type) :: Time_step_ocean ! A time_type version of time_step.
  real :: elapsed_time = 0.0    ! Elapsed time in this run in seconds.
  logical :: elapsed_time_master  ! If true, elapsed time is used to set the
                                ! model's master clock (Time).  This is needed
                                ! if Time_step_ocean is not an exact
                                ! representation of time_step.
  real :: time_step             ! The time step of a call to step_GOLD in seconds.
  real :: dt                    ! The baroclinic dynamics time step, in seconds.
  integer :: ntstep             ! The number of baroclinic dynamics time steps
                                ! within time_step.

  integer :: Restart_control    ! An integer that is bit-tested to determine whether
                                ! incremental restart files are saved and whether they
                                ! have a time stamped name.  +1 (bit 0) for generic
                                ! files and +2 (bit 1) for time-stamped files.  A
                                ! restart file is saved at the end of a run segment
                                ! unless Restart_control is negative.
  real :: Time_unit             ! The time unit in seconds for the following input fields.
  type(time_type) :: restint         ! The time between saves of the restart file.
  type(time_type) :: daymax          ! The final day of the simulation.
  type(time_type) :: energysavedays  ! The interval between writing the energies
                                     ! and other integral quantities of the run.

  integer :: date_init(6)=0     ! The start date of the whole simulation.
  integer :: date(6)=-1         ! Possibly the start date of this run segment.
  integer :: years=0, months=0, days=0     ! These may determine the segment run
  integer :: hours=0, minutes=0, seconds=0 ! length, if read from a namelist.
  integer :: yr, mon, day, hr, min, sec  ! Temp variables for writing the date.
  type(param_file_type) :: param_file ! The structure indicating the file(s)
                                ! containing all run-time parameters.
  character(len=9)  :: month
  character(len=16) :: calendar = 'julian'
  integer :: calendar_type=-1

  integer :: unit, io_status, ierr
  logical :: unit_in_use

  integer :: initClock, mainClock, termClock

  type(GOLD_control_struct), pointer :: GOLD_CSp => NULL()
  type(surface_forcing_CS),  pointer :: surface_forcing_CSp => NULL()
  type(sum_output_CS),       pointer :: sum_output_CSp => NULL()
  type(write_cputime_CS),    pointer :: write_CPU_CSp => NULL()
  type(ice_shelf_CS),        pointer :: ice_shelf_CSp => NULL()
  !-----------------------------------------------------------------------

  character(len=4), parameter :: vers_num = 'v2.0'
  character(len=128) :: version = '$Id: GOLD_driver.F90,v 13.0.2.5.2.17 2011/09/19 16:10:35 Robert.Hallberg Exp $'
  character(len=128) :: tagname = '$Name: GOLD_ogrp $'
  character(len=40)  :: mod = "GOLD_main (GOLD_driver)" ! This module's name.

  namelist /ocean_solo_nml/ date_init, calendar, months, days, hours, minutes, seconds

  !#######################################################################

  call write_cputime_start_clock(write_CPU_CSp)

  call GOLD_infra_init() ; call io_infra_init()

  ! These clocks are on the global pelist.
  initClock = cpu_clock_id( 'Initialization' )
  mainClock = cpu_clock_id( 'Main loop' )
  termClock = cpu_clock_id( 'Termination' )
  call cpu_clock_begin(initClock)

  call GOLD_mesg('======== Model being driven by GOLD_driver ========', 2)

  if (file_exists('input.nml')) then
    ! Provide for namelist specification of the run length and calendar data.
    call open_file(unit, 'input.nml', form=ASCII_FILE, action=READONLY_FILE)
    read(unit, ocean_solo_nml, iostat=io_status)
    call close_file(unit)
    if (years+months+days+hours+minutes+seconds > 0) then
      ierr = check_nml_error(io_status,'ocean_solo_nml')
      if (is_root_pe()) write(*,ocean_solo_nml)
    endif
  endif

  ! Read ocean_solo restart, which can override settings from the namelist.
  if (file_exists(trim(dirs%restart_input_dir)//'ocean_solo.res')) then
    call open_file(unit,trim(dirs%restart_input_dir)//'ocean_solo.res', &
                   form=ASCII_FILE,action=READONLY_FILE)
    read(unit,*) calendar_type
    read(unit,*) date_init
    read(unit,*) date
    call close_file(unit)
  else
    if (calendar(1:6) == 'julian') then ;         calendar_type = JULIAN
    else if (calendar(1:6) == 'NOLEAP') then ;    calendar_type = NOLEAP
    else if (calendar(1:10)=='thirty_day') then ; calendar_type = THIRTY_DAY_MONTHS
    else if (calendar(1:11)=='no_calendar') then; calendar_type = NO_CALENDAR
    else if (calendar(1:1) /= ' ') then
      call GOLD_error(FATAL,'GOLD_driver: Invalid namelist value for calendar')
    else
      call GOLD_error(FATAL,'GOLD_driver: No namelist value for calendar')
    endif
  endif
  call set_calendar_type(calendar_type)

  if (sum(date_init) > 0) then
    Start_time = set_date(date_init(1),date_init(2), date_init(3), &
         date_init(4),date_init(5),date_init(6))
  else
    Start_time = set_time(0,0)
  endif

  if (sum(date) >= 0) then
    ! In this case, the segment starts at a time fixed by ocean_solo.res
    segment_start_time = set_date(date(1),date(2),date(3),date(4),date(5),date(6))
    Time = segment_start_time
    call initialize_GOLD(Time, param_file, dirs, GOLD_CSp, segment_start_time)
  else
    ! In this case, the segment starts at a time read from the GOLD restart file
    ! or left as Start_time by GOLD_initialize.
    Time = Start_time
    call initialize_GOLD(Time, param_file, dirs, GOLD_CSp)
  endif
  Master_Time = Time
  grid => GOLD_CSp%grid
  call calculate_surface_state(state, GOLD_CSp%u(:,:,:,1), GOLD_CSp%v(:,:,:,1), &
           GOLD_CSp%h(:,:,:,1), GOLD_CSp%ave_ssh, grid, GOLD_CSp)

  call read_param(param_file, "ICE_SHELF", use_ice_shelf) ! false by default.
  if (use_ice_shelf) then
! Delete this?  call GOLD_error(FATAL,"The ice-shelf code is not yet available in GOLD.")
    call initialize_ice_shelf(Time, ice_shelf_CSp, fluxes)
  endif

  call surface_forcing_init(Time, grid, param_file, GOLD_CSp%diag, &
                            surface_forcing_CSp, GOLD_CSp%tracer_flow_CSp)
  call GOLD_mesg("Done surface_forcing_init.", 5)

  segment_start_time = Time
  elapsed_time = 0.0

  call read_param(param_file,"DT",dt,.true.)
  time_step = dt ; call read_param(param_file,"DT_FORCING",time_step)

  ntstep = MAX(1,ceiling(time_step/dt - 0.001))

  Time_step_ocean = set_time(int(floor(time_step+0.5)))
  elapsed_time_master = (abs(time_step - time_type_to_real(Time_step_ocean)) > 1.0e-12*time_step)
  if (elapsed_time_master) &
    call GOLD_mesg("Using real elapsed time for the master clock.", 2)

  ! Determine the segment end time, either from the namelist file or parsed input file.
  Time_unit = 86400.0 ;  call read_param(param_file,"TIMEUNIT",Time_unit)
  if (years+months+days+hours+minutes+seconds > 0) then
    Time_end = increment_date(Time, years, months, days, hours, minutes, seconds)
    call GOLD_mesg('Segment run length determied from ocean_solo_nml.', 2)
    daymax = Time_end ; call read_param(param_file,"DAYMAX",daymax,Time_unit)
  else
    call read_param(param_file,"DAYMAX",daymax,Time_unit,.true.)
    Time_end = daymax
  endif

  if (Time >= Time_end) call GOLD_error(FATAL, &
    "GOLD_driver: The run has been started at or after the end time of the run.")

  Restart_control = 1 ; restint = Time_end + Time_step_ocean
  call read_param(param_file,"RESTART_CONTROL",Restart_control)
  call read_param(param_file,"RESTINT",restint,Time_unit)
  energysavedays = set_time(0,int(time_step+0.5))
  call read_param(param_file,"ENERGYSAVEDAYS",energysavedays,Time_unit)

  call GOLD_mesg("Call GOLD_sum_output_init.", 5)
  call GOLD_sum_output_init(grid, param_file, dirs%output_directory, &
                           GOLD_CSp%ntrunc, Start_time, sum_output_CSp)
  call GOLD_write_cputime_init(param_file, dirs%output_directory, Start_time, &
                              write_CPU_CSp)
  call GOLD_mesg("Done GOLD_sum_output_init.", 5)

  ! Write all relevant parameters to the model log.
  call log_version(param_file, mod, version, tagname, "")
  call log_param(param_file, mod, "DT_FORCING", time_step, &
                 "The time step for changing forcing, coupling with other \n"//&
                 "components, or potentially writing certain diagnostics. \n"//&
                 "The default value is given by DT.", units="s", default=dt)
  call log_param(param_file, mod, "TIMEUNIT", Time_unit, &
                 "The time unit for DAYMAX, ENERGYSAVEDAYS, and RESTINT.", &
                 units="s", default=86400.0)
  call log_param(param_file, mod, "DAYMAX", daymax, &
                 "The final time of the whole simulation, in units of \n"//&
                 "TIMEUNIT seconds.  This also sets the potential end \n"//&
                 "time of the present run segment if the end time is \n"//&
                 "not set via ocean_solo_nml in input.nml.", timeunit=Time_unit)
  call log_param(param_file, mod, "RESTART_CONTROL", Restart_control, &
                 "An integer whose bits encode which restart files are \n"//&
                 "written. Add 2 (bit 1) for a time-stamped file, and odd \n"//&
                 "(bit 0) for a non-time-stamped file.  A restart file \n"//&
                 "will be saved at the end of the run segment for any \n"//&
                 "non-negative value.", default=1)
  call log_param(param_file, mod, "RESTINT", restint, &
                 "The interval between saves of the restart file in units \n"//&
                 "of TIMEUNIT.  Use a value that is larger than DAYMAX to \n"//&
                 "not save incremental restart files  within a run.  Use \n"//&
                 "0 not to save restart files at all.  The default is to \n"//&
                 "use a larger value than DAYMAX.", timeunit=Time_unit)
  call log_param(param_file, mod, "ENERGYSAVEDAYS", energysavedays, &
                 "The interval in units of TIMEUNIT between saves of the \n"//&
                 "energies of the run and other globally summed diagnostics.", &
                 default=set_time(0,int(time_step+0.5)), timeunit=Time_unit)
  call log_param(param_file, mod, "ICE_SHELF", use_ice_shelf, &
                 "If true, call the code to apply an ice shelf model over \n"//&
                 "some of the domain.", default=.false.)
  call log_param(param_file, mod, "ELAPSED TIME AS MASTER", elapsed_time_master)

  ! Close the param_file.  No further parsing of input is possible after this.
  call close_param_file(param_file)

  ! Write out a time stamp file.
  call open_file(unit, 'time_stamp.out', form=ASCII_FILE, action=APPEND_FILE, &
                 threading=SINGLE_FILE)
  call get_date(Time, date(1), date(2), date(3), date(4), date(5), date(6))
  month = month_name(date(2))
  if (is_root_pe()) write(unit,'(6i4,2x,a3)') date, month(1:3)
  call get_date(Time_end, date(1), date(2), date(3), date(4), date(5), date(6))
  month = month_name(date(2))
  if (is_root_pe()) write(unit,'(6i4,2x,a3)') date, month(1:3)
  call close_file(unit)

  call write_energy(GOLD_CSp%u(:,:,:,1), GOLD_CSp%v(:,:,:,1), GOLD_CSp%h(:,:,:,1), &
                    GOLD_CSp%tv, Time, 0, grid, sum_output_CSp, GOLD_CSp%tracer_flow_CSp)
  call write_cputime(Time, 0, nmax, write_CPU_CSp)

! restart_time is the next integral multiple of restint.
  restart_time = Start_time + restint * &
      (1 + ((Time + Time_step_ocean) - Start_time) / restint)

  write_energy_time = Start_time + energysavedays * &
      (1 + ((Time + Time_step_ocean) - Start_time) / energysavedays)

  if ((.not.BTEST(Restart_control,1)) .and. &
      (.not.BTEST(Restart_control,0))) permit_restart = .false.

  call cpu_clock_end(initClock) !end initialization

  call cpu_clock_begin(mainClock) !begin main loop

  n = 1 ; m = 1
  do while ((n < nmax) .and. (Time < Time_end))

    ! Set the forcing for the next steps.
    call set_forcing(state, fluxes, Time, Time_step_ocean, grid, &
                     surface_forcing_CSp)
    if (use_ice_shelf) then
      call shelf_calc_flux(state, fluxes, Time, time_step, ice_shelf_CSp)
!###IS     call add_shelf_flux_forcing(fluxes, ice_shelf_CSp)
!###IS  ! With a coupled ice/ocean run, use the following call.
!###IS      call add_shelf_flux_IOB(ice_ocean_bdry_type, ice_shelf_CSp)
    endif


    ! This call steps the model over a time time_step.
    Time1 = Master_Time ; Time = Master_Time
    m = step_GOLD(fluxes, state, Time1, time_step, GOLD_CSp)

!    Time = Time + Time_step_ocean
!  This is here to enable fractional-second time steps.
    elapsed_time = elapsed_time + time_step
    if (elapsed_time > 2e9) then
      ! This is here to ensure that the conversion from a real to an integer
      ! can be accurately represented in long runs (longer than ~63 years).
      ! It will also ensure that elapsed time does not loose resolution of order
      ! the timetype's resolution, provided that the timestep and tick are
      ! larger than 10-5 seconds.  If a clock with a finer resolution is used,
      ! a smaller value would be required.
      segment_start_time = segment_start_time + set_time(int(floor(elapsed_time)))
      elapsed_time = elapsed_time - floor(elapsed_time)
    endif
    if (elapsed_time_master) then
      Master_Time = segment_start_time + set_time(int(floor(elapsed_time+0.5)))
    else
      Master_Time = Master_Time + Time_step_ocean
    endif
    Time = Master_Time

    call enable_averaging(time_step,Time,GOLD_CSp%diag)
    call average_forcing(fluxes, time_step, grid, surface_forcing_CSp)
    call accumulate_net_input(fluxes, state, time_step, grid, sum_output_CSp)
    call disable_averaging(GOLD_CSp%diag)

!  See if it is time to write out the energy.
    if (Time + (Time_step_ocean/2) > write_energy_time) then
      call write_energy(GOLD_CSp%u(:,:,:,m), GOLD_CSp%v(:,:,:,m), GOLD_CSp%h(:,:,:,m), &
                        GOLD_CSp%tv, Time, n+ntstep-1, grid, sum_output_CSp, &
                        GOLD_CSp%tracer_flow_CSp)
      call write_cputime(Time, n+ntstep-1, nmax, write_CPU_CSp)
      write_energy_time = write_energy_time + energysavedays
    endif

!  See if it is time to write out a restart file - timestamped or not.
    if (Time + (Time_step_ocean/2) > restart_time) then
      if (permit_restart) then
        if (BTEST(Restart_control,1)) then
          call save_restart(dirs%restart_output_dir, Time, m, grid, &
                            GOLD_CSp%restart_CSp, .true.)
          call forcing_save_restart(surface_forcing_CSp, grid, Time, &
                              dirs%restart_output_dir, .true.)
          if (use_ice_shelf) call ice_shelf_save_restart(ice_shelf_CSp, Time, &
                                      dirs%restart_output_dir, .true.)
        endif
        if (BTEST(Restart_control,0)) then
          call save_restart(dirs%restart_output_dir, Time, m, grid, &
                            GOLD_CSp%restart_CSp)
          call forcing_save_restart(surface_forcing_CSp, grid, Time, &
                              dirs%restart_output_dir)
          if (use_ice_shelf) call ice_shelf_save_restart(ice_shelf_CSp, Time, &
                                      dirs%restart_output_dir)
        endif
      endif
      restart_time = restart_time + restint
    endif

    n = n + ntstep
  enddo

  call cpu_clock_end(mainClock)
  call cpu_clock_begin(termClock)
  if (Restart_control>=0) then
    call save_restart(dirs%restart_output_dir, Time, m, grid, GOLD_CSp%restart_CSp)
    if (use_ice_shelf) call ice_shelf_save_restart(ice_shelf_CSp, Time, &
                                dirs%restart_output_dir)
    ! Write ocean solo restart file.
    call open_file(unit, trim(dirs%restart_output_dir)//'ocean_solo.res', nohdrs=.true.)
    if (is_root_pe())then
        write(unit, '(i6,8x,a)') calendar_type, &
             '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

        call get_date(Start_time, yr, mon, day, hr, min, sec)
        write(unit, '(6i6,8x,a)') yr, mon, day, hr, min, sec, &
             'Model start time:   year, month, day, hour, minute, second'
        call get_date(Time, yr, mon, day, hr, min, sec)
        write(unit, '(6i6,8x,a)') yr, mon, day, hr, min, sec, &
             'Current model time: year, month, day, hour, minute, second'
    end if
    call close_file(unit)
  endif

  if (is_root_pe()) then
    do unit=10,1967
      INQUIRE(unit,OPENED=unit_in_use)
      if (.not.unit_in_use) exit
    enddo
    open(unit,FILE="exitcode",FORM="FORMATTED",STATUS="REPLACE",action="WRITE")
    if (Time < daymax) then
      write(unit,*) 9
    else
      write(unit,*) 0
    endif
    close(unit)
  endif

  call diag_mediator_end(Time)
  call cpu_clock_end(termClock)

  call io_infra_end ; call GOLD_infra_end

  call GOLD_end(GOLD_CSp)
  if (use_ice_shelf) call ice_shelf_end(ice_shelf_CSp)

end program GOLD_main
