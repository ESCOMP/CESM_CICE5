module ice_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for CICE
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC                 , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                 , only : NUOPC_CompFilterPhaseMap, NUOPC_IsUpdated, NUOPC_IsAtTime
  use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC                 , only : NUOPC_SetAttribute, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC_Model           , only : model_routine_SS           => SetServices
  use NUOPC_Model           , only : model_label_Advance        => label_Advance
  use NUOPC_Model           , only : model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model           , only : model_label_SetRunClock    => label_SetRunClock
  use NUOPC_Model           , only : model_label_Finalize       => label_Finalize
  use NUOPC_Model           , only : NUOPC_ModelGet, SetVM
  use med_constants_mod     , only : R8, CS,CL
  use shr_sys_mod           , only : shr_sys_abort
  use shr_log_mod           , only : shr_log_Unit
  use shr_file_mod          , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_file_mod          , only : shr_file_getloglevel, shr_file_setloglevel
  use shr_file_mod          , only : shr_file_setIO, shr_file_getUnit
  use shr_string_mod        , only : shr_string_listGetNum
  use shr_orb_mod           , only : shr_orb_decl
  use shr_const_mod         , only : shr_const_pi
  use shr_cal_mod           , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use shr_nuopc_scalars_mod , only : flds_scalar_num
  use shr_nuopc_scalars_mod , only : flds_scalar_index_nx
  use shr_nuopc_scalars_mod , only : flds_scalar_index_ny
  use shr_nuopc_scalars_mod , only : flds_scalar_index_nextsw_cday
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_chkerr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_Clock_TimePrint
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_SetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_GetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_Diagnose
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_state_fldDebug
  use shr_nuopc_time_mod    , only : shr_nuopc_time_AlarmInit

  use ice_import_export,      only : ice_import, ice_export
  use ice_import_export,      only : ice_advertise_fields, ice_realize_fields
  use ice_domain_size,        only : nx_global, ny_global
  use ice_domain,             only : nblocks, blocks_ice, distrb_info
  use ice_blocks,             only : block, get_block, nx_block, ny_block, nblocks_x, nblocks_y
  use ice_blocks,             only : nblocks_tot, get_block_parameter
  use ice_distribution,       only : ice_distributiongetblockloc
  use ice_grid,               only : tlon, tlat, hm, tarea, ULON, ULAT
  use ice_constants,          only : rad_to_deg
  use ice_communicate,        only : my_task, master_task, mpi_comm_ice
  use ice_calendar,           only : force_restart_now, write_ic
  use ice_calendar,           only : idate, mday, time, month, daycal, time2sec, year_init
  use ice_calendar,           only : sec, dt, calendar, calendar_type, nextsw_cday, istep
  use ice_orbital,            only : eccen, obliqr, lambm0, mvelpp
  use ice_kinds_mod,          only : dbl_kind, int_kind
  use ice_scam,               only : scmlat, scmlon, single_column
  use ice_fileunits,          only : nu_diag, inst_index, inst_name, inst_suffix, release_all_fileunits
  use ice_ocean,              only : tfrz_option
  use ice_therm_shared,       only : ktherm
  use ice_restart_shared,     only : runid, runtype, restart_dir, restart_file
  use ice_history,            only : accum_hist
  use ice_history_shared,     only : model_doi_url  ! TODO: add this functionality
  use ice_prescribed_mod,     only : ice_prescribed_init
  use ice_atmo,               only : flux_convergence_tolerance, flux_convergence_max_iteration
  use ice_atmo,               only : use_coldair_outbreak_mod
  use CICE_InitMod,           only : CICE_Init
  use CICE_RunMod,            only : CICE_Run
  use perf_mod,               only : t_startf, t_stopf, t_barrierf
  use ice_timers

  implicit none

  public  :: SetServices
  public  :: SetVM

  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelSetRunClock
  private :: ModelFinalize

  integer     , parameter :: dbug = 10
  integer     , parameter :: debug_import = 0 ! internal debug level
  integer     , parameter :: debug_export = 0 ! internal debug level
  character(*), parameter :: modName =  "(ice_comp_nuopc)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

!=======================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    integer :: dbrc
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)

    ! Arguments
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !--------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
         acceptStringList=(/"IPDv01p"/), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local variables
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !--------------------------------

    call ice_advertise_fields(gcomp, importState, exportSTate, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    use shr_nuopc_utils_mod, only: shr_nuopc_set_component_logging
    use shr_nuopc_utils_mod, only: shr_nuopc_get_component_instance

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_DistGrid)     :: distGrid
    type(ESMF_Mesh)         :: Emesh, EmeshTemp
    integer                 :: spatialDim
    integer                 :: numOwnedElements
    real(R8), pointer       :: ownedElemCoords(:)
    real(r8), pointer       :: lat(:), latMesh(:)
    real(r8), pointer       :: lon(:), lonMesh(:)
    integer , allocatable   :: gindex_ice(:)
    integer , allocatable   :: gindex_elim(:)
    integer , allocatable   :: gindex(:)
    integer                 :: globalID
    character(CL)           :: cvalue
    character(ESMF_MAXSTR)  :: convCIM, purpComp
    type(ESMF_VM)           :: vm
    type(ESMF_Time)         :: currTime           ! Current time
    type(ESMF_Time)         :: startTime          ! Start time
    type(ESMF_Time)         :: stopTime           ! Stop time
    type(ESMF_Time)         :: refTime            ! Ref time
    type(ESMF_TimeInterval) :: timeStep           ! Model timestep
    type(ESMF_Calendar)     :: esmf_calendar      ! esmf calendar
    type(ESMF_CalKind_Flag) :: esmf_caltype       ! esmf calendar type
    integer                 :: start_ymd          ! Start date (YYYYMMDD)
    integer                 :: start_tod          ! start time of day (s)
    integer                 :: curr_ymd           ! Current date (YYYYMMDD)
    integer                 :: curr_tod           ! Current time of day (s)
    integer                 :: stop_ymd           ! stop date (YYYYMMDD)
    integer                 :: stop_tod           ! stop time of day (sec)
    integer                 :: ref_ymd            ! Reference date (YYYYMMDD)
    integer                 :: ref_tod            ! reference time of day (s)
    integer                 :: yy,mm,dd           ! Temporaries for time query
    integer                 :: iyear              ! yyyy
    integer                 :: dtime              ! time step
    integer                 :: lmpicom
    integer                 :: shrlogunit         ! original log unit
    integer                 :: shrloglev          ! original log level
    character(len=cs)       :: starttype          ! infodata start type
    integer                 :: lsize              ! local size of coupling array
    character(len=512)      :: diro
    character(len=512)      :: logfile
    logical                 :: isPresent
    integer                 :: localPet
    integer                 :: dbrc
    integer                 :: n,c,g,i,j,m        ! indices
    integer                 :: iblk, jblk         ! indices
    integer                 :: ig, jg             ! indices
    integer                 :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block)             :: this_block         ! block information for current block
    integer                 :: compid             ! component id
    character(len=CL)       :: tempc1,tempc2
    real(R8)                :: diff_lon
    integer                 :: npes
    integer                 :: num_elim_global
    integer                 :: num_elim_local
    integer                 :: num_elim
    integer                 :: num_ice
    integer                 :: num_elim_gcells    ! local number of eliminated gridcells
    integer                 :: num_elim_blocks    ! local number of eliminated blocks
    integer                 :: num_total_blocks
    integer                 :: my_elim_start, my_elim_end
    character(*), parameter     :: F00   = "('(ice_comp_nuopc) ',2a,1x,d21.14)"
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, localPet=localPet, PetCount=npes, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call shr_nuopc_get_component_instance(gcomp, inst_suffix, inst_index)
    inst_name = "ICE"//trim(inst_suffix)

    !----------------------------------------------------------------------------
    ! reset shr logging to my log file
    !----------------------------------------------------------------------------

    ! Note that sets the nu_diag module variable in ice_fileunits
    ! nu_diag in this module is initialized to 0 in the module, and if this reset does not
    ! happen here - then ice_init.F90 will obtain it from the input file ice_modelio.nml

    call shr_nuopc_set_component_logging(gcomp, my_task==master_task, nu_diag, shrlogunit, shrloglev)

    !----------------------------------------------------------------------------
    ! start cice timers
    !----------------------------------------------------------------------------

    call t_startf ('cice_init_total')

    !----------------------------------------------------------------------------
    ! Initialize cice - needed in realize phase to get grid information
    !----------------------------------------------------------------------------

    ! Get orbital values

    call NUOPC_CompAttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) eccen

    call NUOPC_CompAttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) obliqr

    call NUOPC_CompAttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) lambm0

    call NUOPC_CompAttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) mvelpp

    ! Determine start type

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) starttype

    if (trim(starttype) == trim('startup')) then
       runtype = "initial"
    else if (trim(starttype) == trim('continue') ) then
       runtype = "continue"
    else if (trim(starttype) == trim('branch')) then
       runtype = "continue"
    else
       call shr_sys_abort( subname//' ERROR: unknown starttype' )
    end if

    ! Note that in the mct version the atm was initialized first so that nextsw_cday could be passed to the other
    ! components - this assumed that cam or datm was ALWAYS initialized first.
    ! In the nuopc version it will be easier to assume that on startup - nextsw_cday is just what cam was setting
    ! it as the current calendar day
    ! TOOD: need to get the perpetual run working
    ! Set nextsw_cday to -1 (this will skip an orbital calculation on initialization

    if (trim(runtype) /= 'initial') then
       nextsw_cday = -1.0_r8
    else
       call ESMF_ClockGet( clock, currTime=currTime, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Determine single column info

    call NUOPC_CompAttributeGet(gcomp, name='scmlon', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlon

    call NUOPC_CompAttributeGet(gcomp, name='scmlat', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlat

    call NUOPC_CompAttributeGet(gcomp, name='single_column', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) single_column

    ! Determine case name, tfreeze_option, flux convertence before call to cice_init

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) runid

    call NUOPC_CompAttributeGet(gcomp, name="tfreeze_option", value=tfrz_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="flux_convergence", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flux_convergence_tolerance

    call NUOPC_CompAttributeGet(gcomp, name="flux_max_iteration", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flux_convergence_max_iteration

    call NUOPC_CompAttributeGet(gcomp, name="coldair_outbreak_mod", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) use_coldair_outbreak_mod

    ! Get clock information before call to cice_init

    call ESMF_ClockGet( clock, &
         currTime=currTime, startTime=startTime, stopTime=stopTime, refTime=RefTime, &
         timeStep=timeStep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    call ESMF_TimeGet( startTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,start_ymd)

    call ESMF_TimeGet( stopTime, yy=yy, mm=mm, dd=dd, s=stop_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,stop_ymd)

    call ESMF_TimeGet( refTime, yy=yy, mm=mm, dd=dd, s=ref_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,ref_ymd)

    call ESMF_TimeIntervalGet( timeStep, s=dtime, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    dt = real(dtime)

    call ESMF_TimeGet( currTime, calkindflag=esmf_caltype, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar_type = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar_type = shr_cal_gregorian
    else
       call shr_sys_abort( subname//'ERROR:: bad calendar for ESMF' )
    end if

    ! *** Initialize cice ***
    ! Assume that always get atmospheric aerosols to cice (atm_aero no longer needed as flag)-
    ! Note that cice_init also sets time manager info as well as mpi communicator info,
    ! including master_task and my_task

    call t_startf ('cice_init')
    call cice_init( lmpicom )
    call t_stopf ('cice_init')

    ! Now write output to nu_diag - this must happen AFTER call to cice_init
    if (localPet == 0) then
       write(nu_diag,F00) trim(subname),' cice init nextsw_cday = ',nextsw_cday
       write(nu_diag,*) trim(subname),' tfrz_option = ',trim(tfrz_option)
       if (ktherm == 2 .and. trim(tfrz_option) /= 'mushy') then
          write(nu_diag,*) trim(subname),' Warning: Using ktherm = 2 and tfrz_option = ', trim(tfrz_option)
       endif
       write(nu_diag,*) trim(subname),' inst_name   = ',trim(inst_name)
       write(nu_diag,*) trim(subname),' inst_index  = ',inst_index
       write(nu_diag,*) trim(subname),' inst_suffix = ',trim(inst_suffix)
       write(nu_diag,*) trim(subname),' flux_convergence = ', flux_convergence_tolerance
       write(nu_diag,*) trim(subname),' flux_convergence_max_iteration = ', flux_convergence_max_iteration
    endif

    !---------------------------------------------------------------------------
    ! use EClock to reset calendar information on initial start
    !---------------------------------------------------------------------------

    ! - on initial run
    !   - iyear, month and mday obtained from sync clock
    !   - time determined from iyear, month and mday
    !   - istep0 and istep1 are set to 0
    ! - on restart run
    !   - istep0, time and time_forc are read from restart file
    !   - istep1 is set to istep0
    !   - idate is determined from time via the call to calendar (see below)

    if (runtype == 'initial') then
       if (ref_ymd /= start_ymd .or. ref_tod /= start_tod) then
          if (my_task == master_task) then
             write(nu_diag,*) trim(subname),': ref_ymd ',ref_ymd, ' must equal start_ymd ',start_ymd
             write(nu_diag,*) trim(subname),': ref_ymd ',ref_tod, ' must equal start_ymd ',start_tod
          end if
       end if

       if (my_task == master_task) then
          write(nu_diag,*) trim(subname),' idate from sync clock = ', start_ymd
          write(nu_diag,*) trim(subname),'   tod from sync clock = ', start_tod
          write(nu_diag,*) trim(subname),' resetting idate to match sync clock'
       end if
       idate = curr_ymd

       if (idate < 0) then
          if (my_task == master_task) then
             write(nu_diag,*) trim(subname),' ERROR curr_ymd,year_init =',curr_ymd,year_init
             write(nu_diag,*) trim(subname),' ERROR idate lt zero',idate
          end if
          call shr_sys_abort(subname//' :: ERROR idate lt zero')
       endif
       iyear = (idate/10000)                     ! integer year of basedate
       month = (idate-iyear*10000)/100           ! integer month of basedate
       mday  =  idate-iyear*10000-month*100      ! day of month of basedate

       if (my_task == master_task) then
          write(nu_diag,*) trim(subname),' curr_ymd = ',curr_ymd
          write(nu_diag,*) trim(subname),' cice year_init = ',year_init
          write(nu_diag,*) trim(subname),' cice start date = ',idate
          write(nu_diag,*) trim(subname),' cice start ymds = ',iyear,month,mday,start_tod
       endif

       if (calendar_type /= "GREGORIAN") then
          call time2sec(iyear-year_init,month,mday,time)
       else
          call time2sec(iyear-(year_init-1),month,mday,time)
       endif
       time = time+start_tod
    end if

    call calendar(time)     ! update calendar info
    if (write_ic) then
       call accum_hist(dt)  ! write initial conditions
    end if

    !---------------------------------------------------------------------------
    ! Determine the global index space needed for the distgrid
    !---------------------------------------------------------------------------

    ! number the local grid to get allocation size for gindex_ice
    lsize = 0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          do i = ilo, ihi
             lsize = lsize + 1
          enddo
       enddo
    enddo

    ! set global index array
    allocate(gindex_ice(lsize))
    n = 0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
             ig = this_block%i_glob(i)
             jg = this_block%j_glob(j)
             gindex_ice(n) = (jg-1)*nx_global + ig
          enddo
       enddo
    enddo

    ! Determine total number of eliminated blocks globally
    globalID = 0
    num_elim_global = 0  ! number of eliminated blocks
    num_total_blocks = 0
    do jblk=1,nblocks_y
       do iblk=1,nblocks_x
          globalID = globalID + 1
          num_total_blocks = num_total_blocks + 1
          if (distrb_info%blockLocation(globalID) == 0) then
             num_elim_global = num_elim_global + 1
          end if
       end do
    end do

    if (num_elim_global > 0) then

       ! Distribute the eliminated blocks in a round robin fashion amoung processors
       num_elim_local = num_elim_global / npes
       my_elim_start = num_elim_local*localPet + min(localPet, mod(num_elim_global, npes)) + 1
       if (localPet < mod(num_elim_global, npes)) then
          num_elim_local = num_elim_local + 1
       end if
       my_elim_end = my_elim_start + num_elim_local - 1

       ! Determine the number of eliminated gridcells locally
       globalID = 0
       num_elim_blocks = 0  ! local number of eliminated blocks
       num_elim_gcells = 0
       do jblk=1,nblocks_y
          do iblk=1,nblocks_x
             globalID = globalID + 1
             if (distrb_info%blockLocation(globalID) == 0) then
                num_elim_blocks = num_elim_blocks + 1
                if (num_elim_blocks >= my_elim_start .and. num_elim_blocks <= my_elim_end) then
                   this_block = get_block(globalID, globalID)
                   num_elim_gcells = num_elim_gcells + &
                        (this_block%jhi-this_block%jlo+1) * (this_block%ihi-this_block%ilo+1)
                end if
             end if
          end do
       end do

       ! Determine the global index space of the eliminated gridcells
       allocate(gindex_elim(num_elim_gcells))
       globalID = 0
       num_elim_gcells = 0  ! local number of eliminated gridcells
       num_elim_blocks = 0  ! local number of eliminated blocks
       do jblk=1,nblocks_y
          do iblk=1,nblocks_x
             globalID = globalID + 1
             if (distrb_info%blockLocation(globalID) == 0) then
                this_block = get_block(globalID, globalID)
                num_elim_blocks = num_elim_blocks + 1
                if (num_elim_blocks >= my_elim_start .and. num_elim_blocks <= my_elim_end) then
                   do j=this_block%jlo,this_block%jhi
                      do i=this_block%ilo,this_block%ihi
                         num_elim_gcells = num_elim_gcells + 1
                         ig = this_block%i_glob(i)
                         jg = this_block%j_glob(j)
                         gindex_elim(num_elim_gcells) = (jg-1)*nx_global + ig
                      end do
                   end do
                end if
             end if
          end do
       end do

       ! create a global index that includes both active and eliminated gridcells
       num_ice  = size(gindex_ice)
       num_elim = size(gindex_elim)
       allocate(gindex(num_elim + num_ice))
       do n = 1,num_ice
          gindex(n) = gindex_ice(n)
       end do
       do n = num_ice+1,num_ice+num_elim
          gindex(n) = gindex_elim(n-num_ice)
       end do

       deallocate(gindex_elim)

    else

       ! No eliminated land blocks
       num_ice = size(gindex_ice)
       allocate(gindex(num_ice))
       do n = 1,num_ice
          gindex(n) = gindex_ice(n)
       end do

    end if

    !---------------------------------------------------------------------------
    ! Create distGrid from global index array
    !---------------------------------------------------------------------------

    DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------------------------------------------
    ! Create the CICE mesh
    !---------------------------------------------------------------------------

    ! read in the mesh
    call NUOPC_CompAttributeGet(gcomp, name='mesh_ice', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    EMeshTemp = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) then
       write(nu_diag,*)'mesh file for cice domain is ',trim(cvalue)
    end if

    ! recreate the mesh using the above distGrid
    EMesh = ESMF_MeshCreate(EMeshTemp, elementDistgrid=Distgrid, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! obtain mesh lats and lons
    call ESMF_MeshGet(Emesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    allocate(lonMesh(numOwnedElements), latMesh(numOwnedElements))
    call ESMF_MeshGet(Emesh, ownedElemCoords=ownedElemCoords)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1,numOwnedElements
       lonMesh(n) = ownedElemCoords(2*n-1)
       latMesh(n) = ownedElemCoords(2*n)
    end do

    ! obtain internally generated cice lats and lons for error checks
    allocate(lon(lsize))
    allocate(lat(lsize))
    n = 0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
             lon(n) = tlon(i,j,iblk)*rad_to_deg
             lat(n) = tlat(i,j,iblk)*rad_to_deg
          enddo
       enddo
    enddo

    ! error check differences between internally generated lons and those read in
    do n = 1,lsize
       diff_lon = abs(lonMesh(n) - lon(n))
       if ( (diff_lon > 1.e2  .and. abs(diff_lon - 360_r8) > 1.e-1) .or.&
            (diff_lon > 1.e-3 .and. diff_lon < 1._r8) ) then
          !write(6,100)n,lonMesh(n),lon(n), diff_lon
100       format('ERROR: CICE  n, lonmesh(n), lon(n), diff_lon = ',i6,2(f21.13,3x),d21.5)
          !call shr_sys_abort()
       end if
       if (abs(latMesh(n) - lat(n)) > 1.e-1) then
          !write(6,101)n,latMesh(n),lat(n), abs(latMesh(n)-lat(n))
101       format('ERROR: CICE n, latmesh(n), lat(n), diff_lat = ',i6,2(f21.13,3x),d21.5)
          !call shr_sys_abort()
       end if
    end do

    ! deallocate memory
    deallocate(ownedElemCoords)
    deallocate(lon, lonMesh)
    deallocate(lat, latMesh)

    !-----------------------------------------------------------------
    ! Realize the actively coupled fields
    !-----------------------------------------------------------------

    call ice_realize_fields(gcomp, mesh=Emesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------------------------------------------
    ! Prescribed ice initialization - first get compid
    !-----------------------------------------------------------------
    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) compid  ! convert from string to integer

    call ice_prescribed_init(lmpicom, compid, gindex_ice)

    !-----------------------------------------------------------------
    ! Create cice export state
    !-----------------------------------------------------------------

    call ice_export (exportstate, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(nx_global), flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_State_SetScalar(dble(ny_global), flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO (mvertens, 2018-12-21): fill in iceberg_prognostic as .false.

    if (debug_export > 0 .and. my_task==master_task) then
       call shr_nuopc_methods_State_fldDebug(exportState, flds_scalar_name, 'cice_export:', &
            idate, sec, nu_diag, rc=rc)
    end if

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"
    call ESMF_AttributeAdd(comp, convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ShortName", "CICE", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", "CICE Model", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", "CICE5", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "TBD", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Sea Ice",  convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Name", "David Bailey", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "EmailAddress", "dbailey@ucar.edu", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", convention=convCIM, purpose=purpComp, rc=rc)
#endif

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

    !---------------------------------------------------------------------------
    ! Reset shr logging to original values
    !---------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

    call t_stopf ('cice_init_total')

    deallocate(gindex_ice)
    deallocate(gindex)

  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    !---------------------------------------------------------------------------
    ! Run CICE
    !---------------------------------------------------------------------------

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_Clock)           :: clock
    type(ESMF_Alarm)           :: alarm
    type(ESMF_Time)            :: currTime
    type(ESMF_Time)            :: nextTime
    type(ESMF_State)           :: importState, exportState
    character(ESMF_MAXSTR)     :: cvalue
    integer                    :: shrlogunit ! original log unit
    integer                    :: shrloglev  ! original log level
    integer                    :: k,n        ! index
    logical                    :: stop_now   ! .true. ==> stop at the end of this run phase
    integer                    :: ymd        ! Current date (YYYYMMDD)
    integer                    :: tod        ! Current time of day (sec)
    integer                    :: curr_ymd   ! Current date (YYYYMMDD)
    integer                    :: curr_tod   ! Current time of day (s)
    integer                    :: yy,mm,dd   ! year, month, day, time of day
    integer                    :: ymd_sync   ! Sync date (YYYYMMDD)
    integer                    :: yr_sync    ! Sync current year
    integer                    :: mon_sync   ! Sync current month
    integer                    :: day_sync   ! Sync current day
    integer                    :: tod_sync   ! Sync current time of day (sec)
    character(CL)              :: restart_date
    character(CL)              :: restart_filename
    integer                    :: dbrc
    character(*)   , parameter :: F00   = "('(ice_comp_nuopc) ',2a,i8,d21.14)"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !--------------------------------
    ! Turn on timers
    !--------------------------------

    call ice_timer_start(timer_total) ! time entire run
    call t_barrierf('cice_run_total_BARRIER',mpi_comm_ice)
    call t_startf ('cice_run_total')

    !--------------------------------
    ! Reset shr logging to my log file
    !--------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (nu_diag)

    !--------------------------------
    ! Query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 1) then
      call shr_nuopc_methods_Clock_TimePrint(clock,subname//'clock',rc=rc)
    endif

    !--------------------------------
    ! Determine time of next atmospheric shortwave calculation
    !--------------------------------

    call shr_nuopc_methods_State_GetScalar(importState, flds_scalar_index_nextsw_cday, nextsw_cday, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (my_task == master_task) then
       write(nu_diag,F00) trim(subname),' cice istep, nextsw_cday = ',istep, nextsw_cday
    end if

    !--------------------------------
    ! Obtain orbital values
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) eccen

    call NUOPC_CompAttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) obliqr

    call NUOPC_CompAttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) lambm0

    call NUOPC_CompAttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) mvelpp

    !--------------------------------
    ! check that cice internal time is in sync with master clock before timestep update
    !--------------------------------

    ! cice clock
    tod = sec
    ymd = idate

    ! model clock
    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( currTime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr_sync, mon_sync, day_sync, ymd_sync)

    ! error check
    if ( (ymd /= ymd_sync) .and. (tod /= tod_sync) ) then
       if (my_task == master_task) then
          write(nu_diag,*)' cice ymd=',ymd     ,'  cice tod= ',tod
          write(nu_diag,*)' sync ymd=',ymd_sync,'  sync tod= ',tod_sync
       end if
       call ESMF_LogWrite(subname//" CICE clock not in sync with ESMF model clock",ESMF_LOGMSG_ERROR, rc=dbrc)
       rc = ESMF_FAILURE
       return
    end if

    !--------------------------------
    ! Determine if time to write restart
    !--------------------------------

    ! Note this logic triggers off of the component clock rather than the internal cice time
    ! The component clock does not get advanced until the end of the loop - not at the beginning

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       force_restart_now = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockGetNextTime(clock, nextTime=nextTime, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(nexttime, yy=yy, mm=mm, dd=dd, s=tod, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       write(restart_date,"(i4.4,a,i2.2,a,i2.2,a,i5.5)") yy, '-', mm, '-',dd,'-',tod
       write(restart_filename,'(4a)') trim(restart_dir), trim(restart_file), '.', trim(restart_date)
    else
       force_restart_now = .false.
    endif

    !--------------------------------
    ! Unpack import state
    !--------------------------------

    call t_barrierf('cice_run_import_BARRIER',mpi_comm_ice)
    call t_startf ('cice_run_import')
    call ice_timer_start(timer_cplrecv)

    call ice_import(importState, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ice_timer_stop(timer_cplrecv)
    call t_stopf ('cice_run_import')

    ! write Debug output
    if (debug_import  > 0 .and. my_task==master_task) then
       call shr_nuopc_methods_State_fldDebug(importState, flds_scalar_name, 'cice_import:', &
            idate, sec, nu_diag, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !--------------------------------
    ! Advance cice and timestep update
    !--------------------------------

    if (force_restart_now) then
       call CICE_Run(restart_filename=restart_filename)
    else
       call CICE_Run()
    end if

    !--------------------------------
    ! Create export state
    !--------------------------------

    call t_barrierf('cice_run_export_BARRIER',mpi_comm_ice)
    call t_startf ('cice_run_export')
    call ice_timer_start(timer_cplsend)

    call ice_export(exportState, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ice_timer_stop(timer_cplsend)
    call t_stopf ('cice_run_export')

    if (debug_export > 0 .and. my_task==master_task) then
       call shr_nuopc_methods_State_fldDebug(exportState, flds_scalar_name, 'cice_export:', &
            idate, sec, nu_diag, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! reset shr logging to my original values
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

    !--------------------------------
    ! stop timers and print timer info
    !--------------------------------
    ! Need to have this logic here instead of in finalize phase
    ! since the finalize phase will still be called even in aqua-planet mode

    !--------------------------------
    ! Determine if time to stop
    !--------------------------------

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_stop', alarm=alarm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       stop_now = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       stop_now = .false.
    endif

    call t_stopf ('cice_run_total')
    ! Need to stop this at the end of every run phase in a coupled run.
    call ice_timer_stop(timer_total)
    if (stop_now) then
       call ice_timer_print_all(stats=.true.) ! print timing information
       call release_all_fileunits
    endif

  105  format( A, 2i8, A, f10.2, A, f10.2, A)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)
    use shr_nuopc_time_mod, only : shr_nuopc_time_set_component_stop_alarm

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option ! Restart option units
    integer                  :: restart_n      ! Number until restart interval
    integer                  :: restart_ymd    ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    integer                  :: dbrc
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart and stop alarms
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO, rc=dbrc)

       !----------------
       ! Restart alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call shr_nuopc_time_alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------
       ! Stop alarm
       !----------------
       call shr_nuopc_time_set_component_stop_alarm(gcomp, rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(*), parameter :: F00   = "('(ice_comp_nuopc) ',8a)"
    character(*), parameter :: F91   = "('(ice_comp_nuopc) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !--------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    if (my_task == master_task) then
       write(nu_diag,F91)
       write(nu_diag,F00) 'CICE: end of main integration loop'
       write(nu_diag,F91)
    end if

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelFinalize

end module ice_comp_nuopc
