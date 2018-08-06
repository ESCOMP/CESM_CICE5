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
  use NUOPC_Model           , only : NUOPC_ModelGet
  use med_constants_mod     , only : IN, R8, I8, CXX, CS,CL
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
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_Meshinit
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_ArrayToState
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_StateToArray
  use shr_nuopc_time_mod    , only : shr_nuopc_time_AlarmInit

  ! TODO: remove these
  use esmFlds               , only : fldListFr, fldListTo, compice, compname, flds_i2o_per_cat
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Realize
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Concat
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Getnumflds
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Getfldinfo

  use ice_cpl_indices
  use ice_import_export,      only : ice_import, ice_export
  use ice_domain_size,        only : nx_global, ny_global, block_size_x, block_size_y, max_blocks
  use ice_domain,             only : nblocks, blocks_ice
  use ice_blocks,             only : block, get_block, nx_block, ny_block
  use ice_grid,               only : tlon, tlat, tarea, hm, grid_type, gridcpl_file, ocn_gridcell_frac
  use ice_constants,          only : rad_to_deg, radius
  use ice_communicate,        only : my_task, master_task, MPI_COMM_ICE
  use ice_calendar,           only : force_restart_now, write_ic
  use ice_calendar,           only : idate, mday, time, month, daycal, time2sec, year_init
  use ice_calendar,           only : sec, dt, calendar, calendar_type, nextsw_cday, istep
  use ice_orbital,            only : eccen, obliqr, lambm0, mvelpp
  use ice_ocean,              only : tfrz_option
  use ice_kinds_mod,          only : dbl_kind
  use ice_scam,               only : scmlat, scmlon, single_column
  use ice_fileunits,          only : nu_diag, ice_stdout, inst_index, inst_name, inst_suffix, release_all_fileunits
  use ice_therm_shared,       only : ktherm
  use ice_restart_shared,     only : runid, runtype, restart_dir, restart_file
  use ice_history,            only : accum_hist
  use ice_history_shared,     only : history_dir, history_file, model_doi_url
  use ice_prescribed_mod,     only : ice_prescribed_init
  use ice_atmo,               only : flux_convergence_tolerance, flux_convergence_max_iteration
  use ice_atmo,               only : use_coldair_outbreak_mod
  use CICE_InitMod,           only : CICE_Init
  use CICE_RunMod,            only : CICE_Run
  use perf_mod,               only : t_startf, t_stopf, t_barrierf
  use ice_timers
  use mct_mod ! TODO: remove this

  implicit none

  public  :: SetServices

  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelSetRunClock
  private :: ModelFinalize

  private :: ice_set_domain
  private :: ice_set_gsmap

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  character(CXX)             :: flds_i2x = ''
  character(CXX)             :: flds_x2i = ''
  real(r8), allocatable      :: x2i(:,:)
  real(r8), allocatable      :: i2x(:,:)
  integer                    :: nflds_i2x
  integer                    :: nflds_x2i
  integer, parameter         :: dbug = 10
  character(*),parameter     :: modName =  "(ice_comp_nuopc)"
  character(*),parameter     :: u_FILE_u = &
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
    integer       :: n,nflds
    logical       :: activefld
    character(CS) :: stdname, shortname
    integer       :: dbrc
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !--------------------------------
    ! create import and export field list
    !--------------------------------

    call shr_nuopc_fldList_Concat(fldListFr(compice), fldListTo(compice), flds_i2x, flds_x2i, flds_scalar_name)

    !--------------------------------
    ! advertise import and export fields
    !--------------------------------

    nflds = shr_nuopc_fldList_Getnumflds(fldListFr(compice))
    do n = 1,nflds
       call shr_nuopc_fldList_Getfldinfo(fldListFr(compice), n, activefld, stdname, shortname)
       if (activefld) then
          call NUOPC_Advertise(exportState, standardName=stdname, shortname=shortname, name=shortname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       call ESMF_LogWrite(subname//':Fr_'//trim(compname(compice))//': '//trim(shortname), ESMF_LOGMSG_INFO)
    end do

    nflds = shr_nuopc_fldList_Getnumflds(fldListTo(compice))
    do n = 1,nflds
       call shr_nuopc_fldList_Getfldinfo(fldListTo(compice), n, activefld, stdname, shortname)
       if (activefld) then
          call NUOPC_Advertise(importState, standardName=stdname, shortname=shortname, name=shortname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       call ESMF_LogWrite(subname//':To_'//trim(compname(compice))//': '//trim(shortname), ESMF_LOGMSG_INFO)
    end do

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local variables
    real(r8), pointer         :: lat(:)
    real(r8), pointer         :: lon(:)
    real(r8), pointer         :: elemCoords(:,:)
    real(r8), pointer         :: elemCornerCoords(:,:,:)
    integer , pointer         :: gindex(:)
    real(r8)                  :: dx,dy
    character(CL)             :: cvalue
    character(ESMF_MAXSTR)    :: convCIM, purpComp
    type(ESMF_VM)             :: vm
    type(ESMF_Grid)           :: Egrid
    type(ESMF_Mesh)           :: Emesh
    type(ESMF_Time)           :: currTime      ! Current time
    type(ESMF_Time)           :: startTime     ! Start time
    type(ESMF_Time)           :: stopTime      ! Stop time
    type(ESMF_Time)           :: refTime       ! Ref time
    type(ESMF_TimeInterval)   :: timeStep      ! Model timestep
    type(ESMF_Calendar)       :: esmf_calendar ! esmf calendar
    type(ESMF_CalKind_Flag)   :: esmf_caltype  ! esmf calendar type
    integer                   :: start_ymd     ! Start date (YYYYMMDD)
    integer                   :: start_tod     ! start time of day (s)
    integer                   :: curr_ymd      ! Current date (YYYYMMDD)
    integer                   :: curr_tod      ! Current time of day (s)
    integer                   :: stop_ymd      ! stop date (YYYYMMDD)
    integer                   :: stop_tod      ! stop time of day (sec)
    integer                   :: ref_ymd       ! Reference date (YYYYMMDD)
    integer                   :: ref_tod       ! reference time of day (s)
    integer                   :: yy,mm,dd      ! Temporaries for time query
    integer                   :: iyear         ! yyyy
    integer                   :: nyrp          ! yyyy
    integer                   :: dtime         ! time step
    integer                   :: lmpicom
    integer                   :: shrlogunit    ! original log unit
    integer                   :: shrloglev     ! original log level
    logical                   :: connected     ! is field connected?
    integer                   :: n,c,g,i,j,m   ! indices
    character(len=cs)         :: starttype     ! infodata start type
    character(len=cl)         :: caseid        ! case ID
    integer                   :: lsize         ! local size of coupling array
    integer                   :: iam,ierr
    character(len=512)        :: diro
    character(len=512)        :: logfile
    logical                   :: isPresent
    integer                   :: localPet
    integer                   :: dbrc
    integer                   :: iblk, ilon, jlat   ! indices
    integer                   :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block)               :: this_block         ! block information for current block
    integer                   :: compid             ! component id
    character(*), parameter   :: F00   = "('(ice_comp_nuopc) ',2a,1x,d21.14)"
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    ! TODO: remove these
    type(mct_gGrid)           :: dom_ice
    type(mct_gsMap)           :: gsmap_ice
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, localPet=localPet, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! determine instance information - first get compid
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) compid  ! convert from string to integer

    call NUOPC_CompAttributeGet(gcomp, name="inst_name", value=inst_name, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="inst_index", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) inst_index

    call ESMF_AttributeGet(gcomp, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", value=inst_suffix, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       inst_suffix = ''
    end if

    !----------------------------------------------------------------------------
    ! reset shr logging to my log file
    !----------------------------------------------------------------------------

    ! Note that sets the nu_diag module variable in ice_fileunits
    ! nu_diag in this module is initialized to 0 in the module, and if this reset does not
    ! happen here - then ice_init.F90 will obtain it from the input file ice_modelio.nml

    if (localPet == 0) then
       call NUOPC_CompAttributeGet(gcomp, name="diro", value=diro, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(gcomp, name="logfile", value=logfile, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       nu_diag = shr_file_getUnit()
       open(nu_diag,file=trim(diro)//"/"//trim(logfile))
    else
       nu_diag = 6
    endif

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (nu_diag)

    !----------------------------------------------------------------------------
    ! start cice timers
    !----------------------------------------------------------------------------

    call t_startf ('cice_init_total')

    !----------------------
    ! Determine field indices of import/export arrays
    !----------------------

    call ice_cpl_indices_set(flds_x2i, flds_i2x, flds_i2o_per_cat)

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

       call shr_sys_flush(nu_diag)
    end if

    call calendar(time)     ! update calendar info
    if (write_ic) then
       call accum_hist(dt)  ! write initial conditions
    end if

    !---------------------------------------------------------------------------
    ! Initialize MCT attribute vectors and indices
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Generate the EMSF mesh
    !---------------------------------------------------------------------------

    ! number the local grid to get allocation size for gindex
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

    allocate(gindex(lsize))
    allocate(lon(lsize))
    allocate(lat(lsize))
    allocate(elemCoords(2,lsize))          ! (lon+lat) * n_gridcells
    allocate(elemCornerCoords(2,4,lsize))  ! (lon+lat) * n_corners * n_gridcells

    ! set gindex, lon and lat
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
             ilon = this_block%i_glob(i)
             jlat = this_block%j_glob(j)
             gindex(n) = (jlat-1)*nx_global + ilon
             lon(n) = TLON(ilon,jlat,iblk)*rad_to_deg
             lat(n) = TLAT(ilon,jlat,iblk)*rad_to_deg
          enddo
       enddo
    enddo

    do n = 1,lsize
       elemCoords(1,n) = lon(n)
       elemCoords(2,n) = lat(n)
       ! TODO: cice does not define corner values and there is no info about grid sizes (ie. dx, dy)
       ! anywhere so make something up for now.  this has to be fixed if weights are generated on the fly!
       ! corners are defined counterclockwise - is this true?
       do m = 1,4
          if (m == 1 .or. m == 4) dx = -0.05
          if (m == 2 .or. m == 3) dx =  0.05
          if (m == 1 .or. m == 2) dy = -0.05
          if (m == 3 .or. m == 4) dy =  0.05
          elemCornerCoords(1,m,n) = lon(n) + dx
          elemCornerCoords(2,m,n) = lat(n) + dy
       end do
    end do

    Emesh = ESMF_MeshCreate(parametricDim=2, &
                            coordSys=ESMF_COORDSYS_SPH_DEG, &
                            elementIds=gindex, &
                            elementType=ESMF_MESHELEMTYPE_QUAD, &
                            elementCoords=elemCoords, &
                            elementCornerCoords=elemCornerCoords, &
                            rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(gindex)
    deallocate(lon)
    deallocate(lat)
    deallocate(elemCoords)
    deallocate(elemCornerCoords)

    !-----------------------------------------------------------------
    ! Realize the actively coupled fields
    !-----------------------------------------------------------------

    call shr_nuopc_fldList_Realize(importState, fldListTo(compice), flds_scalar_name, flds_scalar_num, &
         mesh=Emesh, tag=subname//':ciceImport', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_fldList_Realize(exportState, fldListFr(compice), flds_scalar_name, flds_scalar_num, &
         mesh=Emesh, tag=subname//':ciceExport', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------------------------------------------
    ! Prescribed ice initialization
    !-----------------------------------------------------------------

    call ice_set_gsmap( lmpicom, compid, gsmap_ice )
    call ice_set_domain( lmpicom, gsmap_ice, dom_ice )
    call ice_prescribed_init(compid, gsmap_ice, dom_ice)

    !-----------------------------------------------------------------
    ! Create cice export state
    !-----------------------------------------------------------------

    ! First allocate memory
    nflds_i2x = shr_string_listGetNum(flds_i2x)
    nflds_x2i = shr_string_listGetNum(flds_x2i)
    allocate(i2x(nflds_i2x,lsize))
    allocate(x2i(nflds_x2i,lsize))

    ! Pack export state -  copy from i2x to exportState and  Set the coupling scalars

    call ice_export (i2x)

    call shr_nuopc_grid_ArrayToState(i2x, flds_i2x, exportState, grid_option='mesh', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(nx_global), flds_scalar_index_nx, exportState, mpi_comm_ice, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(ny_global), flds_scalar_index_ny, exportState, mpi_comm_ice, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO: fill in iceberg_prognostic as .false.

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
    call t_barrierf('cice_run_total_BARRIER',MPI_COMM_ICE)
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

    call shr_nuopc_methods_State_GetScalar(importState, flds_scalar_index_nextsw_cday, nextsw_cday, mpi_comm_ice, &
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
    ! Unpack import state
    !--------------------------------

    call t_barrierf('cice_run_import_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_import')
    call ice_timer_start(timer_cplrecv)

    call shr_nuopc_grid_StateToArray(importState, x2i, flds_x2i, grid_option='mesh', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ice_import( x2i )

    call ice_timer_stop(timer_cplrecv)
    call t_stopf ('cice_run_import')

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
    ! Timestep update
    !--------------------------------

    if (force_restart_now) then
       call CICE_Run(restart_filename=restart_filename)
    else
       call CICE_Run()
    end if

    !--------------------------------
    ! Create export state
    !--------------------------------

    call t_barrierf('cice_run_export_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_export')
    call ice_timer_start(timer_cplsend)

    call ice_export(i2x)

    call shr_nuopc_grid_ArrayToState(i2x, flds_i2x, exportState, grid_option='mesh', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ice_timer_stop(timer_cplsend)
    call t_stopf ('cice_run_export')

    ! reset shr logging to my original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

    !--------------------------------
    ! stop timers and print timer info
    !--------------------------------
    ! Need to have this logic here instead of in finalize phase
    ! since the finalize phase will still be called even in aqua-planet mode

    ! Need to stop this at the end of every run phase in a coupled run.
    call ice_timer_stop(timer_total)

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

       call ice_timer_print_all(stats=.true.) ! print timing information
       call release_all_fileunits
    else
       stop_now = .false.
    endif

    call t_stopf ('cice_run_total')

  105  format( A, 2i8, A, f10.2, A, f10.2, A)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)
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
    character(len=256)       :: stop_option    ! Stop option units
    integer                  :: stop_n         ! Number until stop interval
    integer                  :: stop_ymd       ! Stop date (YYYYMMDD)
    type(ESMF_ALARM)         :: stop_alarm
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
       call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_n

       call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_ymd

       call shr_nuopc_time_alarmInit(mclock, stop_alarm, stop_option, &
            opt_n   = stop_n,           &
            opt_ymd = stop_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_stop', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
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

  !===============================================================================

  subroutine ice_set_domain(  mpicom, gsmap_i, dom_i )

    use shr_flds_mod, only : shr_flds_dom_coord, shr_flds_dom_other
    use mct_mod

    ! Arguments
    integer        , intent(in)    :: mpicom 
    type(mct_gsMap), intent(in)    :: gsMap_i
    type(mct_ggrid), intent(inout) :: dom_i

    ! Local Variables
    integer                 :: lsize
    integer                 :: i, j, iblk, n, gi  ! indices
    integer                 :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    real(dbl_kind), pointer :: data(:)            ! temporary
    integer       , pointer :: idata(:)           ! temporary
    type(block)             :: this_block         ! block information for current block
    !--------------------------------

    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (ocean), 0 (non-ocean)

    lsize = mct_gsMap_lsize(gsmap_i, mpicom)

    call mct_gGrid_init( GGrid=dom_i, &
         CoordChars=trim(shr_flds_dom_coord), OtherChars=trim(shr_flds_dom_other), lsize=lsize )
    call mct_aVect_zero(dom_i%data)

    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT

    call mct_gsMap_orderedPoints(gsMap_i, my_task, idata)
    call mct_gGrid_importIAttr(dom_i,'GlobGridNum',idata,lsize)
    deallocate(idata)

    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value

    allocate(data(lsize))
    data(:) = -9999.0_R8
    call mct_gGrid_importRAttr(dom_i,"lat"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_i,"lon"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_i,"area" ,data,lsize)
    call mct_gGrid_importRAttr(dom_i,"aream",data,lsize)
    data(:) = 0.0_R8
    call mct_gGrid_importRAttr(dom_i,"mask",data,lsize)
    call mct_gGrid_importRAttr(dom_i,"frac",data,lsize)

    ! Fill in correct values for domain components

    data(:) = -9999.0_R8
    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
       do i = ilo, ihi
          n = n+1
          data(n) = TLON(i,j,iblk)*rad_to_deg
       enddo    !i
       enddo    !j
    enddo       !iblk
    call mct_gGrid_importRattr(dom_i,"lon",data,lsize)

    data(:) = -9999.0_R8
    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
       do i = ilo, ihi
          n = n+1
          data(n) = TLAT(i,j,iblk)*rad_to_deg
       enddo   !i
       enddo   !j
    enddo      !iblk
    call mct_gGrid_importRattr(dom_i,"lat",data,lsize)

    data(:) = -9999.0_R8
    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
       do i = ilo, ihi
          n = n+1
          data(n) = tarea(i,j,iblk)/(radius*radius)
       enddo   !i
       enddo   !j
    enddo      !iblk
    call mct_gGrid_importRattr(dom_i,"area",data,lsize)

    data(:) = 0.0_R8
    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
       do i = ilo, ihi
          n = n+1
          data(n) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
       enddo   !i
       enddo   !j
    enddo      !iblk
    call mct_gGrid_importRattr(dom_i,"mask",data,lsize)

    data(:) = 0.0_R8
    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
       do i = ilo, ihi
          n = n+1
         if (trim(grid_type) == 'latlon') then
             data(n) = ocn_gridcell_frac(i,j,iblk)
          else
             data(n) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
          end if
       enddo   !i
       enddo   !j
    enddo      !iblk
    call mct_gGrid_importRattr(dom_i,"frac",data,lsize)

    deallocate(data)

  end subroutine ice_set_domain

  !=======================================================================

  subroutine ice_set_gsmap( mpicom, ID, gsMap_ice)

    use mct_mod

    ! Arguments
    integer        , intent(in)    :: mpicom
    integer        , intent(in)    :: ID
    type(mct_gsMap), intent(inout) :: gsMap_ice

    ! Local variables
    integer     :: lat, lon, i, j, iblk, n
    integer     :: lsize,gsize
    integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block) :: this_block         ! block information for current block
    integer,allocatable :: gindex(:)
    !--------------------------------

    ! Build the CICE grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP

    ! number the local grid to get allocation size for gindex
    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
          enddo
       enddo
    enddo
    lsize = n

    ! set gindex
    allocate(gindex(lsize))
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
             lon = this_block%i_glob(i)
             lat = this_block%j_glob(j)
             gindex(n) = (lat-1)*nx_global + lon
          enddo
       enddo
    enddo

    ! now set gsmap once gindex is available
    gsize = nx_global*ny_global
    call mct_gsMap_init( gsMap_ice, gindex, mpicom, ID, lsize, gsize )

    !  deallocate gindex
    deallocate(gindex)

  end subroutine ice_set_gsmap

end module ice_comp_nuopc
