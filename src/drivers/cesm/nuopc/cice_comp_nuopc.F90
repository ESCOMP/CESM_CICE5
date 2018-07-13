module cice_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for CICE
  !----------------------------------------------------------------------------

  use shr_kind_mod          , only : R8=>SHR_KIND_R8, IN=>SHR_KIND_IN
  use shr_kind_mod          , only : CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, CXX => shr_kind_CXX
  use shr_sys_mod           , only : shr_sys_abort, shr_sys_flush
  use shr_log_mod           , only : shr_log_Unit
  use shr_cal_mod           , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
  use shr_file_mod          , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_file_mod          , only : shr_file_getloglevel, shr_file_setloglevel
  use shr_file_mod          , only : shr_file_setIO, shr_file_getUnit
  use shr_flds_mod          , only : shr_flds_dom_coord, shr_flds_dom_other
  use seq_timemgr_mod       , only : seq_timemgr_EClockGetData, seq_timemgr_EClockDateInSync
  use esmFlds               , only : fldListFr, fldListTo, compice, compname
  use esmFlds               , only : flds_scalar_name, flds_scalar_num
  use esmFlds               , only : flds_scalar_index_nx, flds_scalar_index_ny
  use esmFlds               , only : flds_scalar_index_dead_comps
  use esmFlds               , only : flds_scalar_index_nextsw_cday
  use esmFlds               , only : flds_i2o_per_cat
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Realize
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Concat
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Getnumflds
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Getfldinfo
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_chkerr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_Clock_TimePrint
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_SetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_GetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_Diagnose
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_Meshinit
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_ArrayToState
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_StateToArray

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS           => SetServices,          &
    model_label_Advance        => label_Advance,        &
    model_label_DataInitialize => label_DataInitialize, &
    model_label_SetRunClock    => label_SetRunClock,    &
    model_label_Finalize       => label_Finalize

  use mct_mod
  use ice_cpl_indices
  use ice_import_export,  only : ice_import, ice_export
  use ice_domain_size,    only : nx_global, ny_global, block_size_x, block_size_y, max_blocks
  use ice_domain,         only : nblocks, blocks_ice
  use ice_blocks,         only : block, get_block, nx_block, ny_block
  use ice_grid,           only : tlon, tlat, tarea, hm, grid_type, gridcpl_file, ocn_gridcell_frac
  use ice_constants,      only : rad_to_deg, radius
  use ice_communicate,    only : my_task, master_task, MPI_COMM_ICE
  use ice_calendar,       only : force_restart_now, write_ic
  use ice_calendar,       only : idate, mday, time, month, daycal, time2sec, year_init
  use ice_calendar,       only : sec, dt, calendar, calendar_type, nextsw_cday, istep
  use ice_orbital,        only : eccen, obliqr, lambm0, mvelpp
  use ice_ocean,          only : tfrz_option
  use ice_kinds_mod,      only : dbl_kind
  use ice_scam,           only : scmlat, scmlon, single_column
  use ice_fileunits,      only : nu_diag, ice_stdout, inst_index, inst_name, inst_suffix, release_all_fileunits
  use ice_therm_shared,   only : ktherm
  use ice_restart_shared, only : runid, runtype, restart_dir, restart_file
  use ice_history,        only : accum_hist
  use ice_history_shared, only : history_dir, history_file, model_doi_url
  use ice_prescribed_mod, only : ice_prescribed_init
  use ice_atmo,           only : flux_convergence_tolerance, flux_convergence_max_iteration
  use ice_atmo,           only : use_coldair_outbreak_mod
  use CICE_InitMod,       only : CICE_Init
  use CICE_RunMod,        only : CICE_Run
  use perf_mod,           only : t_startf, t_stopf, t_barrierf
  use ice_timers

  implicit none

  public  :: SetServices

  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelSetRunClock
  private :: ModelFinalize

  private :: ice_domain_mct
  private :: ice_setdef_mct
  private :: ice_coffset_mct
  private :: ice_setcoupling_mct
  private :: ice_SetGSMap_mct

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  character(CXX)             :: flds_i2x = ''
  character(CXX)             :: flds_x2i = ''
  real(r8), allocatable      :: x2i(:,:)
  real(r8), allocatable      :: i2x(:,:)
  integer(IN)                :: nflds_i2x
  integer(IN)                :: nflds_x2i
  character(len=*),parameter :: grid_option = "mesh" ! grid_de, grid_arb, grid_reg, mesh
  integer(IN), parameter     :: dbug = 10
  integer(IN)                :: dbrc
  integer(IN)                :: compid               ! component id

  !----- formats -----
  character(*),parameter :: modName =  "(cice_comp_nuopc)"
  character(*),parameter :: u_FILE_u = __FILE__

  !--- for coupling on other grid from gridcpl_file ---
  type(mct_gsMap) :: gsMap_iloc    ! local gsmaps
  type(mct_gGrid) :: dom_iloc      ! local domain
  type(mct_aVect) :: x2i_ice , i2x_ice
  type(mct_aVect) :: x2i_iloc, i2x_iloc
  type(mct_rearr) :: rearr_ice2iloc
  type(mct_rearr) :: rearr_iloc2ice
  integer         :: nxcpl, nycpl  ! size of coupling grid
  logical         :: other_cplgrid ! using different coupling grid

!=======================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
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
    integer                   :: ncols, ngcols ! number of local and global columsn
    integer                   :: n,c,g,i,j,m   ! indices
    character(len=cs)         :: starttype     ! infodata start type
    character(len=cl)         :: caseid        ! case ID
    integer                   :: lsize         ! local size of coupling array
    integer                   :: lsize_loc
    integer                   :: xoff,yoff
    integer                   :: nxg,nyg
    integer                   :: iam,ierr
    integer         , pointer :: gindex(:)
    type(mct_gsMap)           :: gsmap_extend  ! local gsmaps
    type(mct_gGrid)           :: dom_ice
    type(mct_gsMap)           :: gsMap_ice
    character(len=512)        :: diro
    character(len=512)        :: logfile
    logical                   :: isPresent
    integer                   :: localPet
    character(*), parameter   :: F00   = "('(cice_comp_nuopc) ',2a,1x,d21.14)"
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
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

    ! Initialize ice gsMap

    if (trim(gridcpl_file) == 'unknown_gridcpl_file') then

       ! This is the normal case where the ice coupling grid and model grid are the same

       ! First get the cice gsmap (gsmap_ice)
       call ice_SetGSMap_mct( mpi_comm_ice, compid, gsmap_ice )

       ! Get the local size of the gsmap (lsize)
       lsize = mct_gsMap_lsize(gsMap_ice, mpi_comm_ice)

       ! Create the domain (dom_ice)
       call ice_domain_mct( lsize, gsMap_ice, dom_ice )

       other_cplgrid = .false.
       nxg = nx_global
       nyg = ny_global

    else

       ! This is the case where the input file for cice coupling grid info if coupling grid is different than then
       ! cice computation grid

       ! First get the cice computational grid and associated gsmap (gsmap_iloc)
       call ice_SetGSMap_mct( mpi_comm_ice, compid, gsmap_iloc )

       ! Get the size of the computational grid gsmap (lsize_loc)
       lsize_loc = mct_gsMap_lsize(gsMap_iloc, MPI_COMM_ICE)

       ! Create the domain computational grid (dom_iloc)
       call ice_domain_mct( lsize_loc, gsMap_iloc, dom_iloc )

       ! Create the gsmap of the coupling grid (gsmap_ice) and the domain of the coupling grid (dom_i)
       call ice_setcoupling_mct(MPI_COMM_ICE, compid, gsmap_ice, dom_ice)

       ! Get the size of the coupling grid gsmap
       lsize = mct_gsMap_lsize(gsMap_ice, mpi_comm_ice)

       ! Create offsets: xoff and yoff, given the
       call ice_coffset_mct(xoff, yoff, gsmap_iloc, dom_iloc, gsmap_ice, dom_ice, MPI_COMM_ICE)

       ! Now get the cice coupling grid and its global index space (gsmap_extend)
       call ice_SetGSMap_mct( mpi_comm_ice, compid, gsmap_extend, xoff=xoff, yoff=yoff, nxgin=nxcpl, nygin=nycpl)

       ! Verify that the extended gsmap local size is the same as computational grid local size
       if (lsize_loc /= mct_gsmap_lsize(gsmap_extend, MPI_COMM_ICE)) then
          if (my_task == master_task) then
             write(nu_diag,*) subname,' :: gsmap_extend extended ',lsize_loc, mct_gsmap_lsize(gsmap_extend,MPI_COMM_ICE)
          end if
          call shr_sys_abort(subname//' :: error in gsmap_extend extended')
       endif

       ! Create rearrangers from the computational grid to the coupling grid and vice versa
       call mct_rearr_init(gsmap_ice, gsmap_extend, MPI_COMM_ICE, rearr_ice2iloc)
       call mct_rearr_init(gsmap_extend, gsmap_ice, MPI_COMM_ICE, rearr_iloc2ice)

       ! Initialize attribute vectors for the computational grid
       call mct_aVect_init(x2i_iloc, rList=flds_x2i, lsize=lsize_loc)
       call mct_aVect_zero(x2i_iloc)
       call mct_aVect_init(i2x_iloc, rList=flds_i2x, lsize=lsize_loc)
       call mct_aVect_zero(i2x_iloc)

       ! Clean the
       call mct_gsmap_clean(gsmap_extend)

       other_cplgrid = .true.
       nxg = nxcpl
       nyg = nycpl
    endif

    !---------------------------------------------------------------------------
    ! Generate the EMSF mesh
    !---------------------------------------------------------------------------

    allocate(gindex(lsize))
    allocate(lon(lsize))
    allocate(lat(lsize))
    allocate(elemCoords(2,lsize))          ! (lon+lat) * n_gridcells
    allocate(elemCornerCoords(2,4,lsize))  ! (lon+lat) * n_corners * n_gridcells

    call mct_gGrid_exportRattr(dom_ice, 'lon', lon, lsize)
    call mct_gGrid_exportRattr(dom_ice, 'lat', lat, lsize)
    call mct_gsMap_OrderedPoints(gsMap_ice, my_task, gindex)

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

    if (other_cplgrid) then
       call ice_prescribed_init(compid, gsmap_iloc, dom_iloc)
    else
       call ice_prescribed_init(compid, gsmap_ice, dom_ice)
    endif

    !-----------------------------------------------------------------
    ! Create cice export state
    !-----------------------------------------------------------------

    ! First inialize mct attribute vectors for the coupling grid

    call mct_aVect_init(x2i_ice, rList=flds_x2i, lsize=lsize)
    call mct_aVect_init(i2x_ice, rList=flds_i2x, lsize=lsize)
    call mct_aVect_zero(x2i_ice)
    call mct_aVect_zero(i2x_ice)

    if (other_cplgrid) then
       call ice_export (i2x_iloc%rattr)
       call ice_setdef_mct ( i2x_ice )
       call mct_rearr_rearrange(i2x_iloc, i2x_ice, rearr_iloc2ice)
    else
       call ice_export (i2x_ice%rattr)
    endif

    ! Pack export state -  copy from i2x to exportState and  Set the coupling scalars

    call shr_nuopc_grid_ArrayToState(i2x_ice%rattr, flds_i2x, exportState, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(nxg), flds_scalar_index_nx, exportState, mpi_comm_ice, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(nyg), flds_scalar_index_ny, exportState, mpi_comm_ice, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(0.0_r8, flds_scalar_index_dead_comps, exportState, mpi_comm_ice, &
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
    type(ESMF_Clock) :: clock
    type(ESMF_Alarm) :: alarm
    type(ESMF_Time)  :: time
    type(ESMF_State) :: importState, exportState
    character(CL)    :: cvalue
    integer          :: shrlogunit ! original log unit
    integer          :: shrloglev  ! original log level
    integer          :: k,n        ! index
    logical          :: stop_now   ! .true. ==> stop at the end of this run phase
    integer          :: ymd        ! Current date (YYYYMMDD)
    integer          :: tod        ! Current time of day (sec)
    integer          :: tod_sync   ! Sync current time of day (sec)
    integer          :: ymd_sync   ! Current year of sync clock
    integer          :: curr_ymd   ! Current date (YYYYMMDD)
    integer          :: curr_tod   ! Current time of day (s)
    integer          :: yy,mm,dd   ! year, month, day, time of day
    character(CL)    :: restart_date
    character(CL)    :: restart_filename
    character(*)   , parameter  :: F00   = "('(cice_comp_nuopc) ',2a,i8,d21.14)"
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

    call shr_nuopc_grid_StateToArray(importState, x2i_ice%rattr, flds_x2i, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (other_cplgrid) then
       call mct_rearr_rearrange(x2i_ice, x2i_iloc, rearr_ice2iloc)
       call ice_import( x2i_iloc%rattr )
    else
       call ice_import( x2i_ice%rattr )
    endif

    call ice_timer_stop(timer_cplrecv)
    call t_stopf ('cice_run_import')

    !--------------------------------
    ! check that cice internal time is in sync with master clock before timestep update
    !--------------------------------

    tod = sec
    ymd = idate
    if (.not. seq_timemgr_EClockDateInSync( clock, ymd, tod )) then
       call seq_timemgr_EClockGetData( clock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       if (my_task == master_task) then
          write(nu_diag,*)' cice ymd=',ymd     ,'  cice tod= ',tod
          write(nu_diag,*)' sync ymd=',ymd_sync,'  sync tod= ',tod_sync
       end if
       call shr_sys_abort( SubName// ":: Internal sea-ice clock not in sync with Sync Clock")
    end if

    !--------------------------------
    ! Determine if time to write restart
    !--------------------------------

    ! Note this logic triggers off of the component clock rather than the internal cice time
    ! The component clock does not get advanced until the end of the loop - not at the beginning

    call ESMF_ClockGetAlarm(clock, alarmname='seq_timemgr_alarm_restart', alarm=alarm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       force_restart_now = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call seq_timemgr_EClockGetData( clock, curr_yr=yy, curr_mon=mm, curr_day=dd, curr_tod=tod)
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

    if (other_cplgrid) then
       call ice_export ( i2x_iloc%rattr )
       call ice_setdef_mct ( i2x_ice )
       call mct_rearr_rearrange(i2x_iloc, i2x_ice, rearr_iloc2ice)
    else
       call ice_export ( i2x_ice%rattr )
    endif

    call shr_nuopc_grid_ArrayToState(i2x_ice%rattr, flds_i2x, exportState, grid_option, rc=rc)
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

    call ESMF_ClockGetAlarm(clock, alarmname='seq_timemgr_alarm_stop', alarm=alarm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       stop_now = .true.
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
    character(len=128)       :: mtimestring, dtimestring
    type(ESMF_Alarm),pointer :: alarmList(:)
    type(ESMF_Alarm)         :: dalarm
    integer                  :: alarmcount, n
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! check that the current time in the model and driver are the same
    !--------------------------------

    if (mcurrtime /= dcurrtime) then
      call ESMF_TimeGet(dcurrtime, timeString=dtimestring, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_TimeGet(mcurrtime, timeString=mtimestring, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_LogWrite(subname//" ERROR in time consistency; "//trim(dtimestring)//" ne "//trim(mtimestring),  &
           ESMF_LOGMSG_ERROR, rc=dbrc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      rc=ESMF_Failure
    endif

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! copy alarms from driver to model clock if model clock has no alarms (do this only once!)
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then
      call ESMF_ClockGetAlarmList(dclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      allocate(alarmList(alarmCount))
      call ESMF_ClockGetAlarmList(dclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmList=alarmList, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      do n = 1, alarmCount
         ! call ESMF_AlarmPrint(alarmList(n), rc=rc)
         ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         dalarm = ESMF_AlarmCreate(alarmList(n), rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         call ESMF_AlarmSet(dalarm, clock=mclock, rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      enddo

      deallocate(alarmList)
    endif

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
    character(*), parameter :: F00   = "('(cice_comp_nuopc) ',8a)"
    character(*), parameter :: F91   = "('(cice_comp_nuopc) ',73('-'))"
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

  subroutine ice_domain_mct( lsize, gsMap_i, dom_i )

    ! Arguments
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_i
    type(mct_ggrid), intent(inout) :: dom_i

    ! Local Variables
    integer                 :: i, j, iblk, n, gi  ! indices
    integer                 :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    real(dbl_kind), pointer :: work_dom(:)        ! temporary
    real(dbl_kind), pointer :: data(:)            ! temporary
    integer       , pointer :: idata(:)           ! temporary
    type(block)             :: this_block         ! block information for current block
    !--------------------------------

    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (ocean), 0 (non-ocean)

    call mct_gGrid_init( GGrid=dom_i, &
         CoordChars=trim(shr_flds_dom_coord), OtherChars=trim(shr_flds_dom_other), lsize=lsize )
    call mct_aVect_zero(dom_i%data)

    allocate(data(lsize))

    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT

    call mct_gsMap_orderedPoints(gsMap_i, my_task, idata)
    call mct_gGrid_importIAttr(dom_i,'GlobGridNum',idata,lsize)

    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value

    data(:) = -9999.0_R8
    call mct_gGrid_importRAttr(dom_i,"lat"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_i,"lon"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_i,"area" ,data,lsize)
    call mct_gGrid_importRAttr(dom_i,"aream",data,lsize)
    data(:) = 0.0_R8
    call mct_gGrid_importRAttr(dom_i,"mask",data,lsize)
    call mct_gGrid_importRAttr(dom_i,"frac",data,lsize)

    ! Fill in correct values for domain components

    allocate(work_dom(lsize))
    work_dom(:) = 0.0_dbl_kind

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
    deallocate(idata)
    deallocate(work_dom)

  end subroutine ice_domain_mct

  !=======================================================================

  subroutine ice_setdef_mct( i2x_i )
    type(mct_aVect)   , intent(inout) :: i2x_i

    call mct_aVect_zero(i2x_i)
    ! tcraig : this is where observations could be read in

  end subroutine ice_setdef_mct

  !=======================================================================

  subroutine ice_coffset_mct(xoff,yoff,gsmap_a,dom_a,gsmap_b,dom_b,mpicom_i)

    ! Arguments
    integer        , intent(out)   :: xoff
    integer        , intent(out)   :: yoff
    type(mct_gsmap), intent(in)    :: gsmap_a
    type(mct_ggrid), intent(in)    :: dom_a
    type(mct_gsmap), intent(in)    :: gsmap_b
    type(mct_ggrid), intent(in)    :: dom_b
    integer        , intent(in)    :: mpicom_i

    ! Local variables
    type(mct_aVect)  :: ava
    type(mct_aVect)  :: avag
    integer          :: k1,k2,k
    integer          :: npt
    integer          :: noff,noffg
    real(dbl_kind)   :: x1,y1,x2,y2
    real(dbl_kind)   :: dist,distmin,distming
    integer          :: lsizea,lsizeb
    integer          :: iam,ierr
    integer, pointer :: ipoints(:)
    character(len=*),parameter :: subname = "ice_coffset_mct"
    !--------------------------------

    call mpi_comm_rank(mpicom_i,iam,ierr)

    lsizea = mct_aVect_lsize(dom_a%data)
    lsizeb = mct_aVect_lsize(dom_b%data)

    !--- compute lon/lat at dom_a (local) point (1,1)

    call mct_aVect_init(ava,rList='lon:lat',lsize=lsizea)
    call mct_aVect_copy(dom_a%data,ava,'lon:lat')
    call mct_aVect_gather(ava,avag,gsmap_a,0,mpicom_i)

    if (iam == 0) then
       k1 = mct_aVect_indexRA(avag,'lon',dieWith=subname//'_avag')
       k2 = mct_aVect_indexRA(avag,'lat',dieWith=subname//'_avag')
       npt = 1   ! actual corner points screwed up by U average/wraparound
       npt = nx_global + 2  ! use global point (2,2)
       x1 = mod(avag%rAttr(k1,npt)+360.0_r8,360.0_r8)
       y1 = avag%rAttr(k2,npt)
    endif

    call mct_aVect_clean(avag)
    call mct_aVect_clean(ava)

    call shr_mpi_bcast(x1,mpicom_i)
    call shr_mpi_bcast(y1,mpicom_i)

    !--- find x1,y1 point in dom_b (extended grid)

    noff = -1
    noffg = -1

    call mct_gsMap_orderedPoints(gsMap_b, iam, ipoints)
    if (size(ipoints) /= lsizeb) then
       if (my_task == master_task) then
          write(nu_diag,*) subname,' size ipoints = ',size(ipoints),lsizeb
       end if
       call shr_sys_abort(subname//' :: error size of ipoints')
    endif

    k1 = mct_aVect_indexRA(dom_b%data,'lon',dieWith=subname//'_domb')
    k2 = mct_aVect_indexRA(dom_b%data,'lat',dieWith=subname//'_domb')
    distmin = 1.0e36
    do k = 1,lsizeb
       x2 = mod(dom_b%data%rAttr(k1,k)+360.0_r8,360.0_r8)
       y2 = dom_b%data%rAttr(k2,k)
       dist = abs((x1-x2)*(x1-x2))+abs((y1-y2)*(y1-y2))
       if (dist < distmin) then
          distmin = dist
          noff = ipoints(k)
       endif
       dist = abs((x1-x2-360.0_r8)*(x1-x2-360.0_r8))+abs((y1-y2)*(y1-y2))
       if (dist < distmin) then
          distmin = dist
          noff = ipoints(k)
       endif
       dist = abs((x1-x2+360.0_r8)*(x1-x2+360.0_r8))+abs((y1-y2)*(y1-y2))
       if (dist < distmin) then
          distmin = dist
          noff = ipoints(k)
       endif
    enddo

    deallocate(ipoints)

    call shr_mpi_min(distmin,distming,mpicom_i,'distmin',all=.true.)

    if (distming /= distmin) then
       noff = -1
    endif

    call shr_mpi_max(noff,noffg,mpicom_i,'noffg',all=.true.)

    ! subtract extra -1 and -nxcpl for point (2,2)
    xoff = mod(noffg-1-1,nxcpl) + 1
    yoff = (noffg-1-nxcpl)/nxcpl + 1

    if (iam == 0) then
       write(nu_diag,*) subname,' :: x1,y1  = ',x1,y1
       write(nu_diag,*) subname,' :: offset = ',noffg,xoff,yoff
       call shr_sys_flush(nu_diag)
    endif

    if (noffg < 1) then
       call shr_sys_abort(subname//' :: noffg lt 1')
    endif

  end subroutine ice_coffset_mct

  !=======================================================================

  subroutine ice_setcoupling_mct(mpicom_i, ICEID, gsmap_i, dom_i)

    include 'netcdf.inc'

    ! Arguments
    integer        , intent(in)    :: mpicom_i
    integer        , intent(in)    :: ICEID
    type(mct_gsmap), intent(inout) :: gsmap_i
    type(mct_ggrid), intent(inout) :: dom_i

    ! Local variables
    integer                :: n     ! counter
    integer                :: iam   ! pe rank
    integer                :: npes  ! number of pes
    integer                :: ierr  ! error code
    integer                :: rcode ! error code
    integer                :: nx,ny ! grid size
    integer                :: gsize ! global size
    integer                :: lsize ! local size
    integer, pointer       :: start(:),length(:),pe_loc(:)
    integer, pointer       :: idata(:)
    real(dbl_kind),pointer :: data(:)
    type(mct_avect)        :: avg, av1
    integer                :: fid,did,vid
    character(len=8)       :: avfld,dofld
    character(len=*), parameter  :: SubName = "ice_setcoupling_mct"
    !--------------------------------

    call MPI_comm_rank(mpicom_i,iam,ierr)
    call MPI_comm_size(mpicom_i,npes,ierr)

    allocate(start(npes),length(npes),pe_loc(npes))

    if (iam == 0) then
       rcode = nf_open(gridcpl_file(1:len_trim(gridcpl_file)),NF_NOWRITE,fid)
       rcode = nf_inq_dimid (fid, 'ni', did)
       rcode = nf_inq_dimlen(fid, did, nx)
       rcode = nf_inq_dimid (fid, 'nj', did)
       rcode = nf_inq_dimlen(fid, did, ny)
       gsize = nx*ny
       nxcpl = nx
       nycpl = ny

       length = gsize / npes
       do n = 1,npes
          if (n <= mod(gsize,npes)) length(n) = length(n) + 1
       enddo

       start(1) = 1
       pe_loc(1) = 0
       do n = 2,npes
          pe_loc(n) = n-1
          start(n) = start(n-1) + length(n-1)
       enddo
       if ((start(npes) + length(npes) - 1) /= gsize) then
          write(nu_diag,*) &
            subname,' gsize, start, length = ',gsize,start(npes),length(npes)
          call shr_sys_flush(nu_diag)
          call shr_sys_abort( SubName//":: decomp inconsistent")
       endif

       write(nu_diag,*) subname,' read ',trim(gridcpl_file)
       write(nu_diag,*) subname,' size ',nx,ny,gsize
    endif

    call shr_mpi_bcast(nxcpl,mpicom_i)
    call shr_mpi_bcast(nycpl,mpicom_i)
    call shr_mpi_bcast(gsize,mpicom_i)
    call mct_gsmap_init(gsmap_i,npes,start,length,pe_loc,0,mpicom_i,ICEID,gsize)
    deallocate(start,length,pe_loc)

    lsize = mct_gsmap_lsize(gsmap_i,mpicom_i)
    call mct_gGrid_init( GGrid=dom_i, &
         CoordChars=trim(shr_flds_dom_coord), OtherChars=trim(shr_flds_dom_other), lsize=lsize )
    call mct_aVect_zero(dom_i%data)

    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT

    call mct_gsMap_orderedPoints(gsMap_i, my_task, idata)
    call mct_gGrid_importIAttr(dom_i,'GlobGridNum',idata,lsize)
    deallocate(idata)

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
    deallocate(data)

    ! Read domain arrays

    if (iam == 0) then
       call mct_avect_init(avg,rList='fld',lsize=gsize)
    endif

    do n = 1,5

       if (n == 1) avfld = 'lat'
       if (n == 1) dofld = 'yc'
       if (n == 2) avfld = 'lon'
       if (n == 2) dofld = 'xc'
       if (n == 3) avfld = 'area'
       if (n == 3) dofld = 'area'
       if (n == 4) avfld = 'frac'
       if (n == 4) dofld = 'frac'
       if (n == 5) avfld = 'mask'
       if (n == 5) dofld = 'mask'
       if (iam == 0) then
          rcode = nf_inq_varid(fid,trim(dofld),vid)
          if (n == 5) then
             allocate(idata(gsize))
             rcode = nf_get_var_int(fid,vid,idata)
             avg%rAttr(1,:) = idata
             deallocate(idata)
          else
             rcode = nf_get_var_double(fid,vid,avg%rAttr(1,:))
          endif
       endif

       call mct_aVect_scatter(avg,av1,gsmap_i,0,mpicom_i)
       call mct_aVect_copy(av1, dom_i%data,'fld',avfld)

       if (iam == 0) then
          call mct_avect_clean(av1)
       endif

    enddo

    if (iam == 0) then
       call mct_avect_clean(avg)
    endif

  end subroutine ice_setcoupling_mct

  !=======================================================================

  subroutine ice_SetGSMap_mct( mpicom, ID, gsMap_ice, xoff, yoff, nxgin, nygin )

    ! Arguments
    integer        , intent(in)    :: mpicom
    integer        , intent(in)    :: ID
    type(mct_gsMap), intent(inout) :: gsMap_ice
    integer,optional, intent(in)   :: xoff   ! x offset
    integer,optional, intent(in)   :: yoff   ! y offset
    integer,optional, intent(in)   :: nxgin ! global size
    integer,optional, intent(in)   :: nygin ! global size

    ! Local variables
    integer     :: lat
    integer     :: lon
    integer     :: i, j, iblk, n, gi
    integer     :: lsize,gsize
    integer     :: lxoff,lyoff,nxg,nyg
    integer     :: ier
    integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block) :: this_block         ! block information for current block
    integer,allocatable :: gindex(:)
    !--------------------------------

    ! Build the CICE grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    lxoff = 1
    if (present(xoff)) then
       lxoff = xoff
    endif
    lyoff = 1
    if (present(yoff)) then
       lyoff = yoff
    endif
    nxg = nx_global
    if (present(nxgin)) then
       nxg = nxgin
    endif
    nyg = ny_global
    if (present(nygin)) then
       nyg = nygin
    endif
    gsize = nxg*nyg

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
          enddo !i
       enddo    !j
    enddo       !iblk
    lsize = n

    ! allocate gindex
    allocate(gindex(lsize),stat=ier)

    ! set gindex
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
             lon = this_block%i_glob(i) + lxoff - 1
             lat = this_block%j_glob(j) + lyoff - 1
             gi = (lat-1)*nxg + lon
             gindex(n) = gi
          enddo !i
       enddo    !j
    enddo       !iblk

    ! now set gsmap once gindex is available
    call mct_gsMap_init( gsMap_ice, gindex, mpicom, ID, lsize, gsize )

    ! no longer need gindex - can deallocate it
    deallocate(gindex)

  end subroutine ice_SetGSMap_mct

end module cice_comp_nuopc
