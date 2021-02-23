module ice_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for CICE
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC                  , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                  , only : NUOPC_CompFilterPhaseMap, NUOPC_IsUpdated, NUOPC_IsAtTime
  use NUOPC                  , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC                  , only : NUOPC_SetAttribute, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC_Model            , only : model_routine_SS           => SetServices
  use NUOPC_Model            , only : model_label_Advance        => label_Advance
  use NUOPC_Model            , only : model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model            , only : model_label_SetRunClock    => label_SetRunClock
  use NUOPC_Model            , only : model_label_Finalize       => label_Finalize
  use NUOPC_Model            , only : NUOPC_ModelGet, SetVM
  use shr_kind_mod           , only : r8 => shr_kind_r8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod            , only : shr_sys_abort
  use shr_file_mod           , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_orb_mod            , only : shr_orb_decl, shr_orb_params, SHR_ORB_UNDEF_REAL, SHR_ORB_UNDEF_INT
  use shr_cal_mod            , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
  use ice_shr_methods        , only : chkerr, state_setscalar, state_getscalar, state_diagnose, alarmInit
  use ice_shr_methods        , only : set_component_logging, get_component_instance, state_flddebug
  use ice_import_export      , only : ice_import, ice_export, ice_advertise_fields, ice_realize_fields
  use ice_domain_size        , only : nx_global, ny_global
  use ice_communicate        , only : init_communicate, my_task, master_task, mpi_comm_ice
  use ice_calendar           , only : force_restart_now, write_ic
  use ice_calendar           , only : idate, mday, time, month, daycal, time2sec, year_init
  use ice_calendar           , only : sec, dt, calendar, calendar_type, nextsw_cday, istep
  use ice_orbital            , only : eccen, obliqr, lambm0, mvelpp
  use ice_kinds_mod          , only : dbl_kind, int_kind
  use ice_scam               , only : scmlat, scmlon, single_column
  use ice_fileunits          , only : nu_diag, inst_index, inst_name, inst_suffix, release_all_fileunits
  use ice_ocean              , only : tfrz_option
  use ice_therm_shared       , only : ktherm
  use ice_restart_shared     , only : runid, runtype, restart_dir, restart_file
  use ice_history            , only : accum_hist
  use ice_history_shared     , only : model_doi_url  ! TODO: add this functionality
  use ice_prescribed_mod     , only : ice_prescribed_init
  use ice_atmo               , only : flux_convergence_tolerance, flux_convergence_max_iteration
  use ice_atmo               , only : use_coldair_outbreak_mod
  use ice_mesh_mod           , only : ice_mesh_set_distgrid, ice_mesh_setmask_from_maskfile, ice_mesh_check
  use ice_grid               , only : init_grid2, grid_type
  use CICE_InitMod           , only : CICE_Init1, CICE_Init2
  use CICE_RunMod            , only : CICE_Run
  use perf_mod               , only : t_startf, t_stopf, t_barrierf
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
  private :: cice_orbital_init
  private :: cice_orbital_update

  character(len=CL) :: flds_scalar_name = ''
  integer           :: flds_scalar_num = 0
  integer           :: flds_scalar_index_nx = 0
  integer           :: flds_scalar_index_ny = 0
  integer           :: flds_scalar_index_nextsw_cday = 0

  integer     , parameter :: dbug = 0
  character(*), parameter :: modName =  "(ice_comp_nuopc)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

  character(len=CL)      :: orb_mode        ! attribute - orbital mode
  integer                :: orb_iyear       ! attribute - orbital year
  integer                :: orb_iyear_align ! attribute - associated with model year
  real(R8)               :: orb_obliq       ! attribute - obliquity in degrees
  real(R8)               :: orb_mvelp       ! attribute - moving vernal equinox longitude
  real(R8)               :: orb_eccen       ! attribute and update-  orbital eccentricity

  character(len=*) , parameter :: orb_fixed_year       = 'fixed_year'
  character(len=*) , parameter :: orb_variable_year    = 'variable_year'
  character(len=*) , parameter :: orb_fixed_parameters = 'fixed_parameters'

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
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

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
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local variables
    character(len=CL)  :: cvalue
    character(len=CL)  :: logmsg
    logical            :: isPresent, isSet
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldName')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldCount')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(logmsg,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNX')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(logmsg,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNY')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nextsw_cday
       write(logmsg,*) flds_scalar_index_nextsw_cday
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nextsw_cday = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxNextSwCday')
    endif

    call ice_advertise_fields(gcomp, importState, exportState, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
    integer                 :: lmpicom
    integer                 :: localPet
    integer                 :: npes
    type(ESMF_DistGrid)     :: ice_distgrid
    type(ESMF_Mesh)         :: ice_mesh
    type(ESMF_VM)           :: vm
    type(ESMF_Time)         :: currTime           ! Current time
    type(ESMF_Time)         :: startTime          ! Start time
    type(ESMF_Time)         :: stopTime           ! Stop time
    type(ESMF_Time)         :: refTime            ! Ref time
    type(ESMF_TimeInterval) :: timeStep           ! Model timestep
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
    integer                 :: shrlogunit         ! original log unit
    character(len=cs)       :: starttype          ! infodata start type
    character(CL)           :: cvalue
    logical                 :: isPresent
    character(len=CL)       :: ice_meshfile
    character(len=CL)       :: ice_maskfile
    character(*), parameter     :: F00   = "('(ice_comp_nuopc) ',2a,1x,d21.14)"
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, localPet=localPet, PetCount=npes, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! Initialize cice communicators
    !----------------------------------------------------------------------------

    call init_communicate(lmpicom) ! initial setup for message passing

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    inst_name = "ICE"

    !----------------------------------------------------------------------------
    ! start cice timers
    !----------------------------------------------------------------------------

    call t_startf ('cice_init_total')

    !----------------------------------------------------------------------------
    ! Get orbital values
    !----------------------------------------------------------------------------

    ! Note that these values are obtained in a call to init_orbit in ice_shortwave.F90
    ! if CESMCOUPLED is not defined

    call cice_orbital_init(gcomp, nu_diag, my_task==master_task, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call cice_orbital_update(clock, nu_diag, my_task==master_task, eccen, obliqr, lambm0, mvelpp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! Determine attributes - also needed in realize phase to get grid information
    !----------------------------------------------------------------------------

    ! Determine runtype and possibly nextsw_cday
    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
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
       ! In the nuopc version it will be easier to assume that on startup - nextsw_cday is just the current time

       ! TOOD (mvertens, 2019-03-21): need to get the perpetual run working

       if (trim(runtype) /= 'initial') then
          ! Set nextsw_cday to -1 (this will skip an orbital calculation on initialization
          nextsw_cday = -1.0_r8
       else
          call ESMF_ClockGet( clock, currTime=currTime, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    else
       ! This would be the NEMS branch
       ! Note that in NEMS - nextsw_cday is not needed in ice_orbital.F90 and what is needed is
       ! simply a CPP variable declaratino of NEMSCOUPLED

       runtype = 'initial' ! determined from the namelist in ice_init if CESMCOUPLED is not defined
    end if

    ! Determine single column info
    call NUOPC_CompAttributeGet(gcomp, name='single_column', value=cvalue, isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       read(cvalue,*) single_column
    else
       single_column = .false.
    end if
    if (single_column) then
       ! Must have these attributes present
       call NUOPC_CompAttributeGet(gcomp, name='scmlon', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) scmlon
       call NUOPC_CompAttributeGet(gcomp, name='scmlat', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) scmlat
    end if

    ! Determine runid
    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       read(cvalue,*) runid
    else
       runid = 'unknown'  ! read in from the namelist in ice_init.F90 if CESMCOUPLED is not defined
    end if

    ! Determine tfreeze_option, flux convertence before call to cice_init
    call NUOPC_CompAttributeGet(gcomp, name="tfreeze_option", value=tfrz_option, isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (.not. isPresent) then
       tfrz_option = 'linear_salt'  ! TODO: is this right? This must be the same as mom is using for the calculation.
    end if

    call NUOPC_CompAttributeGet(gcomp, name="flux_convergence", value=cvalue, isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       read(cvalue,*) flux_convergence_tolerance
    else
       flux_convergence_tolerance = 0._r8
    end if

    call NUOPC_CompAttributeGet(gcomp, name="flux_max_iteration", value=cvalue, isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       read(cvalue,*) flux_convergence_max_iteration
    else
       flux_convergence_max_iteration = 5
    end if

    call NUOPC_CompAttributeGet(gcomp, name="coldair_outbreak_mod", value=cvalue, isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       read(cvalue,*) use_coldair_outbreak_mod
    else
       use_coldair_outbreak_mod = .false.
    end if

    ! Get clock information before call to cice_init

    call ESMF_ClockGet( clock, &
         currTime=currTime, startTime=startTime, stopTime=stopTime, refTime=RefTime, &
         timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    call ESMF_TimeGet( startTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,start_ymd)

    call ESMF_TimeGet( stopTime, yy=yy, mm=mm, dd=dd, s=stop_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,stop_ymd)

    call ESMF_TimeGet( refTime, yy=yy, mm=mm, dd=dd, s=ref_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,ref_ymd)

    call ESMF_TimeIntervalGet( timeStep, s=dtime, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    dt = real(dtime)

    call ESMF_TimeGet( currTime, calkindflag=esmf_caltype, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar_type = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar_type = shr_cal_gregorian
    else
       call shr_sys_abort( subname//'ERROR:: bad calendar for ESMF' )
    end if

    !----------------------------------------------------------------------------
    ! Set cice logging
    !----------------------------------------------------------------------------
    ! Note that sets the nu_diag module variable in ice_fileunits
    ! nu_diag in this module is initialized to 0 in the module, and if this reset does not
    ! happen here - then ice_init.F90 will obtain it from the input file ice_modelio.nml

    call set_component_logging(gcomp, my_task==master_task, nu_diag, shrlogunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_file_setLogUnit (shrlogunit)

    !----------------------------------------------------------------------------
    ! Initialize cice1
    !----------------------------------------------------------------------------

    !call t_startf ('cice_init1')
    call cice_init1()
    !call t_stopf ('cice_init1')

    !----------------------------------------------------------------------------
    ! Initialize grid info
    !----------------------------------------------------------------------------

    ! Initialize cice mesh and mask if appropriate
    !call t_startf ('cice_initgrid')

    ! Determine mesh input file
    call NUOPC_CompAttributeGet(gcomp, name='mesh_ice', value=ice_meshfile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine mask input file
    call NUOPC_CompAttributeGet(gcomp, name='mesh_mask', value=ice_maskfile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (my_task == master_task) then
       write(nu_diag,*)'mesh file for cice domain is ',trim(ice_meshfile)
       write(nu_diag,*)'mask file for cice domain is ',trim(ice_maskfile)
    end if

    ! Determine the model distgrid using the decomposition obtained in
    ! call to init_grid1 called from cice_init1
    call ice_mesh_set_distgrid(localpet, npes, ice_distgrid, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read in the ice mesh on the cice distribution
    ice_mesh = ESMF_MeshCreate(filename=trim(ice_meshfile), fileformat=ESMF_FILEFORMAT_ESMFMESH, &
         elementDistGrid=ice_distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize the cice mesh and the cice mask
    if (trim(grid_type) == 'setmask') then
       ! In this case cap code determines the mask file
       call ice_mesh_setmask_from_maskfile(ice_mesh, ice_maskfile, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       ! In this case init_grid2 will initialize tlon, tlat, area and hm
       call init_grid2()
       call ice_mesh_check(ice_mesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !call t_stopf ('cice_initgrid')

    !----------------------------------------------------------------------------
    ! Initialize cice2
    !----------------------------------------------------------------------------

    ! Note that cice_init2 also sets time manager info as well as mpi communicator info,
    ! including master_task and my_task
    !call t_startf ('cice_init2')
    call cice_init2()
    !call t_stopf ('cice_init2')

    !----------------------------------------------------------------------------
    ! reset shr logging to my log file
    !----------------------------------------------------------------------------

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

    !-----------------------------------------------------------------
    ! Realize the actively coupled fields
    !-----------------------------------------------------------------

    call ice_realize_fields(gcomp, mesh=ice_mesh, &
         flds_scalar_name=flds_scalar_name, flds_scalar_num=flds_scalar_num, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------------------------------------------
    ! Prescribed ice initialization
    !-----------------------------------------------------------------

    call ice_prescribed_init(clock, ice_mesh, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------------------------------------------
    ! Create cice export state
    !-----------------------------------------------------------------

    call ice_export (exportstate, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call State_SetScalar(dble(nx_global), flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_SetScalar(dble(ny_global), flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO (mvertens, 2018-12-21): fill in iceberg_prognostic as .false.

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       call State_diagnose(exportState,subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

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
    character(*)   , parameter :: F00   = "('(ice_comp_nuopc) ',2a,i8,d21.14)"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

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
    call shr_file_setLogUnit (nu_diag)

    !--------------------------------
    ! Query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Determine time of next atmospheric shortwave calculation
    !--------------------------------

    call State_GetScalar(importState, flds_scalar_index_nextsw_cday, nextsw_cday, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (my_task == master_task) then
       write(nu_diag,F00) trim(subname),' cice istep, nextsw_cday = ',istep, nextsw_cday
    end if

    !--------------------------------
    ! Obtain orbital values
    !--------------------------------

    call cice_orbital_update(clock, nu_diag, my_task==master_task, eccen, obliqr, lambm0, mvelpp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! check that cice internal time is in sync with master clock before timestep update
    !--------------------------------

    ! cice clock
    tod = sec
    ymd = idate

    ! model clock
    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( currTime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr_sync, mon_sync, day_sync, ymd_sync)

    ! error check
    if ( (ymd /= ymd_sync) .or. (tod /= tod_sync) ) then
       if (my_task == master_task) then
          write(nu_diag,*)' cice ymd=',ymd     ,'  cice tod= ',tod
          write(nu_diag,*)' sync ymd=',ymd_sync,'  sync tod= ',tod_sync
       end if
       call ESMF_LogWrite(subname//" CICE clock not in sync with ESMF model clock",ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    !--------------------------------
    ! Determine if time to write restart
    !--------------------------------

    ! Note this logic triggers off of the component clock rather than the internal cice time
    ! The component clock does not get advanced until the end of the loop - not at the beginning

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       force_restart_now = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockGetNextTime(clock, nextTime=nextTime, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(nexttime, yy=yy, mm=mm, dd=dd, s=tod, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ice_timer_stop(timer_cplrecv)
    call t_stopf ('cice_run_import')

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
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ice_timer_stop(timer_cplsend)
    call t_stopf ('cice_run_export')

    ! reset shr logging to my original values
    call shr_file_setLogUnit (shrlogunit)

    !--------------------------------
    ! stop timers and print timer info
    !--------------------------------
    ! Need to have this logic here instead of in finalize phase
    ! since the finalize phase will still be called even in aqua-planet mode

    !--------------------------------
    ! Determine if time to stop
    !--------------------------------

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_stop', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       stop_now = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
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

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

    ! intput/output variables
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
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart and stop alarms
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       !----------------
       ! Restart alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------
       ! Stop alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_n

       call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_ymd

       call alarmInit(mclock, stop_alarm, stop_option, &
            opt_n   = stop_n,           &
            opt_ymd = stop_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_stop', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: F00   = "('(ice_comp_nuopc) ',8a)"
    character(*), parameter :: F91   = "('(ice_comp_nuopc) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !--------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    if (my_task == master_task) then
       write(nu_diag,F91)
       write(nu_diag,F00) 'CICE: end of main integration loop'
       write(nu_diag,F91)
    end if

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelFinalize

  !===============================================================================

  subroutine cice_orbital_init(gcomp, logunit, mastertask, rc)

    !----------------------------------------------------------
    ! Initialize orbital related values
    !----------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)                 :: gcomp
    integer             , intent(in)    :: logunit
    logical             , intent(in)    :: mastertask
    integer             , intent(out)   :: rc              ! output error

    ! local variables
    character(len=CL) :: msgstr          ! temporary
    character(len=CL) :: cvalue          ! temporary
    character(len=*) , parameter :: subname = "(cice_orbital_init)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine orbital attributes from input
    call NUOPC_CompAttributeGet(gcomp, name="orb_mode", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_mode

    call NUOPC_CompAttributeGet(gcomp, name="orb_iyear", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_iyear

    call NUOPC_CompAttributeGet(gcomp, name="orb_iyear_align", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_iyear_align

    call NUOPC_CompAttributeGet(gcomp, name="orb_obliq", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_obliq

    call NUOPC_CompAttributeGet(gcomp, name="orb_eccen", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_eccen

    call NUOPC_CompAttributeGet(gcomp, name="orb_mvelp", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_mvelp

    ! Error checks
    if (trim(orb_mode) == trim(orb_fixed_year)) then
       orb_obliq = SHR_ORB_UNDEF_REAL
       orb_eccen = SHR_ORB_UNDEF_REAL
       orb_mvelp = SHR_ORB_UNDEF_REAL
       if (orb_iyear == SHR_ORB_UNDEF_INT) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: fixed_year settings = ',orb_iyear
             write (msgstr, *) ' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    elseif (trim(orb_mode) == trim(orb_variable_year)) then
       orb_obliq = SHR_ORB_UNDEF_REAL
       orb_eccen = SHR_ORB_UNDEF_REAL
       orb_mvelp = SHR_ORB_UNDEF_REAL
       if (orb_iyear == SHR_ORB_UNDEF_INT .or. orb_iyear_align == SHR_ORB_UNDEF_INT) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: variable_year settings = ',orb_iyear, orb_iyear_align
             write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    elseif (trim(orb_mode) == trim(orb_fixed_parameters)) then
       !-- force orb_iyear to undef to make sure shr_orb_params works properly
       orb_iyear = SHR_ORB_UNDEF_INT
       orb_iyear_align = SHR_ORB_UNDEF_INT
       if (orb_eccen == SHR_ORB_UNDEF_REAL .or. &
           orb_obliq == SHR_ORB_UNDEF_REAL .or. &
           orb_mvelp == SHR_ORB_UNDEF_REAL) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: orb_eccen = ',orb_eccen
             write(logunit,*) trim(subname),' ERROR: orb_obliq = ',orb_obliq
             write(logunit,*) trim(subname),' ERROR: orb_mvelp = ',orb_mvelp
             write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    else
       write (msgstr, *) subname//' ERROR: invalid orb_mode '//trim(orb_mode)
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       rc = ESMF_FAILURE
       return  ! bail out
    endif

  end subroutine cice_orbital_init

  !===============================================================================

  subroutine cice_orbital_update(clock, logunit,  mastertask, eccen, obliqr, lambm0, mvelpp, rc)

    !----------------------------------------------------------
    ! Update orbital settings
    !----------------------------------------------------------

    ! input/output variables
    type(ESMF_Clock) , intent(in)    :: clock
    integer          , intent(in)    :: logunit
    logical          , intent(in)    :: mastertask
    real(R8)         , intent(inout) :: eccen  ! orbital eccentricity
    real(R8)         , intent(inout) :: obliqr ! Earths obliquity in rad
    real(R8)         , intent(inout) :: lambm0 ! Mean long of perihelion at vernal equinox (radians)
    real(R8)         , intent(inout) :: mvelpp ! moving vernal equinox longitude of perihelion plus pi (radians)
    integer          , intent(out)   :: rc     ! output error

    ! local variables
    type(ESMF_Time)   :: CurrTime ! current time
    integer           :: year     ! model year at current time
    integer           :: orb_year ! orbital year for current orbital computation
    character(len=CL) :: msgstr   ! temporary
    logical           :: lprint
    logical           :: first_time = .true.
    character(len=*) , parameter :: subname = "(cice_orbital_update)"
    !-------------------------------------------

    if (trim(orb_mode) == trim(orb_variable_year)) then
       call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(CurrTime, yy=year, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       orb_year = orb_iyear + (year - orb_iyear_align)
       lprint = mastertask
    else
       orb_year = orb_iyear
       if (first_time) then
          lprint = mastertask
          first_time = .false.
       else
          lprint = .false.
       end if
    end if

    eccen = orb_eccen
    call shr_orb_params(orb_year, eccen, orb_obliq, orb_mvelp, obliqr, lambm0, mvelpp, lprint)

    if ( eccen  == SHR_ORB_UNDEF_REAL .or. obliqr == SHR_ORB_UNDEF_REAL .or. &
         mvelpp == SHR_ORB_UNDEF_REAL .or. lambm0 == SHR_ORB_UNDEF_REAL) then
       write (msgstr, *) subname//' ERROR: orb params incorrect'
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    endif

  end subroutine cice_orbital_update

end module ice_comp_nuopc
