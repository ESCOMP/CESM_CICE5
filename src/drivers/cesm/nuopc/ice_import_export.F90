module ice_import_export

  use ESMF
  use NUOPC
  use NUOPC_Model
  use med_constants_mod     , only : R8, CS,CL
  use shr_sys_mod           , only : shr_sys_abort, shr_sys_flush
  use shr_frz_mod           , only : shr_frz_freezetemp
  use ice_kinds_mod         , only : int_kind, dbl_kind, char_len_long, log_kind
  use ice_constants         , only : c0, c1, puny, tffresh, spval_dbl
  use ice_constants         , only : field_loc_center, field_type_scalar
  use ice_constants         , only : field_type_vector, c100
  use ice_constants         , only : vonkar, zref, iceruf
  use ice_blocks            , only : block, get_block, nx_block, ny_block
  use ice_flux              , only : strairxt, strairyt, strocnxt, strocnyt
  use ice_flux              , only : alvdr, alidr, alvdf, alidf, Tref, Qref, Uref
  use ice_flux              , only : flat, fsens, flwout, evap, fswabs, fhocn, fswthru
  use ice_flux              , only : fresh, fsalt, zlvl, uatm, vatm, potT, Tair, Qa
  use ice_flux              , only : rhoa, swvdr, swvdf, swidr, swidf, flw, frain
  use ice_flux              , only : fsnow, uocn, vocn, sst, ss_tltx, ss_tlty, frzmlt
  use ice_flux              , only : sss, tf, wind, fsw, init_flux_atm, init_flux_ocn
  use ice_flux              , only : faero_atm, faero_ocn
  use ice_flux              , only : fiso_atm, fiso_ocn, fiso_rain, fiso_evap
  use ice_flux              , only : Qa_iso, Qref_iso, HDO_ocn, H2_18O_ocn, H2_16O_ocn
  use ice_flux              , only : send_i2x_per_cat, fswthrun_ai
  use ice_ocean             , only : tfrz_option
  use ice_atmo              , only : Cdn_atm
  use ice_state             , only : vice, vsno, aice, aicen_init, trcr
  use ice_state             , only : tr_aero, tr_iso, tr_iage, tr_FY, tr_pond, tr_lvl
  use ice_domain            , only : nblocks, blocks_ice, halo_info, distrb_info
  use ice_domain_size       , only : nx_global, ny_global, block_size_x, block_size_y, max_blocks, ncat
  use ice_grid              , only : tlon, tlat, tarea, tmask, anglet, hm, ocn_gridcell_frac
  use ice_grid              , only : grid_type, t2ugrid_vector
  use ice_boundary          , only : ice_HaloUpdate
  use ice_fileunits         , only : nu_diag
  use ice_communicate       , only : my_task, master_task, MPI_COMM_ICE
  use ice_calendar          , only : istep, istep1
  use ice_prescribed_mod    , only : prescribed_ice
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use shr_nuopc_scalars_mod , only : flds_scalar_num
  use shr_nuopc_scalars_mod , only : flds_scalar_index_nx
  use shr_nuopc_scalars_mod , only : flds_scalar_index_ny
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_chkerr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_state_reset 
  use perf_mod              , only : t_startf, t_stopf, t_barrierf

  implicit none
  public

  public  :: ice_advertise_fields
  public  :: ice_realize_fields
  public  :: ice_import
  public  :: ice_export

  private :: fldlist_add
  private :: fldlist_realize
  private :: state_fldchk

  interface state_getfldptr 
     module procedure state_getfldptr_1d
     module procedure state_getfldptr_2d
  end interface state_getfldptr
  private :: state_getfldptr

  ! Private module data

  type fld_list_type
    character(len=128) :: stdname
    integer :: ungridded_lbound = 0
    integer :: ungridded_ubound = 0
  end type fld_list_type

  integer, parameter     :: fldsMax = 100
  integer                :: fldsToIce_num = 0
  integer                :: fldsFrIce_num = 0
  type (fld_list_type)   :: fldsToIce(fldsMax)
  type (fld_list_type)   :: fldsFrIce(fldsMax)

#ifdef RASM_MODS
  ! (1)  Andrew Roberts:  Added artificial correction to snow and rain division
  !      This is to be consistent with VIC in the Regional Arctic System Model
  logical, parameter :: rasm_snowrain_split = .true.
#else
  logical, parameter :: rasm_snowrain_split = .false.
#endif
   integer     ,parameter :: dbug = 10 ! i/o debug messages
   integer     ,parameter :: debug = 0 ! internal debug level
   character(*),parameter :: F01 = "('(ice_import_export) ',a, i8,2x,i8,2x,d21.14)"

  character(*), parameter :: u_FILE_u = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine ice_advertise_fields(gcomp, importState, exportState, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    integer, intent(out) :: rc

    ! local variables
    integer       :: n 
    character(CS) :: stdname
    character(CS) :: cvalue
    logical       :: flds_wiso         ! use case
    logical       :: flds_i2o_per_cat  ! .true. => select per ice thickness category
    integer       :: dbrc
    character(len=*), parameter :: subname='(ice_import_export:ice_advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='flds_wiso', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_wiso
    call ESMF_LogWrite('flds_wiso = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='flds_i2o_per_cat', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_i2o_per_cat
    call ESMF_LogWrite('flds_i2o_per_cat = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    !-----------------
    ! advertise import fields
    !-----------------

    call fldlist_add(fldsToIce_num, fldsToIce, trim(flds_scalar_name))

    ! from ocean
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_dhdx')
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_dhdy')
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_t'   )
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_s'   )
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_u'   )
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_v'   )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Fioo_q' )
    if (flds_wiso) then
       call fldlist_add(fldsToIce_num, fldsToIce, 'So_roce_HDO')
       call fldlist_add(fldsToIce_num, fldsToIce, 'So_roce_16O')
       call fldlist_add(fldsToIce_num, fldsToIce, 'So_roce_18O')
    end if

    ! from atmosphere
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_z'         )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_u'         )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_v'         )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_ptem'      )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_shum'      )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_dens'      )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_tbot'      )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_swvdr'   )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_swndr'   )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_swvdf'   )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_swndf'   )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_lwdn'    )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_rain'    )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_snow'    )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_bcphodry')
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_bcphidry')
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_bcphiwet')
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_dstdry1' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_dstdry2' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_dstdry3' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_dstdry4' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_dstwet1' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_dstwet2' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_dstwet3' )
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_dstwet4' )
    if (flds_wiso) then
       call fldlist_add(fldsToIce_num, fldsToIce, 'So_roce_HDO')
       call fldlist_add(fldsToIce_num, fldsToIce, 'So_roce_16O')
       call fldlist_add(fldsToIce_num, fldsToIce, 'So_roce_18O')
    end if

    do n = 1,fldsToIce_num
       call NUOPC_Advertise(importState, standardName=fldsToIce(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    !-----------------
    ! advertise export fields
    !-----------------

    call fldlist_add(fldsFrIce_num, fldsFrIce, trim(flds_scalar_name))

    ! ice states
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_imask')
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_ifrac')
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_t'    )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_tref' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_qref' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_snowh')
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_vice' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_vsno' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_u10'  )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_avsdr')
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_anidr')
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_avsdf')
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_anidf')
    if (send_i2x_per_cat) then
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_ifrac_n', &
            ungridded_lbound=1, ungridded_ubound=ncat)
    end if

    ! ice/atm fluxes computed by ice
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_taux' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_tauy' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_lat'  )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_sen'  )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_lwup' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_evap' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_swnet')

    ! ice/ocn fluxes computed by ice
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_melth' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_swpen' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_meltw' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_salt'  )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_taux'  )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_tauy'  )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_bcpho' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_bcphi' )
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_flxdst')
    if (send_i2x_per_cat) then
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_swpen_ifrac_n', &
            ungridded_lbound=1, ungridded_ubound=ncat)
    end if
    if (flds_wiso) then
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_HDO'     )
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_16O'     )
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_18O'     )
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_evap_HDO')
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_evap_16O')
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_evap_18O')
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_qref_HDO'  )
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_qref_16O'  )
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_qref_18O'  )
    end if

    do n = 1,fldsFrIce_num
       call NUOPC_Advertise(exportState, standardName=fldsFrIce(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ice_advertise_fields

!==============================================================================

  subroutine ice_realize_fields(gcomp, Emesh, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_Mesh)      :: Emesh
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    character(len=*), parameter :: subname='(ice_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrIce, &
         numflds=fldsFrIce_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':CICE_Export',&
         mesh=Emesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToIce, &
         numflds=fldsToIce_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':CICE_Import',&
         mesh=Emesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine ice_realize_fields

  !==============================================================================

  subroutine ice_import( importState, rc )

    ! input/output variables
    type(ESMF_State) , intent(in)  :: importState
    integer          , intent(out) :: rc

    ! local variables
    integer,parameter                :: nflds=15
    integer,parameter                :: nfldv=6
    integer                          :: i, j, iblk, n
    integer                          :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    type(block)                      :: this_block         ! block information for current block
    real (kind=dbl_kind),allocatable :: aflds(:,:,:,:)
    real (kind=dbl_kind)             :: workx, worky
    real (kind=dbl_kind)             :: MIN_RAIN_TEMP, MAX_SNOW_TEMP
    real (kind=dbl_kind), pointer    :: dataPtr_Sa_tbot(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Sa_shum(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Sa_z(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Sa_dens(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Sa_u(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Sa_v(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Sa_ptem(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_lwdn(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_swvdr(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_swvdf(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_swndr(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_swndf(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_rain(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_snow(:)
    real (kind=dbl_kind), pointer    :: dataPtr_So_t(:)
    real (kind=dbl_kind), pointer    :: dataPtr_So_s(:)
    real (kind=dbl_kind), pointer    :: dataPtr_So_dhdx(:)
    real (kind=dbl_kind), pointer    :: dataPtr_So_dhdy(:)
    real (kind=dbl_kind), pointer    :: dataPtr_So_u(:)
    real (kind=dbl_kind), pointer    :: dataPtr_So_v(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Fioo_q(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_bcphodry(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_bcphidry(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_bcphiwet(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_dstdry1(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_dstdry2(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_dstdry3(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_dstdry4(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_dstwet1(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_dstwet2(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_dstwet3(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_dstwet4(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_rain_HDO(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_rain_16O(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_rain_18O(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_snow_HDO(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_snow_16O(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Faxa_snow_18O(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Sa_shum_HDO(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Sa_shum_16O(:)
    real (kind=dbl_kind), pointer    :: dataPtr_Sa_shum_18O(:)
    real (kind=dbl_kind), pointer    :: dataPtr_So_roce_HDO(:)
    real (kind=dbl_kind), pointer    :: dataPtr_So_roce_16O(:)
    real (kind=dbl_kind), pointer    :: dataPtr_So_roce_18O(:)
    character(len=*),   parameter    :: subname = 'ice_import'
    !-----------------------------------------------------

    call State_getFldPtr(importState,'Sa_tbot'       ,dataPtr_Sa_tbot       ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Sa_shum'       ,dataPtr_Sa_shum       ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Sa_z'          ,dataPtr_Sa_z          ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Sa_dens'       ,dataPtr_Sa_dens       ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Sa_u'          ,dataPtr_Sa_u          ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Sa_v'          ,dataPtr_Sa_v          ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Sa_ptem'       ,dataPtr_Sa_ptem       ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_lwdn'     ,dataPtr_Faxa_lwdn       ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_swvdr'    ,dataPtr_Faxa_swvdr    ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_swvdf'    ,dataPtr_Faxa_swvdf    ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_swndr'    ,dataPtr_Faxa_swndr    ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_swndf'    ,dataPtr_Faxa_swndf    ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_rain'     ,dataPtr_Faxa_rain     ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_snow'     ,dataPtr_Faxa_snow     ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'So_t'          ,dataPtr_So_t          ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'So_s'          ,dataPtr_So_s          ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'So_dhdx'       ,dataPtr_So_dhdx       ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'So_dhdy'       ,dataPtr_So_dhdy       ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'So_u'          ,dataPtr_So_u          ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'So_v'          ,dataPtr_So_v          ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Fioo_q'        ,dataPtr_Fioo_q        ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_bcphodry' ,dataPtr_Faxa_bcphodry ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_bcphidry' ,dataPtr_Faxa_bcphidry ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_bcphiwet' ,dataPtr_Faxa_bcphiwet ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_dstdry1'  ,dataPtr_Faxa_dstdry1  ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_dstdry2'  ,dataPtr_Faxa_dstdry2  ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_dstdry3'  ,dataPtr_Faxa_dstdry3  ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_dstdry4'  ,dataPtr_Faxa_dstdry4  ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_dstwet1'  ,dataPtr_Faxa_dstwet1  ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_dstwet2'  ,dataPtr_Faxa_dstwet2  ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_dstwet3'  ,dataPtr_Faxa_dstwet3  ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(importState,'Faxa_dstwet4'  ,dataPtr_Faxa_dstwet4  ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (State_FldChk(importState, 'Sa_shum_HDO')) then
       call State_getFldPtr(importState,'Sa_shum_HDO'  ,dataPtr_Sa_shum_HDO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(importState,'Sa_shum_18O'  ,dataPtr_Sa_shum_18O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(importState,'Sa_shum_16O'  ,dataPtr_Sa_shum_16O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call State_getFldPtr(importState,'So_roce_HDO'  ,dataPtr_So_roce_HDO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(importState,'So_roce_HDO'  ,dataPtr_So_roce_18O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(importState,'So_roce_HDO'  ,dataPtr_So_roce_16O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call State_getFldPtr(importState,'Faxa_rain_HDO'  ,dataPtr_Faxa_rain_HDO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(importState,'Faxa_rain_HDO'  ,dataPtr_Faxa_rain_18O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(importState,'Faxa_rain_HDO'  ,dataPtr_Faxa_rain_16O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call State_getFldPtr(importState,'Faxa_snow_HDO'  ,dataPtr_Faxa_snow_HDO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(importState,'Faxa_snow_HDO'  ,dataPtr_Faxa_snow_18O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(importState,'Faxa_snow_HDO'  ,dataPtr_Faxa_snow_16O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Note that the precipitation fluxes received from the mediator
    ! are in units of kg/s/m^2 which is what CICE requires.
    ! Note also that the read in below includes only values needed
    ! by the thermodynamic component of CICE.  Variables uocn, vocn,
    ! ss_tltx, and ss_tlty are excluded. Also, because the SOM and
    ! DOM don't  compute SSS.   SSS is not read in and is left at
    ! the initilized value (see ice_flux.F init_coupler_flux) of
    ! 34 ppt

    ! Use aflds to gather the halo updates of multiple fields
    ! Need to separate the scalar from the vector halo updates

    allocate(aflds(nx_block,ny_block,nflds,nblocks))
    aflds = c0

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
             aflds(i,j, 1,iblk)   = dataPtr_So_t(n)
             aflds(i,j, 2,iblk)   = dataPtr_So_s(n)
             aflds(i,j, 3,iblk)   = dataPtr_Sa_z(n)
             aflds(i,j, 4,iblk)   = dataPtr_Sa_ptem(n)
             aflds(i,j, 5,iblk)   = dataPtr_Sa_tbot(n)
             aflds(i,j, 6,iblk)   = dataPtr_Sa_shum(n)
             aflds(i,j, 7,iblk)   = dataPtr_Sa_dens(n)
             aflds(i,j, 8,iblk)   = dataPtr_Fioo_q(n)
             aflds(i,j, 9,iblk)   = dataPtr_Faxa_swvdr(n)
             aflds(i,j,10,iblk)   = dataPtr_Faxa_swndr(n)
             aflds(i,j,11,iblk)   = dataPtr_Faxa_swvdf(n)
             aflds(i,j,12,iblk)   = dataPtr_Faxa_swndf(n)
             aflds(i,j,13,iblk)   = dataPtr_Faxa_lwdn(n)
             aflds(i,j,14,iblk)   = dataPtr_Faxa_rain(n)
             aflds(i,j,15,iblk)   = dataPtr_Faxa_snow(n)
          enddo  !i
       enddo     !j
    enddo        !iblk

    if (.not.prescribed_ice) then
       call t_startf ('cice_imp_halo')
       call ice_HaloUpdate(aflds, halo_info, field_loc_center, field_type_scalar)
       call t_stopf ('cice_imp_halo')
    endif

    if (rasm_snowrain_split) then
       MIN_RAIN_TEMP = Tffresh-c1
       MAX_SNOW_TEMP = Tffresh+c0
    endif

    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1,ny_block
          do i = 1,nx_block
             sst  (i,j,iblk)   = aflds(i,j, 1,iblk)
             sss  (i,j,iblk)   = aflds(i,j, 2,iblk)
             zlvl (i,j,iblk)   = aflds(i,j, 3,iblk)
             potT (i,j,iblk)   = aflds(i,j, 4,iblk)
             Tair (i,j,iblk)   = aflds(i,j, 5,iblk)
             Qa   (i,j,iblk)   = aflds(i,j, 6,iblk)
             rhoa (i,j,iblk)   = aflds(i,j, 7,iblk)
             frzmlt (i,j,iblk) = aflds(i,j, 8,iblk)
             swvdr(i,j,iblk)   = aflds(i,j, 9,iblk)
             swidr(i,j,iblk)   = aflds(i,j,10,iblk)
             swvdf(i,j,iblk)   = aflds(i,j,11,iblk)
             swidf(i,j,iblk)   = aflds(i,j,12,iblk)
             flw  (i,j,iblk)   = aflds(i,j,13,iblk)
             frain(i,j,iblk)   = aflds(i,j,14,iblk)
             fsnow(i,j,iblk)   = aflds(i,j,15,iblk)
          enddo    !i
       enddo    !j
    enddo        !iblk
    !$OMP END PARALLEL DO

    if (rasm_snowrain_split) then
       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
       do iblk = 1, nblocks
          do j = 1,ny_block
             do i = 1,nx_block
                !--- Artificial correction to snow and rain for RASM
                if (Tair(i,j,iblk)<MIN_RAIN_TEMP) then
                   fsnow(i,j,iblk)=fsnow(i,j,iblk)+frain(i,j,iblk)
                   frain(i,j,iblk)=0
                elseif (Tair(i,j,iblk)>MAX_SNOW_TEMP) then
                   frain(i,j,iblk)=fsnow(i,j,iblk)+frain(i,j,iblk)
                   fsnow(i,j,iblk)=0
                else
                   frain(i,j,iblk)=fsnow(i,j,iblk)+frain(i,j,iblk)
                   fsnow(i,j,iblk)=frain(i,j,iblk)
                   frain(i,j,iblk)=frain(i,j,iblk)*(Tair(i,j,iblk)-MIN_RAIN_TEMP) / &
                                                   (MAX_SNOW_TEMP-MIN_RAIN_TEMP)
                   fsnow(i,j,iblk)=fsnow(i,j,iblk)-frain(i,j,iblk)
                endif
                !--- end artificial RASM correction
             enddo    !i
          enddo    !j
       enddo        !iblk
       !$OMP END PARALLEL DO
    endif  ! rasm_snowrain_split

    deallocate(aflds)
    allocate(aflds(nx_block,ny_block,nfldv,nblocks))
    aflds = c0

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
             aflds(i,j, 1,iblk)   = dataPtr_So_u(n)
             aflds(i,j, 2,iblk)   = dataPtr_So_v(n)
             aflds(i,j, 3,iblk)   = dataPtr_Sa_u(n)
             aflds(i,j, 4,iblk)   = dataPtr_Sa_v(n)
             aflds(i,j, 5,iblk)   = dataPtr_So_dhdx(n)
             aflds(i,j, 6,iblk)   = dataPtr_So_dhdy(n)
          enddo
       enddo
    enddo

    if (.not.prescribed_ice) then
       call t_startf ('cice_imp_halo')
       call ice_HaloUpdate(aflds, halo_info, field_loc_center, field_type_vector)
       call t_stopf ('cice_imp_halo')
    endif

    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1,ny_block
          do i = 1,nx_block
             uocn (i,j,iblk)   = aflds(i,j, 1,iblk)
             vocn (i,j,iblk)   = aflds(i,j, 2,iblk)
             uatm (i,j,iblk)   = aflds(i,j, 3,iblk)
             vatm (i,j,iblk)   = aflds(i,j, 4,iblk)
             ss_tltx(i,j,iblk) = aflds(i,j, 5,iblk)
             ss_tlty(i,j,iblk) = aflds(i,j, 6,iblk)
          enddo  !i
       enddo     !j
    enddo        !iblk
    !$OMP END PARALLEL DO

    deallocate(aflds)

    !-------------------------------------------------------
    ! Set aerosols from coupler
    !-------------------------------------------------------

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
             faero_atm(i,j,1,iblk) = dataPtr_Faxa_bcphodry(n)
             faero_atm(i,j,2,iblk) = dataPtr_Faxa_bcphidry(n) + dataPtr_Faxa_bcphiwet(n)

             ! Combine all of the dust into one category
             faero_atm(i,j,3,iblk) = dataPtr_Faxa_dstwet1(n) + dataPtr_Faxa_dstdry1(n) &
                                   + dataPtr_Faxa_dstwet2(n) + dataPtr_Faxa_dstdry2(n) &
                                   + dataPtr_Faxa_dstwet3(n) + dataPtr_Faxa_dstdry3(n) &
                                   + dataPtr_Faxa_dstwet4(n) + dataPtr_Faxa_dstdry4(n)
          end do
       end do
    end do

    !-------------------------------------------------------
    ! Water isotopes form the mediator
    !-------------------------------------------------------

    if (State_FldChk(importState, 'shum_HDO')) then
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
                Qa_iso(i,j,1,iblk)  = dataPtr_Sa_shum_HDO(n)
                Qa_iso(i,j,2,iblk)  = dataPtr_Sa_shum_16O(n)
                Qa_iso(i,j,3,iblk)  = dataPtr_Sa_shum_18O(n)

                fiso_rain(i,j,1,iblk) = dataPtr_Faxa_rain_HDO(n)
                fiso_rain(i,j,2,iblk) = dataPtr_Faxa_rain_16O(n)
                fiso_rain(i,j,3,iblk) = dataPtr_Faxa_rain_18O(n)

                fiso_atm(i,j,1,iblk) = dataPtr_Faxa_snow_HDO(n)
                fiso_atm(i,j,2,iblk) = dataPtr_Faxa_snow_16O(n)
                fiso_atm(i,j,3,iblk) = dataPtr_Faxa_snow_18O(n)

                HDO_ocn(i,j,iblk)    = dataPtr_So_roce_HDO(n)
                H2_16O_ocn(i,j,iblk) = dataPtr_So_roce_16O(n)
                H2_18O_ocn(i,j,iblk) = dataPtr_So_roce_18O(n)
             enddo  !i
          enddo     !j
       enddo        !iblk
    end if

    !-----------------------------------------------------------------
    ! rotate zonal/meridional vectors to local coordinates
    ! compute data derived quantities
    !-----------------------------------------------------------------

    ! Vector fields come in on T grid, but are oriented geographically
    ! need to rotate to pop-grid FIRST using ANGLET
    ! then interpolate to the U-cell centers  (otherwise we
    ! interpolate across the pole)
    ! use ANGLET which is on the T grid !

    call t_startf ('cice_imp_ocn')
    !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
    do iblk = 1, nblocks

       do j = 1,ny_block
          do i = 1,nx_block

             ! ocean
             workx      = uocn  (i,j,iblk) ! currents, m/s
             worky      = vocn  (i,j,iblk)
             uocn(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                            + worky*sin(ANGLET(i,j,iblk))
             vocn(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                            - workx*sin(ANGLET(i,j,iblk))

             workx      = ss_tltx  (i,j,iblk)           ! sea sfc tilt, m/m
             worky      = ss_tlty  (i,j,iblk)
             ss_tltx(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                               + worky*sin(ANGLET(i,j,iblk))
             ss_tlty(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                               - workx*sin(ANGLET(i,j,iblk))

             sst(i,j,iblk) = sst(i,j,iblk) - Tffresh       ! sea sfc temp (C)

             sss(i,j,iblk) = max(sss(i,j,iblk),c0)

          enddo
       enddo

       ! Use shr_frz_mod for this
       Tf(:,:,iblk) = shr_frz_freezetemp(sss(:,:,iblk))

    enddo
    !$OMP END PARALLEL DO
    call t_stopf ('cice_imp_ocn')

    ! Interpolate ocean dynamics variables from T-cell centers to
    ! U-cell centers.

    if (.not.prescribed_ice) then
       call t_startf ('cice_imp_t2u')
       call t2ugrid_vector(uocn)
       call t2ugrid_vector(vocn)
       call t2ugrid_vector(ss_tltx)
       call t2ugrid_vector(ss_tlty)
       call t_stopf ('cice_imp_t2u')
    end if

    ! Atmosphere variables are needed in T cell centers in
    ! subroutine stability and are interpolated to the U grid
    ! later as necessary.

    call t_startf ('cice_imp_atm')
    !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
    do iblk = 1, nblocks
       do j = 1, ny_block
          do i = 1, nx_block

             ! atmosphere
             workx      = uatm(i,j,iblk) ! wind velocity, m/s
             worky      = vatm(i,j,iblk)
             uatm (i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                             + worky*sin(ANGLET(i,j,iblk))   ! note uatm, vatm, wind
             vatm (i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) & ! are on the T-grid here
                             - workx*sin(ANGLET(i,j,iblk))

             wind (i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
             fsw  (i,j,iblk) = swvdr(i,j,iblk) + swvdf(i,j,iblk) &
                             + swidr(i,j,iblk) + swidf(i,j,iblk)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call t_stopf ('cice_imp_atm')

    !-----------------------------------------------------------------
    ! debug output
    !-----------------------------------------------------------------

    if (debug > 0 .and. my_task==master_task) then
       n = 0
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo; ihi = this_block%ihi
          jlo = this_block%jlo; jhi = this_block%jhi
          do j = jlo, jhi
             do i = ilo, ihi
                n = n+1
                write(nu_diag,F01)'import: istep, n, So_dhdx       = ',istep1,n,dataPtr_So_dhdx(n)
                write(nu_diag,F01)'import: istep, n, So_dhdy       = ',istep1,n,dataPtr_So_dhdy(n)
                write(nu_diag,F01)'import: istep, n, So_t          = ',istep1,n,dataPtr_So_t(n)
                write(nu_diag,F01)'import: istep, n, So_s          = ',istep1,n,dataPtr_So_s(n)
                write(nu_diag,F01)'import: istep, n, So_u          = ',istep1,n,dataPtr_so_u(n)
                write(nu_diag,F01)'import: istep, n, So_v          = ',istep1,n,dataPtr_So_v(n)
                write(nu_diag,F01)'import: istep, n, Sa_u          = ',istep1,n,dataPtr_Sa_u(n)
                write(nu_diag,F01)'import: istep, n, Sa_v          = ',istep1,n,dataPtr_Sa_v(n)
                write(nu_diag,F01)'import: istep, n, Sa_z          = ',istep1,n,dataPtr_Sa_z(n)
                write(nu_diag,F01)'import: istep, n, Sa_ptem       = ',istep1,n,dataPtr_Sa_ptem(n)
                write(nu_diag,F01)'import: istep, n, Sa_tbot       = ',istep1,n,dataPtr_Sa_tbot(n)
                write(nu_diag,F01)'import: istep, n, Sa_shum       = ',istep1,n,dataPtr_Sa_shum(n)
                write(nu_diag,F01)'import: istep, n, Sa_dens       = ',istep1,n,dataPtr_Sa_dens(n)
                write(nu_diag,F01)'import: istep, n, Fioo_q        = ',istep1,n,dataPtr_Fioo_q(n)
                write(nu_diag,F01)'import: istep, n, Faxa_swvdr    = ',istep1,n,dataPtr_Faxa_swvdr(n)
                write(nu_diag,F01)'import: istep, n, Faxa_swndr    = ',istep1,n,dataPtr_Faxa_swndr(n)
                write(nu_diag,F01)'import: istep, n, Faxa_swvdf    = ',istep1,n,dataPtr_Faxa_swvdf(n)
                write(nu_diag,F01)'import: istep, n, Faxa_swndf    = ',istep1,n,dataPtr_Faxa_swndf(n)
                write(nu_diag,F01)'import: istep, n, Faxa_lwdn     = ',istep1,n,dataPtr_Faxa_lwdn(n)
                write(nu_diag,F01)'import: istep, n, Faxa_rain     = ',istep1,n,dataPtr_Faxa_rain(n)
                write(nu_diag,F01)'import: istep, n, Faxa_snow     = ',istep1,n,dataPtr_Faxa_snow(n)
                write(nu_diag,F01)'import: istep, n, Faxa_bcphodry = ',istep1,n,dataPtr_Faxa_bcphodry(n)
                write(nu_diag,F01)'import: istep, n, Faxa_bcphidry = ',istep1,n,dataPtr_Faxa_bcphidry(n)
                write(nu_diag,F01)'import: istep, n, Faxa_bcphiwet = ',istep1,n,dataPtr_Faxa_bcphiwet(n)
                write(nu_diag,F01)'import: istep, n, Faxa_dstwet1  = ',istep1,n,dataPtr_Faxa_dstwet1(n)
                write(nu_diag,F01)'import: istep, n, Faxa_dstdry1  = ',istep1,n,dataPtr_Faxa_dstdry1(n)
                write(nu_diag,F01)'import: istep, n, Faxa_dstwet2  = ',istep1,n,dataPtr_Faxa_dstwet2(n)
                write(nu_diag,F01)'import: istep, n, Faxa_dstdry2  = ',istep1,n,dataPtr_Faxa_dstdry2(n)
                write(nu_diag,F01)'import: istep, n, Faxa_dstwet3  = ',istep1,n,dataPtr_Faxa_dstwet3(n)
                write(nu_diag,F01)'import: istep, n, Faxa_dstdry3  = ',istep1,n,dataPtr_Faxa_dstdry3(n)
                write(nu_diag,F01)'import: istep, n, Faxa_dstwet4  = ',istep1,n,dataPtr_Faxa_dstwet4(n)
                write(nu_diag,F01)'import: istep, n, Faxa_dstdry4  = ',istep1,n,dataPtr_Faxa_dstdry4(n)
                if (State_FldChk(importState, 'shum_HDO')) then
                   write(nu_diag,F01)'import: istep, n, Faxa_rain_HDO = ',istep1,n,dataPtr_Faxa_rain_HDO(n)
                   write(nu_diag,F01)'import: istep, n, Faxa_rain_16O = ',istep1,n,dataPtr_Faxa_rain_16O(n)
                   write(nu_diag,F01)'import: istep, n, Faxa_rain_18O = ',istep1,n,dataPtr_Faxa_rain_18O(n)
                   write(nu_diag,F01)'import: istep, n, Faxa_snow_HDO = ',istep1,n,dataPtr_Faxa_snow_HDO(n)
                   write(nu_diag,F01)'import: istep, n, Faxa_snow_16O = ',istep1,n,dataPtr_Faxa_snow_16O(n)
                   write(nu_diag,F01)'import: istep, n, Faxa_snow_18O = ',istep1,n,dataPtr_Faxa_snow_18O(n)
                   write(nu_diag,F01)'import: istep, n, Sa_shum_HDO   = ',istep1,n,dataPtr_Sa_shum_HDO(n)
                   write(nu_diag,F01)'import: istep, n, Sa_shum_16O   = ',istep1,n,dataPtr_Sa_shum_16O(n)
                   write(nu_diag,F01)'import: istep, n, Sa_shum_18O   = ',istep1,n,dataPtr_Sa_shum_18O(n)
                   write(nu_diag,F01)'import: istep, n, So_roce_HDO   = ',istep1,n,dataPtr_So_roce_HDO(n)
                   write(nu_diag,F01)'import: istep, n, So_roce_16O   = ',istep1,n,dataPtr_So_roce_16O(n)
                   write(nu_diag,F01)'import: istep, n, So_roce_18O   = ',istep1,n,dataPtr_So_roce_18O(n)
                end if
             end do
          end do
       end do
    end if

  end subroutine ice_import

  !===============================================================================

  subroutine ice_export( exportState, rc )

    ! input/output variables
    type(ESMF_State), intent(inout) :: exportState
    integer         , intent(out)   :: rc

    ! local variables
    real(kind=dbl_kind), pointer :: dataPtr_Si_imask(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_ifrac(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_ifrac_n(:,:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_t(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_avsdr(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_anidr(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_avsdf(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_anidf(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_u10(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_tref(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_qref(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_snowh(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_vsno(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_vice(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_logz0(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_qref_HDO(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_qref_16O(:)
    real(kind=dbl_kind), pointer :: dataPtr_Si_qref_18O(:)
    !
    real(kind=dbl_kind), pointer :: dataPtr_Faii_taux(:)
    real(kind=dbl_kind), pointer :: dataPtr_Faii_tauy(:)
    real(kind=dbl_kind), pointer :: dataPtr_Faii_lat(:)
    real(kind=dbl_kind), pointer :: dataPtr_Faii_sen(:)
    real(kind=dbl_kind), pointer :: dataPtr_Faii_lwup(:)
    real(kind=dbl_kind), pointer :: dataPtr_Faii_evap(:)
    real(kind=dbl_kind), pointer :: dataPtr_Faii_swnet(:)
    real(kind=dbl_kind), pointer :: dataPtr_Fioi_melth(:)
    real(kind=dbl_kind), pointer :: dataPtr_Fioi_swpen(:)
    real(kind=dbl_kind), pointer :: dataPtr_Fioi_swpen_ifrac_n(:,:)
    real(kind=dbl_kind), pointer :: dataPtr_Fioi_meltw(:)
    real(kind=dbl_kind), pointer :: dataPtr_Fioi_salt(:)
    real(kind=dbl_kind), pointer :: dataPtr_Fioi_taux(:)
    real(kind=dbl_kind), pointer :: dataPtr_Fioi_tauy(:)
    real(kind=dbl_kind), pointer :: dataPtr_Fioi_bcpho(:)
    real(kind=dbl_kind), pointer :: dataPtr_Fioi_bcphi(:)
    real(kind=dbl_kind), pointer :: dataPtr_Fioi_flxdst(:)
    real(kind=dbl_kind), pointer :: dataPtr_Fioi_meltw_HDO(:)
    real(kind=dbl_kind), pointer :: dataPtr_Fioi_meltw_16O(:)
    real(kind=dbl_kind), pointer :: dataPtr_Fioi_meltw_18O(:)
    real(kind=dbl_kind), pointer :: dataPtr_Faii_evap_HDO(:)
    real(kind=dbl_kind), pointer :: dataPtr_Faii_evap_16O(:)
    real(kind=dbl_kind), pointer :: dataPtr_Faii_evap_18O(:)
    !
    type(block)             :: this_block         ! block information for current block
    integer                 :: i, j, iblk, n      ! incides
    integer                 :: n2                 ! thickness category index
    integer                 :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    real    (kind=dbl_kind) :: workx, worky       ! tmps for converting grid
    integer (kind=int_kind) :: icells             ! number of ocean/ice cells
    logical                 :: flag
    integer                 :: dbrc
    !
    integer (kind=int_kind), dimension (nx_block*ny_block)           :: indxi  ! compressed indices in i
    integer (kind=int_kind), dimension (nx_block*ny_block)           :: indxj  ! compressed indices in i
    real    (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: Tsrf   ! surface temperature
    real    (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: tauxa  ! atmo/ice stress
    real    (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: tauya  ! atm/ice stress
    real    (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: tauxo  ! ice/ocean stress
    real    (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: tauyo  ! ice/ocean stress
    real    (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: ailohi ! fractional ice area
    character(len=*),parameter :: subname = 'ice_export'
    !-----------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    call State_getFldPtr(exportState, 'Si_imask',  dataPtr_Si_imask ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Si_ifrac',  dataPtr_Si_ifrac ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (State_FldChk(exportState, 'Si_ifrac_n')) then 
       call State_getFldPtr(exportState, 'Si_ifrac_n', dataPtr_Si_ifrac_n, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call State_getFldPtr(exportState, 'Si_t',	  dataPtr_Si_t	   ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Si_avsdr', dataPtr_Si_avsdr ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Si_anidr', dataPtr_Si_anidr ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Si_avsdf', dataPtr_Si_avsdf ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Si_anidf', dataPtr_Si_anidf ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Si_u10',	  dataPtr_Si_u10   ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Si_tref',  dataPtr_Si_tref  ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Si_qref',  dataPtr_Si_qref  ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Si_snowh', dataPtr_Si_snowh ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Si_vsno',  dataPtr_Si_vsno  ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Si_vice',  dataPtr_Si_vice  ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (State_FldChk(exportState, 'Si_logz0')) then 
       call State_getFldPtr(exportState, 'Si_logz0', dataPtr_Si_logz0 ,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if (State_FldChk(exportState, 'Si_qref_HDO')) then 
       call State_getFldPtr(exportState, 'Si_qref_HDO',	dataPtr_Si_qref_HDO ,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(exportState, 'Si_qref_16O',	dataPtr_Si_qref_16O ,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(exportState, 'Si_qref_18O',	dataPtr_Si_qref_18O ,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call State_getFldPtr(exportState, 'Faii_taux',	dataPtr_Faii_taux      ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Faii_tauy',	dataPtr_Faii_tauy      ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Faii_lat',	dataPtr_Faii_lat       ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Faii_sen',	dataPtr_Faii_sen       ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Faii_lwup',	dataPtr_Faii_lwup      ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Faii_evap',	dataPtr_Faii_evap      ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Faii_swnet',	dataPtr_Faii_swnet     ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Fioi_melth',	dataPtr_Fioi_melth     ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Fioi_swpen',	dataPtr_Fioi_swpen     ,rc=rc) 
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Fioi_meltw',	dataPtr_Fioi_meltw     ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Fioi_salt',	dataPtr_Fioi_salt      ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Fioi_taux',	dataPtr_Fioi_taux      ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Fioi_tauy',	dataPtr_Fioi_tauy      ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Fioi_bcpho',	dataPtr_Fioi_bcpho     ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Fioi_bcphi',	dataPtr_Fioi_bcphi     ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_getFldPtr(exportState, 'Fioi_flxdst',	dataPtr_Fioi_flxdst    ,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (State_FldChk(exportState, 'Fioi_meltw_HDO')) then 
       call State_getFldPtr(exportState, 'Fioi_meltw_HDO',	dataPtr_Fioi_meltw_HDO ,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(exportState, 'Fioi_meltw_16O',	dataPtr_Fioi_meltw_16O ,rc=rc) 
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(exportState, 'Fioi_meltw_18O',	dataPtr_Fioi_meltw_18O ,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(exportState, 'Faii_evap_HDO',	dataPtr_Faii_evap_HDO  ,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(exportState, 'Faii_evap_16O',	dataPtr_Faii_evap_16O  ,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_getFldPtr(exportState, 'Faii_evap_18O',	dataPtr_Faii_evap_18O  ,rc=rc)  
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if (State_FldChk(exportState, 'Fioi_swpen_ifrac_n')) then 
       call State_getFldPtr(exportState, 'Fioi_swpen_ifrac_n', dataPtr_Fioi_swpen_ifrac_n ,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !calculate ice thickness from aice and vice. Also
    !create Tsrf from the first tracer (trcr) in ice_state.F

    ailohi(:,:,:) = c0
    Tsrf(:,:,:)  = c0
    tauxa(:,:,:) = c0
    tauya(:,:,:) = c0
    tauxo(:,:,:) = c0
    tauyo(:,:,:) = c0

    !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky, this_block, ilo, ihi, jlo, jhi)
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo,jhi
          do i = ilo,ihi

             ! ice fraction
             ailohi(i,j,iblk) = min(aice(i,j,iblk), c1)

             ! surface temperature
             Tsrf(i,j,iblk)  = Tffresh + trcr(i,j,1,iblk)     !Kelvin (original ???)

             ! wind stress  (on POP T-grid:  convert to lat-lon)
             workx = strairxT(i,j,iblk)                             ! N/m^2
             worky = strairyT(i,j,iblk)                             ! N/m^2
             tauxa(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) - worky*sin(ANGLET(i,j,iblk))
             tauya(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) + workx*sin(ANGLET(i,j,iblk))

             ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
             workx = -strocnxT(i,j,iblk)                            ! N/m^2
             worky = -strocnyT(i,j,iblk)                            ! N/m^2
             tauxo(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) - worky*sin(ANGLET(i,j,iblk))
             tauyo(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) + workx*sin(ANGLET(i,j,iblk))

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    flag=.false.
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo,jhi
          do i = ilo,ihi
             if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
                flag = .true.
             endif
          end do
       end do
    end do
    if (flag) then
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          do j = jlo,jhi
             do i = ilo,ihi
                if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
                   write(nu_diag,*) &
                        ' (ice) send: ERROR ailohi < 0.0 ',i,j,ailohi(i,j,iblk)
                   call shr_sys_flush(nu_diag)
                endif
             end do
          end do
       end do
    endif

    ! Zero out fields with tmask for proper coupler accumulation in ice free areas
    call shr_nuopc_methods_State_reset(exportState, value=c0, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return


    ! Create the export state
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

             dataPtr_Si_ifrac(n) = ailohi(i,j,iblk)

             if (trim(grid_type) == 'latlon') then
                dataPtr_Si_imask(n) = ocn_gridcell_frac(i,j,iblk)
             else
                dataPtr_Si_imask(n) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
             end if

             if ( tmask(i,j,iblk) .and. ailohi(i,j,iblk) > c0 ) then

                !-------states--------------------
                dataPtr_Si_t(n)     = Tsrf(i,j,iblk)
                dataPtr_Si_avsdr(n) = alvdr(i,j,iblk)
                dataPtr_Si_anidr(n) = alidr(i,j,iblk)
                dataPtr_Si_avsdf(n) = alvdf(i,j,iblk)
                dataPtr_Si_anidf(n) = alidf(i,j,iblk)
                dataPtr_Si_u10(n)   = Uref(i,j,iblk)
                dataPtr_Si_tref(n)  = Tref(i,j,iblk)
                dataPtr_Si_qref(n)  = Qref(i,j,iblk)
                dataPtr_Si_snowh(n) = vsno(i,j,iblk) / ailohi(i,j,iblk)
                dataPtr_Si_vsno(n)  = vsno(i,j,iblk)
                dataPtr_Si_vice(n)  = vice(i,j,iblk)

                if (State_FldChk(exportState, 'Si_logz0')) then
                   if (Cdn_atm(i,j,iblk) > c0) then
                      dataPtr_Si_logz0(n) = log(zref)-(vonkar/sqrt(Cdn_atm(i,j,iblk)))
                   else
                      ! this should not happen but if it does, continue gracefully
                      write(nu_diag,*) trim(subname), ' WARNING: Cdn_atm error ',Cdn_atm(i,j,iblk),i,j,iblk
                      dataPtr_Si_logz0(n) = log(iceruf)
                   endif
                endif

                !--- a/i fluxes computed by ice
                dataPtr_Faii_taux(n)  = tauxa(i,j,iblk)
                dataPtr_Faii_tauy(n)  = tauya(i,j,iblk)
                dataPtr_Faii_lat(n)   = flat(i,j,iblk)
                dataPtr_Faii_sen(n)   = fsens(i,j,iblk)
                dataPtr_Faii_lwup(n)  = flwout(i,j,iblk)
                dataPtr_Faii_evap(n)  = evap(i,j,iblk)
                dataPtr_Faii_swnet(n) = fswabs(i,j,iblk)

                !--- i/o fluxes computed by ice
                dataPtr_Fioi_melth(n) = fhocn(i,j,iblk)
                dataPtr_Fioi_swpen(n) = fswthru(i,j,iblk) ! hf from melting
                dataPtr_Fioi_meltw(n) = fresh(i,j,iblk)   ! h2o flux from melting
                dataPtr_Fioi_salt(n)  = fsalt(i,j,iblk)   ! salt flux from melting
                dataPtr_Fioi_taux(n)  = tauxo(i,j,iblk)   ! stress n i/o zonal
                dataPtr_Fioi_tauy(n)  = tauyo(i,j,iblk)   ! stress n i/o meridional

                if (State_FldChk(exportState, 'Fioi_bcpho')) then
                   dataPtr_Fioi_bcpho(n)  = faero_ocn(i,j,1,iblk)  ! hydrophobic bc
                end if
                if (State_FldChk(exportState, 'Fioi_bcphi')) then
                   dataPtr_Fioi_bcphi(n)  = faero_ocn(i,j,2,iblk)  ! hydrophilic bc
                end if
                if (State_FldChk(exportState, 'Fioi_flxdst')) then
                   dataPtr_Fioi_flxdst(n)  = faero_ocn(i,j,3,iblk)  ! dust
                end if
                if (State_FldChk(exportState, 'Fioi_meltw_HDO')) then
                   dataPtr_Fioi_meltw_HDO(n) = fiso_ocn (i,j,1,iblk) !  Isotopes to ocean
                   dataPtr_Fioi_meltw_16O(n) = fiso_ocn (i,j,2,iblk) !  Isotopes to ocean
                   dataPtr_Fioi_meltw_18O(n) = fiso_ocn (i,j,3,iblk) !  Isotopes to ocean
                   dataPtr_Faii_evap_HDO(n)  = fiso_evap(i,j,1,iblk) !  Isotope evap to atm
                   dataPtr_Faii_evap_16O(n)  = fiso_evap(i,j,2,iblk) !  Isotope evap to atm
                   dataPtr_Faii_evap_18O(n)  = fiso_evap(i,j,3,iblk) !  Isotope evap to atm
                   dataPtr_Si_qref_HDO(n)    = Qref_iso(i,j,1,iblk)  !  Isotope qref to atm
                   dataPtr_Si_qref_16O(n)    = Qref_iso(i,j,2,iblk)  !  Isotope qref to atm
                   dataPtr_Si_qref_18O(n)    = Qref_iso(i,j,3,iblk)  !  Isotope qref to atm
                endif
             end if
          enddo    !i
       enddo    !j
    enddo        !iblk

    if (State_FldChk(exportState, 'Fioi_swpen_ifrac_n')) then
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

                ! ice fraction
                do n2 = 1, ncat
                   dataPtr_Si_ifrac_n(n,n2) = aicen_init(i,j,n2,iblk)
                enddo

                ! penetrative shortwave
                if ( tmask(i,j,iblk) .and. ailohi(i,j,iblk) > c0 ) then
                   do n2 = 1, ncat
                      dataPtr_Fioi_swpen_ifrac_n(n,n2) = fswthrun_ai(i,j,n2,iblk)
                   enddo
                else
                   !--- zero out pass-through fields over land for benefit of x2oacc fields in cpl hist files
                   do n2 = 1, ncat
                      dataPtr_Fioi_swpen_ifrac_n(n,n2) = c0
                   enddo
                end if

             enddo    !i
          enddo       !j
       enddo          !iblk
    end if ! send_i2x_per_cat

    !-----------------------------------------------------------------
    ! Debug output
    !-----------------------------------------------------------------

    if (debug > 0 .and. my_task==master_task) then
       n=0
       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo; ihi = this_block%ihi
          jlo = this_block%jlo; jhi = this_block%jhi
          do j = jlo, jhi
             do i = ilo, ihi
                n = n+1

                !--- ice states
                write(nu_diag,F01)'export: istep, n, Si_imask   = ',istep1,n,dataPtr_Si_imask(n)
                write(nu_diag,F01)'export: istep, n, Si_ifrac   = ',istep1,n,dataPtr_Si_ifrac(n)
                write(nu_diag,F01)'export: istep, n, Si_t       = ',istep1,n,dataPtr_Si_t(n)
                write(nu_diag,F01)'export: istep, n, Si_avsdr   = ',istep1,n,dataPtr_Si_avsdr(n)
                write(nu_diag,F01)'export: istep, n, Si_anidr   = ',istep1,n,dataPtr_Si_anidr(n)
                write(nu_diag,F01)'export: istep, n, Si_avsdf   = ',istep1,n,dataPtr_Si_avsdf(n)
                write(nu_diag,F01)'export: istep, n, Si_anidf   = ',istep1,n,dataPtr_Si_anidf(n)
                write(nu_diag,F01)'export: istep, n, Si_u10     = ',istep1,n,dataPtr_Si_u10(n)
                write(nu_diag,F01)'export: istep, n, Si_tref    = ',istep1,n,dataPtr_Si_tref(n)
                write(nu_diag,F01)'export: istep, n, Si_qref    = ',istep1,n,dataPtr_Si_qref(n)
                write(nu_diag,F01)'export: istep, n, Si_snowh   = ',istep1,n,dataPtr_Si_snowh(n)

                !--- a/i fluxes computed by ice
                write(nu_diag,F01)'export: istep, n, Faii_taux  = ',istep1,n,dataPtr_Faii_taux(n)
                write(nu_diag,F01)'export: istep, n, Faii_tauy  = ',istep1,n,dataPtr_Faii_tauy(n)
                write(nu_diag,F01)'export: istep, n, Faii_lat   = ',istep1,n,dataPtr_Faii_lat(n)
                write(nu_diag,F01)'export: istep, n, Faii_sen   = ',istep1,n,dataPtr_Faii_sen(n)
                write(nu_diag,F01)'export: istep, n, Faii_lwup  = ',istep1,n,dataPtr_Faii_lwup(n)
                write(nu_diag,F01)'export: istep, n, Faii_evap  = ',istep1,n,dataPtr_Faii_evap(n)
                write(nu_diag,F01)'export: istep, n, Faii_swnet = ',istep1,n,dataPtr_Faii_swnet(n)

                !--- i/o fluxes computed by ice
                write(nu_diag,F01)'export: istep, n, Fioi_melth = ',istep1,n,dataPtr_Fioi_melth(n)
                write(nu_diag,F01)'export: istep, n, Fioi_swpen = ',istep1,n,dataPtr_Fioi_swpen(n)
                write(nu_diag,F01)'export: istep, n, Fioi_meltw = ',istep1,n,dataPtr_Fioi_meltw(n)
                write(nu_diag,F01)'export: istep, n, Fioi_salt  = ',istep1,n,dataPtr_Fioi_salt(n)
                write(nu_diag,F01)'export: istep, n, Fioi_taux  = ',istep1,n,dataPtr_Fioi_taux(n)
                write(nu_diag,F01)'export: istep, n, Fioi_tauy  = ',istep1,n,dataPtr_Fioi_tauy(n)
                if (State_FldChk(exportState, 'Fioi_bchpo')) then
                   write(nu_diag,F01)'export: istep, n, Fioi_bcpho  = ',istep1,n,dataPtr_Fioi_bcpho(n)
                end if
                if (State_FldChk(exportState, 'Fioi_bchpi')) then
                   write(nu_diag,F01)'export: istep, n, Fioi_bcphi  = ',istep1,n,dataPtr_Fioi_bcpho(n)
                end if
                if (State_FldChk(exportState, 'Fioi_flxdst')) then
                   write(nu_diag,F01)'export: istep, n, Fioi_flxdst = ',istep1,n,dataPtr_Fioi_flxdst(n)
                end if
                if (State_FldChk(exportState, 'Fioi_HDO')) then
                   write(nu_diag,F01)'export: istep, n, Fioi_HDO      = ',istep1,n,dataPtr_Fioi_meltw_HDO(n)
                   write(nu_diag,F01)'export: istep, n, Fioi_16O      = ',istep1,n,dataPtr_Fioi_meltw_16O(n)
                   write(nu_diag,F01)'export: istep, n, Fioi_18O      = ',istep1,n,dataPtr_Fioi_meltw_18O(n)
                   write(nu_diag,F01)'export: istep, n, Faii_evap_HDO = ',istep1,n,dataPtr_Faii_evap_HDO(n)
                   write(nu_diag,F01)'export: istep, n, Faii_evap_16O = ',istep1,n,dataPtr_Faii_evap_16O(n)
                   write(nu_diag,F01)'export: istep, n, Faii_evap_18O = ',istep1,n,dataPtr_Faii_evap_18O(n)
                   write(nu_diag,F01)'export: istep, n, Si_qref_HDO   = ',istep1,n,dataPtr_Si_qref_HDO(n)
                   write(nu_diag,F01)'export: istep, n, Si_qref_16O   = ',istep1,n,dataPtr_Si_qref_16O(n)
                   write(nu_diag,F01)'export: istep, n, Si_qref_18O   = ',istep1,n,dataPtr_Si_qref_18O(n)
                end if
             end do
          end do
       end do
    end if

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ice_export

  !===============================================================================

  subroutine fldlist_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer             , intent(inout) :: num
    type(fld_list_type) , intent(inout) :: fldlist(:)
    character(len=*)    , intent(in)    :: stdname
    integer, optional   , intent(in)    :: ungridded_lbound
    integer, optional   , intent(in)    :: ungridded_ubound

    ! local variables
    character(len=*), parameter :: subname='(fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call shr_sys_abort(trim(subname)//": ERROR num > fldsMax "//trim(stdname))
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine fldlist_add

  !===============================================================================

  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU
    use ESMF  , only : ESMF_VM

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: dbrc
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(dshr_nuopc_mod:fld_list_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                     ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                     ungriddedUbound=(/fldlist(n)%ungridded_ubound/), rc=rc)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(dshr_nuopc_mod:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !===============================================================================

  subroutine State_GetFldPtr_1d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)             , intent(in)     :: State
    character(len=*)             , intent(in)     :: fldname
    real(ESMF_KIND_R8) , pointer , intent(inout)  :: fldptr(:)
    integer, optional            , intent(out)    :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    character(len=*),parameter :: subname='(ice_import_export:State_GetFldPtr_1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
  end subroutine State_GetFldPtr_1d

  subroutine State_GetFldPtr_2d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)             , intent(in)     :: State
    character(len=*)             , intent(in)     :: fldname
    real(ESMF_KIND_R8) , pointer , intent(inout)  :: fldptr(:,:)
    integer, optional            , intent(out)    :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    character(len=*),parameter :: subname='(ice_import_export:State_GetFldPtr_2d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
  end subroutine State_GetFldPtr_2d

  !===============================================================================

  logical function State_FldChk(State, fldname)
    ! ----------------------------------------------
    ! Determine if field is in state
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State) , intent(in)  :: State
    character(len=*) , intent(in)  :: fldname

    ! local variables
    type(ESMF_StateItem_Flag) :: itemType
    ! ----------------------------------------------

    call ESMF_StateGet(State, trim(fldname), itemType)

    State_FldChk = (itemType /= ESMF_STATEITEM_NOTFOUND)

  end function State_FldChk

end module ice_import_export
