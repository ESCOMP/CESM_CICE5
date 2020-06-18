module ice_prescribed_mod

  ! !DESCRIPTION:
  ! The prescribed ice model reads in ice concentration data from a netCDF
  ! file.  Ice thickness, temperature, the ice temperature profile are
  ! prescribed.  Air/ice fluxes are computed to get surface temperature,
  ! Ice/ocean fluxes are set to zero, and ice dynamics are not calculated.
  ! Regridding and data cycling capabilities are included.

  use ESMF
  use shr_nl_mod       , only : shr_nl_find_group_name
  use dshr_strdata_mod , only : shr_strdata_type, shr_strdata_print
  use dshr_strdata_mod , only : shr_strdata_init_from_inline, shr_strdata_advance
  use dshr_methods_mod , only : dshr_fldbun_getfldptr
  use ice_communicate  , only : my_task, master_task, MPI_COMM_ICE
  use ice_kinds_mod    , only : dbl_kind, int_kind, log_kind, char_len_long
  use ice_fileunits    , only : nu_diag, nml_filename
  use ice_exit         , only : abort_ice
  use ice_domain_size  , only : ncat, nilyr, nslyr, max_blocks
  use ice_blocks       , only : nx_block, ny_block, block, get_block
  use ice_domain       , only : nblocks, distrb_info, blocks_ice
  use ice_grid         , only : TLAT, TLON, hm, tmask
  use ice_itd          , only : hin_max
  use ice_shr_methods  , only : chkerr
  use ice_broadcast    , only : broadcast_scalar
  use ice_constants

  implicit none
  private ! except

  ! public member functions:
  public :: ice_prescribed_init      ! initialize input data stream
  public :: ice_prescribed_run       ! get time slices and time interp
  public :: ice_prescribed_phys      ! set prescribed ice state and fluxes

  ! public data members:
  logical(kind=log_kind), public :: prescribed_ice      ! true if prescribed ice

  ! private data members:
  type(shr_strdata_type)         :: sdat                ! prescribed data stream
  real(kind=dbl_kind)            :: ice_cov(nx_block,ny_block,max_blocks) ! ice cover

  ! private parameters:
  ! note: heat capacity of sea ice, rhoi*C=rcpi+rLfidepressT*salinity/T^2
  real (kind=dbl_kind), parameter :: cp_sno = 0.0_dbl_kind        ! specific heat of snow         (J/kg/K)
  real (kind=dbl_kind), parameter :: rLfi = Lfresh*rhoi           ! latent heat of fusion ice     (J/m^3)
  real (kind=dbl_kind), parameter :: rLfs = Lfresh*rhos           ! latent heat of fusion snow    (J/m^3)
  real (kind=dbl_kind), parameter :: rLvi = Lvap*rhoi             ! latent heat of vapor*rhoice   (J/m^3)
  real (kind=dbl_kind), parameter :: rLvs = Lvap*rhos             ! latent heat of vapor*rhosno   (J/m^3)
  real (kind=dbl_kind), parameter :: rcpi = cp_ice*rhoi           ! heat capacity of fresh ice    (J/m^3)
  real (kind=dbl_kind), parameter :: rcps = cp_sno*rhos           ! heat capacity of snow         (J/m^3)
  real (kind=dbl_kind), parameter :: rcpidepressT = rcpi*depressT ! param for finding T(z) from q (J/m^3)
  real (kind=dbl_kind), parameter :: rLfidepressT = rLfi*depressT ! param for heat capacity       (J deg/m^3)

  character(*), parameter :: u_FILE_u = &
       __FILE__

!=======================================================================
contains
!===============================================================================

  subroutine ice_prescribed_init(compid, clock, mesh, rc)

    ! Prescribed ice initialization

    include 'mpif.h'

    ! input/output parameters
    integer(kind=int_kind) , intent(in)  :: compid
    type(ESMF_Clock)       , intent(in)  :: clock
    type(ESMF_Mesh)        , intent(in)  :: mesh
    integer                , intent(out) :: rc

    ! local parameters
    integer(kind=int_kind),parameter :: nFilesMaximum = 400 ! max number of files
    integer(kind=int_kind)           :: n, nFile, ierr
    integer(kind=int_kind)           :: nml_error ! namelist i/o error flag
    character(len=char_len_long)     :: stream_meshFile
    character(len=char_len_long)     :: stream_dataFiles(nFilesMaximum)
    character(len=char_len_long)     :: stream_varname
    integer(kind=int_kind)           :: stream_yearfirst   ! first year in stream to use
    integer(kind=int_kind)           :: stream_yearlast    ! last year in stream to use
    integer(kind=int_kind)           :: stream_yearalign   ! align stream_year_first
    integer(kind=int_kind)           :: nu_nml
    logical                          :: prescribed_ice_mode
    character(*),parameter           :: subName = "('ice_prescribed_init')"
    character(*),parameter           :: F00 = "('(ice_prescribed_init) ',4a)"
    character(*),parameter           :: F01 = "('(ice_prescribed_init) ',a,i0)"
    character(*),parameter           :: F02 = "('(ice_prescribed_init) ',2a,i0,)"
    !--------------------------------

    namelist /ice_prescribed_nuopc_nml/ &
         prescribed_ice_mode,           &
         stream_varname ,               &
         stream_meshfile,               &
         stream_datafiles,              &
         stream_yearalign,              &
         stream_yearfirst ,             &
         stream_yearlast                

    rc = ESMF_SUCCESS

    ! default values for namelist
    prescribed_ice_mode = .false. ! if true, prescribe ice
    stream_yearfirst    = 1       ! first year in  pice stream to use
    stream_yearlast     = 1       ! last  year in  pice stream to use
    stream_yearalign    = 1       ! align stream_year_first with this model year
    stream_varname      = 'ice_cov'
    stream_meshfile     = ' '
    stream_datafiles(:) = ' '

    ! read namelist on master task
    if (my_task == master_task) then
       open (newunit=nu_nml, file=nml_filename, status='old',iostat=nml_error)
       call shr_nl_find_group_name(nu_nml, 'ice_prescribed_nuopc_nml', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, ice_prescribed_nuopc_nml, iostat=nml_error)
          if (nml_error > 0) then
             write(nu_diag,F00) "ERROR: problem on read of ice_prescribed namelist in ice_prescribed_mod"
             call abort_ice(subName)
          endif
       endif
       close(nu_nml)
    end if

    ! --------------------------------------------------
    ! only do the following if prescribed ice mode is on
    ! --------------------------------------------------

    ! broadcast namelist input
    call broadcast_scalar(prescribed_ice_mode, master_task)

    ! set module variable
    prescribed_ice = prescribed_ice_mode

    if (prescribed_ice_mode) then

       call broadcast_scalar(stream_yearalign , master_task)
       call broadcast_scalar(stream_yearfirst , master_task)
       call broadcast_scalar(stream_yearlast  , master_task)
       call broadcast_scalar(stream_meshfile  , master_task)
       call broadcast_scalar(stream_varname   , master_task)
       call mpi_bcast(stream_dataFiles, len(stream_datafiles(1))*NFilesMaximum, MPI_CHARACTER, 0, MPI_COMM_ICE, ierr)

       nFile = 0
       do n = 1,nFilesMaximum
          if (stream_datafiles(n) /= ' ') nFile = nFile + 1
       end do

       if (my_task == master_task) then
          write(nu_diag,*) ' '
          write(nu_diag,F00) 'This is the prescribed ice coverage option.'
          write(nu_diag,F01) '  stream_yearfirst   = ',stream_yearfirst
          write(nu_diag,F01) '  stream_yearlast    = ',stream_yearlast
          write(nu_diag,F01) '  stream_yearalign   = ',stream_yearalign
          write(nu_diag,F00) '  stream_meshfile    = ',trim(stream_meshfile)
          write(nu_diag,F00) '  stream_varname     = ',trim(stream_varname)
          do n = 1,nFile
             write(nu_diag,F00) '  stream_datafiles   = ',trim(stream_dataFiles(n))
          end do
          write(nu_diag,*) ' '
       endif

       ! initialize sdat
       call shr_strdata_init_from_inline(sdat,       &
            my_task             = my_task,           &
            logunit             = nu_diag,           &
            compid              = compid ,           &
            model_clock         = clock ,            &
            model_mesh          = mesh,              &
            stream_meshfile     = stream_meshfile,   &
            stream_filenames    = stream_datafiles(1:nfile),  &
            stream_fldlistFile  = (/'ice_cov'/),     &
            stream_fldListModel = (/'ice_cov'/),     &
            stream_yearFirst    = stream_yearFirst,  &
            stream_yearLast     = stream_yearLast,   &
            stream_yearAlign    = stream_yearAlign , &
            stream_offset       = 0,                 &
            stream_taxmode      = 'cycle',           &
            rc                  = rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! print out sdat info
       if (my_task == master_task) then
          call shr_strdata_print(sdat,'ice coverage prescribed data')
       endif

       ! For one ice category, set hin_max(1) to something big
       if (ncat == 1) then
          hin_max(1) = 999._dbl_kind
       end if

    end if  ! end of if prescribed ice mode

  end subroutine ice_prescribed_init

  !=======================================================================
  subroutine ice_prescribed_run(mDateIn, secIn)

    ! Finds two time slices bounding current model time, remaps if necessary
    ! Interpolate to new ice coverage

    ! input/output parameters:
    integer(kind=int_kind), intent(in)  :: mDateIn  ! Current model date (yyyymmdd)
    integer(kind=int_kind), intent(in)  :: secIn    ! Elapsed seconds on model date

    ! local variables
    integer(kind=int_kind)       :: i,j,n,iblk       ! loop indices and counter
    integer(kind=int_kind)       :: ilo,ihi,jlo,jhi  ! beginning and end of physical domain
    type (block)                 :: this_block
    real(kind=dbl_kind)          :: aice_max         ! maximun ice concentration
    real(kind=dbl_kind), pointer :: dataptr(:)
    integer                      :: rc               ! ESMF return code
    character(*),parameter       :: subName = "('ice_prescribed_run')"
    character(*),parameter       :: F00 = "('(ice_prescribed_run) ',a,2g20.13)"
    logical                      :: first_time = .true.
    !------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Advance sdat stream
    call shr_strdata_advance(sdat, ymd=mDateIn, tod=SecIn, logunit=nu_diag, istr='cice_pice', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Get pointer for stream data that is time and spatially interpolate to model time and grid
    call dshr_fldbun_getFldPtr(sdat%pstrm(1)%fldbun_model, 'ice_cov', dataptr,  rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Fill in module ice_cov array
    ice_cov(:,:,:) = c0  ! This initializes ghost cells as well
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
             ice_cov(i,j,iblk) = dataptr(n)
          end do
       end do
    end do

    ! Check to see that ice concentration is in fraction, not percent
    if (first_time) then
       aice_max = maxval(ice_cov)
       if (aice_max > c10) then
          write(nu_diag,F00) "ERROR: Ice conc data must be in fraction, aice_max= ", aice_max
          rc = ESMF_FAILURE
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if
       end if
       first_time = .false.
    end if

    ! Set prescribed ice state and fluxes
    call ice_prescribed_phys()

  end subroutine ice_prescribed_run

  !=======================================================================
  subroutine ice_prescribed_phys()

    ! Set prescribed ice state using input ice concentration;
    ! set surface ice temperature to atmospheric value; use
    ! linear temperature gradient in ice to ocean temperature.

    use ice_flux
    use ice_state
    use ice_itd, only  : aggregate
    use ice_dyn_evp

    !----- Local ------
    integer(kind=int_kind) :: layer    ! level index
    integer(kind=int_kind) :: nc       ! ice category index
    integer(kind=int_kind) :: i,j,k    ! longitude, latitude and level indices
    integer(kind=int_kind) :: iblk

    real(kind=dbl_kind) :: slope     ! diff in underlying ocean tmp and ice surface tmp
    real(kind=dbl_kind) :: Ti        ! ice level temperature
    real(kind=dbl_kind) :: Tmlt      ! ice level melt temperature
    real(kind=dbl_kind) :: qin_save(nilyr)
    real(kind=dbl_kind) :: qsn_save(nslyr)
    real(kind=dbl_kind) :: hi        ! ice prescribed (hemispheric) ice thickness
    real(kind=dbl_kind) :: hs        ! snow thickness
    real(kind=dbl_kind) :: zn        ! normalized ice thickness
    real(kind=dbl_kind) :: salin(nilyr)  ! salinity (ppt)

    real(kind=dbl_kind), parameter :: nsal    = 0.407_dbl_kind
    real(kind=dbl_kind), parameter :: msal    = 0.573_dbl_kind
    real(kind=dbl_kind), parameter :: saltmax = 3.2_dbl_kind   ! max salinity at ice base (ppm)
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! Set ice cover over land to zero, not sure if this should be
    ! be done earier, before time/spatial interp??????
    !-----------------------------------------------------------------
    do iblk = 1,nblocks
       do j = 1,ny_block
          do i = 1,nx_block
             if (tmask(i,j,iblk)) then
                if (ice_cov(i,j,iblk) .lt. eps04) ice_cov(i,j,iblk) = c0
                if (ice_cov(i,j,iblk) .gt. c1)    ice_cov(i,j,iblk) = c1
             else
                ice_cov(i,j,iblk) = c0
             end if
          enddo
       enddo
    enddo

    do iblk = 1,nblocks
       do j = 1,ny_block
          do i = 1,nx_block

             if (tmask(i,j,iblk)) then   ! Over ocean points

                !--------------------------------------------------------------
                ! Place ice where ice concentration > .0001
                !--------------------------------------------------------------

                if (ice_cov(i,j,iblk) >= eps04) then

                   hi = 0.0_dbl_kind
                   !----------------------------------------------------------
                   ! Set ice thickness in each hemisphere
                   !----------------------------------------------------------
                   if(TLAT(i,j,iblk)*rad_to_deg > 40.0_dbl_kind) then
                      hi  = 2.0_dbl_kind
                   else if(TLAT(i,j,iblk)*rad_to_deg < -40.0_dbl_kind) then
                      hi  = 1.0_dbl_kind
                   end if

                   !----------------------------------------------------------
                   ! All ice in appropriate thickness category
                   !----------------------------------------------------------
                   do nc = 1,ncat

                      if(hin_max(nc-1) < hi .and. hi < hin_max(nc)) then

                         if (aicen(i,j,nc,iblk) > c0) then
                            hs = vsnon(i,j,nc,iblk) / aicen(i,j,nc,iblk)
                         else
                            hs = c0
                         endif

                         aicen(i,j,nc,iblk) = ice_cov(i,j,iblk)
                         vicen(i,j,nc,iblk) = hi*aicen(i,j,nc,iblk)
                         vsnon(i,j,nc,iblk) = hs*aicen(i,j,nc,iblk)

                         !---------------------------------------------------------
                         ! make linear temp profile and compute enthalpy
                         !---------------------------------------------------------

                         if (abs(trcrn(i,j,nt_qice,nc,iblk)) < puny) then

                            if (aice(i,j,iblk) < puny) &
                                 trcrn(i,j,nt_Tsfc,nc,iblk) = Tf(i,j,iblk)

                            slope = Tf(i,j,iblk) - trcrn(i,j,nt_Tsfc,nc,iblk)
                            do k = 1, nilyr
                               zn = (real(k,kind=dbl_kind)-p5) / real(nilyr,kind=dbl_kind)
                               Ti = trcrn(i,j,nt_Tsfc,nc,iblk) + slope*zn
                               salin(k) = (saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
                               Tmlt = -salin(k)*depressT
                               trcrn(i,j,nt_sice+k-1,nc,iblk) = salin(k)
                               trcrn(i,j,nt_qice+k-1,nc,iblk) =                      &
                                    -(rhoi * (cp_ice*(Tmlt-Ti) &
                                    + Lfresh*(c1-Tmlt/Ti) - cp_ocn*Tmlt))
                            enddo

                            do k=1,nslyr
                               trcrn(i,j,nt_qsno+k-1,nc,iblk) =                      &
                                    -rhos*(Lfresh - cp_ice*trcrn(i,j,nt_Tsfc,nc,iblk))
                            enddo

                         endif  ! aice < puny
                      end if    ! hin_max
                   enddo        ! ncat
                else
                   trcrn(i,j,nt_Tsfc,:,iblk) = Tf(i,j,iblk)
                   aicen(i,j,:,iblk) = c0
                   vicen(i,j,:,iblk) = c0
                   vsnon(i,j,:,iblk) = c0
                   trcrn(i,j,nt_sice:nt_sice+nilyr-1,:,iblk) = c0
                   trcrn(i,j,nt_qice:nt_qice+nilyr-1,:,iblk) = c0
                   trcrn(i,j,nt_qsno:nt_qsno+nslyr-1,:,iblk) = c0
                end if          ! ice_cov >= eps04
             end if             ! tmask
          enddo                 ! i
       enddo                 ! j

       !--------------------------------------------------------------------
       ! compute aggregate ice state and open water area
       !--------------------------------------------------------------------
       call aggregate (nx_block,          ny_block,             &
            aicen(:,:,:,iblk),                       &
            trcrn(:,:,1:ntrcr,:,iblk),               &
            vicen(:,:,:,iblk), vsnon(:,:,  :,iblk),  &
            aice (:,:,  iblk),                       &
            trcr (:,:,1:ntrcr,  iblk),               &
            vice (:,:,  iblk), vsno (:,:,    iblk),  &
            aice0(:,:,  iblk), tmask(:,:,    iblk),  &
            ntrcr, trcr_depend(1:ntrcr))

    enddo                 ! iblk

    do iblk = 1, nblocks
       do j = 1, ny_block
          do i = 1, nx_block
             aice_init(i,j,iblk) = aice(i,j,iblk)
          enddo
       enddo
    enddo

    !--------------------------------------------------------------------
    ! set non-computed fluxes, ice velocities, ice-ocn stresses to zero
    !--------------------------------------------------------------------

    frzmlt    (:,:,:) = c0
    uvel      (:,:,:) = c0
    vvel      (:,:,:) = c0
    strocnxT  (:,:,:) = c0
    strocnyT  (:,:,:) = c0

    !-----------------------------------------------------------------
    ! other atm and ocn fluxes
    !-----------------------------------------------------------------
    call init_flux_atm
    call init_flux_ocn

  end subroutine ice_prescribed_phys

end module ice_prescribed_mod
