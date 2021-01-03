module ice_mesh_mod

  use ESMF
  use ice_kinds_mod    , only : dbl_kind, int_kind
  use ice_constants    , only : rad_to_deg
  use ice_grid         , only : tlon, tlat, hm, tarea, ULON, ULAT, grid_type
  use ice_domain       , only : nblocks, blocks_ice, distrb_info
  use ice_blocks       , only : block, get_block, nx_block, ny_block, nblocks_x, nblocks_y
  use ice_blocks       , only : nblocks_tot, get_block_parameter
  use ice_scam         , only : scmlat, scmlon, single_column
  use ice_shr_methods  , only : chkerr

  implicit none
  private

  public  :: ice_mesh_init

  private :: ice_mesh_set_distgrid 
  private :: ice_mesh_set_mask

  ! Only relevant for lat-lon grids gridcell value of [1 - (land fraction)] (T-cell)
  real (kind=dbl_kind), allocatable, public :: ocn_gridcell_frac(:,:,:)

!=======================================================================
contains
!=======================================================================

  subroutine ice_mesh_init(ice_meshfile, ice_maskfile, ice_mesh, rc)

    ! Create the CICE mesh

    ! input/output parameters
    character(len=*) , intent(in)  :: ice_meshfile
    character(len=*) , intent(in)  :: ice_maskfile
    type(ESMF_Mesh)  , intent(out) :: ice_mesh
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_DistGrid) :: distGrid
    integer             :: n,c,g,i,j,m        ! indices
    integer             :: iblk, jblk         ! indices
    integer             :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    integer             :: spatialDim
    integer             :: numOwnedElements
    real(R8), pointer   :: ownedElemCoords(:)
    real(r8), pointer   :: lat(:), latmesh(:)
    real(r8), pointer   :: lon(:), lonmesh(:)
    !---------------------------------------------------

    ! Allocate module variable ocn_gridcell_frac
    allocate(ocn_gridcell_frac(nx_block,ny_block,max_blocks))

    ! read in the ice mesh and create the mesh on the cice distribution
    ice_mesh = ESMF_MeshCreate(filename=trim(ice_meshfile), fileformat=ESMF_FILEFORMAT_ESMFMESH, &
         elementDistGrid=distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) then
       write(nu_diag,*)'mesh file for cice domain is ',trim(cvalue)
    end if

    ! Determine the model distgrid using the decomposition obtained in
    ! call to init_grid1 called from cice_init1
    call ice_mesh_set_distgrid(distgrid, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (trim(grid_type) == 'latlon') then

       call ice_set_latlon_grid(ice_mesh, ice_maskfile, 

    else

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
100          format('ERROR: CICE  n, lonmesh(n), lon(n), diff_lon = ',i6,2(f21.13,3x),d21.5)
             !call shr_sys_abort()
          end if
          if (abs(latMesh(n) - lat(n)) > 1.e-1) then
             !write(6,101)n,latMesh(n),lat(n), abs(latMesh(n)-lat(n))
101          format('ERROR: CICE n, latmesh(n), lat(n), diff_lat = ',i6,2(f21.13,3x),d21.5)
             !call shr_sys_abort()
          end if
       end do

    end if

    ! deallocate memory
    deallocate(ownedElemCoords)
    deallocate(lon, lonMesh)
    deallocate(lat, latMesh)

  end subroutine ice_mesh_init

  !=======================================================================
  subroutine ice_mesh_set_distgrid(distgrid, rc)

    ! Determine the global index space needed for the distgrid

    ! input/output variables
    type(ESMF_DistGrid) , intent(out) :: distgrid
    integer             , intent(out) :: rc

    ! local variables
    integer     :: n,c,g,i,j,m        ! indices
    integer     :: iblk, jblk         ! indices
    integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    integer     :: lsize              ! local size of coupling array
    integer     :: ig, jg             ! indices
    integer     :: localPet
    type(block) :: this_block         ! block information for current block
    real(R8)    :: diff_lon
    integer     :: num_elim_global
    integer     :: num_elim_local
    integer     :: num_elim
    integer     :: num_ice
    integer     :: num_elim_gcells    ! local number of eliminated gridcells
    integer     :: num_elim_blocks    ! local number of eliminated blocks
    integer     :: num_total_blocks
    integer     :: my_elim_start, my_elim_end
    integer , allocatable   :: gindex(:)
    integer , allocatable   :: gindex_ice(:)
    integer , allocatable   :: gindex_elim(:)
    character(len=*), parameter :: subname = ' ice_mesh_set_distgrid: '
    !----------------------------------------------------------------

    rc = ESMF_SUCCESS

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
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(gindex_ice)
    deallocate(gindex)

  end subroutine ice_mesh_set_distgrid

  !===============================================================================
  subroutine ice_mesh_set_mask(ice_mesh, ice_maskfile, ice_mask, ice_frac, rc)

    ! input/out variables
    type(ESMF_Mesh)     , intent(in)  :: ice_mesh
    character(len=*)    , intent(in)  :: ice_maskfile
    integer , pointer   , intent(out) :: ice_mask(:)
    real(r8), pointer   , intent(out) :: ice_frac(:)
    integer             , intent(out) :: rc

    ! local variables:
    type(ESMF_Mesh)        :: mesh_mask
    type(ESMF_Field)       :: field_mask
    type(ESMF_Field)       :: field_dst
    type(ESMF_RouteHandle) :: rhandle
    integer                :: srcMaskValue = 0
    integer                :: dstMaskValue = -987987 ! spval for RH mask values
    integer                :: srcTermProcessing_Value = 0
    logical                :: checkflag = .false.
    real(r8) , pointer     :: mask_src(:) ! on mesh created from ice_maskfile
    real(r8) , pointer     :: dataptr1d(:)
    type(ESMF_DistGrid)    :: distgrid_mask
    type(ESMF_Array)       :: elemMaskArray
    integer                :: lsize_mask, lsize_dst
    integer                :: n, spatialDim
    real(r8)               :: fminval = 0.001_r8
    real(r8)               :: fmaxval = 1._r8
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    mesh_mask = ESMF_MeshCreate(trim(ice_maskfile), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MeshGet(ice_mesh, spatialDim=spatialDim, numOwnedElements=lsize_dst, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ice_mask(lsize_dst))
    allocate(ice_frac(lsize_dst))

    ! create fields on source and destination meshes
    field_mask = ESMF_FieldCreate(mesh_mask, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    field_dst = ESMF_FieldCreate(ice_mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create route handle to map source mask (assume ocean) to destination mesh (assume atm/lnd)
    call ESMF_FieldRegridStore(field_mask, field_dst, routehandle=rhandle, &
         srcMaskValues=(/srcMaskValue/), dstMaskValues=(/dstMaskValue/), &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, normType=ESMF_NORMTYPE_DSTAREA, &
         srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! fill in values for field_mask with mask on source mesh
    call ESMF_MeshGet(mesh_mask, elementdistGrid=distgrid_mask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distgrid_mask, localDe=0, elementCount=lsize_mask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(mask_src(lsize_mask))
    elemMaskArray = ESMF_ArrayCreate(distgrid_mask, mask_src, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! The following call fills in the values of mask_src
    call ESMF_MeshGet(mesh_mask, elemMaskArray=elemMaskArray, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! The following call fills in the values of field_mask
    call ESMF_FieldGet(field_mask, farrayptr=dataptr1d, rc=rc)
    dataptr1d(:) = mask_src(:)

    ! map source mask to destination mesh - to obtain destination mask and frac
    call ESMF_FieldRegrid(field_mask, field_dst, routehandle=rhandle, &
         termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! now determine ice_mask and ice_frac 
    call ESMF_MeshGet(ice_mesh, spatialDim=spatialDim, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_dst, farrayptr=dataptr1d, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    ice_frac(:) = dataptr1d(:)
    do n = 1,lsize_dst
       if (ice_frac(n) > fmaxval) then
          ice_frac(n) = 1._r8
       end if
       if (ice_frac(n) < fminval) then
          ice_frac(n) = 0._r8
       end if
       if (ice_frac(n) /= 0._r8) then
          ice_mask(n) = 1
       else
          ice_mask(n) = 0
       end if
    enddo

    ! reset the model mesh mask
    call ESMF_MeshSet(ice_mesh, elementMask=ice_mask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! deallocate memory
    call ESMF_RouteHandleDestroy(rhandle, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldDestroy(field_mask, rc=rc) 
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldDestroy(field_dst, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    deallocate(mask_src)

  end subroutine ice_mesh_set_mask

!=======================================================================
  subroutine ice_set_latlon_grid(ice_mesh, ice_maskfile, rc)

    ! obtain the model mask by mapping the mesh created by reading in the model_maskfile to the
    ! model mesh and then reset the model mesh mask
    call ice_mesh_set_mask(ice_mesh, ice_maskfile, ice_mask, ice_frac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set hm and ocn_gridcell_frac
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
             hm(i,j,iblk) = ice_mask(n)
             ocn_gridcell_frac(i,j,iblk) = ice_frac(n)
          enddo
       enddo
    enddo

    ! obtain mesh lats, lons, areas
    call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    allocate(lonmesh(numOwnedElements))
    allocate(latmesh(numOwnedElements))
    call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO: Set area
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
          do i = ilo, ihi
             ! Convert from radians^2 to m^2 (area in domain file is in radians^2 and tarea is in m^2)
             tarea(i,j,iblk) = tarea(i,j,iblk) * (radius*radius)
          end do
       end do
    end do

    ! set internally generated cice lats and lons
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
             tlon(i,j,iblk) = lonmesh(n)/rad_to_deg 
             tlat(i,j,iblk) = latmesh(n)/rad_to_deg 
          enddo
       enddo
    enddo

    call ice_ HaloUpdate (TLON  , halo_info, field_loc_center, field_type_scalar, fillValue=c1)
    call ice_ HaloUpdate (TLAT  , halo_info, field_loc_center, field_type_scalar, fillValue=c1)
    call ice_ HaloUpdate (tarea , halo_info, field_loc_center, field_type_scalar, fillValue=c1)
    call ice_ HaloUpdate (hm    , halo_info, field_loc_center, field_type_scalar, fillValue=c1)

    !-----------------------------------------------------------------
    ! CALCULATE various geometric 2d arrays
    ! The U grid (velocity) is not used when run with sequential CAM
    ! because we only use thermodynamic sea ice.  However, ULAT is used
    ! in the default initialization of CICE so we calculate it here as
    ! a "dummy" so that CICE will initialize with ice.  If a no ice
    ! initialization is OK (or desired) this can be commented out and
    ! ULAT will remain 0 as specified above.  ULAT is located at the
    ! NE corner of the grid cell, TLAT at the center, so here ULAT is
    ! hacked by adding half the latitudinal spacing (in radians) to TLAT.
    !-----------------------------------------------------------------

    ANGLET(:,:,:) = c0

    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
          do i = ilo, ihi

             if (ny_global == 1) then
                uarea(i,j,iblk)  = tarea(i,j,  iblk)
             else
                uarea(i,j,iblk)  = p25*  &
                     (tarea(i,j,  iblk) + tarea(i+1,j,  iblk) &
                     + tarea(i,j+1,iblk) + tarea(i+1,j+1,iblk))
             endif
             tarear(i,j,iblk)   = c1/tarea(i,j,iblk)
             uarear(i,j,iblk)   = c1/uarea(i,j,iblk)
             tinyarea(i,j,iblk) = puny*tarea(i,j,iblk)

             if (single_column) then
                ULAT  (i,j,iblk) = TLAT(i,j,iblk)+(pi/nj)
             else
                if (ny_global == 1) then
                   ULAT  (i,j,iblk) = TLAT(i,j,iblk)
                else
                   ULAT  (i,j,iblk) = TLAT(i,j,iblk)+(pi/ny_global)
                endif
             endif
             ULON  (i,j,iblk) = c0
             ANGLE (i,j,iblk) = c0

             HTN   (i,j,iblk) = 1.e36_dbl_kind
             HTE   (i,j,iblk) = 1.e36_dbl_kind
             dxt   (i,j,iblk) = 1.e36_dbl_kind
             dyt   (i,j,iblk) = 1.e36_dbl_kind
             dxu   (i,j,iblk) = 1.e36_dbl_kind
             dyu   (i,j,iblk) = 1.e36_dbl_kind
             dxhy  (i,j,iblk) = 1.e36_dbl_kind
             dyhx  (i,j,iblk) = 1.e36_dbl_kind
             cyp   (i,j,iblk) = 1.e36_dbl_kind
             cxp   (i,j,iblk) = 1.e36_dbl_kind
             cym   (i,j,iblk) = 1.e36_dbl_kind
             cxm   (i,j,iblk) = 1.e36_dbl_kind
          enddo
       enddo
    enddo

    call ice_HaloUpdate (ULAT, halo_info, field_loc_center, field_type_scalar, fillValue=c1)

    call makemask

  end subroutine ice_set_latlon_grid

end module ice_mesh_mod
