module ice_latlon_grid

  use ice_kinds_mod

  implicit none
  private

  public :: latlon_grid

  ! Only relevant for lat-lon grids gridcell value of [1 - (land fraction)] (T-cell)
  real (kind=dbl_kind), allocatable, public :: ocn_gridcell_frac(:,:,:)

!=======================================================================
contains
!=======================================================================

  subroutine latlon_grid()

    ! Read in kmt file that matches CAM lat-lon grid and has single column functionality

    use ice_domain_size
    use ice_scam        , only : scmlat, scmlon, single_column
    use ice_constants   , only : c0, c1, pi, pi2, rad_to_deg, puny, p5, p25
    use ice_constants   , only : field_loc_center, field_type_scalar, radius
    use ice_exit        , only : abort_ice
    use ice_boundary    , only : ice_HaloUpdate, ice_HaloExtrapolate
    use ice_communicate , only : my_task, master_task
    use ice_blocks      , only : block, get_block, nx_block, ny_block
    use ice_domain      , only : blocks_ice, nblocks, halo_info, distrb_info
    use ice_fileunits   , only : nu_diag
    use ice_read_write  , only : ice_read_nc, ice_open_nc, ice_close_nc
    use ice_timers      , only : timer_bound, ice_timer_start, ice_timer_stop
    use ice_grid        , only : tlon, tlat, hm, tarea, ULON, ULAT, HTN, HTE, ANGLE, ANGLET
    use ice_grid        , only : uarea, uarear, tarear, tinyarea
    use ice_grid        , only : dxt, dyt, dxu, dyu, dyhx, dxhy, cyp, cxp, cym, cxm
    use ice_grid        , only : kmt_file,  makemask
    use netcdf

    integer (kind=int_kind) :: &
         i, j, iblk

    integer (kind=int_kind) :: &
         ni, nj, ncid, dimid, varid, ier

    character (len=char_len) :: &
         subname='latlongrid' ! subroutine name

    type (block) :: &
         this_block           ! block information for current block

    integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

    real (kind=dbl_kind) :: &
         closelat, &        ! Single-column latitude value
         closelon, &        ! Single-column longitude value
         closelatidx, &     ! Single-column latitude index to retrieve
         closelonidx        ! Single-column longitude index to retrieve

    integer (kind=int_kind) :: &
         start(2), &        ! Start index to read in
         count(2)           ! Number of points to read in

    integer (kind=int_kind) :: &
         start3(3), &        ! Start index to read in
         count3(3)           ! Number of points to read in

    integer (kind=int_kind) :: &
         status                ! status flag

    real (kind=dbl_kind), allocatable :: &
         lats(:),lons(:),pos_lons(:), glob_grid(:,:)  ! temporaries

    real (kind=dbl_kind) :: &
         pos_scmlon,&         ! temporary
         scamdata             ! temporary

    !-----------------------------------------------------------------
    ! - kmt file is actually domain fractional land file
    ! - Determine consistency of dimensions
    ! - Read in lon/lat centers in degrees from kmt file
    ! - Read in ocean from "kmt" file (1 for ocean, 0 for land)
    !-----------------------------------------------------------------

    ! Allocate ocn_gridcell_frac
    allocate(ocn_gridcell_frac(nx_block,ny_block,max_blocks))

    ! Determine dimension of domain file and check for consistency

    if (my_task == master_task) then
       call ice_open_nc(kmt_file, ncid)

       status = nf90_inq_dimid (ncid, 'ni', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=ni)
       status = nf90_inq_dimid (ncid, 'nj', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len=nj)
    end if

    ! Determine start/count to read in for either single column or global lat-lon grid
    ! If single_column, then assume that only master_task is used since there is only one task

    if (single_column) then

       ! Check for consistency
       if (my_task == master_task) then
          if ((nx_global /= 1).or. (ny_global /= 1)) then
             write(nu_diag,*) 'Because you have selected the column model flag'
             write(nu_diag,*) 'Please set nx_global=ny_global=1 in file'
             write(nu_diag,*) 'ice_domain_size.F and recompile'
             call abort_ice ('latlongrid: check nx_global, ny_global')
          endif
       end if

       ! Read in domain file for single column
       allocate(lats(nj))
       allocate(lons(ni))
       allocate(pos_lons(ni))
       allocate(glob_grid(ni,nj))

       start3=(/1,1,1/)
       count3=(/ni,nj,1/)
       status = nf90_inq_varid(ncid, 'xc' , varid)
       if (status /= nf90_noerr) call abort_ice (subname//' inq_varid xc')
       status = nf90_get_var(ncid, varid, glob_grid, start3, count3)
       if (status /= nf90_noerr) call abort_ice (subname//' get_var xc')
       do i = 1,ni
          lons(i) = glob_grid(i,1)
       end do

       status = nf90_inq_varid(ncid, 'yc' , varid)
       if (status /= nf90_noerr) call abort_ice (subname//' inq_varid yc')
       status = nf90_get_var(ncid, varid, glob_grid, start3, count3)
       if (status /= nf90_noerr) call abort_ice (subname//' get_var yc')
       do j = 1,nj
          lats(j) = glob_grid(1,j)
       end do

       ! convert lons array and scmlon to 0,360 and find index of value closest to 0
       ! and obtain single-column longitude/latitude indices to retrieve

       pos_lons(:)= mod(lons(:) + 360._dbl_kind,360._dbl_kind)
       pos_scmlon = mod(scmlon  + 360._dbl_kind,360._dbl_kind)
       start(1) = (MINLOC(abs(pos_lons-pos_scmlon),dim=1))
       start(2) = (MINLOC(abs(lats    -scmlat    ),dim=1))

       deallocate(lats)
       deallocate(lons)
       deallocate(pos_lons)
       deallocate(glob_grid)

       status = nf90_inq_varid(ncid, 'xc' , varid)
       if (status /= nf90_noerr) call abort_ice (subname//' inq_varid xc')
       status = nf90_get_var(ncid, varid, scamdata, start)
       if (status /= nf90_noerr) call abort_ice (subname//' get_var xc')
       TLON = scamdata
       status = nf90_inq_varid(ncid, 'yc' , varid)
       if (status /= nf90_noerr) call abort_ice (subname//' inq_varid yc')
       status = nf90_get_var(ncid, varid, scamdata, start)
       if (status /= nf90_noerr) call abort_ice (subname//' get_var yc')
       TLAT = scamdata
       status = nf90_inq_varid(ncid, 'area' , varid)
       if (status /= nf90_noerr) call abort_ice (subname//' inq_varid area')
       status = nf90_get_var(ncid, varid, scamdata, start)
       if (status /= nf90_noerr) call abort_ice (subname//' get_var are')
       tarea = scamdata
       status = nf90_inq_varid(ncid, 'mask' , varid)
       if (status /= nf90_noerr) call abort_ice (subname//' inq_varid mask')
       status = nf90_get_var(ncid, varid, scamdata, start)
       if (status /= nf90_noerr) call abort_ice (subname//' get_var mask')
       hm = scamdata
       status = nf90_inq_varid(ncid, 'frac' , varid)
       if (status /= nf90_noerr) call abort_ice (subname//' inq_varid frac')
       status = nf90_get_var(ncid, varid, scamdata, start)
       if (status /= nf90_noerr) call abort_ice (subname//' get_var frac')
       ocn_gridcell_frac = scamdata

    else

       ! Check for consistency
       if (my_task == master_task) then
          if (nx_global /= ni .and. ny_global /= nj) then
             call abort_ice ('latlongrid: ni,nj not equal to nx_global,ny_global')
          end if
       end if

       ! Read in domain file for global lat-lon grid
       call ice_read_nc(ncid, 1, 'xc'  , TLON             , diag=.true.)
       call ice_read_nc(ncid, 1, 'yc'  , TLAT             , diag=.true.)
       call ice_read_nc(ncid, 1, 'area', tarea            , diag=.true., &
            field_loc=field_loc_center,field_type=field_type_scalar)
       call ice_read_nc(ncid, 1, 'mask', hm               , diag=.true.)
       call ice_read_nc(ncid, 1, 'frac', ocn_gridcell_frac, diag=.true.)

    end if

    if (my_task == master_task) then
       call ice_close_nc(ncid)
    end if

    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
          do i = ilo, ihi
             ! Convert from degrees to radians
             TLON(i,j,iblk) = pi*TLON(i,j,iblk)/180._dbl_kind

             ! Convert from degrees to radians
             TLAT(i,j,iblk) = pi*TLAT(i,j,iblk)/180._dbl_kind

             ! Convert from radians^2 to m^2
             ! (area in domain file is in radians^2 and tarea is in m^2)
             tarea(i,j,iblk) = tarea(i,j,iblk) * (radius*radius)
          end do
       end do
    end do

    call ice_HaloUpdate (TLON  , halo_info, field_loc_center, field_type_scalar, fillValue=c1)
    call ice_HaloUpdate (TLAT  , halo_info, field_loc_center, field_type_scalar, fillValue=c1)
    call ice_HaloUpdate (tarea , halo_info, field_loc_center, field_type_scalar, fillValue=c1)
    call ice_HaloUpdate (hm    , halo_info, field_loc_center, field_type_scalar, fillValue=c1)

    !-----------------------------------------------------------------
    ! Calculate various geometric 2d arrays
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

  end subroutine latlon_grid

end module ice_latlon_grid
