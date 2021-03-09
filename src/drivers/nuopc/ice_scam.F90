module ice_scam

  use ice_kinds_mod

  implicit none

  ! single column control variables (only used for prescribed mode)

  logical              :: single_column ! true => single column mode
  logical              :: scol_valid    ! true => single column mask is 1
  real (kind=dbl_kind) :: scmlat        ! single column latitude (degrees)
  real (kind=dbl_kind) :: scmlon        ! single column longitude (degrees)
  real (kind=dbl_kind) :: scol_frac     ! single column ice fraction
  real (kind=dbl_kind) :: scol_mask     ! single column ice mask
  real (kind=dbl_kind) :: scol_area     ! single column ice area
  integer              :: scol_nj       ! nj size of single column domain file

end module ice_scam
